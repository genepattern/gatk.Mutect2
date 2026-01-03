#!/usr/bin/env python3

import os
import sys
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ---------------------------
# tolerant option parsing
# ---------------------------

def _next_is_value(argv, i):
    return (i + 1 < len(argv)) and (not argv[i + 1].startswith("-"))

def parse_args_tolerant(argv):
    
    spec = {
        # Reference
        "-R": ("ref", True), "--ref": ("ref", True),

        # Tumor BAM + sample
        "-T": ("tumor_bam", True), "--tumor-bam": ("tumor_bam", True),
        "-I": ("tumor_bam", True),  # GP style
        "-Ts": ("tumor_sample", True),
        "--tumor-sample": ("tumor_sample", True),
        "--tumor_sample": ("tumor_sample", True),  # GP style

        # Normal BAM + sample
        "-N": ("normal_bam", True), "--normal-bam": ("normal_bam", True),
        "-Ns": ("normal_sample", True),
        "--normal-sample": ("normal_sample", True),
        "--normal_sample": ("normal_sample", True),  # GP style

        # Resources
        "-Pon": ("pon_vcf", True), "--pon-vcf": ("pon_vcf", True),
        "-G": ("germline_resource", True),
        "--germline-resource": ("germline_resource", True),
        "--germline_resource": ("germline_resource", True),  # GP style

        # Intervals
        "-L": ("intervals", True), "--intervals": ("intervals", True),

        # Output
        "-O": ("out", True), "--out": ("out", True), "--out-prefix": ("out", True),

        # GATK / execution
        "-g": ("gatk", True), "--gatk": ("gatk", True),
        "-P": ("pairhmm_threads", True), "--pairhmm-threads": ("pairhmm_threads", True),

        "-C": ("contigs", True), "--contigs": ("contigs", True),

        "--tmp-dir": ("tmp_dir", True),
        "--dry-run": ("dry_run", False),

        # Scatter / parallel
        "-S": ("split5", True), "--split5": ("split5", True),
        "--max-parallel": ("max_parallel", True),
        "-J": ("cpu_count", True), "--cpu-count": ("cpu_count", True),
    }
    opt = {
        "ref": None, "tumor_bam": None, "tumor_sample": None,
        "normal_bam": None, "normal_sample": None,
        "pon_vcf": None, "germline_resource": None,
        "intervals": None, "out": None,
        "gatk": "gatk", "pairhmm_threads": "4",
        "contigs": None, "tmp_dir": None, "dry_run": False,
        "split5": "false",
        "max_parallel": None,
        "cpu_count": None,
    }
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in spec:
            key, expects_val = spec[tok]
            if expects_val:
                if _next_is_value(argv, i):
                    opt[key] = argv[i + 1]
                    i += 2
                else:
                    # flag present but no value
                    if key == "split5":
                        opt[key] = "true"
                    else:
                        opt[key] = None
                    i += 1
            else:
                opt[key] = True
                i += 1
        else:
            i += 1
    return opt

# ---------------------------
# helpers
# ---------------------------

def is_bgzf(path: str) -> bool:
    """
    Detect BGZF vs regular gzip by inspecting the gzip extra subfields.
    BGZF is gzip with FEXTRA and subfield 'BC'.
    """
    p = Path(path)
    if not p.exists() or p.is_dir():
        return False

    try:
        with open(p, "rb") as fh:
            hdr = fh.read(64)
        if len(hdr) < 18:
            return False

        # gzip magic
        if hdr[0:2] != b"\x1f\x8b":
            return False
        # compression method must be deflate
        if hdr[2] != 8:
            return False

        flg = hdr[3]
        FEXTRA = 0x04
        if (flg & FEXTRA) == 0:
            return False

        # skip basic header: 10 bytes
        xlen = int.from_bytes(hdr[10:12], "little")
        # Extra field begins at byte 12, length xlen
        extra = hdr[12:12 + xlen]
        i = 0
        while i + 4 <= len(extra):
            si1 = extra[i:i+1]
            si2 = extra[i+1:i+2]
            slen = int.from_bytes(extra[i+2:i+4], "little")
            sub = extra[i+4:i+4+slen]
            if si1 == b"B" and si2 == b"C":
                return True
            i += 4 + slen
        return False
    except Exception:
        return False


def ensure_uncompressed_fasta(ref: str, dry: bool = False) -> str:
    """
    If ref is a regular gzip FASTA (not BGZF), decompress to .fa and return that path.
    If ref is BGZF or uncompressed, return ref unchanged.
    """
    if not ref:
        return ref

    ref_path = Path(ref)

    # If already uncompressed, keep it
    if ref_path.suffix != ".gz":
        return ref

    # If BGZF, keep as-is
    if is_bgzf(ref):
        print(f"[info] Reference is BGZF; keeping as-is: {ref}")
        return ref

    # Regular gzip -> decompress to .fa next to it
    out_path = Path(str(ref_path)[:-3])  # strip ".gz" -> ".fa"
    if out_path.exists():
        print(f"[info] Using existing uncompressed reference: {out_path}")
        return str(out_path)

    print(f"[info] Reference is plain gzip; decompressing to: {out_path}")
    cmd = shlex.join(["bash", "-lc", f"gzip -dc {shlex.quote(str(ref_path))} > {shlex.quote(str(out_path))}"])
    rc = run(cmd, dry=dry)
    if rc != 0 and not dry:
        fail(f"Failed to decompress reference FASTA: {ref_path} -> {out_path}")

    return str(out_path)

def fail(msg, code=2):
    print(f"[error] {msg}", file=sys.stderr)
    sys.exit(code)

def must(path, label):
    """For file-path like arguments only."""
    if not path:
        fail(f"Missing required {label}.")
    if not Path(path).exists():
        fail(f"{label} not found: {path}")

def run(cmd, dry=False):
    print(f"[cmd] {cmd}", flush=True)
    if dry:
        return 0
    p = subprocess.run(cmd, shell=True)
    return p.returncode

def ensure_vcf_index(vcffile, gatk, dry=False):
    """
    Ensure a VCF.gz has a .tbi index using GATK IndexFeatureFile.
    """
    if not vcffile:
        return 0
    tbi = vcffile + ".tbi"
    if Path(tbi).exists() or dry:
        return 0
    cmd = shlex.join([gatk, "IndexFeatureFile", "-I", vcffile])
    rc = run(cmd, dry)
    if rc != 0 and not dry:
        fail(f"Failed to create VCF index (.tbi) for: {vcffile}")
    return rc

def ensure_reference_indices(ref, gatk, dry=False):
    """
    Ensure reference FASTA has .fai and .dict.
    - .fai via: samtools faidx ref
    - .dict via: gatk CreateSequenceDictionary -R ref -O ref_base.dict
    """
    if not ref:
        return
    ref_path = Path(ref)
    fai_path = Path(str(ref) + ".fai")
    dict_path = ref_path.with_suffix(".dict")

    # .fai
    if not fai_path.exists() and not dry:
        cmd = shlex.join(["samtools", "faidx", str(ref_path)])
        rc = run(cmd, dry)
        if rc != 0:
            fail(
                f"Failed to create FASTA index (.fai) for: {ref_path}. "
                f"Make sure samtools is installed and on PATH."
            )
    else:
        print(f"[info] Using existing FASTA index: {fai_path}")

    # .dict
    if not dict_path.exists() and not dry:
        cmd = shlex.join([gatk, "CreateSequenceDictionary", "-R", str(ref_path), "-O", str(dict_path)])
        rc = run(cmd, dry)
        if rc != 0:
            fail(f"Failed to create sequence dictionary (.dict) for: {ref_path}")
    else:
        print(f"[info] Using existing sequence dictionary: {dict_path}")

def ensure_bam_index(bam, dry=False):
    """
    Ensure BAM has a .bai index using samtools index.
    """
    if not bam:
        return
    bam_path = Path(bam)
    bai_path = Path(str(bam_path) + ".bai")
    if bai_path.exists() or dry:
        print(f"[info] Using existing BAM index: {bai_path}")
        return
    cmd = shlex.join(["samtools", "index", str(bam_path)])
    rc = run(cmd, dry)
    if rc != 0:
        fail(
            f"Failed to create BAM index (.bai) for: {bam_path}. "
            f"Make sure samtools is installed and on PATH."
        )

def normalized_contigs(contigs_opt):
    if not contigs_opt:
        # default human contigs
        return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    txt = contigs_opt.replace(",", " ")
    return [c for c in txt.split() if c]

def truthy(x):
    return str(x).strip().lower() in {"1", "true", "yes", "y", "on"}

def derive_prefix_and_merged_path(out):
    if not out:
        fail("Missing output path/prefix (-O/--out/--out-prefix).")
    s = str(out)
    if s.endswith(".vcf.gz") or s.endswith(".vcf"):
        merged = out
        base = s.replace(".vcf.gz", "").replace(".vcf", "")
        prefix = base.replace(".unfiltered", "")
        return prefix, merged
    return out, f"{out}.unfiltered.vcf.gz"

# ---------------------------
# command builders
# ---------------------------

def m2_cmd_for_contig(
    gatk,
    ref,
    tumor_bam,
    tumor_sample,
    germline_resource,
    pairhmm_threads,
    out_prefix_shard,
    contig,
    normal_bam=None,
    normal_sample=None,
    pon_vcf=None,
    intervals=None,
    tmp_dir=None,
):
    vcf = f"{out_prefix_shard}.unfiltered.vcf.gz"
    f1r2 = f"{out_prefix_shard}.f1r2.tar.gz"
    cmd = [
        gatk,
        "Mutect2",
        "-R",
        ref,
        "-I",
        tumor_bam,
        "-tumor",
        tumor_sample,
        "--f1r2-tar-gz",
        f1r2,
        "--native-pair-hmm-threads",
        str(pairhmm_threads),
        "-O",
        vcf,
        "-L",
        contig,
    ]
    if germline_resource:
        cmd += ["--germline-resource", germline_resource]
    if normal_bam and normal_sample:
        cmd += ["-I", normal_bam, "-normal", normal_sample]
    if pon_vcf:
        cmd += ["--panel-of-normals", pon_vcf]
    # If user gave intervals file, pass it too (GATK will union multiple -L).
    if intervals:
        cmd += ["-L", intervals]
    if tmp_dir:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        cmd += ["--tmp-dir", tmp_dir]
    return vcf, f1r2, shlex.join(cmd)

def m2_cmd_single(
    gatk,
    ref,
    tumor_bam,
    tumor_sample,
    germline_resource,
    pairhmm_threads,
    out_prefix,
    normal_bam=None,
    normal_sample=None,
    pon_vcf=None,
    intervals=None,
    tmp_dir=None,
):
    vcf = f"{out_prefix}.unfiltered.vcf.gz"
    f1r2 = f"{out_prefix}.f1r2.tar.gz"
    cmd = [
        gatk,
        "Mutect2",
        "-R",
        ref,
        "-I",
        tumor_bam,
        "-tumor",
        tumor_sample,
        "--f1r2-tar-gz",
        f1r2,
        "--native-pair-hmm-threads",
        str(pairhmm_threads),
        "-O",
        vcf,
    ]
    if germline_resource:
        cmd += ["--germline-resource", germline_resource]
    if normal_bam and normal_sample:
        cmd += ["-I", normal_bam, "-normal", normal_sample]
    if pon_vcf:
        cmd += ["--panel-of-normals", pon_vcf]
    if intervals:
        cmd += ["-L", intervals]
    if tmp_dir:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        cmd += ["--tmp-dir", tmp_dir]
    return vcf, f1r2, shlex.join(cmd)

# ---------------------------
# main
# ---------------------------

def main():
    opt = parse_args_tolerant(sys.argv[1:])

    # Required file-like inputs
    must(opt["ref"], "reference FASTA (-R/--ref)")
    must(opt["tumor_bam"], "tumor BAM (-T/--tumor-bam/-I)")

    # Required *string* input: tumor sample name
    if not opt["tumor_sample"]:
        fail("Missing required tumor sample name (-Ts/--tumor-sample/--tumor_sample).")

    # Normal logic: if normal BAM provided but sample missing → ignore normal
    if opt["normal_bam"] and not opt["normal_sample"]:
        print("[warn] -N provided without -Ns/--normal-sample; running tumor-only.", file=sys.stderr)
        opt["normal_bam"] = None

    gatk = opt["gatk"] or "gatk"
    pairhmm_threads = int(opt["pairhmm_threads"] or "4")
    dry = bool(opt["dry_run"])

    # If reference is plain gzip, decompress first
    opt["ref"] = ensure_uncompressed_fasta(opt["ref"], dry=dry)

    # Ensure indices/dicts for all key inputs
    print("[info] Ensuring reference indices (.fai/.dict)...")
    ensure_reference_indices(opt["ref"], gatk, dry=dry)

    print("[info] Ensuring BAM indices (.bai)...")
    ensure_bam_index(opt["tumor_bam"], dry=dry)
    if opt["normal_bam"]:
        ensure_bam_index(opt["normal_bam"], dry=dry)

    print("[info] Ensuring VCF indices (.tbi) for resources (if provided)...")
    ensure_vcf_index(opt["germline_resource"], gatk, dry=dry)
    if opt["pon_vcf"]:
        ensure_vcf_index(opt["pon_vcf"], gatk, dry=dry)

    prefix, merged_vcf = derive_prefix_and_merged_path(opt["out"])
    out_dir = Path(prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = opt["tmp_dir"]
    split5 = truthy(opt["split5"])
    contigs = normalized_contigs(opt["contigs"])

    # Decide max_parallel using either explicit --max-parallel or cpu_count
    if opt["max_parallel"]:
        try:
            max_parallel = int(opt["max_parallel"])
        except ValueError:
            max_parallel = 5
        max_parallel = max(1, min(5, max_parallel))
    else:
        cpu_count = None
        if opt["cpu_count"]:
            try:
                cpu_count = int(opt["cpu_count"])
            except ValueError:
                cpu_count = None
        if cpu_count and cpu_count > 0:
            max_parallel = max(1, min(5, cpu_count))
        else:
            max_parallel = 5

    # SINGLE RUN (default)
    if not split5:
        print("[info] split5 disabled → single Mutect2 run.")
        vcf, f1r2, cmd = m2_cmd_single(
            gatk=gatk,
            ref=opt["ref"],
            tumor_bam=opt["tumor_bam"],
            tumor_sample=opt["tumor_sample"],
            germline_resource=opt["germline_resource"],
            pairhmm_threads=pairhmm_threads,
            out_prefix=prefix,
            normal_bam=opt["normal_bam"],
            normal_sample=opt["normal_sample"],
            pon_vcf=opt["pon_vcf"],
            intervals=opt["intervals"],
            tmp_dir=tmp_dir,
        )
        if run(cmd, dry) != 0:
            fail("Mutect2 failed (single run).")

        # Link or copy to requested merged path
        if Path(vcf).resolve() != Path(merged_vcf).resolve() and not dry:
            Path(merged_vcf).unlink(missing_ok=True)
            try:
                os.link(vcf, merged_vcf)
            except Exception:
                Path(merged_vcf).write_bytes(Path(vcf).read_bytes())

        if ensure_vcf_index(merged_vcf, gatk, dry) != 0:
            fail("IndexFeatureFile failed on merged VCF.")

        manifest = f"{prefix}.f1r2_manifest.txt"
        if not dry:
            with open(manifest, "w") as fh:
                fh.write(f1r2 + "\n")

        print("\n[done] Outputs:")
        print(f"  Merged unfiltered VCF: {merged_vcf}")
        print(f"  F1R2 manifest:         {manifest}")
        print("  (single-run mode)")
        return

    # SPLIT5 MODE
    print(f"[info] split5 enabled → up to {max_parallel} parallel chromosomes.")
    print(f"[info] Contigs: {', '.join(contigs)}")

    shard_vcfs, shard_f1r2 = [], []

    def run_one(idx, contig):
        shard_prefix = f"{prefix}.shard{idx:03d}"
        vcf, f1r2, cmd = m2_cmd_for_contig(
            gatk=gatk,
            ref=opt["ref"],
            tumor_bam=opt["tumor_bam"],
            tumor_sample=opt["tumor_sample"],
            germline_resource=opt["germline_resource"],
            pairhmm_threads=pairhmm_threads,
            out_prefix_shard=shard_prefix,
            contig=contig,
            normal_bam=opt["normal_bam"],
            normal_sample=opt["normal_sample"],
            pon_vcf=opt["pon_vcf"],
            intervals=opt["intervals"],  # union with -L contig
            tmp_dir=tmp_dir,
        )
        rc = run(cmd, dry)
        return (rc, vcf, f1r2, contig)

    with ThreadPoolExecutor(max_workers=max_parallel) as ex:
        jobs = {ex.submit(run_one, i + 1, c): (i + 1, c) for i, c in enumerate(contigs)}
        for fut in as_completed(jobs):
            idx, contig = jobs[fut]
            try:
                rc, vcf, f1r2, c = fut.result()
                if rc != 0:
                    fail(f"Mutect2 failed for {c} (shard {idx:03d}) with exit {rc}")
                shard_vcfs.append(vcf)
                shard_f1r2.append(f1r2)
                print(f"[ok] {c} finished → {vcf}")
            except Exception as e:
                fail(f"Shard {idx:03d} raised exception: {e}")

    # Merge per-shard VCFs
    if len(shard_vcfs) == 1:
        if not dry:
            Path(merged_vcf).unlink(missing_ok=True)
            try:
                os.link(shard_vcfs[0], merged_vcf)
            except Exception:
                Path(merged_vcf).write_bytes(Path(shard_vcfs[0]).read_bytes())
        print(f"[info] Single-shard result → {merged_vcf}")
    else:
        merge_cmd = [gatk, "MergeVcfs"]
        for v in shard_vcfs:
            merge_cmd += ["-I", v]
        merge_cmd += ["-O", merged_vcf]
        if run(shlex.join(merge_cmd), dry) != 0:
            fail("MergeVcfs failed.")

    if ensure_vcf_index(merged_vcf, gatk, dry) != 0:
        fail("IndexFeatureFile failed on merged VCF.")

    # F1R2 manifest
    manifest = f"{prefix}.f1r2_manifest.txt"
    if not dry:
        with open(manifest, "w") as fh:
            for f in shard_f1r2:
                fh.write(f + "\n")

    print("\n[done] Outputs:")
    print(f"  Merged unfiltered VCF: {merged_vcf}")
    print(f"  F1R2 manifest:         {manifest}")
    print(f"  Shard VCFs:            {len(shard_vcfs)}")
    print(f"  Shard F1R2 tars:       {len(shard_f1r2)}")

if __name__ == "__main__":
    main()
