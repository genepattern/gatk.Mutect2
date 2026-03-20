#!/usr/bin/env bash
# =============================================================================
# mutect2_wrapper.sh  --  GenePattern wrapper for GATK Mutect2
#
# Supports:
#   - Tumor-normal paired calling and tumor-only calling
#   - Automatic index generation for BAM, FASTA, and VCF inputs
#   - Optional FilterMutectCalls post-processing (--run.filter)
#   - Optional orientation bias artifact correction
#     (--run.orientation.bias.filter)
#
# Validation test samples: tumor=HG008-T  normal=HG008-N-D
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
log_info()  { echo "[INFO]  $*"; }
log_warn()  { echo "[WARN]  $*"; }
log_error() { echo "[ERROR] $*" >&2; }

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Required:
  --tumor.bam FILE                  Sorted tumor BAM/CRAM file
  --reference FILE                  Reference genome FASTA

Optional inputs:
  --normal.bam FILE                 Matched normal BAM/CRAM (tumor-normal mode)
  --germline.resource FILE          Population AF VCF.gz (e.g. gnomAD af-only)
  --panel.of.normals FILE           Panel of Normals VCF.gz
  --intervals FILE                  Genomic intervals file

Sample names:
  --tumor.sample.name TEXT          SM tag for tumor (auto-detected if omitted)
  --normal.sample.name TEXT         SM tag for normal (auto-detected if omitted)

Output:
  --output.prefix TEXT              Prefix for output files (default: mutect2_output)

Pipeline control:
  --run.filter CHOICE               Run FilterMutectCalls [yes|no] (default: yes)
  --run.orientation.bias.filter CHOICE
                                    Run LearnReadOrientationModel [yes|no] (default: no)

Performance:
  --java.heap.size TEXT             JVM heap (e.g. 8g, 16g) (default: 8g)
  --pairhmm.threads INT             PairHMM threads (default: 4)

Advanced Mutect2 flags:
  --max.reads.per.alignment.start INT
  --initial.tumor.lod FLOAT
  --tumor.lod.to.emit FLOAT
  --normal.lod FLOAT
  --af.of.alleles.not.in.resource FLOAT
  --mitochondria.mode CHOICE        [yes|no]
  --pcr.indel.model CHOICE          [NONE|HOSTILE|AGGRESSIVE|CONSERVATIVE]
  --max.mnp.distance INT
  --base.quality.score.threshold INT
  --dont.use.soft.clipped.bases CHOICE   [yes|no]
  --genotype.germline.sites CHOICE       [yes|no]
  --genotype.pon.sites CHOICE            [yes|no]
EOF
    exit 1
}

# ---------------------------------------------------------------------------
# Default values
# ---------------------------------------------------------------------------
TUMOR_BAM=""
NORMAL_BAM=""
REFERENCE=""
GERMLINE_RESOURCE=""
PANEL_OF_NORMALS=""
INTERVALS=""
TUMOR_SAMPLE_NAME=""
NORMAL_SAMPLE_NAME=""
OUTPUT_PREFIX="mutect2_output"
RUN_FILTER="yes"
RUN_ORIENTATION_BIAS_FILTER="no"
JAVA_HEAP_SIZE="8g"
PAIRHMM_THREADS="4"
MAX_READS_PER_ALIGNMENT_START=""
INITIAL_TUMOR_LOD=""
TUMOR_LOD_TO_EMIT=""
NORMAL_LOD=""
AF_OF_ALLELES_NOT_IN_RESOURCE=""
MITOCHONDRIA_MODE=""
PCR_INDEL_MODEL=""
MAX_MNP_DISTANCE=""
BASE_QUALITY_SCORE_THRESHOLD=""
DONT_USE_SOFT_CLIPPED_BASES=""
GENOTYPE_GERMLINE_SITES=""
GENOTYPE_PON_SITES=""

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parse_arguments() {
    if [[ $# -eq 0 ]]; then
        usage
    fi

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --tumor.bam)                        TUMOR_BAM="$2";                        shift 2 ;;
            --normal.bam)                       NORMAL_BAM="$2";                       shift 2 ;;
            --reference)                        REFERENCE="$2";                        shift 2 ;;
            --germline.resource)                GERMLINE_RESOURCE="$2";                shift 2 ;;
            --panel.of.normals)                 PANEL_OF_NORMALS="$2";                 shift 2 ;;
            --intervals)                        INTERVALS="$2";                        shift 2 ;;
            --tumor.sample.name)                TUMOR_SAMPLE_NAME="$2";                shift 2 ;;
            --normal.sample.name)               NORMAL_SAMPLE_NAME="$2";               shift 2 ;;
            --output.prefix)                    OUTPUT_PREFIX="$2";                    shift 2 ;;
            --run.filter)                       RUN_FILTER="$2";                       shift 2 ;;
            --run.orientation.bias.filter)      RUN_ORIENTATION_BIAS_FILTER="$2";      shift 2 ;;
            --java.heap.size)                   JAVA_HEAP_SIZE="$2";                   shift 2 ;;
            --pairhmm.threads)                  PAIRHMM_THREADS="$2";                  shift 2 ;;
            --max.reads.per.alignment.start)    MAX_READS_PER_ALIGNMENT_START="$2";    shift 2 ;;
            --initial.tumor.lod)                INITIAL_TUMOR_LOD="$2";                shift 2 ;;
            --tumor.lod.to.emit)                TUMOR_LOD_TO_EMIT="$2";                shift 2 ;;
            --normal.lod)                       NORMAL_LOD="$2";                       shift 2 ;;
            --af.of.alleles.not.in.resource)    AF_OF_ALLELES_NOT_IN_RESOURCE="$2";    shift 2 ;;
            --mitochondria.mode)                MITOCHONDRIA_MODE="$2";                shift 2 ;;
            --pcr.indel.model)                  PCR_INDEL_MODEL="$2";                  shift 2 ;;
            --max.mnp.distance)                 MAX_MNP_DISTANCE="$2";                 shift 2 ;;
            --base.quality.score.threshold)     BASE_QUALITY_SCORE_THRESHOLD="$2";     shift 2 ;;
            --dont.use.soft.clipped.bases)      DONT_USE_SOFT_CLIPPED_BASES="$2";      shift 2 ;;
            --genotype.germline.sites)          GENOTYPE_GERMLINE_SITES="$2";          shift 2 ;;
            --genotype.pon.sites)               GENOTYPE_PON_SITES="$2";               shift 2 ;;
            -h|--help)                          usage ;;
            *)
                log_error "Unknown option: $1"
                usage
                ;;
        esac
    done
}

# ---------------------------------------------------------------------------
# Helper: normalize boolean-style choices to true/false
# Returns 0 (true) if value is yes/true/1, 1 (false) otherwise
# ---------------------------------------------------------------------------
is_yes() {
    local val
    val=$(echo "$1" | tr '[:upper:]' '[:lower:]')
    [[ "$val" == "yes" || "$val" == "true" || "$val" == "1" ]]
}

# ---------------------------------------------------------------------------
# Index helpers
# ---------------------------------------------------------------------------

# Ensure a BAM or CRAM file has its companion index (.bai / .crai).
ensure_bam_index() {
    local bam="$1"
    local bai_inline="${bam}.bai"
    local bai_sidecar="${bam%.bam}.bai"
    local crai_inline="${bam}.crai"

    if [[ -f "$bai_inline" || -f "$bai_sidecar" || -f "$crai_inline" ]]; then
        log_info "BAM/CRAM index found for ${bam}"
        return 0
    fi

    log_info "BAM/CRAM index not found for ${bam} -- running gatk BuildBamIndex"
    gatk --java-options "-Xmx${JAVA_HEAP_SIZE}" BuildBamIndex \
        -I "${bam}" \
        -O "${bai_sidecar}"
    log_info "Index created: ${bai_sidecar}"
}

# Ensure a VCF.gz file has a companion tabix index (.tbi).
ensure_vcf_index() {
    local vcf="$1"
    local tbi="${vcf}.tbi"

    if [[ -f "$tbi" ]]; then
        log_info "VCF index found for ${vcf}"
        return 0
    fi

    log_info "VCF index (.tbi) not found for ${vcf} -- running gatk IndexFeatureFile"
    gatk --java-options "-Xmx${JAVA_HEAP_SIZE}" IndexFeatureFile \
        -I "${vcf}"
    log_info "VCF index created: ${tbi}"
}

# Ensure a FASTA file has both .fai and .dict companion files.
ensure_fasta_index() {
    local fa="$1"
    local fai="${fa}.fai"
    # Convention: reference.fa -> reference.dict  OR  reference.fasta -> reference.dict
    local dict="${fa%.*}.dict"

    if [[ ! -f "$fai" ]]; then
        log_info "FASTA index (.fai) not found for ${fa} -- running samtools faidx"
        samtools faidx "${fa}"
        log_info "FASTA index created: ${fai}"
    else
        log_info "FASTA index (.fai) found for ${fa}"
    fi

    if [[ ! -f "$dict" ]]; then
        log_info "Sequence dictionary (.dict) not found -- running gatk CreateSequenceDictionary"
        gatk --java-options "-Xmx${JAVA_HEAP_SIZE}" CreateSequenceDictionary \
            -R "${fa}"
        log_info "Sequence dictionary created: ${dict}"
    else
        log_info "Sequence dictionary found: ${dict}"
    fi
}

# Auto-detect the SM tag from a BAM header.
get_sm_tag() {
    local bam="$1"
    local sm
    sm=$(samtools view -H "${bam}" 2>/dev/null \
         | grep '^@RG' \
         | grep -oP '(?<=\tSM:)[^\t]+' \
         | head -1) || true
    echo "${sm}"
}

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------
validate_inputs() {
    local errors=0

    # --- Required parameters ---
    if [[ -z "${TUMOR_BAM}" ]]; then
        log_error "Required parameter --tumor.bam is missing."
        errors=$((errors + 1))
    elif [[ ! -f "${TUMOR_BAM}" ]]; then
        log_error "Tumor BAM/CRAM file not found: ${TUMOR_BAM}"
        errors=$((errors + 1))
    fi

    if [[ -z "${REFERENCE}" ]]; then
        log_error "Required parameter --reference is missing."
        errors=$((errors + 1))
    elif [[ ! -f "${REFERENCE}" ]]; then
        log_error "Reference FASTA file not found: ${REFERENCE}"
        errors=$((errors + 1))
    fi

    # --- Optional file parameters ---
    if [[ -n "${NORMAL_BAM}" && ! -f "${NORMAL_BAM}" ]]; then
        log_error "Normal BAM/CRAM file not found: ${NORMAL_BAM}"
        errors=$((errors + 1))
    fi

    if [[ -n "${GERMLINE_RESOURCE}" && ! -f "${GERMLINE_RESOURCE}" ]]; then
        log_error "Germline resource VCF not found: ${GERMLINE_RESOURCE}"
        errors=$((errors + 1))
    fi

    if [[ -n "${PANEL_OF_NORMALS}" && ! -f "${PANEL_OF_NORMALS}" ]]; then
        log_error "Panel of Normals VCF not found: ${PANEL_OF_NORMALS}"
        errors=$((errors + 1))
    fi

    if [[ -n "${INTERVALS}" && ! -f "${INTERVALS}" ]]; then
        log_error "Intervals file not found: ${INTERVALS}"
        errors=$((errors + 1))
    fi

    if [[ ${errors} -gt 0 ]]; then
        log_error "Input validation failed with ${errors} error(s). Exiting."
        exit 1
    fi

    log_info "Input validation passed."
}

# ---------------------------------------------------------------------------
# Ensure all companion index files exist
# ---------------------------------------------------------------------------
prepare_indexes() {
    log_info "Checking/generating companion index files..."

    ensure_fasta_index "${REFERENCE}"
    ensure_bam_index   "${TUMOR_BAM}"

    if [[ -n "${NORMAL_BAM}" ]]; then
        ensure_bam_index "${NORMAL_BAM}"
    fi

    if [[ -n "${GERMLINE_RESOURCE}" ]]; then
        ensure_vcf_index "${GERMLINE_RESOURCE}"
    fi

    if [[ -n "${PANEL_OF_NORMALS}" ]]; then
        ensure_vcf_index "${PANEL_OF_NORMALS}"
    fi
}

# ---------------------------------------------------------------------------
# Resolve sample names
# ---------------------------------------------------------------------------
resolve_sample_names() {
    if [[ -z "${TUMOR_SAMPLE_NAME}" ]]; then
        log_info "Auto-detecting tumor sample name from BAM header..."
        TUMOR_SAMPLE_NAME=$(get_sm_tag "${TUMOR_BAM}")
        if [[ -z "${TUMOR_SAMPLE_NAME}" ]]; then
            log_error "Could not auto-detect tumor sample name from BAM header."
            log_error "Please provide --tumor.sample.name explicitly."
            exit 1
        fi
        log_info "Detected tumor sample name: ${TUMOR_SAMPLE_NAME}"
    fi

    if [[ -n "${NORMAL_BAM}" && -z "${NORMAL_SAMPLE_NAME}" ]]; then
        log_info "Auto-detecting normal sample name from BAM header..."
        NORMAL_SAMPLE_NAME=$(get_sm_tag "${NORMAL_BAM}")
        if [[ -z "${NORMAL_SAMPLE_NAME}" ]]; then
            log_error "Could not auto-detect normal sample name from BAM header."
            log_error "Please provide --normal.sample.name explicitly."
            exit 1
        fi
        log_info "Detected normal sample name: ${NORMAL_SAMPLE_NAME}"
    fi
}

# ---------------------------------------------------------------------------
# Run Mutect2
# ---------------------------------------------------------------------------
run_mutect2() {
    local unfiltered_vcf="${OUTPUT_PREFIX}_unfiltered.vcf.gz"
    local f1r2_tar="${OUTPUT_PREFIX}_f1r2.tar.gz"
    local stats_file="${unfiltered_vcf}.stats"

    log_info "============================================================"
    log_info "Starting Mutect2"
    log_info "  Tumor BAM        : ${TUMOR_BAM}"
    log_info "  Tumor sample     : ${TUMOR_SAMPLE_NAME}"
    log_info "  Normal BAM       : ${NORMAL_BAM:-<tumor-only mode>}"
    log_info "  Normal sample    : ${NORMAL_SAMPLE_NAME:-N/A}"
    log_info "  Reference        : ${REFERENCE}"
    log_info "  Output prefix    : ${OUTPUT_PREFIX}"
    log_info "============================================================"

    local cmd=(
        gatk
        --java-options "-Xmx${JAVA_HEAP_SIZE}"
        Mutect2
        -R "${REFERENCE}"
        -I "${TUMOR_BAM}"
        -tumor "${TUMOR_SAMPLE_NAME}"
        -O "${unfiltered_vcf}"
        --native-pair-hmm-threads "${PAIRHMM_THREADS}"
    )

    # Normal BAM (tumor-normal mode)
    if [[ -n "${NORMAL_BAM}" ]]; then
        cmd+=( -I "${NORMAL_BAM}" )
        cmd+=( -normal "${NORMAL_SAMPLE_NAME}" )
    fi

    # Optional resource files
    if [[ -n "${GERMLINE_RESOURCE}" ]]; then
        cmd+=( --germline-resource "${GERMLINE_RESOURCE}" )
    fi

    if [[ -n "${PANEL_OF_NORMALS}" ]]; then
        cmd+=( --panel-of-normals "${PANEL_OF_NORMALS}" )
    fi

    if [[ -n "${INTERVALS}" ]]; then
        cmd+=( -L "${INTERVALS}" )
    fi

    # Orientation bias: collect F1R2 data
    if is_yes "${RUN_ORIENTATION_BIAS_FILTER}"; then
        cmd+=( --f1r2-tar-gz "${f1r2_tar}" )
    fi

    # LOD thresholds
    if [[ -n "${INITIAL_TUMOR_LOD}" ]]; then
        cmd+=( --initial-tumor-lod "${INITIAL_TUMOR_LOD}" )
    fi
    if [[ -n "${TUMOR_LOD_TO_EMIT}" ]]; then
        cmd+=( --tumor-lod-to-emit "${TUMOR_LOD_TO_EMIT}" )
    fi
    if [[ -n "${NORMAL_LOD}" ]]; then
        cmd+=( --normal-lod "${NORMAL_LOD}" )
    fi

    # AF prior
    if [[ -n "${AF_OF_ALLELES_NOT_IN_RESOURCE}" ]]; then
        cmd+=( --af-of-alleles-not-in-resource "${AF_OF_ALLELES_NOT_IN_RESOURCE}" )
    fi

    # Mitochondria mode
    if [[ -n "${MITOCHONDRIA_MODE}" ]] && is_yes "${MITOCHONDRIA_MODE}"; then
        cmd+=( --mitochondria-mode )
    fi

    # PCR indel model
    if [[ -n "${PCR_INDEL_MODEL}" ]]; then
        cmd+=( --pcr-indel-model "${PCR_INDEL_MODEL}" )
    fi

    # MNP distance
    if [[ -n "${MAX_MNP_DISTANCE}" ]]; then
        cmd+=( --max-mnp-distance "${MAX_MNP_DISTANCE}" )
    fi

    # Base quality threshold
    if [[ -n "${BASE_QUALITY_SCORE_THRESHOLD}" ]]; then
        cmd+=( --base-quality-score-threshold "${BASE_QUALITY_SCORE_THRESHOLD}" )
    fi

    # Max reads per alignment start
    if [[ -n "${MAX_READS_PER_ALIGNMENT_START}" ]]; then
        cmd+=( --max-reads-per-alignment-start "${MAX_READS_PER_ALIGNMENT_START}" )
    fi

    # Soft-clipped bases
    if [[ -n "${DONT_USE_SOFT_CLIPPED_BASES}" ]] && is_yes "${DONT_USE_SOFT_CLIPPED_BASES}"; then
        cmd+=( --dont-use-soft-clipped-bases )
    fi

    # Germline / PoN site genotyping
    if [[ -n "${GENOTYPE_GERMLINE_SITES}" ]] && is_yes "${GENOTYPE_GERMLINE_SITES}"; then
        cmd+=( --genotype-germline-sites )
    fi
    if [[ -n "${GENOTYPE_PON_SITES}" ]] && is_yes "${GENOTYPE_PON_SITES}"; then
        cmd+=( --genotype-pon-sites )
    fi

    log_info "Mutect2 command: ${cmd[*]}"
    "${cmd[@]}"

    if [[ ! -f "${unfiltered_vcf}" ]]; then
        log_error "Mutect2 did not produce expected output: ${unfiltered_vcf}"
        exit 2
    fi

    log_info "Mutect2 complete. Unfiltered VCF: ${unfiltered_vcf}"
}

# ---------------------------------------------------------------------------
# Optionally run LearnReadOrientationModel
# ---------------------------------------------------------------------------
learn_orientation_model() {
    local f1r2_tar="${OUTPUT_PREFIX}_f1r2.tar.gz"
    local rom_tar="${OUTPUT_PREFIX}_read-orientation-model.tar.gz"

    if [[ ! -f "${f1r2_tar}" ]]; then
        log_error "F1R2 counts file not found: ${f1r2_tar}"
        exit 2
    fi

    log_info "Running LearnReadOrientationModel..."
    gatk --java-options "-Xmx${JAVA_HEAP_SIZE}" LearnReadOrientationModel \
        -I "${f1r2_tar}" \
        -O "${rom_tar}"

    if [[ ! -f "${rom_tar}" ]]; then
        log_error "LearnReadOrientationModel did not produce expected output: ${rom_tar}"
        exit 2
    fi

    log_info "Read orientation model created: ${rom_tar}"
}

# ---------------------------------------------------------------------------
# Optionally run FilterMutectCalls
# ---------------------------------------------------------------------------
run_filter_mutect_calls() {
    local unfiltered_vcf="${OUTPUT_PREFIX}_unfiltered.vcf.gz"
    local filtered_vcf="${OUTPUT_PREFIX}_filtered.vcf.gz"
    local rom_tar="${OUTPUT_PREFIX}_read-orientation-model.tar.gz"

    log_info "Running FilterMutectCalls..."

    local filter_cmd=(
        gatk
        --java-options "-Xmx${JAVA_HEAP_SIZE}"
        FilterMutectCalls
        -R "${REFERENCE}"
        -V "${unfiltered_vcf}"
        -O "${filtered_vcf}"
    )

    # Add orientation bias priors if available
    if is_yes "${RUN_ORIENTATION_BIAS_FILTER}" && [[ -f "${rom_tar}" ]]; then
        filter_cmd+=( --orientation-bias-artifact-priors "${rom_tar}" )
    fi

    log_info "FilterMutectCalls command: ${filter_cmd[*]}"
    "${filter_cmd[@]}"

    if [[ ! -f "${filtered_vcf}" ]]; then
        log_error "FilterMutectCalls did not produce expected output: ${filtered_vcf}"
        exit 2
    fi

    log_info "Filtering complete. Filtered VCF: ${filtered_vcf}"
}

# ---------------------------------------------------------------------------
# Print summary of outputs
# ---------------------------------------------------------------------------
print_summary() {
    log_info "============================================================"
    log_info "Mutect2 wrapper completed successfully."
    log_info "Output files:"

    local unfiltered_vcf="${OUTPUT_PREFIX}_unfiltered.vcf.gz"
    local filtered_vcf="${OUTPUT_PREFIX}_filtered.vcf.gz"
    local f1r2_tar="${OUTPUT_PREFIX}_f1r2.tar.gz"
    local rom_tar="${OUTPUT_PREFIX}_read-orientation-model.tar.gz"

    [[ -f "${unfiltered_vcf}" ]] && log_info "  Unfiltered VCF   : ${unfiltered_vcf}"
    [[ -f "${filtered_vcf}"   ]] && log_info "  Filtered VCF     : ${filtered_vcf}"
    [[ -f "${f1r2_tar}"       ]] && log_info "  F1R2 counts      : ${f1r2_tar}"
    [[ -f "${rom_tar}"        ]] && log_info "  Orientation model: ${rom_tar}"
    log_info "============================================================"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main() {
    parse_arguments "$@"
    validate_inputs
    prepare_indexes
    resolve_sample_names
    run_mutect2

    if is_yes "${RUN_ORIENTATION_BIAS_FILTER}"; then
        learn_orientation_model
    fi

    if is_yes "${RUN_FILTER}"; then
        run_filter_mutect_calls
    fi

    print_summary
    exit 0
}

main "$@"
