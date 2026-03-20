# Dockerfile for Mutect2 GenePattern Module
# Generated from planning data
# Resource requirements: 8 CPU cores, 16GB memory
# Supported input formats: bam, sam, cram, fa, fasta, vcf, vcf.gz, intervals, list, bed, interval_list
# Wrapper: mutect2_wrapper.sh (bash) — uses gatk and samtools (no pip/R packages)
FROM broadinstitute/gatk:4.6.2.0

# Metadata labels
LABEL maintainer="GenePattern"
LABEL module.name="Mutect2"
LABEL module.version="latest"
LABEL module.language="java"

# Set working directory
WORKDIR /module

# Install system dependencies verified against broadinstitute/gatk:4.6.2.0:
#   samtools  — BAM/FASTA indexing (samtools faidx, samtools view)
#   bcftools  — VCF/BCF utilities
#   bedtools  — genome arithmetic (intervals support)
#   tabix     — VCF .tbi index support
#   wget/curl — general download utilities
#   git, ca-certificates — baseline tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        curl \
        git \
        ca-certificates \
        samtools \
        bcftools \
        bedtools \
        tabix \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy required module files
# The wrapper script handles:
#   - Auto-generation of BAM indexes (.bai) via gatk BuildBamIndex
#   - Auto-generation of FASTA indexes (.fai / .dict) via samtools faidx
#     and gatk CreateSequenceDictionary
#   - Auto-generation of VCF tabix indexes (.tbi) via gatk IndexFeatureFile
#   - Tumor-normal paired calling and tumor-only mode
#   - Optional FilterMutectCalls and LearnReadOrientationModel post-processing
COPY mutect2_wrapper.sh /module/
RUN chmod +x /module/mutect2_wrapper.sh

# Environment variables for resource management
ENV MODULE_CPU_CORES=8
ENV MODULE_MEMORY=16GB
ENV MODULE_NAME=Mutect2

# Default entrypoint: bash (GenePattern will invoke the wrapper directly)
CMD ["/bin/bash"]
