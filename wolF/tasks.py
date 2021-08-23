import wolf

def absolute(seg_file, maf, skew, pairName):

    return wolf.Task(
      name = "ABSOLUTE",
      inputs = {
        "seg_file" : seg_file,
        "maf" : maf,
        "skew" : skew,
        "pairName" : pairName
      },
      outputs = {
        "absolute_highres_plot" : "*.ABSOLUTE_plot.pdf",
        # An R file containing an object seg.dat which provides all of the information used to generate the plot.
        "absolute_rdata" : "*.PP-modes.data.RData"
      },
      script = [
        "set -euxo pipefail",

        """SNV_MAF="${pairName}.snv.maf"
        INDEL_MAF="${pairName}.indel.maf"
        python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
        python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $INDEL_MAF -f Variant_Type -v "INS|DEL" """,

        'grep -v "NA" ${seg_file} > no_nan_segs.tsv',

        """
        awk 'BEGIN{FS=OFS="\t"} {gsub(/^chr/, "", $5)} 1' $SNV_MAF > reformat_snv.maf
        awk 'BEGIN{FS=OFS="\t"} {gsub(/^chr/, "", $5)} 1' $INDEL_MAF > reformat_indel.maf
        awk 'BEGIN{FS=OFS="\t"} {gsub(/^chr/, "", $1)} 1' no_nan_segs.tsv > reformat_seg.tsv
        """,
        

        """Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
        --seg_dat_fn reformat_seg.tsv \
        --maf_fn reformat_snv.maf \
        --indelmaf_fn reformat_indel.maf \
        --sample_name ${pairName} \
        --results_dir . \
        --ssnv_skew ${skew} \
        --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5"""
      ],
      resources = {
        "mem" : "8G"
      },
      docker = "gcr.io/broad-getzlab-workflows/absolute_wolf:v2"
    )