import wolf

#def absolute(seg_file, maf, skew, pairName):
class absolute(wolf.Task): 
  name = "ABSOLUTE" 
  inputs = {
    "seg_file",
    "maf",
    "skew",
    "pairName"
  } 
  output_patterns = {
    "absolute_highres_plot" : "*.ABSOLUTE_plot.pdf",
    # An R file containing an object seg.dat which provides all of the information used to generate the plot.
    "absolute_rdata" : "*.PP-modes.data.RData"
  } 
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
  ] 
  resources = {
    "mem" : "8G"
  } 
  docker = "gcr.io/broad-getzlab-workflows/absolute_wolf:v6"

class absolute_extract(wolf.Task):
    name = "ABSOLUTE_extract"
    inputs = {
      "abs_solution_number" : None,
      "absolute_rdata" : None,
      "sample_name" : None,
      "analyst_id" : "wolf"
    }
    script = """
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_extract_cli_start.R \
      --solution_num ${abs_solution_number} \
      --analyst_id ${analyst_id} \
      --rdata_modes_fn ${absolute_rdata} \
      --sample_name ${sample_name} \
      --results_dir . \
      --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

    cp reviewed/SEG_MAF/${sample_name}_ABS_MAF.txt .
    cp reviewed/SEG_MAF/${sample_name}.segtab.txt .
    cp reviewed/samples/${sample_name}.ABSOLUTE.${analyst_id}.called.RData .
    cp reviewed/${sample_name}.${analyst_id}.ABSOLUTE.table.txt .
    
    cut -f4 ${sample_name}.${analyst_id}.ABSOLUTE.table.txt | tail -n1 > purity
    cut -f5 ${sample_name}.${analyst_id}.ABSOLUTE.table.txt | tail -n1 > ploidy
    """
    output_patterns = {
        "absolute_annotated_maf_capture" : "*ABS_MAF.txt",
        "absolute_seg_file" : "*.segtab.txt",
        "absolute_segdat_file" : "*.called.RData",
        "absolute_table" : "*.ABSOLUTE.table.txt",
        "purity" : ("purity", wolf.read_file),
        "ploidy" : ("ploidy", wolf.read_file)
    }
    resources = {
    "mem" : "8G"
    } 
    docker = "gcr.io/broad-getzlab-workflows/absolute_wolf:v6"

class absolute_forcecall(wolf.Task):
    name = "ABSOLUTE_forcecall"
    inputs = {
      "seg_file",
      "maf",
      "ploidy",
      "purity",
      "skew",
      "sample_name"
    }
    script = """
    SNV_MAF="${pairName}.snv.maf"
    INDEL_MAF="${pairName}.indel.maf"
    python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP"
    python /usr/local/bin/split_maf_indel_snp.py -i ${maf} -o $INDEL_MAF -f Variant_Type -v "INS|DEL"

    #python /opt/re_center.py no_nan_segs.tsv no_nan
    grep -v "NA" ${seg_file} > no_nan_segs.tsv

    awk 'BEGIN{FS=OFS="\t"} {gsub(/^chr/, "", $5)} 1' $SNV_MAF > reformat_snv.maf
    awk 'BEGIN{FS=OFS="\t"} {gsub(/^chr/, "", $5)} 1' $INDEL_MAF > reformat_indel.maf
    awk 'BEGIN{FS=OFS="\t"} {gsub(/^chr/, "", $1)} 1' no_nan_segs.tsv > reformat_seg.tsv
    
    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
      --seg_dat_fn reformat_seg.tsv \
      --maf_fn reformat_snv.maf \
      --indelmaf_fn reformat_indel.maf \
      --sample_name ${sample_name} \
      --results_dir . \
      --ssnv_skew ${skew} \
      --force_alpha ${purity} \
      --force_tau ${ploidy} \
      --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_extract_cli_start.R \
      --solution_num 1 \
      --analyst_id force_called \
      --rdata_modes_fn ${sample_name}.PP-modes.data.RData \
      --sample_name ${sample_name} \
      --results_dir . \
      --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

    cp reviewed/samples/${sample_name}.ABSOLUTE.force_called.called.RData .
    cp reviewed/SEG_MAF/${sample_name}_ABS_MAF.txt .
    cp reviewed/SEG_MAF/${sample_name}.segtab.txt .
    """
    output_patterns = {
        "absolute_annotated_maf_capture" : "*ABS_MAF.txt",
        "absolute_seg_file" : "*.segtab.txt",
        "absolute_segdat_file" : "*.called.RData",
        "absolute_highres_plot" : "*.ABSOLUTE_plot.pdf",
        "absolute_rdata" : "*.PP-modes.data.RData"
    }
    resources = {
    "mem" : "8G"
    } 
    docker = "gcr.io/broad-getzlab-workflows/absolute_wolf:no_indel_filter_v6"

class absolute_to_cn_profile(wolf.Task):
    name = "absolute_to_cn_profile"
    inputs = {
            "sample_id",
            "rdata_fn"
            }
    script = """
    Rscript /usr/local/bin/get_CN_Absolute.Phylogic_SinglePatientTiming.R \
            --rdata_fn ${rdata_fn} \
            --pair_id ${sample_id} \
            --pcawg TRUE \
            --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/
    """
    output_patterns = {
            "CN_Profile_ccf":"*segments.txt",
            "CN_Profile_tsv":"*.tsv"
            }
    resources = {
            "mem": "8G"
            }
    docker = "gcr.io/broad-getzlab-workflows/absolute_wolf:v19"
