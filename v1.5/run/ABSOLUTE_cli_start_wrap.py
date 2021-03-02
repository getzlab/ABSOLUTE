import os
import sys
import subprocess
print sys.argv

abs_cmd="""Rscript /xchip/cga_home/igleshch/soft/local/absolute/run/ABSOLUTE_cli_start.R \
 --seg_dat_fn {seg_file} --maf_fn {snv_maf} --indelmaf_fn {indel_maf} --sample_name {sample_name} --results_dir results \
--ssnv_skew {skew} """

abs_cmd=abs_cmd.format(sample_name=sys.argv[1],seg_file=sys.argv[3],snv_maf=sys.argv[4],indel_maf=sys.argv[5],skew=sys.argv[2])

log_file=open(sys.argv[1]+".log",'w')
print "starting absolute"
print >>log_file,subprocess.check_output(abs_cmd,shell=True)

log_file.close()
