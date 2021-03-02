## FIXME: Point lib.loc to a more appropriate location
library(ABSOLUTE, lib.loc="/xchip/cga2/jgentry/firehoseRLibrary")

firehose_extract_reviewed_results = function(call_override, analyst.id, seg.dat, 
                                             out.dir.base, obj.name, copy_num_type, 
				             verbose=FALSE) {
  if (copy_num_type == "total") {
    ABSOLUTE:::set_total_funcs()
  } else if (copy_num_type == "allelic") {
    ABSOLUTE:::set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  
  ## provides segobj.list
  segobj.list = list(seg.dat)
  names(segobj.list) = seg.dat$array.name
  
  pp.calls = ABSOLUTE:::GetPpTabVals(seg.dat, seg.dat$mode.res$call.status)
  
  called.segobj.list = ABSOLUTE:::override_absolute_calls(segobj.list, call_override)

  ABSOLUTE:::process_extract_reviewed_results(out.dir.base, called.segobj.list, pp.calls, obj.name, analyst.id)
}

library(optparse)
option_list = list(
  make_option("--sample", dest="sample"),
  make_option("--purity", dest="purity"),
  make_option("--ploidy", dest="ploidy"),
  make_option("--solution", dest="solution"),
  make_option("--abs_fn", dest="abs_fn"))

opt = parse_args(OptionParser(option_list=option_list))
save(opt, file="debug.RData")

tmp_analyst = "temp"
cwd = getwd()

load(opt$abs_fn)

if (opt$solution == -1) {
  ## the '-1' is a hack to flag the case where we have a legacy call but not a *real*
  ## called solution, it isn't used here so just pass it through.
  opt$solution = 1
}

firehose_extract_reviewed_results(opt$solution, tmp_analyst, seg.dat, cwd, opt$sample, "allelic")
pp_fn = paste(opt$sample, tmp_analyst, "ABSOLUTE", "table", "txt", sep=".")
pp_call = read.delim(file.path(cwd, "reviewed", pp_fn))

get_value = function(val, val_name, pp_call, mode.flag) {
  if ((is.na(val)) || (is.null(val)) || (val %in% c("null", "NULL", "NA"))) {
    out = pp_call[1, val_name]
  } else {
    print(paste("Using previously called value of", val, "for", val_name))
    out = val
  }

  return(out)
}

out_purity = get_value(opt$purity, "purity", pp_call, seg.dat$mode.res$mode.flag)
out_ploidy = get_value(opt$ploidy, "ploidy", pp_call, seg.dat$mode.res$mode.flag)

out_mat = cbind(opt$sample, out_purity, out_ploidy)
colnames(out_mat) = c("sample_id", "absolute_extract_purity", "absolute_extract_ploidy")

write.table(out_mat, file="import_file.txt", quote=FALSE, row.names=FALSE, sep="\t")
