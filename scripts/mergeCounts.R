####################################
## Combine counts in 1 big matrix ##
####################################

suppressPackageStartupMessages({
  library("optparse", quietly = TRUE)
})


#command line options
option_list = list(
  make_option(c("-d","--countdir"), type="character", default=NULL, 
              help="directory containing the gene level counts per sample", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#make sure all required options are filled in
if (is.null(opt$countdir)){
  print_help(opt_parser)
  stop("Please provide the count directory", call.=FALSE)
}

countdir <- opt$countdir
outdir <- opt$countdir   #output in same directory as input

#get all files
files = list.files(path = countdir, pattern = "\\.counts$")

#reads counted by HTseq-count
counts <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c('gene_ID'))
for (f in files){
  #empty bam files are ignored
  if (file.size(paste0(countdir, "/", f)) == 0L){
    next
  }
  count_list <- read.table(paste0(countdir, "/", f), header = F, sep = "\t")
  sample <- sub('\\.counts$', '', f)
  
  if (ncol(count_list) == 2){
    colnames(count_list) <- c('gene_ID', sample)
  }else{
    colnames(count_list) <- c('gene_ID', 'gene_name', sample)
    count_list$gene_name <- NULL
  }
  
  counts <- merge(x = counts, y = count_list, by = c('gene_ID'), all= TRUE)
}

#set NAs to 0
counts[is.na(counts)] <- 0

#remove columns
rownames(counts) <- counts$gene_ID
counts$gene_ID <- NULL

# #write count matrix to file
write.table(counts, paste0(outdir, "/scCounts.tsv"), sep = '\t', row.names = T, col.names = T, quote = F)
