run_liftover <- function(bim_df, output_path, from_build = "hg19") {
  gz_path    <- tempfile(fileext = ".chain.gz")
  chain_path <- sub("\\.gz$", "", gz_path)
  
  chain_url <- switch(from_build,
                      "hg19" = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                      "hg18" = "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz",
                      stop("❌ Genoma d'origen no suportat: ", from_build)
  )
  
  download.file(chain_url, destfile = gz_path, mode = "wb")
  R.utils::gunzip(gz_path, destname = chain_path, remove = FALSE)
  chain <- rtracklayer::import.chain(chain_path)
  
  bim_df <- bim_df[bim_df$CHR %in% 1:22, ]
  
  gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", bim_df$CHR),
    ranges   = IRanges::IRanges(start = bim_df$POS, width = 1),
    strand   = "*",
    SNP = bim_df$SNP,
    A1  = bim_df$A1,
    A2  = bim_df$A2
  )
  
  lifted <- rtracklayer::liftOver(gr, chain)
  mapped <- unlist(lifted[S4Vectors::elementNROWS(lifted) == 1])
  
  message("🔎 SNPs mapejats LiftOver: ", length(mapped))
  
  bim_hg38 <- data.table::data.table(
    CHR = gsub("chr", "", as.character(GenomicRanges::seqnames(mapped))),
    SNP = mapped$SNP,
    CM  = 0,
    POS = GenomicRanges::start(mapped),
    A1  = mapped$A1,
    A2  = mapped$A2
  )
  
  data.table::fwrite(bim_hg38, file = paste0(output_path, ".bim"),
                     sep = "\t", col.names = FALSE)
  invisible(nrow(bim_hg38))
}
