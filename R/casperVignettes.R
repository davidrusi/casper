casperDesign <- function(view=TRUE) {
  f <- system.file("doc","DesignRNASeq.pdf",package="casper")
  if (view) {
    if (.Platform$OS.type == "windows") {
       shell.exec(f)
   } else {
       system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
   }
   return(f)
}
