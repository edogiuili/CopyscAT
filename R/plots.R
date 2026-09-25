#' Plot the distribution of collapsed signal
#'
#' Writes a violin plot of per-arm signal across cells to the session output
#' path. Useful as a QC check before and after scaling.
#'
#' @param inputMatrix Collapsed matrix from [collapseChrom3N()].
#' @param outputSuffix Suffix for the PDF file name.
#'
#' @return Invisibly `NULL`. Called for the file it writes.
#' @family reporting
#' @export
#' @examples
#' \dontrun{
#' graphCNVDistribution(scData_collapse, outputSuffix = "violins")
#' }
graphCNVDistribution <- function(inputMatrix,outputSuffix="_all")
{
  pdf(file=str_c(scCNVCaller$locPrefix,scCNVCaller$outPrefix,outputSuffix,"_graph.pdf"),width=8,height=6)
  print(ggplot(inputMatrix  %>% gather(Cell,Density,dplyr::ends_with(scCNVCaller$cellSuffix)),aes(chrom,Density))+geom_violin(scale="width",trim=FALSE) +
          theme(axis.text.x=element_text(angle=-90,vjust = 0.5,hjust = 0, color = "#000000"),
                axis.text.y=element_text(color="#000000"),
                axis.line.x.bottom = element_line(colour="#000000"),
                axis.line.y.left =  element_line(colour="#000000"),
                panel.background = element_rect(fill="#ffffff"),
                plot.background = element_rect(color="#ffffff",fill="#ffffff")) +
          xlab("Segment"))
  dev.off()
}
