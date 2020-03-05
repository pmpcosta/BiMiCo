#'Pricipal coordinate analysis (PCoA)
#'
#'Principal coordinate analysis using sample-wise Bray-Curtis dissimilarity index.
#'
#' @param data (Required) The microbiota data on which you want to perform ordination needs to represent phyloseq-class.i.e. the data should include the components of otu_table_class, taxonomyTable_class and sample_data_class.otu_table_class should be organised so that rows represent samples and columns represent OTUs.
#' @param dotcol (Optional) column name in sample_data of the phyloseq object specified to argument "data". Uses phyloseq default color palette (default:black)
#' @return Returns principal coordinate analysis object.
#'
#' @examples
#' pcoa_ph()
#'
#' @export
pcoa_ph <- function(data, dotcol = "black") {

  pcoa <- phyloseq::ordinate(data,
                    method = "PCoA",
                    distance = "bray"
                    )

  plot1 <- phyloseq::plot_ordination(data,
                                     pcoa,
                                     color = dotcol)

  return(plot1)

}
