## TODO: define optional args ?

#' Summarize ASV table, taxonomy and phenodata in a Phyloseq object, with ASVs renamed to more convenient shorter form, and actual ASV sequences stored as a DNAStringSet (refseq slot)
#'
#' Takes as input an ASV table and taxonomy table created in the BiMiCo pipeline, with optional sample phenotypes imported e.g. from a csv table, and creates a Phyloseq object.
#' @param asvtab (Required) ASV table generated in BimiCo (ASVs in rows)
#' @param taxtab (Required) Taxonomy table generated in BiMiCo
#' @param phenodata (Required) sample phenotype matrix, if available; sample names must match with ASV table sample names.
#' @keywords read processing dada2
#' @export
#' @examples
#' create_phylo()

create_phylo <- function(asvtb, taxtb, phenodata){

    phs <- phyloseq::phyloseq(phyloseq::otu_table(asvtb, taxa_are_rows=TRUE),
                                phyloseq::tax_table(taxtb),
                                phyloseq::sample_data((phenodata)))
    seqs <- Biostrings::DNAStringSet(taxa_names(phs))
    names(seqs) <- phyloseq::taxa_names(phs)
    phs <- phyloseq::merge_phyloseq(phs, seqs)
    phyloseq::taxa_names(phs) <- paste0("ASV", seq(phyloseq::ntaxa(phs)))

  return(phs)

}

