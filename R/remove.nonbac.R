remove.nonbac = function(x){
  if (any(is.na(tax_table(x)))){
    tax_table(x)[is.na(tax_table(x))]<-"unidentified"
  } else {
    tax_table(x) <- x
  }
  if (any(tax_table(x)[,1]!="Bacteria" )) {
    nonbac = c(taxa_names(subset_taxa(x, Kingdom!="Bacteria")))
  } else {
    nonbac = as.character()
  }
  if (any(tax_table(x)[,3]=="Chloroplast", na.rm=TRUE)) {
    chl = c(taxa_names(subset_taxa(x, Class=="Chloroplast")))
  }else {
    chl = as.character()
  }
  if (any(tax_table(x)[,6]=="Mitochondria", na.rm=TRUE)) {
    mito = c(taxa_names(subset_taxa(x, Family=="Mitochondria")))
  }else {
    mito = as.character()
  }
  removetaxa = c(nonbac, chl, mito)
  allTaxa = taxa_names(x)
  myTaxa <- allTaxa[!(allTaxa %in% removetaxa)]
  new_phy = prune_taxa(myTaxa, x)
  return(new_phy)
}
