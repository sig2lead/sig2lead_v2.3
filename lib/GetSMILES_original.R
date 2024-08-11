GetSMILES<- function()
    {
  Compounds<<-read.csv(file="LINCSCompounds.csv", stringsAsFactors = FALSE)
  
  compound_rows <- list("integer")
  rows_to_add <- list("integer")

#if (input$Signature == "Input a Gene"){
  for (i in 1:length(con_comma))
  {
    if (length(con_comma[[i]][[11]])!=0){
      rows_to_add <- which(Compounds$SM_LINCS_ID %in% con_comma[[i]][[11]])
      compound_rows[[i]] <- rows_to_add
    }
    else {compound_rows[[i]] <- 0}
    }
  compound_rows <- compound_rows
#}

  
#######con_comma currently gets LINCSCP instead of LSM. We need LSM in the download from LINCS (ilincsresult)  
#else if (input$Signature == "Upload a Signature"){
#  for (i in 1:length(con_comma))
#  {
#    if (length(con_comma[[i]])!=0){
#      rows_to_add <- which(Compounds$SM_LINCS_ID %in% con_comma[[i]])
#      compound_rows[[i]] <- rows_to_add
#    }
#    else {compound_rows[[i]] <- 0}
#  }
#  compound_rows <- compound_rows
#}
  
##################################
  
  
if (length(which(unlist(compound_rows)==0))!=0){
  removeemptyrows <- unlist(compound_rows)[-which(unlist(compound_rows)==0)]
  compound_rows <<- removeemptyrows
}
  else {
    compound_rows <<- compound_rows
  }
  
selected_compounds <- list("character")
for (i in 1:length(compound_rows))
  {
the_chosen_ones <- Compounds[compound_rows[[i]], c("SM_LINCS_ID", "SM_SMILES_Parent")]
if (length(which(the_chosen_ones$SM_SMILES_Parent==""))!=0){
selected_compounds[[i]] <- the_chosen_ones[-which(the_chosen_ones$SM_SMILES_Parent==""),]
}
else{
  selected_compounds[[i]] <- the_chosen_ones
}
}
selected_compounds <<- selected_compounds

########Get Unique
all_compounds <- ldply(selected_compounds, data.frame)
all_compounds <<- all_compounds
#dups <- which(duplicated(all_compounds[,1]))
#final_compounds <<- all_compounds[-dups,]
final_compounds <<- unique(all_compounds)
}
