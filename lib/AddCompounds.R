AddCompounds <- function(SMILESFile)
  {
  AddedSMILES <- read.SMIset(SMILESFile)
  AddedSDF<-smiles2sdf(AddedSMILES)
}