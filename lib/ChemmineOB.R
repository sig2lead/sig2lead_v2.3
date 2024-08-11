#library(ChemmineOB)
#library(ChemmineR)

Get_SDF<-function()
{
var<-as(final_compounds$SM_SMILES_Parent, "SMIset")
cid(var)<-final_compounds$SM_LINCS_ID

###Remove any illegal aromatic elements
illegal_aromatic <- grep("[i+]", var)
var <- var[-illegal_aromatic,]
##########################

#write.SMI(var, file="smiles_smi.smi")
#sdf_smiles <<- smiles2sdf(var)
sdf_smiles <<- var
}

#rm(sdf_smiles)
