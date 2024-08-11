ClueIO <- function(Concordant_signature)
  {
    url3 <<- paste("https://api.clue.io/api/perts?filter={%22fields%22:[%22canonical_smiles%22,%22pert_iname%22],%22where%22:{%22pert_id%22:{%22inq%22:[%22",Concordant_signature,"%22]}}}&user_key=e201f96e629b5741375fd1c58a46e3d6", sep="")
    raw.result3 <- GET(url = url3)
    raw.result3$status_code
    x3<<-content(raw.result3)
    }