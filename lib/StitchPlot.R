####Not Working Currently

PlotSTITCH <- function(clusternumber){
  stitchinput <- ""
  ClusterMembers <- as.data.frame(ClusterMembers)
  ClusterMembers$labels <- rownames(ClusterMembers)
  for (i in 1:length(ClusterMembers)){
    if (ClusterMembers$clusterID[[i]]==clusternumber){
      stitchinput <- paste(stitchinput, rownames(ClusterMembers$labels[[i]]), sep="%0D")
      print(ClusterMembers$labels[[i]])
    }
  }
  print(stitchinput)
}

#http://stitch.embl.de/api/image/networkList?identifiers=nilotinib%0Dmycn&species=9606&required_score=950