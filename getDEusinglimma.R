#get DE using limma
getDEusinglimma<-function(expre,group_list,outFile){
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(expre)
  contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
  ##step1
  fit <- lmFit(expre,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  ##step3
  DE_BH = topTable(fit2, coef=1,p.value=1.1, n=Inf,adjust.method="BH")
  nr_DE_BH = na.omit(DE_BH)
  DE_BH.005<-nr_DE_BH[nr_DE_BH$adj.P.Val<0.05,]
  outFile_BH<-paste(outFile,".BH.005.csv",sep="")
  write.csv(DE_BH.005,outFile_BH,quote = F)
}