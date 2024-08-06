library("pROC")
library("tidyverse")

ROCStatFunc <- function(dat, group, var,retype = c("threshold", "specificity", "sensitivity"),
                        auc = T,youden = T, digit = 3){
  subgroup <- levels(as.factor(dat[[group]]))
  subgroup1 <- paste0(subgroup[2], " vs ", subgroup[1])
  rocmodel <- roc(dat[[group]], dat[[var]])
  other <- coords(rocmodel, "b", ret = retype)
  other <- round(other, digit)
  if(auc == T){
    auc <- round(ci.auc(rocmodel),digit)
    auc <- paste0(auc[2],"(",auc[1],"-",auc[3],")")
    if(youden == T){
      abc <- coords(rocmodel, "b", ret = c("specificity", "sensitivity"))
      youdenres <- abc[1] + abc[2] - 1
      youdenres <- round(youdenres, digit)
      result <- c(group, subgroup1, auc, other, youdenres)
      names(result) <- c("group", "subgroup","auc(95%CI)", retype, "youden")
    }else{
      result <- c(group, subgroup1, auc, other)
      names(result) <- c("group", "subgroup", "auc(95%CI)", retype)
    }
  }else{
    if(youden == T){
      abc <- coords(rocmodel, "b", ret = c("specificity", "sensitivity"))
      youdenres <- abc[1] + abc[2] - 1
      youdenres <- round(youdenres, digit)
      result <- c(group, subgroup1, other, youdenres)
      names(result) <- c("group","subgroup", retype, "youden")
    }else{
      result <- c(group, subgroup1,other)
      names(result) <- c("group", "subgroup",retype)
    }
  }
  return(result)
}
quiteROCFunc <- quietly(ROCStatFunc)

DATAset1 = read.table("Input_ROC_Kim.txt",header = TRUE)
DATAset2 = read.table("Input_ROC_Gide.txt",header = TRUE)
DATAset3 = read.table("Input_ROC_Liu.txt",header = TRUE)

roc1 <- roc(DATAset1$Response, DATAset1$NETscore,smooth = F)
roc2 <- roc(DATAset2$Response, DATAset2$NETscore,smooth = F)
roc3 <- roc(DATAset3$Response, DATAset3$NETscore,smooth = F)
roc.test(roc2, roc3)
coords(roc1, "b", ret="t")
ggroc(list(Kim=roc1,Gide=roc2,Liu=roc3))

AUC_Noh= quiteROCFunc(DATAset1, group = "Response",var = "NETscore")$result
AUC_IMvigor210 = quiteROCFunc(DATAset2, group = "Response", var = "NETscore")$result
AUC_Liu = quiteROCFunc(DATAset3, group = "Response", var = "NETscore")$result

Matrix= rbind(AUC_Noh,AUC_IMvigor210,AUC_Liu)
write.table(Matrix,file = "AUC.txt")

