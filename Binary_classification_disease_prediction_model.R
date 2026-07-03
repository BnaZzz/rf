###########################################################
# Random Forest 疾病预测模型
###########################################################

#=============================
# 1. 加载R包
#=============================
library(readxl)
library(randomForest)
library(caret)
library(pROC)

#=============================
# 2. 读取数据
#=============================
data <- read_excel("D:/r studio 4.6.0/template-0608 修改版.xlsx")

# 删除第一行（如果第一行为说明信息）
data <- data[-1, ]

#=============================
# 3. 删除无关列（根据实际情况修改）
#=============================
remove_cols <- c(
  "Sample number",
  "Sample name",
  "Other phenotypic information",
  "…(Other phenotypic information)"
)

remove_cols <- intersect(remove_cols, names(data))

data2 <- data[, !(names(data) %in% remove_cols)]

#=============================
# 4. 数据类型转换
#=============================
data2[] <- lapply(data2, function(x){
  suppressWarnings(as.numeric(as.character(x)))
})

# 删除缺失值
data2 <- na.omit(data2)

#=============================
# 5. 设置标签(Group)
#=============================
data2$Group <- factor(data2$Group)

cat("样本数：", nrow(data2), "\n")
print(table(data2$Group))

#=============================
# 6. 找到代谢物列
#=============================
met_cols <- grep("^Metabolite", names(data2))

# 如果你的列名不是Metabolite开头，可以注释掉下面两步
if(length(met_cols)>0){
  
  # log2转换
  data2[,met_cols] <- log2(data2[,met_cols]+1)
  
  # 标准化
  data2[,met_cols] <- scale(data2[,met_cols])
  
}

#=============================
# 7. 划分训练集和测试集
#=============================
set.seed(123)

trainIndex <- createDataPartition(
  data2$Group,
  p=0.8,
  list=FALSE
)

train <- data2[trainIndex,]

test <- data2[-trainIndex,]

cat("\n训练集：\n")
print(table(train$Group))

cat("\n测试集：\n")
print(table(test$Group))

#=============================
# 8. 建立随机森林
#=============================
set.seed(123)

rf_model <- randomForest(
  Group~.,
  data=train,
  ntree=1000,
  importance=TRUE
)

cat("\n=========================\n")
print(rf_model)

#=============================
# 9. 预测
#=============================
pred <- predict(rf_model,test)

prob <- predict(
  rf_model,
  test,
  type="prob"
)[,2]

#=============================
# 10. 混淆矩阵
#=============================
cm <- confusionMatrix(
  pred,
  test$Group
)

cat("\n=========================\n")
print(cm)

#=============================
# 11. Accuracy
#=============================
acc <- mean(pred==test$Group)

cat("\nAccuracy =",acc,"\n")

#=============================
# 12. ROC
#=============================
roc_rf <- roc(
  response=test$Group,
  predictor=prob
)

auc_rf <- auc(roc_rf)

cat("\nAUC =",auc_rf,"\n")

plot(
  roc_rf,
  print.auc=TRUE,
  main="Random Forest ROC"
)

#=============================
# 13. 最佳Cutoff
#=============================
best <- coords(
  roc_rf,
  "best",
  ret=c(
    "threshold",
    "sensitivity",
    "specificity"
  )
)

cat("\n最佳Cutoff：\n")
print(best)

#=============================
# 14. 变量重要性
#=============================
imp <- importance(rf_model)

imp <- data.frame(imp)

imp <- imp[
  order(
    imp$MeanDecreaseGini,
    decreasing=TRUE
  ),
]

cat("\nTop20重要变量：\n")
print(head(imp,20))

#=============================
# 15. 变量重要性图
#=============================
varImpPlot(
  rf_model,
  n.var=20,
  main="Top20 Important Features"
)

#=============================
# 16. Top10变量ROC
#=============================
top10 <- rownames(head(imp,10))

auc_result <- data.frame()

for(i in top10){
  
  roc_obj <- roc(
    data2$Group,
    data2[[i]]
  )
  
  auc_result <- rbind(
    auc_result,
    data.frame(
      Feature=i,
      AUC=as.numeric(auc(roc_obj))
    )
  )
  
}

auc_result <- auc_result[
  order(
    auc_result$AUC,
    decreasing=TRUE
  ),
]

cat("\nTop10变量AUC：\n")
print(auc_result)

#=============================
# 17. Top10变量t检验
#=============================
cat("\n=========================\n")
cat("Top10变量t检验\n")

for(i in top10){
  
  cat("\n-----------------------\n")
  cat(i,"\n")
  
  print(
    t.test(
      data2[data2$Group==0,i],
      data2[data2$Group==1,i]
    )
  )
  
}

###########################################################
# End
###########################################################
