#=============================
# 图片保存路径
#=============================

save_path <- "D:/r studio 4.6.0/RF_Result/"

# 不存在这个文件夹，则自动创建
if(!dir.exists(save_path)){
  dir.create(save_path)
}

###########################################################
# Random Forest 二分类疾病预测模型
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

# 列名
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

# 保存roc曲线
roc_df <- data.frame(
  FPR = 1 - roc_rf$specificities,
  TPR = roc_rf$sensitivities
)

p1 <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(linewidth=1.5, colour="#D55E00") +
  geom_abline(intercept=0, slope=1,
              linetype=2,
              colour="grey60") +
  theme_classic(base_size=16) +
  labs(title="Random Forest ROC Curve",
       x="False Positive Rate",
       y="True Positive Rate") +
  annotate("text",
           x=0.65,
           y=0.15,
           label=paste("AUC =", round(as.numeric(auc_rf),3)),
           size=6)

ggsave(
  filename=paste0(save_path,"ROC_RF.png"),
  plot=p1,
  width=6,
  height=5,
  dpi=600
)

library(reshape2)

cm <- table(True=test$Group,
            Predicted=pred)

cm_df <- melt(cm)

p2 <- ggplot(cm_df,
             aes(Predicted,True))+
  geom_tile(aes(fill=value),
            colour="white")+
  geom_text(aes(label=value),
            size=8,
            fontface="bold")+
  scale_fill_viridis_c()+
  theme_classic(base_size=16)+
  labs(title="Confusion Matrix")

ggsave(
  filename=paste0(save_path,"Confusion_Matrix.png"),
  plot=p2,
  width=5,
  height=5,
  dpi=600
)

# 保存混淆矩阵
library(reshape2)
library(ggplot2)

#=========================
# 混淆矩阵
#=========================
cm <- table(
  True = test$Group,
  Predicted = pred
)

cm_df <- melt(cm)

# 将0和1改成文字标签
cm_df$True <- factor(
  cm_df$True,
  levels = c(0,1),
  labels = c("Group 0","Group 1")
)

cm_df$Predicted <- factor(
  cm_df$Predicted,
  levels = c(0,1),
  labels = c("Group 0","Group 1")
)

# 绘图
p2 <- ggplot(
  cm_df,
  aes(x = Predicted,
      y = True,
      fill = value)
) +
  
  geom_tile(
    color = "white",
    linewidth = 1
  ) +
  
  geom_text(
    aes(label = value),
    size = 10,
    fontface = "bold"
  ) +
  
  scale_fill_gradient(
    low = "#EAF2F8",
    high = "#2166AC"
  ) +
  
  labs(
    title = "Confusion Matrix",
    x = "Predicted Class",
    y = "True Class"
  ) +
  
  theme_classic(base_size = 18) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 15
    ),
    
    legend.position = "none"
    
  )

# 保存图片
ggsave(
  filename = paste0(save_path,"Confusion_Matrix.png"),
  plot = p2,
  width = 6,
  height = 5,
  dpi = 600
)
# 保存top20重要性
imp <- importance(rf_model)

imp <- data.frame(imp)

imp$Feature <- rownames(imp)

imp <- imp[order(imp$MeanDecreaseGini,
                 decreasing=TRUE),]

imp20 <- head(imp,20)

p3 <- ggplot(
  imp20,
  aes(reorder(Feature,
              MeanDecreaseGini),
      MeanDecreaseGini)
)+
  geom_col(fill="#4DBBD5")+
  coord_flip()+
  theme_classic(base_size=16)+
  labs(title="Top20 Important Features",
       x="",
       y="Mean Decrease Gini")

ggsave(
  filename=paste0(save_path,"Top20_Importance.png"),
  plot=p3,
  width=8,
  height=7,
  dpi=600
)

#保存性能模型图
library(caret)

cm <- confusionMatrix(pred,test$Group)

acc <- mean(pred==test$Group)

sen <- as.numeric(cm$byClass["Sensitivity"])

spe <- as.numeric(cm$byClass["Specificity"])

result <- data.frame(
  
  Metric=c(
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "AUC"),
  
  Value=c(
    acc,
    sen,
    spe,
    as.numeric(auc_rf))
)

p4 <- ggplot(result,
             aes(Metric,
                 Value,
                 fill=Metric))+
  geom_col(width=0.6)+
  geom_text(aes(label=round(Value,3)),
            vjust=-0.5,
            size=6)+
  ylim(0,1.1)+
  theme_classic(base_size=16)+
  theme(legend.position="none")

ggsave(
  filename=paste0(save_path,"Model_Performance.png"),
  plot=p4,
  width=6,
  height=5,
  dpi=600
)

cat("所有图片已保存至：\n")
