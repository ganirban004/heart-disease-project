# REQUIRED LIBRARIES
library(dplyr)
library(ggcorrplot)
library(car)
library(MASS)
library(caret)
library(yardstick)
library(e1071)
library(randomForest)
library(naivebayes)
library(kernlab)
library(klaR)


######## DATA PRE-PROCESSING ########
# IMPORTING HEART DISEASE DATA
heart_df <- read.csv('heart.csv')

# OMIITING MISSING VALUES
sum(is.na(heart_df))

# REPLACING THE 0 VALUES OF THE RestingBP BY ITS MEDIAN
heart_df <- heart_df %>%
  mutate(RestingBP=replace(RestingBP, RestingBP==0, median(heart_df$RestingBP)))
heart_df$RestingBP <- as.integer(heart_df$RestingBP)

# REPLACING THE 0 VALUES OF THE CHOLESTEROL BY ITS MEDIAN
heart_df <- heart_df %>%
  mutate(Cholesterol=replace(Cholesterol, Cholesterol==0, median(heart_df$Cholesterol)))
heart_df$Cholesterol <- as.integer(heart_df$Cholesterol)

# CHANGING NEGETIVE VALUES OF OLDPEAK TO POSITIVE
heart_df$Oldpeak <- abs(heart_df$Oldpeak)

# DROPPING FastingBS COLUMN
heart_df$FastingBS <- NULL

# OUTLIER DETECTION & REMOVAL
MD <- mahalanobis(heart_df[, c(1, 4, 5, 7, 9)], colMeans(heart_df[, c(1, 4, 5, 7, 9)]), cov(heart_df[, c(1, 4, 5, 7, 9)]))
heart_df$MD <- round(MD, 3)
heart_df$is_outlier <- FALSE
heart_df$is_outlier[heart_df$MD > 12] <- TRUE
data.frame(table(heart_df$is_outlier))
heart_df <- heart_df[!(heart_df$is_outlier == TRUE), ]
heart_df$MD <- NULL
heart_df$is_outlier <- NULL

# CODING ATTRIBUTES TO VARIABLES
# Sex
heart_df_coded <- heart_df %>%
  mutate(Sex=replace(Sex, Sex=='M', 1))
heart_df_coded <- heart_df_coded %>%
  mutate(Sex=replace(Sex, Sex=='F', 2))
heart_df_coded$Sex <- as.integer(heart_df_coded$Sex)
# ChestPainType
heart_df_coded <- heart_df_coded %>%
  mutate(ChestPainType=replace(ChestPainType, ChestPainType=='ATA', 1))
heart_df_coded <- heart_df_coded %>%
  mutate(ChestPainType=replace(ChestPainType, ChestPainType=='NAP', 2))
heart_df_coded <- heart_df_coded %>%
  mutate(ChestPainType=replace(ChestPainType, ChestPainType=='ASY', 3))
heart_df_coded <- heart_df_coded %>%
  mutate(ChestPainType=replace(ChestPainType, ChestPainType=='TA', 4))
heart_df_coded$ChestPainType <- as.integer(heart_df_coded$ChestPainType)
# RestingECG
heart_df_coded <- heart_df_coded %>%
  mutate(RestingECG=replace(RestingECG, RestingECG=='Normal', 1))
heart_df_coded <- heart_df_coded %>%
  mutate(RestingECG=replace(RestingECG, RestingECG=='ST', 2))
heart_df_coded <- heart_df_coded %>%
  mutate(RestingECG=replace(RestingECG, RestingECG=='LVH', 3))
heart_df_coded$RestingECG <- as.integer(heart_df_coded$RestingECG)
# ExerciseAngina
heart_df_coded <- heart_df_coded %>%
  mutate(ExerciseAngina=replace(ExerciseAngina, ExerciseAngina=='Y', 1))
heart_df_coded <- heart_df_coded %>%
  mutate(ExerciseAngina=replace(ExerciseAngina, ExerciseAngina=='N', 0))
heart_df_coded$ExerciseAngina <- as.integer(heart_df_coded$ExerciseAngina)
# ST_Slope
heart_df_coded <- heart_df_coded %>%
  mutate(ST_Slope=replace(ST_Slope, ST_Slope=='Down', -1))
heart_df_coded <- heart_df_coded %>%
  mutate(ST_Slope=replace(ST_Slope, ST_Slope=='Flat', 0))
heart_df_coded <- heart_df_coded %>%
  mutate(ST_Slope=replace(ST_Slope, ST_Slope=='Up', 1))
heart_df_coded$ST_Slope <- as.integer(heart_df_coded$ST_Slope)


######## VISUALIZATION ########
# CORRELATION MATRIX
cor_mat <- round(cor(heart_df_coded), 2)
head(cor_mat)
ggcorrplot(cor_mat)

# BAR DIAGRAMS AND PIE CHARTS
# Age
hist(heart_df_coded$Age, main='Age', xlab = 'Age', col='green', border='black')
# Sex
pie(table(heart_df_coded$Sex), labels = paste(round(prop.table(table(heart_df_coded$Sex))*100, 1), '%', sep = ""), col = heat.colors(2), main = 'Sex')
legend('topright', c('M','F'), fill = heat.colors(2))
# ChestPainType
pie(table(heart_df_coded$ChestPainType), labels = paste(round(prop.table(table(heart_df_coded$ChestPainType))*100, 1), '%', sep = ""), col = heat.colors(4), main = 'ChestPainType')
legend("topright", c('ATA', 'NAP', 'ASY', 'TA'), fill = heat.colors(4))
# RestingBP
hist(heart_df_coded$RestingBP, main='RestingBP', xlab = 'RestingBP', col='green', border='black')
# Cholesterol
hist(heart_df_coded$Cholesterol, main='Cholesterol', xlab = 'Cholesterol', col='green', border='black')
# RestingECG
pie(table(heart_df_coded$RestingECG), labels = paste(round(prop.table(table(heart_df_coded$RestingECG))*100, 1), '%', sep = ""), col = heat.colors(3), main = 'RestingECG')
legend("topright", c('Normal', 'ST', 'LVH'), fill = heat.colors(3))
# MaxHR
hist(heart_df_coded$MaxHR, main='MaxHR', xlab = 'MaxHR', col='green', border='black')
# ExerciseAngina
pie(table(heart_df_coded$ExerciseAngina), labels = paste(round(prop.table(table(heart_df_coded$ExerciseAngina))*100, 1), '%', sep = ""), col = heat.colors(2), main = 'ExerciseAngina')
legend("topright", c('N', 'Y'), fill = heat.colors(2))
# Oldpeak
hist(heart_df_coded$Oldpeak, main='Oldpeak', xlab = 'Oldpeak', col='green', border='black')
# ST_Slope
pie(table(heart_df_coded$ST_Slope), labels = paste(round(prop.table(table(heart_df_coded$ST_Slope))*100, 1), '%', sep = ""), col = heat.colors(3), main = 'ST_Slope')
legend("topright", c('Down', 'Flat', 'Up'), fill = heat.colors(3))
#HeartDisease
pie(table(heart_df_coded$HeartDisease), labels = paste(round(prop.table(table(heart_df_coded$HeartDisease))*100, 1), '%', sep = ""), col = heat.colors(2), main = 'HeartDisease')
legend("topright", c('0', '1'), fill = heat.colors(2))


######## VIF CALCULATION ########
# FITTING LINEAR REGRESSION
LiR_model <- lm(HeartDisease~., data = heart_df_coded)
summary(LiR_model)

# CALCULATING VIF
vif(LiR_model)


######## PRINCIPAL COMPONENT ANALYSIS ########
HeartDisease_pc <- prcomp(heart_df_coded, center = TRUE, scale = TRUE)
summary(HeartDisease_pc)
pc_summary <- summary(HeartDisease_pc)
pc_summary$importance


######## TRAIN-TEST SPLITTING ########
# SPLITTING DATA INTO TRAIN AND TEST SET
set.seed(2345)
rows <- sample(nrow(heart_df_coded))
heart_df_shuffled <- heart_df_coded[rows, ]
size_train_set <- round(nrow(heart_df_coded)*0.8, 0)
train_set <- heart_df_shuffled[1:size_train_set, ]
test_set <- heart_df_shuffled[(size_train_set+1):nrow(heart_df_coded), ]

# DESCRIPTION OF TRAIN AND TEST SET
nrow(train_set)
nrow(test_set)
nrow(train_set)/nrow(heart_df_coded)*100


######## LOGISTIC REGRESSION ########
# FITTING LOGISTIC REGRESSION WITH TRAIN DATA
model_LR <- glm(HeartDisease ~., data = train_set, family = 'binomial')
summary(model_LR)

# MODIFICATION OF MODEL USING stepAIC FUNCTION
modified_model_LR <- stepAIC(model_LR)
summary(modified_model_LR)

# PREDICTING TEST DATA USING THE LOGISTIC REGRESSION MODEL
pred_value_test <- predict(modified_model_LR, test_set, type = 'response')
pred_value_test
pred_class_test <- ifelse(pred_value_test>=0.5, 1, 0)
pred_class_test

# CREATING CONFUSION MATRIX FOR TEST DATA
conf_mat_test <- confusionMatrix(as.factor(test_set$HeartDisease), as.factor(pred_class_test), positive='1')
conf_mat_test
fourfoldplot(conf_mat_test$table, color = c('cyan', 'pink'), conf.level = 0, margin = 1, main = 'Confusion Matrix')

# CALCULATING RESIDUAL SUM OF SQUARES FOR TEST DATA
rss_test <- sum((test_set$HeartDisease-pred_class_test)**2)
rss_test

# PERFORMANCE METRIC FOR TEST DATA
performance <- metric_set(accuracy, precision, f_meas, recall)
performance(data = data.frame(test_set$HeartDisease, pred_class_test),
            estimate = as.factor(pred_class_test), truth = as.factor(test_set$HeartDisease))


######## SUPPORT VECTOR MACHINE ########
# FITTING SUPPORT VECTOR MACHINE WITH TRAIN DATA
model_SVM <- svm(HeartDisease~., data = train_set, type = 'C')
summary(model_SVM)

# PREDICTING TEST DATA USING THE SUPPORT VECTOR MACHINE MODEL
pred_value_test <- predict(model_SVM, test_set)  #gives values as factor 1 and 2
pred_value_test
pred_class_test <- ifelse(as.integer(pred_value_test)==1, 0, 1) #so we need to change it as number to find rss
pred_class_test

# CREATING CONFUSION MATRIX FOR TEST DATA
conf_mat_test <- confusionMatrix(as.factor(test_set$HeartDisease), as.factor(pred_class_test), positive='1')
conf_mat_test
fourfoldplot(conf_mat_test$table, color = c('cyan', 'pink'), conf.level = 0, margin = 1, main = 'Confusion Matrix')

# CALCULATING RESIDUAL SUM OF SQUARES FOR TEST DATA
rss_test <- sum((test_set$HeartDisease-pred_class_test)**2)
rss_test

# PERFORMANCE METRIC FOR TEST DATA
performance <- metric_set(accuracy, precision, f_meas, recall)
performance(data = data.frame(test_set$HeartDisease, pred_class_test),
            estimate = as.factor(pred_class_test), truth = as.factor(test_set$HeartDisease))


######## RANDOM FOREST MODEL ########
# FITTING RANDOM FOREST MODEL WITH TRAIN DATA
set.seed(1111)
model_RF <- randomForest(as.factor(HeartDisease)~., data = train_set)  #we need to change the target variable
summary(model_RF)                       #as factor to perform classification using randomForest function

# PREDICTING TEST DATA USING THE RANDOM FOREST MODEL
pred_value_test <- predict(model_RF, test_set)  #gives values as factor 1 and 2
pred_value_test
pred_class_test <- ifelse(as.integer(pred_value_test)==1, 0, 1) #so we need to change it as number to find rss
pred_class_test


# VARIABLE IMPORTANCE MEASUREMENT
importance(model_RF, type=2)

# CREATING CONFUSION MATRIX FOR TEST DATA
conf_mat_test <- confusionMatrix(as.factor(test_set$HeartDisease), as.factor(pred_class_test), positive='1')
conf_mat_test
fourfoldplot(conf_mat_test$table, color = c('cyan', 'pink'), conf.level = 0, margin = 1, main = 'Confusion Matrix')

# CALCULATING RESIDUAL SUM OF SQUARES FOR TEST DATA
rss_test <- sum((test_set$HeartDisease-pred_class_test)**2)
rss_test

# PERFORMANCE METRIC FOR TEST DATA
performance <- metric_set(accuracy, precision, f_meas, recall)
performance(data = data.frame(test_set$HeartDisease, pred_class_test),
            estimate = as.factor(pred_class_test), truth = as.factor(test_set$HeartDisease))


######## NAIVE BAYES CLASSIFIER ########
# FITTING NAIVE BAYES MODEL WITH TRAIN DATA
model_NB <- naive_bayes(as.factor(HeartDisease)~., data = train_set)  #we need to change the target variable
summary(model_NB)                       #as factor to perform naive_bayes function

# PREDICTING TEST DATA USING THE NAIVE BAYES CLASSIFIER
pred_value_test <- predict(model_NB, test_set, type = 'prob')  #gives values as factor 1 and 2
pred_value_test
pred_class_test <- pred_class_test <- ifelse(pred_value_test[ ,2]>=0.5, 1, 0)
pred_class_test

# CREATING CONFUSION MATRIX FOR TEST DATA
conf_mat_test <- confusionMatrix(as.factor(test_set$HeartDisease), as.factor(pred_class_test), positive='1')
conf_mat_test
fourfoldplot(conf_mat_test$table, color = c('cyan', 'pink'), conf.level = 0, margin = 1, main = 'Confusion Matrix')

# CALCULATING RESIDUAL SUM OF SQUARES FOR TEST DATA
rss_test <- sum((test_set$HeartDisease-pred_class_test)**2)
rss_test

# PERFORMANCE METRIC FOR TEST DATA
performance <- metric_set(accuracy, precision, f_meas, recall)
performance(data = data.frame(test_set$HeartDisease, pred_class_test),
            estimate = as.factor(pred_class_test), truth = as.factor(test_set$HeartDisease))


######## LINEAR DISCRIMINANT ANALYSIS ########
# FITTING LINEAR DISCRIMINANT ANALYSIS MODEL WITH TRAIN DATA
model_LDA <- lda(HeartDisease~., data = train_set)  
summary(model_LDA)                      

# PREDICTING TEST DATA USING THE LINEAR DISCRIMINANT ANALYSIS MODEL
pred_value_test <- predict(model_LDA, test_set)  #gives values as factor 1 and 2
pred_value_test
pred_class_test <- ifelse(as.integer(pred_value_test$class)==1, 0, 1) #so we need to change it as number to find rss
pred_class_test

# CREATING CONFUSION MATRIX FOR TEST DATA
conf_mat_test <- confusionMatrix(as.factor(test_set$HeartDisease), as.factor(pred_class_test), positive='1')
conf_mat_test
fourfoldplot(conf_mat_test$table, color = c('cyan', 'pink'), conf.level = 0, margin = 1, main = 'Confusion Matrix')


# CALCULATING RESIDUAL SUM OF SQUARES FOR TEST DATA
rss_test <- sum((test_set$HeartDisease-pred_class_test)**2)
rss_test


# PERFORMANCE METRIC FOR TEST DATA
performance <- metric_set(accuracy, precision, f_meas, recall)
performance(data = data.frame(test_set$HeartDisease, pred_class_test),
            estimate = as.factor(pred_class_test), truth = as.factor(test_set$HeartDisease))


######## CROSS-VALIDATION ########
# FITTING 10-FOLD CROSS-VALIDATION
control <- trainControl(method = 'repeatedcv', number = 10, repeats = 3)

# CROSS-VALIDATED LOGISTIC REGRESSION
set.seed(345)
model_LR_CV <- train(as.factor(HeartDisease) ~., data = heart_df_coded, trControl = control, method = 'glm', family = 'binomial')

# CROSS-VALIDATED SUPPORT VECTOR MACHINE
set.seed(345)
model_SVM_CV <- train(as.factor(HeartDisease) ~., data = heart_df_coded, trControl = control, method = 'svmRadial')

# CROSS-VALIDATED RANDOM FOREST MODEL
set.seed(345)

model_RF_CV <- train(as.factor(HeartDisease) ~., data = heart_df_coded, trControl = control, method = 'rf')

# CROSS-VALIDATED NAIVE BAYES MODEL
model_NB_CV <- train(as.factor(HeartDisease) ~., data = heart_df_coded, trControl = control, method = 'nb')

# CROSS-VALIDATED LINEAR DISCRIMINANT ANALYSIS MODEL
model_LDA_CV <- train(as.factor(HeartDisease) ~., data = heart_df_coded, trControl = control, method = 'lda')

# COMPARISON
results <- resamples(list(LR=model_LR_CV, SVM=model_SVM_CV, RF=model_RF_CV, LDA=model_LDA_CV))
summary(results)
bwplot(results)
dotplot(results)
