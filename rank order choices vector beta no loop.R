# sample code for rank order choices

if (!require("vctrs")) install.packages("vctrs"); library(vctrs)
if (!require("pracma")) install.packages("pracma"); library(pracma)
if (!require("rje")) install.packages("rje"); library(rje)
if (!require("combinat")) install.packages("combinat"); library(combinat)
if (!require("EnvStats")) install.packages("EnvStats"); library(EnvStats)
if (!require("matrixStats")) install.packages("matrixStats"); library(matrixStats)


n                            <- 1000                    #number of observations
size_of_choice_set           <- 3                       #number of alternatives
parameter_of_logit_attention <- 5                       #parameter of logit attention model

# what we want to estimate
true_beta                    <- c(1, 2)                 #parameter of interest, true value

empty_list <- list()

for (o in 1:1){

value_added                 <- rep(c(3, 0, 0, 1, 0, 0), n)

### monte carlo starts here
number_of_simulation <- 1
for (m in 1:number_of_simulation){
cat("Current simulation:", m, "\n")

tic()
  
# what we can observe
# every size_of_choice_set * length(true_beta) elements are attributes of a single individual
# for example, if length(true_beta) = 2, size_of_choice_set = 3
# the first 6 elements are c(2, 3, 5, 1, 0, 2) means that for the first subject, his attributes of school 1 is (2, 3)
# (5, 1) for school 2, (0, 2) for school 3
parameter_of_choice_set     <- sample.int(5,size=n * size_of_choice_set * length(true_beta),replace=TRUE)
# add some preference to the first dimension of school 1 and the second dimension of school 2

parameter_of_choice_set     <- parameter_of_choice_set + value_added

# generate a matrix recording all possible consideration set
for (i in 2:size_of_choice_set) {
  temp_all_menu <- t(combn(size_of_choice_set, i))
  # for set size i, find all combinations
  
  temp_menu <- matrix(0, nrow = nrow(temp_all_menu), ncol = size_of_choice_set)
  # initialize temp_menu to be appended
  
  for (j in 1:nrow(temp_all_menu)) {
      temp_menu[j, temp_all_menu[j, ]] <- 1
  }
  if (i == 2) {
    menu <- temp_menu
  }
  else {
    menu <- rbind(menu, temp_menu) 
  }
}

# given all possible consideration set, calculate the weight of attention
length_of_menu <- rowSums(menu)
weight_of_attention <- length_of_menu ** parameter_of_logit_attention
weight_of_attention <- weight_of_attention / sum(weight_of_attention)

# recording the realization of consideration set for each individual
index_of_consideration_set <- rmultinom(n, 1, weight_of_attention)

# map index to real consideration set
consideration_set <- t(menu) %*% index_of_consideration_set
consideration_set <- t(consideration_set)

# now we calculate individual's utility of each alternative
# the unobservable part follows extreme value type 1 distribution
unobservable <- revd(n * size_of_choice_set, location = 0, scale = 1)
unobservable <- matrix(unobservable, nrow = n, byrow = TRUE)
# base utility, linear combination of attributes according to true beta
base_utility <- parameter_of_choice_set * true_beta
base_utility <- colSums(matrix(base_utility, nrow=length(true_beta)))
base_utility <- matrix(base_utility, nrow = n, byrow = TRUE)
# sum them up
utility <- base_utility + unobservable

# map utilities to ranks
rank_order <- rowRanks(-utility)
# drop those not in the consideration set
rank_order <- rank_order * consideration_set
# replace 0 with NA
rank_order[rank_order == 0] <- NA
# reorder those ranks
rank_order <- rowRanks(rank_order)

################################################################################
################################################################################
################################################################################
###  Above are data generation process, where we can only observe
###  rank_order, parameter_of_choice_set



# find unique rank orders and their frequencies
unique_rank_order <- unique(rank_order)
data_prob_rank_order <- matrix(0, nrow(unique_rank_order), 1)
for (i in 1:nrow(unique_rank_order)) {
  for (j in 1:nrow(rank_order)){
    data_prob_rank_order[i,] = data_prob_rank_order[i,] + identical(unique_rank_order[i,], rank_order[j,])
  }
}

data_prob_rank_order <- data_prob_rank_order / n

# generate a matrix recording compatibility between rank orders and full rankings
# each row of full_rank_order record ranks of choices, for example, (3, 1, 2) means that school 1 is the 3rd choice
# school 2 is the best choice and school 3 is the second choice
full_rank_order <- permn(size_of_choice_set)
full_rank_order <- do.call(rbind, full_rank_order)
compatibility_matrix <- matrix(0, nrow(unique_rank_order), nrow(full_rank_order))
for (i in 1:nrow(compatibility_matrix)) {
  for (j in 1:ncol(compatibility_matrix)){
    temp_matrix <- rbind(unique_rank_order[i,], full_rank_order[j,])
    temp_matrix <- temp_matrix[ , colSums(is.na(temp_matrix))==0]
    temp_matrix <- matrix(temp_matrix,nrow=2)
    temp_matrix[2,] <- rank(temp_matrix[2,])
    compatibility_matrix[i,j] <- identical(temp_matrix[1,], temp_matrix[2,])
  }
}

### matrix product to calculate data prob of inequality

full_rank_order_RHS <- powerSetMat(nrow(full_rank_order))
full_rank_order_RHS <- full_rank_order_RHS[2:nrow(full_rank_order_RHS),]
num_of_inequality <- nrow(full_rank_order_RHS)

rank_order_LHS <- compatibility_matrix %*% t(full_rank_order_RHS)
rank_order_LHS <- (rank_order_LHS >= 1) * 1

data_prob_inequality <- t(rank_order_LHS) %*% data_prob_rank_order
data_prob_inequality <- c(data_prob_inequality)


### matrix product to calculate model prob of inequality

# generate a matrix to reorder utilities according to all possible full ranks, for example, base utility is 5,3,2 for
# three schools, and a full rank order is 3,1,2 then this matrix is used to reorder 5,3,2 as 3,2,5 (exhuasting all full ranks)
# it is almost a diagonal matrix

matrix_reorder_utility <- matrix(0, nrow = size_of_choice_set * nrow(full_rank_order), ncol = size_of_choice_set * nrow(full_rank_order))

for (i in 1:nrow(full_rank_order)){
  for (j in 1:size_of_choice_set){
    matrix_reorder_utility[size_of_choice_set * (i-1) + j, size_of_choice_set * (i-1) + full_rank_order[i,j]] <- 1
  }
}

# generate a matrix to turn e^u_1, e^u_2, e^u_3 into 
#e^u_1+e^u_2+e^u_3, e^u_2+e^u_3, e^u_3

matrix_calculate_denominator <- diag(size_of_choice_set * nrow(full_rank_order))

for (i in 1:nrow(full_rank_order)){
  matrix_calculate_denominator[(size_of_choice_set * (i-1) + 1):(size_of_choice_set * i), (size_of_choice_set * (i-1) + 1):(size_of_choice_set * i)] <- upper.tri(diag(size_of_choice_set), diag = TRUE)
}



#calculate model implied probability
beta_1 <- linspace(-5, 10, 151)
beta_2 <- linspace(-5, 10, 151)
temp_beta <- t(data.matrix(expand.grid(beta_1, beta_2)))
num_of_beta <- ncol(temp_beta)
    

#step 1: calculate numerator by base utility, which is n * (size_of_choice_set * nrow(full_rank_order))

numerator <- matrix(parameter_of_choice_set, ncol = nrow(temp_beta), byrow = TRUE) %*% temp_beta
numerator <- exp(numerator)
numerator <- matrix(numerator, ncol = size_of_choice_set, byrow = TRUE)
numerator <- matrix(rep(numerator, nrow(full_rank_order)), nrow = nrow(numerator))

#step 2: reorder those base utilities according to reported ranks
#numerator finished, every size_of_choice_set column are e^u_choice1, e^u_choice2, etc for a specific full rank
#every n row record results for a specific candidate of beta
#the whole matrix records those for all possible full ranks for all betas

numerator <- numerator %*% matrix_reorder_utility

#step 3: calculate denominator based on numerator, which is turning e^u_1, e^u_2, e^u_3 into 
#e^u_1+e^u_2+e^u_3, e^u_2+e^u_3, e^u_3

denominator <- t(matrix_calculate_denominator %*% t(numerator))

#step 4, take product to calculate choice probability
model_prob_inequality <- numerator / denominator
#calculate e^u_1 / (e^u_1+e^u_2+e^u_3) * e^u_2 / (e^u_2+e^u_3)
model_prob_inequality <- c(t(model_prob_inequality))
    
#each row is e^u_1 / (e^u_1+e^u_2+e^u_3), e^u_2 / (e^u_2+e^u_3), e^u_3 / (e^u_3), for a specific individual, a specific beta and a specific rank
#every nrow(full_rank_order) row record the above for all full ranks
model_prob_inequality <- matrix(model_prob_inequality, ncol = size_of_choice_set, byrow = TRUE)
    
#take row product of the above matrix and put results for all rank orders in 1 row

for (i in 1:size_of_choice_set){
  if (i == 1){
    temp_model_prob_inequality <- model_prob_inequality[,1]
  }
  if (i > 1){
    temp_model_prob_inequality <- temp_model_prob_inequality * model_prob_inequality[,i]
  }
}

model_prob_inequality <- temp_model_prob_inequality


# an alternative way, which is slower for some unclear reason
# model_prob_inequality <- rowProds(model_prob_inequality)



model_prob_inequality <- matrix(model_prob_inequality, ncol = nrow(full_rank_order), byrow = TRUE)
#so each row has nrow(full_rank_order) elements, every element stand for prob of a full rank for a specific individual and a specific beta
#current size (n*num_of_beta) * nrow(full_rank_order)


#take mean across n individuals
model_prob_inequality <- rowsum(model_prob_inequality, rep(1:num_of_beta, each=n)) / n
model_prob_inequality <- model_prob_inequality %*% t(full_rank_order_RHS)



#avoid rounding issues
model_prob_inequality[,ncol(model_prob_inequality)] <- 1
data_prob_inequality[length(data_prob_inequality)] <- 1

data_prob_inequality <- t(replicate(num_of_beta, data_prob_inequality))
violation_inequality <- (data_prob_inequality < model_prob_inequality) * 1

#find inequalities never violated

never_violated_inequality <- (colSums(violation_inequality) == 0) * 1

index_of_never_violated_inequality <- setdiff(unique(never_violated_inequality * c(1:num_of_inequality)), c(0))


# now we need to see if some columns are always larger than others
# column i is always weakly larger than column j means that
# every time inequality j is violated, i is also violated, so it suffices to check i only

matrix_redundant_inequality <- matrix(0, num_of_inequality, num_of_inequality)
for (i in 1:num_of_inequality){
  temp_matrix <- replicate(num_of_inequality, violation_inequality[,i])
  temp_vector <- colSums(temp_matrix >= violation_inequality)
  matrix_redundant_inequality[i,] <- (temp_vector == num_of_beta) - never_violated_inequality
}

matrix_redundant_inequality <- matrix_redundant_inequality * (1 - diag(num_of_inequality))

index_of_never_violated_inequality

inequality_informative <- matrix(0, num_of_inequality, 2)
inequality_informative[, 1] <- c(1:num_of_inequality)
inequality_informative[, 2] <- rowSums(matrix_redundant_inequality)

if (m == 1){
  matrix_redundant_inequality_intersection <- matrix_redundant_inequality
}

if (m > 1){
  matrix_redundant_inequality_intersection <- matrix_redundant_inequality_intersection * matrix_redundant_inequality
}

toc()

} #close loop for m

empty_list[[o]] <- which(matrix_redundant_inequality_intersection !=0, arr.ind = T)

} #close loop for o

empty_list

containment_RHS <- 1 - full_rank_order_RHS
containment_LHS <- 1 - rank_order_LHS

# which ineq imply other
# uniquely violated
# do the exercise for fixed value_added and beta (core determining set as a function of beta)

# check inequality 1 first, giving a set of beta
# check 2, 3 etc in the identified set

# plot those beta, with different colors and draw the convex hull



index_of_violation_inequality <- list()
for (i in 1:num_of_beta){
  index_of_violation_inequality[[i]] <- unique(violation_inequality[i,] * c(1:num_of_inequality))
}

num_of_violation_inequality <- matrix(rowSums(violation_inequality), nrow = length(beta_2), ncol = length(beta_1))


for (i in (65*151+1):(66*151)){
  print(index_of_violation_inequality[[i]])
}