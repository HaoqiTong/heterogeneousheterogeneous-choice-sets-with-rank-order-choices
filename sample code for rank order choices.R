# sample code for rank order choices

# packages
# install.packages('EnvStats')
# install.packages('matrixStats')
# install.packages('combinat')
# install.packages('rje')
# install.packages('vctrs')
library(vctrs)
library(pracma)
library(rje)
library(combinat)
library(EnvStats)
library(matrixStats)

n                            <- 1000                    #number of observations
size_of_choice_set           <- 3                       #number of alternatives
parameter_of_logit_attention <- 5                       #parameter of logit attention model

# what we want to estimate
true_beta                    <- 2                       #parameter of interest, true value

### monte carlo starts here
number_of_simulation <- 20
results <- vector(mode = "list", 1)

for (m in 1:number_of_simulation){
  
cat("Current simulation:", m, "\n")

# what we can observe
parameter_of_choice_set      <- matrix(sample.int(5,size=n * size_of_choice_set,replace=TRUE), n, size_of_choice_set)
parameter_of_choice_set[,1] <- parameter_of_choice_set[,1] + 3

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
unobservable <- revd(n * size_of_choice_set, location = 0, scale = 1)
unobservable <- matrix(unobservable, nrow = n, byrow = TRUE)
base_utility <- true_beta * parameter_of_choice_set
utility <- base_utility + unobservable

# map utilities to ranks
rank_order <- rowRanks(-utility)
# drop those not in the consideration set
rank_order <- rank_order * consideration_set
# replace 0 with NA
rank_order[rank_order == 0] <- NA
# reorder those ranks
rank_order <- rowRanks(rank_order)



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

# in this example, it suffices to check 63 inequalities (the power set of 6 full rankings)
# consider Artstein's inequality (containment version)
# we first calculate the LHS of it, which is the data probability
# to find LHS, we need to first specify RHS, which is the combination of full rank
# first find all the possible combination of full rank
index_of_full_rank_order_RHS <- powerSet(1:nrow(full_rank_order))
index_of_full_rank_order_RHS <- index_of_full_rank_order_RHS[2:length(index_of_full_rank_order_RHS)]


# for each inequality, find the compatible rank order in data (which is the smallest LHS)
index_of_rank_order_LHS <- vector(mode = "list", length(index_of_full_rank_order_RHS))

for (i in 1:length(index_of_rank_order_LHS)){
  # generate a temp variable recording index of rank order
  temp <- rep(0, nrow(unique_rank_order))
  for (j in index_of_full_rank_order_RHS[[i]]){
    temp <- temp + compatibility_matrix[,j]
  }
  index_of_rank_order_LHS[[i]] <- which(temp > 0)
}


# given the index of rank orders, we now calculate the data implied probability of LHS
data_prob_inequality <- rep(0, length(index_of_full_rank_order_RHS))
for (i in 1:length(index_of_full_rank_order_RHS)){
  for (j in index_of_rank_order_LHS[[i]]){
    data_prob_inequality[i] <- data_prob_inequality[i] + data_prob_rank_order[j]
  }
}

# initialize model implied probability (RHS of Artstein)
model_prob_inequality <- matrix(0, nrow(full_rank_order), 1)

beta_grid <- linspace(-5, 5, 101)
beta_identified <- vector(mode = "list", 1)
violation_inequality <- rep(0, length(beta_grid))

# for any candidate of beta, we need to calculate the 
# model implied probability of a full rank with the following function

# function input are observables, beta to be tested and full rank order, 
# currently they are scalars
predict_prob <- function(individual_char, temp_beta, full_rank_order) {
  # calculate base utility for all choices in default order
  base_utility <- individual_char * temp_beta
  # reorder utilities of choices by rank orders
  # initialize ranked base utility
  base_utility_rank <- matrix(0, nrow(individual_char) , length(full_rank_order))
  for (i in 1:length(full_rank_order)) {
    base_utility_rank[, full_rank_order[i]] <- base_utility[, i]
  }
  # calculate probability of the full rank order
  # initialize numerator and denominator
  numerator <- matrix(0, nrow(individual_char) , length(full_rank_order))
  denominator <- matrix(0, nrow(individual_char) , length(full_rank_order))
  
  numerator <- exp(base_utility_rank)
  for (i in 1:length(full_rank_order)) {
    if (i < length(full_rank_order)){
      denominator[, i] <- rowSums(numerator[, i:length(full_rank_order)])
    }    
    if (i == length(full_rank_order)){
      denominator[, i] <- numerator[, i]
    }
  }
  predict_prob <- numerator / denominator
  predict_prob <- apply(predict_prob, 1, prod)
  predict_prob <- sum(predict_prob) / nrow(individual_char)
  return (predict_prob)
}

for (b in 1:length(beta_grid)){
  temp_beta <- beta_grid[b]
  # print(temp_beta)
  # initialize model implied probability
  model_prob_inequality <- rep(0, length(index_of_full_rank_order_RHS))
  for (i in 1:length(model_prob_inequality)){
    for (j in index_of_full_rank_order_RHS[[i]]){
      model_prob_inequality[i] <- model_prob_inequality[i] + predict_prob(parameter_of_choice_set, temp_beta, full_rank_order[j,])
    }
  }
  violation_inequality[b] <- length(data_prob_inequality) - sum(model_prob_inequality <= data_prob_inequality)
  if (sum(model_prob_inequality <= data_prob_inequality) == length(data_prob_inequality)){
    beta_identified <- append(beta_identified, temp_beta)
    #print(temp_beta)
  }
}

beta_identified

violation_inequality

# record results for this simulation
results <- append(results, violation_inequality)

}

results <- results[2:length(results)]
results <- unlist(results)
results <- matrix(results, nrow = number_of_simulation, byrow = TRUE)
results
