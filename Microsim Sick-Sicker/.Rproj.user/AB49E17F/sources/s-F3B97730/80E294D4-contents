
###############  Simple 3-state microsimulation model              #############
# Includes: 
# individual characteristics: sex
################################################################################

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (2)	
# M.G. Myriam Hunink, MD, PhD (3,4)
# Hawre J. Jalal, MD, PhD (5) 
# Eline M. Krijkamp, MSc (3)
# Petros Pechlivanoglou, PhD (6) 

# In collaboration of: 		
# 1 Center for Research and Teaching in Economics (CIDE), Drug Policy Program, Mexico
# 2 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 3 Erasmus MC, Rotterdam, The Netherlands
# 4 Harvard T.H. Chan School of Public Health, Boston, USA
# 5 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 6 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada


################################################################################
# Please cite our publications when using this code. Publications can be found at # http://darthworkgroup.com/publications/

################################################################################
# Copyright 2017, 
# THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual 
# property are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, 
# distributed or adapted in any way without appropriate citation.

################################################################################



#### 01 Load packages ####
if (!require('devtools')) install.packages('devtools'); library(devtools)
if (!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if (!require('dampack')) install_github('DARTH-git/dampack'); library(dampack)
if (!require('dplyr')) install.packages('dplyr'); library(dplyr)


#### 02 Load Functions ####
source("Functions.R")


#### 03 Input Model Parameters ####

# Model structure
v_n   <- c("healthy", "sick", "dead")          # vector with state names
n_s   <- length(v_n)                           # number of states
n_t   <- 60                                    # number of cycles
n_i   <- 10000                                 # number of individuals
d_r   <- 0.03                                  # discount rate of 3% per cycle
v_dwe <- v_dwc <- 1 / ((1 + d_r) ^ (0:n_t))    # discount weight 

# Transition probabilities
p_HS  <- 0.05          # probability healthy -> sick
p_HD_female <- 0.0382  # probability health -> dead when female
p_HD_male   <- 0.0463  # probability health -> dead when male
m_p_HD      <- data.frame(Sex = c("Female", "Male"), p_HD = c(p_HD_female, p_HD_male))

p_SD        <- 0.1     # probability sick -> dead

# Costs inputs
c_H   <- 1500      # cost of one cycle in healthy state
c_S   <- 5000      # cost of one cycle in sick state
c_D   <- 0         # cost of one cycle in dead state

# utility inputs
u_H   <- 1         # utility when healthy 
u_S   <- 0.85      # utility when sick 
u_D   <- 0         # utility when dead 


#### 04 Sample individual level characteristics ####
#### 04.1 Static characteristics ####

set.seed(1)  # set the seed  

v_sex    <- sample(x = c("Female", "Male"), prob = c(0.5, 0.5), size = n_i, replace = TRUE) # randomly sample the sex of an individual (50% female)
df_X     <- data.frame(ID = 1:n_i, Sex = v_sex)

#### 04.2 Dynamic characteristics 
# Specify the initial health state of the individuals 
# everyone begins in the healthy state (in this example)
v_M_init <- rep("healthy", times = n_i)   


#### 05 Define Simulation Functions ####

#### 05.1 Probability function ####
# The Probs function that updates the transition probabilities of every cycle is shown below.

Probs <- function(M_t, df_X) { 
  # Arguments:
  # M_t:  health state occupied at cycle t (character variable)
  # df_X: data frame with individual characteristics data  
    ## Sex:      sex of the individuals 
  # Returns: 
  #   transition probabilities for that cycle
  
  m_p_t           <- matrix(0, nrow = n_s, ncol = n_i)  # create matrix of state transition probabilities
  rownames(m_p_t) <-  v_n                               # give the state names to the rows
  
  # lookup baseline probability and rate of dying based on individual characteristics
  p_HD_all <- inner_join(df_X, m_p_HD, by = c("Sex"))
  p_HD     <- p_HD_all[M_t == "healthy","p_HD"]
  
  # update m_p_t with the appropriate probabilities   
  m_p_t[, M_t == "healthy"] <- rbind(1 - p_HD - p_HS, p_HS, p_HD)  # transition probabilities when healthy 
  m_p_t[, M_t == "sick"]    <- rbind(0, 1 - p_SD, p_SD)            # transition probabilities when sick 
  m_p_t[, M_t == "dead"]    <- rbind(0, 0, 1)                      # transition probabilities when dead     
  return(t(m_p_t))
}       

#### 05.2 Cost function ####
# The Costs function estimates the costs at every cycle.
Costs <- function (M_t) {
  # M_t: current health state
  c_t <- c()
  c_t[M_t == "dead"]    <- c_D     # costs at dead state
  c_t[M_t == "healthy"] <- c_H     # costs accrued by being healthy this cycle
  c_t[M_t == "sick"]    <- c_S     # costs accrued by being sick this cycle
  
  return(c_t)  # return costs accrued this cycle
}

#### 05.3 Health outcome function ####
# The Effs function to update the utilities at every cycle.

Effs <- function (M_t) {
  # M_t: current health state
  q_t <- c() 
  q_t[M_t == "healthy"] <- u_H     # QALYs accrued by being healthy this cycle
  q_t[M_t == "dead"]    <- u_D     # QALYs at dead state
  q_t[M_t == "sick"]    <- u_S     # QALYs accrued by being sick this cycle
  
  return(q_t)  # return the QALYs accrued this cycle
}


#### 06 Run Microsimulation ####
MicroSim <- function(n_i, df_X, seed = 1) {
# Arguments:  
  # n_i:     number of individuals
  # df_X:    data frame with individual data 
  # seed:    defauls is 1
# Returns
  # a list with information about the individuals transitions, associated costs and effects and total costs and rewards 
  
  set.seed(seed) # set the seed
  
  # create three matrices called m_M, m_C and m_E
  # number of rows is equal to the n_i, the number of columns is equal to n_t (the initial state and all the n_t cycles)
  # m_M is used to store the health state information over time for every individual
  # m_C is used to store the costs information over time for every individual
  # m_E is used to store the effects information over time for every individual
  
  m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                       dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                       paste("cycle", 0:n_t, sep = " ")))  
 
  m_M[, 1] <- v_M_init          # initial health state 
  m_C[, 1] <- Costs(m_M[, 1])   # costs accrued during cycle 0
  m_E[, 1] <- Effs(m_M[, 1])    # QALYs accrued during cycle 0
  
  # open a loop for time running cycles 1 to n_t 
  for (t in 1:n_t) {
    m_P <- Probs(m_M[, t], df_X)        # calculate the transition probabilities for the cycle based on health state t
    m_M[, t + 1]  <- samplev(m_P, 1)       # sample the current health state and store that state in matrix m_M 
    m_C[, t + 1]  <- Costs(m_M[, t + 1])   # calculate costs 
    m_E[, t + 1]  <- Effs (m_M[, t + 1])   # calculate QALYs 
    
    # Display simulation progress
    if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_t * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # calculate  
  tc <- m_C %*% v_dwc    # total (discounted) cost per individual
  te <- m_E %*% v_dwe    # total (discounted) QALYs per individual 
  tc_hat <- mean(tc)     # average (discounted) cost 
  te_hat <- mean(te)     # average (discounted) QALYs
  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, tc_hat = tc_hat, te_hat = te_hat)   
  return(results)  # return the results
} # end of the MicroSim function  

# By specifying all the arguments in the `MicroSim()` the simulation can be started
# In this example the outcomes are of the simulation are stored in the variables `outcomes_no_tr` and `outcomes_trt`.

# Run the simulation for both no treatment and treatment options
outcomes <- MicroSim(n_i = n_i, df_X = df_X, seed = 1)


# Show results
results  <- data.frame("Total Cost" = outcomes$tc_hat, "Total QALYs" = outcomes$te_hat)
results


#### 07 Visualize results ####
options(scipen = 999)

plot(density(outcomes$tc), main = paste("Total cost per person"), xlab = "Cost ($)")
plot(density(outcomes$te), main = paste("Total QALYs per person"), xlab = "QALYs")
plot_m_TR(outcomes$m_M)    # health state trace


