
############  Microsimulation Sick-Sicker model ##########
# Includes: 
# individual characteristicss: age
# age dependent mortality probabilities 
# individual treatment effect modifyer 
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
# Please cite our publications when using this code:
# darthworkgroup.com
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400-422.
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial
################################################################################
# Copyright 2017, 
# THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property
# are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the collaborating 
# institutions and may not be used, reproduced, modified, distributed or adapted 
# in any way without appropriate citation.
################################################################################

#### 01 Load packages ####
if (!require('dampack')) install_github('DARTH-git/dampack'); library(dampack)
if (!require('dplyr')) install.packages('dplyr'); library(dplyr)

#### 02 Load Functions ####
source("Functions.R")


#### 03 Input Model Parameters ####
set.seed(1)                    # set the seed  

# Model structure 
n_t   <- 30                    # time horizon, 30 cycles
n_i   <- 100000                # number of simulated individuals
v_n   <- c("H","S1","S2","D")  # the model states names
n_s   <- length(v_n)           # the number of health states
d_r   <- 0.03                  # discount rate of 3% per cycle
v_dwe <- v_dwc <- 1 / ((1 + d_r) ^ (0:n_t))   # discount weight 
v_names_str <- c("no treatment", "treatment") # strategy names

# Event probabilities (per cycle)
## Annual transition probabilities
p_HS1   <- 0.15                # probability of becoming sick when healthy
p_S1H   <- 0.5                 # probability of recovering to healthy when sick
p_S1S2  <- 0.105               # probability of becoming sicker when sick

## Annual probabilities of death
p_mort   <- read.csv("mortProb_age.csv")        # load age dependent probability
dist_Age <- read.csv("MyPopulation-AgeDistribution.csv") # load age distribution

p_S1D    <- 0.0149          # probability to die in S1 by cycle 
p_S2D    <- 0.048           # probability to die in S2

# Cost inputs
c_H     <- 2000             # cost of one cycle in the healthy state
c_S1    <- 4000             # cost of one cycle in the sick state
c_S2    <- 15000            # cost of one cycle in the sicker state
c_D     <- 0                # cost of one cycle in the dead state
c_Trt   <- 12000            # cost of treatment (per cycle)

# Utility inputs
u_H     <- 1                # utility when healthy 
u_S1    <- 0.75             # utility when sick 
u_S2    <- 0.5              # utility when sicker
u_D     <- 0                # utility when dead
u_Trt   <- 0.95             # utility when sick(er) and being treated


#### 04 Sample individual level characteristics ####
#### 04.1 Static characteristics ####


#### 04.2 Dynamic characteristics 



#### 05 Define Simulation Functions ####
#### 05.1 Probability function ####


#### 05.2 Cost function ####
# The Costs function estimates the costs at every cycle.


#### 05.3 Health outcome function ####
# The Effs function to update the utilities at every cycle.


#### 06 Run Microsimulation ####


#### 07 Visualize results ####


#### 08 Cost Effectiveness Analysis ####
