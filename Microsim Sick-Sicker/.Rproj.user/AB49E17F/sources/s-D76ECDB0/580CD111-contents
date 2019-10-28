#---------------------------------------------------------------------------#
#### R functions calculate the CE outcomes                               ####
#---------------------------------------------------------------------------#

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




# used for sensitivity analysis

calculate_ce_out <- function (l_params_all, n_wtp = 100000) {
  with(as.list(l_params_all), {
    
    ### update sample individual level characteristics 
    # Static characteristics
    v_sex    <- sample(x = c("Female", "Male"), prob = c(0.5, 0.5), size = n_i, replace = TRUE) # randomly sample the sex of an individual (50% female)
    df_X     <- data.frame(ID = 1:n_i, Sex = v_sex)
    
    #### 05.1 Probability function ####
    # The Probs function that updates the transition probabilities of every cycle is shown below.
    
    Probs <- function(M_it, df_X, v_Ts, t) { 
      # Arguments:
      # M_it: health state occupied by individual i at cycle t (character variable)
      # t:     current cycle 
      # Returns: 
      #   transition probabilities for that cycle
      
      m_p_it           <- matrix(0, nrow = n_s, ncol = n_i)  # create matrix of state transition probabilities
      rownames(m_p_it) <-  v_n                               # give the state names to the rows
      
      # lookup baseline probability and rate of dying based on individual characteristics
      p_HD_all <- inner_join(df_X, m_p_HD, by = c("Sex"))
      
      p_HD     <- p_HD_all[M_it == "healthy","p_HD"]
      # p_HD <- case_when(df_X$Sex == "Male"   & M_it == "healthy" ~ p_HD_male,
      #                   df_X$Sex == "Female" & M_it == "healthy" ~ p_HD_female)
      # p_HD <- na.omit(p_HD)
      # update m_p_it with the appropriate probabilities   
      m_p_it[, M_it == "healthy"] <- rbind(1 - p_HD - p_HS, p_HS, p_HD)    # transition probabilities when healthy 
      m_p_it[, M_it == "sick"]    <- rbind(0, 1 - p_SD[v_Ts], p_SD[v_Ts])  # transition probabilities when sick 
      m_p_it[, M_it == "dead"]    <- c(0, 0, 1)                            # transition probabilities when dead     
      return(t(m_p_it))
    }       
    
    #### 05.2 Cost function ####
    # The Costs function estimates the costs at every cycle.
    Costs <- function (M_it) {
      # M_it: current health state
      c_it <- c()
      c_it[M_it == "dead"]    <- c_D     # costs at dead state
      c_it[M_it == "healthy"] <- c_H     # costs accrued by being healthy this cycle
      c_it[M_it == "sick"]    <- c_S     # costs accrued by being sick this cycle
      
      return(c_it)  # return costs accrued this cycle
    }
    
    #### 05.3 Health outcome function ####
    # The Effs function to update the utilities at every cycle.
    
    Effs <- function (M_it) {
      # M_it: current health state
      q_it <- c() 
      q_it[M_it == "dead"]    <- u_D     # QALYs at dead state
      q_it[M_it == "healthy"] <- u_H     # QALYs accrued by being healthy this cycle
      q_it[M_it == "sick"]    <- u_S     # QALYs accrued by being sick this cycle
      
      return(q_it)  # return the QALYs accrued this cycle
    }
    
    
    #### 06 Run Microsimulation ####
    MicroSim <- function(n_i, df_X, seed = 1) {
      # Arguments:  
      # n_i:     number of individuals
      # df_X     data frame with individual data 
      # Sex      sex of the individuals 
      # seed:    default is 1
      
      set.seed(seed) # set the seed
      
      n_s <- length(v_n) # the number of health states
      
      # create three matrices called m_M, m_C and m_E
      # number of rows is equal to the n_i, the number of columns is equal to n_t (the initial state and all the n_t cycles)
      # m_M is used to store the health state information over time for every individual
      # m_C is used to store the costs information over time for every individual
      # m_E is used to store the effects information over time for every individual
      
      m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = n_t + 1, 
                                   dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                                   paste("cycle", 0:n_t, sep = " ")))  
      
      m_M[, 1] <- v_M_init          # initial health state for individual i
      v_Ts     <- v_Ts_init         # initialize time since illnes onset for individual i
      
      m_C[, 1] <- Costs(m_M[, 1])   # costs accrued individual i during cycle 0
      m_E[, 1] <- Effs(m_M[, 1])    # QALYs accrued individual i during cycle 0
      
      # open a loop for time running cycles 1 to n_t 
      for (t in 1:n_t) {
        v_p <- Probs(m_M[, t], df_X, v_Ts, t)  # calculate the transition probabilities for the cycle based on health state t
        m_M[, t + 1]  <- samplev(v_p, 1)       # sample the current health state and store that state in matrix m_M 
        m_C[, t + 1]  <- Costs(m_M[, t + 1])   # calculate costs per individual during cycle t + 1
        m_E[, t + 1]  <- Effs (m_M[, t + 1])   # calculate QALYs per individual during cycle t + 1
        
        v_Ts <- if_else(m_M[, t + 1] == "sick", v_Ts + 1, 0) # update time since illness onset for t + 1 
        
        # # Display simulation progress
        # if(t/(n_t/10) == round(t/(n_t/10), 0)) { # display progress every 10%
        #   cat('\r', paste(t/n_t * 100, "% done", sep = " "))
        # }
        
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

    # Run the simulation for both no treatment and treatment options
    outcomes <- MicroSim(n_i, df_X, seed = 1)
    
    ## Dataframe with discounted cost and effectiveness 
    df_ce <- data.frame(Cost     = outcomes$tc_hat,
                        Effect   = outcomes$te_hat)
    return(df_ce)
  }
  )
}
