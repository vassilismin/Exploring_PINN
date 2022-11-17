setwd("C:/Users/vassi/Documents/GitHub/PBPK_Genetic_Algorithm/Kreyling/NLOPTR")


#####################################
### Function to create Parameters ###
#####################################
create.params <- function(comp_names, w){
  
  # List with names of all possible compartments
  all_comps <- list("RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"="Brain", "Spleen"="Spleen",
                    "Lungs"="Lungs", "Liver"="Liver", "Uterus"="Uterus", "Bone"="Bone", "Adipose"="Adipose", "Skin"="Skin", "Muscles"="Muscles",
                    "GIT"="GIT") # List with names of all possible compartments
  
  ### Density of tissues/organs
  d_tissue <- 1 #g/ml
  d_skeleton <- 1.92 #g/ml
  d_adipose <- 0.940 #g/ml
  
  Q_total <- (1.54*w^0.75)*60 # Total Cardiac Output (ml/h)
  
  Total_Blood <- 0.06*w+0.77 # Total blood volume (ml)
  
  fr_ad <- 0.0199*w + 1.644 # w in g,  Brown et al.1997 p.420. This equation gives the  adipose % of body weight 
  
  #read data from excel
  fractions <- openxlsx::read.xlsx("Rat physiological parameters.xlsx", sheet = 1, colNames = T, rowNames = T)
  fractions <- as.matrix(sapply(fractions, as.numeric))
  rownames(fractions) <- all_comps
  
  #Tissue weight fraction 
  Tissue_fractions <- fractions[,1]/100 # % of BW. Na values refers to the volume of the rest organs(RoB)
  Tissue_fractions[10] <- fr_ad/100
  #Regional blood flow fraction
  Regional_flow_fractions <- fractions[,2]/100 # % of total cardiac output
  #Capillary volume fractions (fractions of tissue volume)
  Capillary_fractions <- fractions[,3] # of tissue volume
  
  W_tis <- rep(0,length(comp_names))
  V_tis <- rep(0,length(comp_names))
  V_cap <- rep(0,length(comp_names))
  Q <- rep(0,length(comp_names))
  
  # The following values were calculated by dividing the %ID/ g tissue with the %ID w/o free 48 from Table 2 of Kreyling et al. (2017)
  # Thus, they represent the average mass, in grams, of the respective tissues in each time group.
  liver_expw <- mean(c(8.57, 8.92, 9.30, 8.61, 9.20))
  spleen_expw <- mean(c(0.93, 0.75, 0.97, 0.68, 0.71))
  kidneys_expw <- mean(c(2.27, 2.36, 2.44, 2.11, 2.26))
  lungs_expw <- mean(c(1.87, 1.60, 1.80, 1.48, 1.31))
  heart_expw <- mean(c(0.89, 1.00, 1.00, 1.00, 0.88))
  blood_expw <- mean(c(16.52, 17.45, 15.33, 18.50, 18.00))
  carcass_expw <- mean(c(206.00, 203.33, 184.00, 202.00, 203.75))
  skeleton_expw <- mean(c(26.15, 27.50, 25.56, 25.79, 25.26))
  soft_tissues <- mean(c(228.57, 253.85, 214.29, 225.93, 231.04))
  
  ### Calculation of tissue weights  
  W_tis[2] <- heart_expw
  W_tis[3] <- kidneys_expw
  W_tis[5] <- spleen_expw
  W_tis[6] <- lungs_expw
  W_tis[7] <- liver_expw
  W_tis[9] <- skeleton_expw
  W_tis[13] <- Tissue_fractions[13]*w
  
  for (i in 1:length(comp_names)) {
    control <- comp_names[i]
    
    Regional_flow_fractions[i] <- ifelse(is.na(control), NA, Regional_flow_fractions[i])
    Capillary_fractions[i] <- ifelse(is.na(control), NA, Capillary_fractions[i])
    
    
    ###Calculation of tissue volumes
    
    if (i==9){
      V_tis[i] <- W_tis[i]/d_skeleton
    } else if(i==10){
      V_tis[i] <- W_tis[i]/d_adipose
    } else{
      V_tis[i] <- W_tis[i]/d_tissue 
    }
    
    ###Calculation of capillary volumes
    V_cap[i] <- V_tis[i]*Capillary_fractions[i]
    
    
    ###Calculation of regional blood flows
    Q[i] <- Q_total*Regional_flow_fractions[i]
  }
  
  
  ### Calculations for "Soft tissue" compartment
  W_tis[1] <- w - sum(W_tis[2:length(W_tis)], na.rm = TRUE)-Total_Blood
  V_tis[1] <- W_tis[1]/d_adipose     
  Q[1] <- Q_total - sum(Q[2:length(Q)],na.rm = TRUE) + Q[6]
  V_cap[1] <- V_tis[1]*Capillary_fractions[1] #Total_Blood - Vven - Vart - sum(V_cap[2:length(V_cap)], na.rm = TRUE)
  
  parameters <- matrix(c(W_tis[],V_tis[],V_cap[],Q[]), ncol = 4)
  colnames(parameters) <- c("W_tis", "V_tis", "V_cap", "Q")
  rownames(parameters) <- all_comps
  
  Vven=0.64*Total_Blood
  Vart=0.15*Total_Blood
  Wm_ven=0.01*Vven
  Wm_art=0.01*Vart
  
  return(c(
    "Q_total"=Q_total, "V_blood"=Total_Blood, "V_ven"=Vven, "V_art"=Vart,
    
    "w_rob"=parameters[1,1], "w_ht"=parameters[2,1], "w_ki"=parameters[3,1], "w_spl"=parameters[5,1], "w_lu"=parameters[6,1], "w_li"=parameters[7,1], "w_bone"=parameters[9,1], "w_git"=parameters[13,1],
    
    "V_tis_rob"=parameters[1,2], "V_tis_ht"=parameters[2,2], "V_tis_ki"=parameters[3,2], "V_tis_spl"=parameters[5,2], "V_tis_lu"=parameters[6,2], "V_tis_li"=parameters[7,2], "V_tis_bone"=parameters[9,2], "V_tis_git"=parameters[13,2], 
    
    "V_cap_rob"=parameters[1,3], "V_cap_ht"=parameters[2,3], "V_cap_ki"=parameters[3,3], "V_cap_spl"=parameters[5,3], "V_cap_lu"=parameters[6,3], "V_cap_li"=parameters[7,3], "V_cap_bone"=parameters[9,3], "V_cap_git"=parameters[13,3],
    
    "Q_rob"=parameters[1,4], "Q_ht"=parameters[2,4], "Q_ki"=parameters[3,4], "Q_spl"=parameters[5,4], "Q_lu"=parameters[6,4], "Q_li"=parameters[7,4], "Q_bone"=parameters[9,4], "Q_git"=parameters[13,4]
    
  ))
}

# Physiological parameters units
# V_blood, V_ven, V_art (ml): Volume of total blood, venous blood and arterial blood
# w_i (g):                    mass of tissue or organ "i"
# V_tis_i (ml):                volume of tissue or organ "i"
# V_cap_i (ml):                volume of capillary blood in tissue "i"
# Q_i, Q_total (ml/h):        regional blood flow of tissue or organ "i"
#===============================================
#2. Function to create initial values for ODEs 
#===============================================

create.inits <- function(parameters, dose){
  with(as.list(parameters),{
    M_ht<-0; M_lu<-0; M_li<-0; M_spl<-0; M_ki<-0; M_git<-0; M_bone<-0; M_rob<-0;
    
    M_cap_ht<-0; M_cap_lu<-0; M_cap_li<-0; M_cap_spl<-0; M_cap_ki<-0; M_cap_git<-0; M_cap_bone<-0; M_cap_rob<-0;
    
    M_lumen <- 0;
    M_ven <- dose; M_art<-0
    M_feces<-0; M_urine<-0 
    
    return(c("M_ht" = M_ht, "M_lu" = M_lu, 
             "M_li" = M_li, "M_spl" = M_spl, 
             "M_ki" = M_ki, "M_git" = M_git, 
             "M_bone" = M_bone,"M_rob"=M_rob,
             
             "M_cap_ht" = M_cap_ht, "M_cap_lu" = M_cap_lu, 
             "M_cap_li" = M_cap_li, "M_cap_spl" = M_cap_spl, 
             "M_cap_ki" = M_cap_ki, "M_cap_git" = M_cap_git, 
             "M_cap_bone" = M_cap_bone,"M_cap_rob"=M_cap_rob,
             
             "M_lumen" = M_lumen,
             "M_ven" = M_ven, "M_art" = M_art, "M_feces" = M_feces, "M_urine" = M_urine))
    
  })
}

#==============
#3. ODEs System
#==============
ode.func <- function(time, inits, params){
  # position <- params[37:(37+15)]
  # fit_pars <- params[(37+16):length(params)]
  
  with(as.list(c(inits, params)),{
    
    P_ht <- P1
    P_lu <- P2
    P_li <- P3
    P_spl <- P4
    P_ki <- P5
    P_git <- P6
    P_bone <- P7
    P_rob <- P8
    
    x_ht <- X1
    x_lu <- X2
    x_li <- X3
    x_spl <- X4
    x_ki <- X5
    x_git <- X6
    x_bone <-X7
    x_rob <- X8
    
    CLE_f <- CLE_f        
    CLE_h <- CLE_h 
    CLE_u <- 0
    
    
    # Concentrations (mg of NPs)/(g of wet tissue)
    C_ht <- M_ht/w_ht
    C_cap_ht <- M_cap_ht/V_cap_ht
    C_lu <- M_lu/w_lu
    C_cap_lu <- M_cap_lu/V_cap_lu
    C_li <- M_li/w_li
    C_cap_li <- M_cap_li/V_cap_li
    C_spl <- M_spl/w_spl
    C_cap_spl <- M_cap_spl/V_cap_spl
    C_ki <- M_ki/w_ki
    C_cap_ki <- M_cap_ki/V_cap_ki
    C_git <- M_git/w_git
    C_cap_git <- M_cap_git/V_cap_git
    C_bone <- M_bone/w_bone
    C_cap_bone <- M_cap_bone/V_cap_bone
    C_rob <- M_rob/w_rob
    C_cap_rob <- M_cap_rob/V_cap_rob
    
    C_ven <- M_ven/V_ven
    C_art <- M_art/V_art
    
    # Heart
    dM_cap_ht <- Q_ht*(C_art - C_cap_ht) - x_ht*Q_ht*(C_cap_ht - C_ht/P_ht)
    dM_ht <- x_ht*Q_ht*(C_cap_ht - C_ht/P_ht) 
    
    # Lungs
    dM_cap_lu <- Q_total*(C_ven - C_cap_lu) - x_lu*Q_total*(C_cap_lu - C_lu/P_lu)
    dM_lu <-  x_lu*Q_total*(C_cap_lu - C_lu/P_lu)
    
    # Liver 
    dM_cap_li <- Q_li*(C_art - C_cap_li) + Q_spl*(C_cap_spl - C_cap_li) + Q_git*(C_cap_git - C_cap_li) -
      x_li*(Q_li)*(C_cap_li - C_li/P_li)
    dM_li <- x_li*Q_li*(C_cap_li - C_li/P_li) - CLE_h*M_li
    
    # Spleen
    dM_cap_spl <- Q_spl*(C_art - C_cap_spl) - x_spl*Q_spl*(C_cap_spl - C_spl/P_spl)
    dM_spl <- x_spl*Q_spl*(C_cap_spl - C_spl/P_spl) 
    
    # Kidneys
    dM_cap_ki <- Q_ki*(C_art - C_cap_ki) - x_ki*Q_ki*(C_cap_ki - C_ki/P_ki)- CLE_u*M_cap_ki
    dM_ki <- x_ki*Q_ki*(C_cap_ki - C_ki/P_ki) 
    
    # GIT - Gastrointestinal Track
    dM_cap_git <- Q_git*(C_art - C_cap_git) - x_git*Q_git*(C_cap_git - C_git/P_git)
    dM_git <- x_git*Q_git*(C_cap_git - C_git/P_git) 
    dM_lumen <- CLE_h*M_li - CLE_f *M_lumen 
    
    # Bone
    dM_cap_bone <- Q_bone*(C_art - C_cap_bone) - x_bone*Q_bone*(C_cap_bone - C_bone/P_bone)
    dM_bone <- x_bone*Q_bone*(C_cap_bone - C_bone/P_bone) 
    
    
    # RoB - Rest of Body
    dM_cap_rob <- Q_rob*(C_art - C_cap_rob) - x_rob*Q_rob*(C_cap_rob - C_rob/P_rob)
    dM_rob <- x_rob*Q_rob*(C_cap_rob - C_rob/P_rob) 
    
    # Urine
    dM_urine <- CLE_u*M_cap_ki
    
    # Feces
    dM_feces <- CLE_f*M_lumen
    
    # Venous Blood
    dM_ven <- Q_ht*C_cap_ht + (Q_li + Q_spl+Q_git)*C_cap_li + Q_ki*C_cap_ki +
      Q_bone*C_cap_bone + Q_rob*C_cap_rob - Q_total*C_ven
    
    # Arterial Blood
    dM_art <-  Q_total*C_cap_lu - Q_total*C_art
    
    Blood_total <- M_ven + M_art + M_cap_ht + M_cap_lu +M_cap_li+M_cap_spl+
      M_cap_ki+ M_cap_git+M_cap_bone+M_cap_rob
    Blood <- Blood_total/(V_blood)
    
    C_soft <- (M_git+M_lumen+M_rob)/(w_git + w_rob)
    
    list(c("dM_ht" = dM_ht, "dM_lu" = dM_lu, 
           "dM_li" = dM_li, "dM_spl" = dM_spl, 
           "dM_ki" = dM_ki, "dM_git" = dM_git, 
           "dM_bone" = dM_bone,"dM_rob"=dM_rob,
           
           "dM_cap_ht" = dM_cap_ht, "dM_cap_lu" = dM_cap_lu, 
           "dM_cap_li" = dM_cap_li, "dM_cap_spl" = dM_cap_spl, 
           "dM_cap_ki" = dM_cap_ki, "dM_cap_git" = dM_cap_git, 
           "dM_cap_bone" = dM_cap_bone,"dM_cap_rob"=dM_cap_rob,
           
           "dM_lumen" = dM_lumen,
           "dM_ven" = dM_ven, "dM_art" = dM_art, "dM_feces" = dM_feces, "dM_urine" = dM_urine),
         
         "Blood"=Blood,
         "C_ht"=C_ht, "C_lu"=C_lu, "C_li"=C_li, "C_spl"=C_spl,
         "C_ki"=C_ki,  "C_bone"=C_bone, "C_soft"=C_soft,
         "Feces"=M_feces)
  })
}




decode_ga <- function(real_num){ 
  # Partition coefficient grouping
  P1 <- floor(real_num[1])
  P2 <- floor(real_num[2])
  P3 <- floor(real_num[3])
  P4 <- floor(real_num[4])
  P5 <- floor(real_num[5])
  P6 <- floor(real_num[6])
  P7 <- floor(real_num[7])
  P8 <- floor(real_num[8])
  
  
  # Permeability coefficient grouping
  X1 <- floor(real_num[9])
  X2 <- floor(real_num[10])
  X3 <- floor(real_num[11])
  X4 <- floor(real_num[12])
  X5 <- floor(real_num[13])
  X6 <- floor(real_num[14])
  X7 <- floor(real_num[15])
  X8 <- floor(real_num[16])
  
  out <- structure(c(P1,P2,P3,P4,P5,P6,P7,P8, X1,X2,X3,X4,X5,
                     X6,X7,X8), names = c("P1","P2","P3","P4",
                                          "P5","P6", "P7", "P8", "X1",
                                          "X2", "X3", "X4", "X5", "X6", "X7", "X8"))
  return(out)
}


decode_ga_sppcg <- function(real_num){ 
  # Partition coefficient grouping
  P1 <- floor(real_num[1])
  P2 <- floor(real_num[2])
  P3 <- floor(real_num[3])
  P4 <- floor(real_num[4])
  P5 <- floor(real_num[5])
  P6 <- floor(real_num[6])
  P7 <- floor(real_num[7])
  P8 <- floor(real_num[8])
  
  
  out <- structure(c(P1,P2,P3,P4,P5,P6,P7,P8), names = c("P1","P2","P3","P4",
                                                         "P5","P6", "P7", "P8"))
  return(out)
}

#=============================
#8. Create position  
#=============================  
# Function for creating the position from which to draw each param from the fitted params vector
create.position <- function(grouping){
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  N_p <-8 #   Number of partition coefficients
  N_x <- 8#   Number of permeability coefficients  
  # Define size of P and X groups
  P_groups <- length(unique(grouping[1:N_p]))  # sample size
  X_groups <- length(unique(grouping[(N_p+1):(N_p+N_x)]))  # sample size
  # set.seed(0)
  # Initilise parameter values
  fitted <- rep(NA, P_groups+X_groups+2)
  # Initialise naming vectors
  pnames <- rep(NA, P_groups)
  xnames <- rep(NA, X_groups)
  
  #Define names for P and X groups
  for (i in 1:P_groups){
    pnames[i] <- paste0("P", as.character(unique(grouping[1:N_p])[i]))
  }
  for (j in 1:X_groups){
    xnames[j] <- paste0("X", as.character(unique(grouping[(N_p+1):(N_p+N_x)])[j]))
  }
  # Define the total parameter vector names
  names(fitted) <- c(pnames, xnames,"CLE_f",  "CLE_h")
  # Variable for keeping which value in the fitted params vector corresponds to each coefficient
  position = rep(NA, length(grouping))
  for (i in 1:(length(position))){
    if(i<=8){
      position[i] <- which(names(fitted) == paste0("P", as.character(grouping[i])))
    }else{
      position[i] <- which(names(fitted) == paste0("X", as.character(grouping[i])))
    }
  }
  fitted[] <- c(log(exp(runif(P_groups, 2,3))),log(exp(runif(X_groups, -4,-3))), log(exp(runif(1, -2,-1))), log(exp(runif(1, -11,-9))))
  
  return(list("position"=position,"fitted"=fitted, 'P_groups' = P_groups, X_groups = X_groups))
}


create.position_constrained <- function(grouping){
  #---------------------------
  # Define fitting parameters 
  #---------------------------
  N_group <- 8 #   Number of groups fro partition and permeability coefficients
  # Define size of P and X groups
  P_groups <- length(unique(grouping))  # sample size
  X_groups <- P_groups
  # set.seed(0)
  # Initilise parameter values
  fitted <- rep(NA, P_groups+X_groups+2)
  # Initialise naming vectors
  pnames <- rep(NA, P_groups)
  xnames <- rep(NA, X_groups)
  
  #Define names for P and X groups
  for (i in 1:P_groups){
    pnames[i] <- paste0("P", as.character(unique(grouping[1:N_group])[i]))
  }
  for (j in 1:X_groups){
    xnames[j] <- paste0("X", as.character(unique(grouping[1:N_group])[j]))
  }
  # Define the total parameter vector names
  names(fitted) <- c(pnames, xnames,"CLE_f",  "CLE_h")
  # Variable for keeping which value in the fitted params vector corresponds to each coefficient
  position = rep(NA, 2*length(grouping))
  for (i in 1:(length(position)/2)){
    position[i] <- which(names(fitted) == paste0("P", as.character(grouping[i])))
    position[i+8] <- position[i] + P_groups
  }
  fitted[] <- c(log(exp(runif(P_groups, 2,3))),log(exp(runif(X_groups, -4,-3))), log(exp(runif(1, -2,-1))), log(exp(runif(1, -11,-9))))
  
  return(list("position"=position,"fitted"=fitted, 'P_groups' = P_groups, X_groups = X_groups))
}

# Create the parameter grouping for the max and ga problems
grouping_MANG <- c(1:8, 1:8)

# Create the position vector to match the ODE parameters with the fitted parameter values
get_position_MANG <- create.position(grouping_MANG)
position_MANG <- get_position_MANG$position
fitted_MANG <-  get_position_MANG$fitted

x0= fitted_MANG

sample_time <- seq(0, 28*24, 1)
sample_time <- seq(0, 9, 1)
# Initialise vector of physiological parameters
dose <- 18.15 # ug # Since results are in %ID we used a random dose that is within the dose range given to the rats
mass <- 263 #g, female Wistar Kyoto rats
compartments <- list( "RoB"="RoB","Heart"="Heart", "Kidneys"="Kidneys", "Brain"= NA, "Spleen"="Spleen",
                      "Lungs"="Lungs", "Liver"="Liver", "Uterus"= NA, "Bone"="Bone", "Adipose"=NA, "Skin"=NA, "Muscles"=NA, "GIT"="GIT") #used as input in function, compartments that are used in pbpk

phys_pars <- create.params(compartments,mass)
inits <- create.inits(phys_pars, dose)

params <- c(phys_pars, exp(x0))
solution <- data.frame(deSolve::ode(times = sample_time,  func = ode.func,
                                    y = inits, parms = params, 
                                    method="lsodes",rtol = 1e-3, atol = 1e-3))