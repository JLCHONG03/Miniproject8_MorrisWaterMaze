run_trial <- function(weights0, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun)
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  
  pool_diameter <- 1.4 #Maze diameter in metres (m)
  platform_radius <- 0.06 #Platform radius
  
  N_pc <- 211 #Population of place cells
  N_ac <- 36 #Population of action cells
  which <- 0
  
  dist <- 0
  wall_zone <- 0
  quadrants <- c(0,0,0,0) #Percentage spent on each quadrant
  
  weights <- weights0 #Initialize modifiable weights
  
  el_tr <- matrix(rep(0, N_pc*N_ac), nrow = N_pc) #Initialize eligibility traces matrix
  
  #Initialize trajectories
  track_x <- starting_x #Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0
  
  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2)
  {
    weights <- weights*(1-noise) + matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult*noise
    
    #Calculate PC activation
    PC_activation <- rep(0, N_pc)
    for (i in 1:N_pc){
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 + (track_y[length(track_y)] - PC_y[i])^2)/(2*sigma_pc^2))
    }
    
    #Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1){
      prevQ <- AC_activation[which] #Displays the Q value before movement
    }
    
    AC_activation <- PC_activation %*% weights
    
    #Make an action
    ACsel <- pmax(0, as.numeric(AC_activation))^beta
    sel_sum <- sum(ACsel)
    
    if (!is.finite(sel_sum) || sel_sum <= 0) {
      ACsel <- rep(1 / N_ac, N_ac)
    } else {
      ACsel <- ACsel / sel_sum
    }
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    while (which < N_ac && ASsum < ASrand){
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }
    
    #Eligibility traces
    el_tr <- el_tr * etdecay
    
    for (j in 1:N_ac){
      itmp <- min(abs(j-which), N_ac-abs(j-which))
      actgaus <- exp(-(itmp*itmp)/(2*sigma_ac*sigma_ac))
      el_tr[,j] <- el_tr[,j] + actgaus*AC_activation[j]*t(t(PC_activation))
    }
    
    vel_x = c(vel_x, (vel_x[length(vel_x)]+ac_const*cos(which/N_ac*2*pi))*Vdecay)
    vel_y = c(vel_y, (vel_y[length(vel_y)]+ac_const*sin(which/N_ac*2*pi))*Vdecay)
    #velocity per time step (not second)
    track_x = c(track_x, track_x[length(track_x)]+vel_x[length(vel_x)])
    track_y = c(track_y, track_y[length(track_y)]+vel_y[length(vel_y)])
    
    #Check if not out of bounds, reset location & speed if so
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > (pool_diameter/2)^2)
    {
      ratio = (track_x[length(track_x)]^2 + track_y[length(track_y)]^2)/((pool_diameter/2)^2)
      track_x[length(track_x)] = track_x[length(track_x)]/sqrt(ratio)
      track_y[length(track_y)] = track_y[length(track_y)]/sqrt(ratio)
      vel_x[length(vel_x)] = track_x[length(track_x)] - track_x[length(track_x)-1]
      vel_y[length(vel_y)] = track_y[length(track_y)] - track_y[length(track_y)-1]
    }
    
    
    if (length(track_x) > 2)
    { if ((track_x[length(track_x)]  - platform_x)^2 + (track_y[length(track_y)]  - platform_y)^2 < platform_radius^2)
    { rew = 10 } #found platform - reward
      else if (track_x[length(track_x)]^2+track_y[length(track_y)]^2 > (0.99*pool_diameter/2)^2)
      { rew = -wall_pun } #hit wall - punishment
      else
      { rew = 0 } #didn't find - no reward
      
      currQ = AC_activation[which]
      tderr = rew + discf*currQ - prevQ #temporal difference error
      weights = pmax(weights + lrate*tderr*el_tr, 0)
    }
    
    laststep = sqrt((track_x[length(track_x)]-track_x[length(track_x)-1])^2 + (track_y[length(track_y)]-track_y[length(track_y)-1])^2)
    dist = dist + laststep
    
    if (track_x[length(track_x)]^2 + track_y[length(track_y)]^2 > 0.8*(pool_diameter/2)^2)
    { wall_zone = wall_zone + 1 }
    else if (track_x[length(track_x)] > 0 && track_y[length(track_y)] > 0)
    { quadrants[1] = quadrants[1] + 1 }
    else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] > 0)
    { quadrants[2] = quadrants[2] + 1 }
    else if (track_x[length(track_x)] < 0 && track_y[length(track_y)] < 0)
    { quadrants[3] = quadrants[3] + 1 }
    else
    { quadrants[4] = quadrants[4] + 1 }
    
    if (length(track_x) > 100) # evaluate latency only after 100+ steps to be accurate
    { speed_ts = mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
    latency = (length(track_x)-1) * speed_ts / speed # convert to seconds
    if (latency > 60) # if more than a minute, stop
    { break }
    }
    
  }
  
  latency <- length(track_x)-1 # latency in time steps
  wall_zone <- wall_zone/latency
  quadrants <- quadrants/latency
  speed_ts <- mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
  
  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency))
}

N_pc <- 211 #Population of place cells [100..300]
N_ac <- 36 #Population of action cells [25..50]

plot_trajectories <- 0 #yes - 1, no - 0
plot_cognitive_maps <- 0 #yes - 1, no - 0
pln <- plot_trajectories + plot_cognitive_maps
Nruns <- 50 #how many runs to run if not plotting anything
## VARIABLE PLATFORM
variable_platform <- 1 

pool_diameter <- 1.4 #Maze diameter (\phi) in metres (m)
platform_radius <- 0.06 #Platform radius (m)
sigma_pc <- 0.1 #place cell sigma (standard deviation), in meters [0.05..0.2]
sigma_ac <- 2 #action cell sigma (standard deviation), in action cells [1..3]

etdecay <- 0.83 #Eligibility trace decay (lambda) [0.75..0.95] LESS THAN GAMMA!
beta <- 6 #Exploration-exploitation factor (\beta) [0.5..12]
alpha <- 0.01 #Learning rate (\alpha) [0.005..0.02]
gamma <- 0.85 #Discount factor (\gamma) [0.75..0.95]

Vdecay <- 0.82 #velocity decay [0.75..0.95]
ac_const <- 0.02 #acceleration const [0.01..0.03]
Wnoise <- 0.0004 #Weight noise [0.0001..0.0007]
Wmult <- 0.1 #Weight multiplier [0.05..0.15]
hitwall <- 0.5 #punishment for hitting the wall [0..1]
speed <- 0.175 #mouse speed (m/s) [0.1..0.25]

Ntrials <- 4 #number of trials per day
Ndays <- 8 #number of days

#performance measures to compute: latency, distance, time in target quadrant, opposite quadrant, and wall zone
if (pln > 0.5) #if any plots
{ PMs <- array(rep(0,5*Ndays*Ntrials), c(5,Ndays,Ntrials))
} else {
  PMs <- array(rep(0,5*Ndays*Ntrials*Nruns), c(5,Ndays,Ntrials,Nruns)) }
#multiple runs


#Platform coordinates:
base_platform_x <- cos(-pi/4) * pool_diameter / 4
base_platform_y <- sin(-pi/4) * pool_diameter / 4

platform_x <- base_platform_x
platform_y <- base_platform_y

#Starting locations of the modeled animal (4 different ones):
strad <- pool_diameter/2*0.85 #15% of maze radius to the wall
starting_xs <- strad * c(cos(pi/6), cos(pi/3), cos(7*pi/6), cos(4*pi/3)) #x coordinates
starting_ys <- strad * c(sin(pi/6), sin(pi/3), sin(7*pi/6), sin(4*pi/3)) #y coordinates

th <- (0:100)/50*pi #for plotting circles :)

if (pln > 0.5) {
  
  #Generate initial weights
  weights <- matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult
  
  #Generate place cells
  PC_x <- rep(0,N_pc) #1xN_pc matrix containing the x = 0 coordinate for each place cell
  PC_y <- rep(0,N_pc) #1xN_pc matrix containing the y = 0 coordinate for each place cell
  for (i in 1:N_pc) {
    #For each place cell:
    PC_x[i] <- (runif(1) - 0.5)*pool_diameter#Random positions of place cells
    PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    while ((PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2)){
      #Checks for out of bounds
      PC_x[i] <- (runif(1) - 0.5)*pool_diameter
      PC_y[i] <- (runif(1) - 0.5)*pool_diameter
    }
  }
  par(mfrow=c(2,2))
  
  for (day in 1:Ndays) {
    idxs = sample(4) #randomly choose 4 starting locations
    for (trial in 1:Ntrials){
      
      ## RANDOMLY FLIP THE PLATFORMS BETWEEN TWO POSSIBLE QUADRANTS
      if (variable_platform == 1) {
        whichplatform <- sample(c(1, -1), 1)
        platform_x <- base_platform_x * whichplatform
        platform_y <- base_platform_y * whichplatform
      } else {
        platform_x <- base_platform_x
        platform_y <- base_platform_y
      }
      
      idx <- idxs[trial] #take each location
      starting_x <- starting_xs[idx]
      starting_y <- starting_ys[idx]
      
      modresults <- run_trial (weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y, Vdecay, ac_const, beta, etdecay, alpha, gamma, Wnoise, platform_x, platform_y, starting_x, starting_y, speed, hitwall)
      #run trial
      weights <- modresults[[1]]
      track_x <- modresults[[2]]
      track_y <- modresults[[3]]
      vel_x <- modresults[[4]]
      vel_y <- modresults[[5]]
      
      #        weights <- wres
      
      PMs[1,day,trial] <- modresults[[9]] #latency
      PMs[2,day,trial] <- modresults[[6]] #dist
      PMs[3,day,trial] <- modresults[[8]][4]*100 #target quadrant
      PMs[4,day,trial] <- modresults[[8]][2]*100 #opposite quadrant
      PMs[5,day,trial] <- modresults[[7]]*100 #wall zone
      #record performance measures
      
      if (plot_trajectories & ((day==1 & trial==1) | (day==8 & trial==4)))
      {
        #plot the maze
        plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l", xlab = paste("day ",day,", trial ",trial), ylab = "trajectory")
        #plot the trajectory
        lines(track_x, track_y, type = "l")
        #plot the platform
        lines(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),type = "l")
      }
      
      if (plot_cognitive_maps & ((day==1 & trial==1) | (day==8 & trial==4)))
      {
        #plot the maze
        plot(pool_diameter/2*cos(th),pool_diameter/2*sin(th),type = "l",xlab = paste("day ",day,", trial ",trial), ylab = "cognitive map")
        #plot the cognitive map
        for (x in (-3:3)*(pool_diameter/6)){
          for (y in (-3:3)*(pool_diameter/6)){
            if (x^2 + y^2 <= (pool_diameter/2)^2){
              x2 = x
              y2 = y
              for (k in 1:N_ac){
                PC_activation <- rep(0,N_pc)
                for (i in 1:N_pc){
                  PC_activation[i] <- exp(-((x - PC_x[i])^2 + (y - PC_y[i])^2)/(2*sigma_pc^2))
                }
                #Calculate AC activation (i.e. value of the action)
                AC_activation <- rep(0,N_ac)
                for (i in 1:N_ac){
                  for (j in 1:N_pc){
                    AC_activation[i] <- AC_activation[i] + PC_activation[j]*weights[j,i]
                  }
                }
                x2 <- c(x2, x + (AC_activation[k]/10)*cos(k/N_ac*2*pi))
                y2 <- c(y2, y + (AC_activation[k]/10)*sin(k/N_ac*2*pi))
                #                 line([x x2],[y y2],'Color',[k/N_ac 0 1-k/N_ac])
              }
              lines(x2,y2,type = "l",col = "blue")
            }
          }
        }
        # plot the platform
        lines(platform_x+platform_radius*cos(th),platform_y+platform_radius*sin(th),type = "l")
      }
    }
  }
  
} else {
  # run multiple times without plotting!
  
  for (reps in 1:Nruns){
    
    # Generate initial weights for each run
    weights <- matrix(runif(N_pc*N_ac), nrow = N_pc) * Wmult
    
    # Generate place cells for each run
    PC_x <- rep(0, N_pc)
    PC_y <- rep(0, N_pc)
    
    for (i in 1:N_pc) {
      PC_x[i] <- (runif(1) - 0.5) * pool_diameter
      PC_y[i] <- (runif(1) - 0.5) * pool_diameter
      while (PC_x[i]^2 + PC_y[i]^2 > (pool_diameter/2)^2) {
        PC_x[i] <- (runif(1) - 0.5) * pool_diameter
        PC_y[i] <- (runif(1) - 0.5) * pool_diameter
      }
    }
    
    for (day in 1:Ndays) {
      idxs <- sample(4)
      
      for (trial in 1:Ntrials) {
        
        if (variable_platform == 1) {
          whichplatform <- sample(c(1, -1), 1)
          platform_x <- base_platform_x * whichplatform
          platform_y <- base_platform_y * whichplatform
        } else {
          platform_x <- base_platform_x
          platform_y <- base_platform_y
        }
        
        idx <- idxs[trial]
        starting_x <- starting_xs[idx]
        starting_y <- starting_ys[idx]
        
        modresults <- run_trial(
          weights, Wmult, sigma_pc, sigma_ac, PC_x, PC_y,
          Vdecay, ac_const, beta, etdecay, alpha, gamma,
          Wnoise, platform_x, platform_y,
          starting_x, starting_y, speed, hitwall
        )
        
        weights <- modresults[[1]]
        
        PMs[1, day, trial, reps] <- modresults[[9]]
        PMs[2, day, trial, reps] <- modresults[[6]]
        PMs[3, day, trial, reps] <- modresults[[8]][4] * 100
        PMs[4, day, trial, reps] <- modresults[[8]][2] * 100
        PMs[5, day, trial, reps] <- modresults[[7]] * 100
      }
    }
  }
  
  ## AVERAGE 50 RUNS FOR EACH TRIAL AND COMPUTE SE
  rows_list <- list()
  k <- 1
  
  for (day in 1:Ndays) {
    for (trial in 1:Ntrials) {
      latency_vals <- PMs[1, day, trial, ]
      dist_vals <- PMs[2, day, trial, ]
      target_vals <- PMs[3, day, trial, ]
      opposite_vals <- PMs[4, day, trial, ]
      wall_vals <- PMs[5, day, trial, ]
      
      rows_list[[k]] <- data.frame(
        day = day,
        trial = trial,
        trial_global = (day - 1) * Ntrials + trial,
        
        latency_pc = mean(latency_vals, na.rm = TRUE),
        latency_pc_se = sd(latency_vals, na.rm = TRUE) / sqrt(sum(!is.na(latency_vals))),
        
        dist_pc = mean(dist_vals, na.rm = TRUE),
        dist_pc_se = sd(dist_vals, na.rm = TRUE) / sqrt(sum(!is.na(dist_vals))),
        
        target_quadrant_pc = mean(target_vals, na.rm = TRUE),
        target_quadrant_pc_se = sd(target_vals, na.rm = TRUE) / sqrt(sum(!is.na(target_vals))),
        
        opposite_quadrant_pc = mean(opposite_vals, na.rm = TRUE),
        opposite_quadrant_pc_se = sd(opposite_vals, na.rm = TRUE) / sqrt(sum(!is.na(opposite_vals))),
        
        wall_zone_pc = mean(wall_vals, na.rm = TRUE),
        wall_zone_pc_se = sd(wall_vals, na.rm = TRUE) / sqrt(sum(!is.na(wall_vals)))
      )
      
      k <- k + 1
    }
  }
  
  ## STORE THE RAW RUN LEVELS DATA TOO
  rows_run <- list()
  k_run <- 1
  
  for (r in 1:Nruns) {
    for (day in 1:Ndays) {
      for (trial in 1:Ntrials) {
        rows_run[[k_run]] <- data.frame(
          run = r,
          day = day,
          trial = trial,
          trial_global = (day - 1) * Ntrials + trial,
          latency = PMs[1, day, trial, r],
          distance = PMs[2, day, trial, r],
          target_quadrant = PMs[3, day, trial, r],
          opposite_quadrant = PMs[4, day, trial, r],
          wall_zone = PMs[5, day, trial, r]
        )
        k_run <- k_run + 1
      }
    }
  }
  
  ## OUTPUT THE RAW RUNS CSV, MEAN + SE CSV, AND THE RAW PM ARRAY
  run_level_pc <- do.call(rbind, rows_run)
  write.csv(run_level_pc, "place_cell_variable_run_level.csv", row.names = FALSE)
  
  trial_summary_pc <- do.call(rbind, rows_list)
  write.csv(trial_summary_pc, "place_cell_variable_trial_mean_over_runs.csv", row.names = FALSE)
  
  saveRDS(PMs, "PMs_place_variable.rds")
}
