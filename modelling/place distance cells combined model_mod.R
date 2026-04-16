# ======================= run trial function =======================
run_trial <- function(weights0_pc, weights0_dc, Wmult, sigma_pc, sigma_dc, sigma_ac, PC_x, PC_y, DC, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, platform_x, platform_y, starting_x, starting_y, speed, wall_pun, weight_wall, weight_place)# c
{
  # FIXED PARAMETERS OF THE EXPERIMENT
  
  pool_diameter <- 1.4 #Maze diameter in metres (m)
  platform_radius <- 0.06 #Platform radius
  
  N_dc <- 10 #Population of place cells
  N_pc <- 211 #Population of place cells # c
  N_ac <- 36 #Population of action cells
  which_pc <- 0 # c
  which_dc <- 0 # c
  which <- 0 # c
  
  dist <- 0
  wall_zone <- 0
  quadrants <- c(0,0,0,0) #Percentage spent on each quadrant
  
  weights_pc <- weights0_pc #Initialize modifiable weights # c
  weights_dc <- weights0_dc #Initialize modifiable weights # c
  
  
  el_tr_dc <- matrix(rep(0, N_dc*N_ac), nrow = N_dc) #Initialize eligibility traces matrix # c
  el_tr_pc <- matrix(rep(0, N_pc*N_ac), nrow = N_pc) #Initialize eligibility traces matrix # c
  
  #Initialize trajectories
  track_x <- starting_x #Current position of trajectory is equal to the starting location of the animal
  track_y <- starting_y
  vel_x <- 0
  vel_y <- 0
  
  
  # NAVIGATION LOOP
  while ((track_x[length(track_x)] - platform_x)^2 + (track_y[length(track_y)] - platform_y)^2 > platform_radius^2) # while out of platform
  {
    
    weights_dc <- weights_dc*(1-noise) + matrix(runif(N_dc*N_ac), nrow = N_dc)*Wmult*noise # c
    weights_pc <- weights_pc*(1-noise) + matrix(runif(N_pc*N_ac), nrow = N_pc)*Wmult*noise # c
    
    dist.to.wall <- pool_diameter/2 - sqrt(track_x[length(track_x)]^2+track_y[length(track_y)]^2)
    
    #Calculate DC and PC activation
    DC_activation <- rep(0, N_dc)
    for (i in 1:N_dc){
      DC_activation[i] <- exp(-(dist.to.wall -DC[i])^2/(2*sigma_dc^2))
    }
    PC_activation <- rep(0, N_pc) # c
    for (i in 1:N_pc){# c
      PC_activation[i] <- exp(-((track_x[length(track_x)] - PC_x[i])^2 + (track_y[length(track_y)] - PC_y[i])^2)/(2*sigma_pc^2))
    }# c
    
    #Calculate AC activation (i.e. value of the action, Q)
    if (length(track_x) > 1){
      prevQ_pc <- AC_activation_pc[which_pc] #Displays the Q value before movement # c
      prevQ_dc <- AC_activation_dc[which_dc] #Displays the Q value before movement # c
    }
    #print(DC_activation)
    AC_activation_dc <- DC_activation %*% weights_dc # c
    AC_activation_pc <- PC_activation %*% weights_pc # c
    
    
    #Make an action
    ACsel_pc <- pmax(0, as.numeric(AC_activation_pc))^beta
    sum_pc <- sum(ACsel_pc)
    
    if (!is.finite(sum_pc) || sum_pc <= 0) {
      ACsel_pc <- rep(1 / N_ac, N_ac)
    } else {
      ACsel_pc <- ACsel_pc / sum_pc
    }
    
    ACsel_dc <- pmax(0, as.numeric(AC_activation_dc))^beta
    sum_dc <- sum(ACsel_dc)
    
    if (!is.finite(sum_dc) || sum_dc <= 0) {
      ACsel_dc <- rep(1 / N_ac, N_ac)
    } else {
      ACsel_dc <- ACsel_dc / sum_dc
    }
    
    shift <- round((180 - atan2(track_y[length(track_y)], track_x[length(track_x)]) / pi * 180) / 10)
    if (shift == 36) {
      shift <- 0
    }
    
    if (shift < 0) {
      ACsel_dc <- c(ACsel_dc[(shift + 37):36], ACsel_dc[1:(shift + 36)])
    } else if (shift > 0) {
      ACsel_dc <- c(ACsel_dc[(shift + 1):36], ACsel_dc[1:shift])
    }
    
    ACsel <- ACsel_dc * weight_wall + ACsel_pc * weight_place
    sum_comb <- sum(ACsel)
    
    if (!is.finite(sum_comb) || sum_comb <= 0) {
      ACsel <- rep(1 / N_ac, N_ac)
    } else {
      ACsel <- ACsel / sum_comb
    }
    
    ASrand <- runif(1)
    which <- 1
    ASsum <- ACsel[1]
    
    while (which < N_ac && is.finite(ASsum) && ASsum < ASrand) {
      which <- which + 1
      ASsum <- ASsum + ACsel[which]
    }
    
    if (which + shift> 36){
      which_dc = which + shift - 36
    }else if(which + shift <= 0){
      which_dc = which + shift + 36
    }else{
      which_dc = which + shift
    }
    
    
    which_pc = which
    
    
    #Eligibility traces
    el_tr_pc <- el_tr_pc * etdecay # c
    el_tr_dc <- el_tr_dc * etdecay # c
    
    for (j in 1:N_ac){# c
      itmp_pc <- min(abs(j-which_pc), N_ac-abs(j-which_pc))
      itmp_dc <- min(abs(j-which_dc), N_ac-abs(j-which_dc))
      actgaus_pc <- exp(-(itmp_pc*itmp_pc)/(2*sigma_ac*sigma_ac))
      actgaus_dc <- exp(-(itmp_dc*itmp_dc)/(2*sigma_ac*sigma_ac))
      el_tr_pc[,j] <- el_tr_pc[,j] + actgaus_pc*AC_activation_pc[j]*t(t(PC_activation))
      el_tr_dc[,j] <- el_tr_dc[,j] + actgaus_dc*AC_activation_dc[j]*t(t(DC_activation))
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
      
      currQ_pc = AC_activation_pc[which_pc]# c
      currQ_dc = AC_activation_dc[which_dc]# c
      tderr_pc = rew + discf*currQ_pc - prevQ_pc #temporal difference error # c
      tderr_dc = rew + discf*currQ_dc - prevQ_dc #temporal difference error # c
      weights_pc = pmax(weights_pc + lrate*tderr_pc*el_tr_pc, 0) # c
      weights_dc = pmax(weights_dc + lrate*tderr_dc*el_tr_dc, 0) # c
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
    latency = (length(track_x)-1) * speed_ts / speed #convert to seconds
    if (latency > 60) # if more than a minute, stop
    { break }
    }
    
  }
  
  latency <- length(track_x)-1 # latency in time steps
  wall_zone <- wall_zone/latency
  quadrants <- quadrants/latency
  speed_ts <- mean(sqrt((vel_x[-1]^2+vel_y[-1]^2))) # speed in meters/time step
  # speed per action step from Hanbing
  speed_ps = (vel_x[-1]^2+vel_y[-1]^2)^0.5
  
  # time step
  time_step = speed_ts/speed
  
  # mean turning angle 
  vel = as.matrix(data.frame(vel_x,vel_y))
  angle=c()
  for (steps in 2:(length(vel_x)-1)){
    A = vel[steps,]
    B = vel[steps+1,]
    angle = c(angle, acos((A%*%B)[1,1])/norm(as.matrix(A)*norm(as.matrix(B))))
  }
  angle = angle * 180 / pi
  mean_angle = mean(angle)
  
  # speed standard deviation
  speed_std = sd((vel_x[-1]^2+vel_y[-1]^2)^0.5)
  speed_std = speed_std / time_step
  
  latency <- latency * speed_ts / speed # latency in seconds
  return(list(weights_pc, weights_dc, track_x, track_y, vel_x, vel_y, dist, wall_zone, quadrants, latency, speed_std, speed_ps,mean_angle,time_step))
}


# =========================== main ===============================
timestart <- Sys.time()

N_dc <- 10
N_pc <- 211
N_ac <- 36

plot_trajectories <- 0
plot_cognitive_maps <- 0
plot_integrated_cognitive_map <- 0
pln <- plot_trajectories + plot_cognitive_maps + plot_integrated_cognitive_map
Nruns <- 50

pool_diameter <- 1.4
platform_radius <- 0.06
sigma_pc <- 0.1
sigma_dc <- c(0.1, 0.1)
sigma_ac <- c(2, 2)

etdecay <- 0.83
beta <- c(6, 6)
alpha <- c(0.006, 0.006)
gamma <- c(0.85, 0.85)

Vdecay <- c(0.82, 0.82)
ac_const <- c(0.02, 0.02)
Wnoise <- c(0.0004, 0.0004)
Wmult <- c(0.1, 0.1)
hitwall <- c(0, 0)
speed <- c(0.175, 0.175)

Ntrials <- 4
Ndays <- 8
Npsets <- 1

# choose platform condition here
variable_platform <- 1   # 0 = fixed, 1 = variable

## RATIOS TO RUN AUTOMATICALLY
place_ratios <- sort(unique(c(seq(0, 1, by = 0.2), 0.5)))

## SAVE OUTPUTS
out_dir <- if (variable_platform == 1) "combined_variable_outputs" else "combined_fixed_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## STORE ALL RAW PMs ARRAY IN ONE LIST
PMs_all <- list()

params <- array(rep(0, 12 * Npsets), c(12, Npsets))
params[1, ] <- runif(Npsets) * (sigma_dc[2] - sigma_dc[1]) + sigma_dc[1]
params[2, ] <- runif(Npsets) * (sigma_ac[2] - sigma_ac[1]) + sigma_ac[1]
params[3, ] <- runif(Npsets) * (beta[2] - beta[1]) + beta[1]
params[4, ] <- runif(Npsets) * (alpha[2] - alpha[1]) + alpha[1]
params[5, ] <- runif(Npsets) * (gamma[2] - gamma[1]) + gamma[1]
params[6, ] <- runif(Npsets) * (Vdecay[2] - Vdecay[1]) + Vdecay[1]
params[7, ] <- runif(Npsets) * (ac_const[2] - ac_const[1]) + ac_const[1]
params[8, ] <- runif(Npsets) * (Wnoise[2] - Wnoise[1]) + Wnoise[1]
params[9, ] <- runif(Npsets) * (Wmult[2] - Wmult[1]) + Wmult[1]
params[10, ] <- runif(Npsets) * (hitwall[2] - hitwall[1]) + hitwall[1]
params[11, ] <- runif(Npsets) * (speed[2] - speed[1]) + speed[1]
params[12, ] <- etdecay

for (place_ratio in place_ratios) {
  
  weight_place <- place_ratio
  weight_wall <- 1 - place_ratio
  
  ratio_tag <- gsub("\\.", "p", sprintf("%.1f", place_ratio))
  wall_tag  <- gsub("\\.", "p", sprintf("%.1f", weight_wall))
  
  message("Running ratio: place = ", place_ratio, ", wall = ", weight_wall)
  
  track_x_sum <- vector(mode = "list", length = Ntrials * Ndays)
  track_y_sum <- vector(mode = "list", length = Ntrials * Ndays)
  
  if (pln > 0.5) {
    PMs <- array(rep(0, 9 * Ndays * Ntrials), c(9, Ndays, Ntrials))
    AMs <- array(rep(0, Ndays * Ntrials), c(Ndays, Ntrials))
  } else {
    PMs <- array(rep(0, 9 * Ndays * Ntrials * Nruns), c(9, Ndays, Ntrials, Nruns, Npsets))
    AMs <- array(rep(0, Ndays * Ntrials * Nruns), c(Ndays, Ntrials, Nruns))
  }
  
  for (pset in 1:Npsets) {
    
    message("  Parameter set ", pset)
    
    sigma_dc_now <- params[1, pset]
    sigma_ac_now <- params[2, pset]
    beta_now <- params[3, pset]
    alpha_now <- params[4, pset]
    gamma_now <- params[5, pset]
    Vdecay_now <- params[6, pset]
    ac_const_now <- params[7, pset]
    Wnoise_now <- params[8, pset]
    Wmult_now <- params[9, pset]
    hitwall_now <- params[10, pset]
    speed_now <- params[11, pset]
    etdecay_now <- params[12, pset]
    
    platform_x <- cos(-pi / 4) * pool_diameter / 4
    platform_y <- sin(-pi / 4) * pool_diameter / 4
    
    strad <- pool_diameter / 2 * 0.85
    starting_xs <- strad * c(cos(pi / 6), cos(pi / 3), cos(7 * pi / 6), cos(4 * pi / 3))
    starting_ys <- strad * c(sin(pi / 6), sin(pi / 3), sin(7 * pi / 6), sin(4 * pi / 3))
    
    th <- (0:100) / 50 * pi
    
    if (pln > 0.5) {
      
      weights_pc <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult_now
      weights_dc <- matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult_now
      
      DC <- rep(0, N_dc)
      PC_x <- rep(0, N_pc)
      PC_y <- rep(0, N_pc)
      
      for (i in 1:N_dc) {
        DC[i] <- runif(1) * (pool_diameter / 2)
      }
      
      for (i in 1:N_pc) {
        PC_x[i] <- (runif(1) - 0.5) * pool_diameter
        PC_y[i] <- (runif(1) - 0.5) * pool_diameter
        while ((PC_x[i]^2 + PC_y[i]^2) > (pool_diameter / 2)^2) {
          PC_x[i] <- (runif(1) - 0.5) * pool_diameter
          PC_y[i] <- (runif(1) - 0.5) * pool_diameter
        }
      }
      
      par(mfrow = c(4, 8))
      
      for (day in 1:Ndays) {
        idxs <- sample(4)
        for (trial in 1:Ntrials) {
          
          whichplatform <- if (variable_platform == 1) sample(c(1, -1), 1) else 1
          
          idx <- idxs[trial]
          starting_x <- starting_xs[idx]
          starting_y <- starting_ys[idx]
          
          modresults <- run_trial(
            weights_pc, weights_dc, Wmult_now, sigma_pc, sigma_dc_now, sigma_ac_now,
            PC_x, PC_y, DC, Vdecay_now, ac_const_now, beta_now, etdecay_now,
            alpha_now, gamma_now, Wnoise_now,
            whichplatform * platform_x, whichplatform * platform_y,
            starting_x, starting_y, speed_now, hitwall_now,
            weight_wall, weight_place
          )
          
          weights_pc <- modresults[[1]]
          weights_dc <- modresults[[2]]
          track_x <- modresults[[3]]
          track_y <- modresults[[4]]
          
          track_x_sum[[(day - 1) * 4 + trial]] <- track_x
          track_y_sum[[(day - 1) * 4 + trial]] <- track_y
          
          PMs[1, day, trial] <- modresults[[10]]
          PMs[2, day, trial] <- modresults[[7]]
          
          if (whichplatform == 1) {
            PMs[3, day, trial] <- modresults[[9]][4] * 100
            PMs[4, day, trial] <- modresults[[9]][2] * 100
          } else {
            PMs[3, day, trial] <- modresults[[9]][2] * 100
            PMs[4, day, trial] <- modresults[[9]][4] * 100
          }
          
          PMs[5, day, trial] <- modresults[[8]] * 100
          PMs[6, day, trial] <- modresults[[11]] * 100
          PMs[7, day, trial] <- modresults[[13]]
          PMs[8, day, trial] <- modresults[[14]]
        }
      }
      
    } else {
      
      rows_run <- list()
      k_run <- 1
      
      for (reps in 1:Nruns) {
        message("    Run ", reps, "/", Nruns)
        
        weights_pc <- matrix(runif(N_pc * N_ac), nrow = N_pc) * Wmult_now
        weights_dc <- matrix(runif(N_dc * N_ac), nrow = N_dc) * Wmult_now
        
        DC <- rep(0, N_dc)
        PC_x <- rep(0, N_pc)
        PC_y <- rep(0, N_pc)
        
        for (i in 1:N_dc) {
          DC[i] <- runif(1) * (pool_diameter / 2)
        }
        
        for (i in 1:N_pc) {
          PC_x[i] <- (runif(1) - 0.5) * pool_diameter
          PC_y[i] <- (runif(1) - 0.5) * pool_diameter
          while ((PC_x[i]^2 + PC_y[i]^2) > (pool_diameter / 2)^2) {
            PC_x[i] <- (runif(1) - 0.5) * pool_diameter
            PC_y[i] <- (runif(1) - 0.5) * pool_diameter
          }
        }
        
        for (day in 1:Ndays) {
          idxs <- sample(4)
          
          for (trial in 1:Ntrials) {
            
            whichplatform <- if (variable_platform == 1) sample(c(1, -1), 1) else 1
            
            idx <- idxs[trial]
            starting_x <- starting_xs[idx]
            starting_y <- starting_ys[idx]
            
            modresults <- run_trial(
              weights_pc, weights_dc, Wmult_now, sigma_pc, sigma_dc_now, sigma_ac_now,
              PC_x, PC_y, DC, Vdecay_now, ac_const_now, beta_now, etdecay_now,
              alpha_now, gamma_now, Wnoise_now,
              whichplatform * platform_x, whichplatform * platform_y,
              starting_x, starting_y, speed_now, hitwall_now,
              weight_wall, weight_place
            )
            
            weights_pc <- modresults[[1]]
            weights_dc <- modresults[[2]]
            
            PMs[1, day, trial, reps, pset] <- modresults[[10]]
            PMs[2, day, trial, reps, pset] <- modresults[[7]]
            
            if (whichplatform == 1) {
              PMs[3, day, trial, reps, pset] <- modresults[[9]][4] * 100
              PMs[4, day, trial, reps, pset] <- modresults[[9]][2] * 100
            } else {
              PMs[3, day, trial, reps, pset] <- modresults[[9]][2] * 100
              PMs[4, day, trial, reps, pset] <- modresults[[9]][4] * 100
            }
            
            PMs[5, day, trial, reps, pset] <- modresults[[8]] * 100
            PMs[6, day, trial, reps, pset] <- modresults[[11]] * 100
            PMs[7, day, trial, reps, pset] <- modresults[[13]]
            PMs[8, day, trial, reps, pset] <- modresults[[14]]
            
            # -----------------------------
            # SAVE RUN-LEVEL RAW DATA
            # -----------------------------
            rows_run[[k_run]] <- data.frame(
              run = reps,
              pset = pset,
              day = day,
              trial = trial,
              trial_global = (day - 1) * Ntrials + trial,
              latency = PMs[1, day, trial, reps, pset],
              distance = PMs[2, day, trial, reps, pset],
              target_quadrant = PMs[3, day, trial, reps, pset],
              opposite_quadrant = PMs[4, day, trial, reps, pset],
              wall_zone = PMs[5, day, trial, reps, pset],
              speed_sd = PMs[6, day, trial, reps, pset],
              mean_angle = PMs[7, day, trial, reps, pset],
              time_step = PMs[8, day, trial, reps, pset],
              place_ratio = weight_place,
              wall_ratio = weight_wall,
              condition = ifelse(variable_platform == 1, "variable", "fixed")
            )
            k_run <- k_run + 1
          }
        }
      }
    }
  }
  
  ## SAVE RAW PMs ARRAY FOR THIS RATIO
  PMs_all[[paste0("place_", ratio_tag, "_wall_", wall_tag)]] <- PMs
  
  saveRDS(
    PMs,
    file = file.path(out_dir, paste0("PMs_place", ratio_tag, "_wall", wall_tag, ".rds"))
  )
  
  ## AVERAGE 50 RUNS FOR EACH TRIAL AND COMPUTE SE 
  rows_list <- list()
  k <- 1
  
  for (day in 1:Ndays) {
    for (trial in 1:Ntrials) {
      latency_vals <- PMs[1, day, trial, , 1]
      dist_vals <- PMs[2, day, trial, , 1]
      target_vals <- PMs[3, day, trial, , 1]
      opposite_vals <- PMs[4, day, trial, , 1]
      wall_vals <- PMs[5, day, trial, , 1]
      
      rows_list[[k]] <- data.frame(
        day = day,
        trial = trial,
        trial_global = (day - 1) * Ntrials + trial,
        
        latency = mean(latency_vals, na.rm = TRUE),
        latency_se = sd(latency_vals, na.rm = TRUE) / sqrt(sum(!is.na(latency_vals))),
        
        dist = mean(dist_vals, na.rm = TRUE),
        dist_se = sd(dist_vals, na.rm = TRUE) / sqrt(sum(!is.na(dist_vals))),
        
        target_quadrant = mean(target_vals, na.rm = TRUE),
        target_quadrant_se = sd(target_vals, na.rm = TRUE) / sqrt(sum(!is.na(target_vals))),
        
        opposite_quadrant = mean(opposite_vals, na.rm = TRUE),
        opposite_quadrant_se = sd(opposite_vals, na.rm = TRUE) / sqrt(sum(!is.na(opposite_vals))),
        
        wall_zone = mean(wall_vals, na.rm = TRUE),
        wall_zone_se = sd(wall_vals, na.rm = TRUE) / sqrt(sum(!is.na(wall_vals))),
        
        place_ratio = weight_place,
        wall_ratio = weight_wall,
        condition = ifelse(variable_platform == 1, "variable", "fixed")
      )
      
      k <- k + 1
    }
  }
  
  trial_summary <- do.call(rbind, rows_list)
  
  write.csv(
    trial_summary,
    file = file.path(
      out_dir,
      paste0(
        "combined_place", ratio_tag,
        "_wall", wall_tag,
        ifelse(variable_platform == 1, "_var", "_fixed"),
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  ## STORE THE RAW RUN LEVELS DATA TOO
  run_level <- do.call(rbind, rows_run)
  
  write.csv(
    run_level,
    file = file.path(
      out_dir,
      paste0(
        "combined_place", ratio_tag,
        "_wall", wall_tag,
        ifelse(variable_platform == 1, "_var_run_level", "_fixed_run_level"),
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# SAVE ALL PMs ARRAY TOGETHER
saveRDS(
  PMs_all,
  file = file.path(
    out_dir,
    paste0("PMs_all_ratios", ifelse(variable_platform == 1, "_var", "_fixed"), ".rds")
  )
)

timeend <- Sys.time()
print(timeend - timestart)

# p1 = read.csv("fixed_start-wall_p-0.004_1.csv")[,2:7]
# p2 = read.csv("fixed_start-wall_p-0.004_2.csv")[,2:7]
# p3 = read.csv("fixed_start-wall_p-0.004_3.csv")[,2:7]
# p4 = read.csv("fixed_start-wall_p-0.004_4.csv")[,2:7]
# p5 = read.csv("fixed_start-wall_p-0.004_5.csv")[,2:7]
# p6 = read.csv("fixed_start-wall_p-0.004_6.csv")[,2:7]
# p7 = read.csv("fixed_start-wall_p-0.004_7.csv")[,2:7]
# p8 = read.csv("fixed_start-wall_p-0.004_8.csv")[,2:7]
# p9 = read.csv("fixed_start-wall_p-0.004_9.csv")[,2:7]
# p10 = read.csv("fixed_start-wall_p-0.004_10.csv")[,2:7]
# 
# all = rbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
# 
# write.csv(all,"fixed_start-wall_p-0.004.csv")

