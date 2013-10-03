 
### Options with arguments (upper case)
 
 short |     long             |     default                           | algorithm                                                  |                   description
------ | -------------------- | --------------------------------------|------------------------------------------------------------|-------------------------------------------------
  D    |   dt                 |      0.0                              |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  integration time step
  I    |   id                 |       0                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  general id (unique integer identifier that will be appended to the output)
  P    |   path               |      '/'                              |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  path where the outputs will be stored
  N    |   n_thread           |      1                                |  smc, pmcmc, mif, simul                                    |  number of threads to be used
  J    |   n_parts            |      100                              |  smc, pmcmc, pmcmc, mif, simul                             |  number of particles
  O    |   n_obs              |      -1                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif              |  number of observations to be fitted (for tempering) 
  A    |   cooling            |      0.98 (kmcmc) / 0.98 (mif)        |  kmcmc, pmcmc, mif                                         |  cooling factor for sampling covariance live tuning
  C    |   cov_switch         |      -1 (kmcmc) / 5 (mif)             |  kmcmc, pmcmc, mif                                         |  select switching iteration from initial covariance to empirical one (mcmc) or to update formula introduced in Ionides et al. 2006 (mif)
  E    |   eps_switch         |      50                               |  kmcmc, pmcmc                                              |  select number of burnin iterations before tuning epsilon
  T    |   n_traj             |      1000                             |  kmcmc, pmcmc                                              |  number of trajectories stored
  M    |   iter               |      10                               |  kmcmc, pmcmc, ksimplex, simplex, mif                      |  number of pmcmc iterations
  B    |   t0                 |      0.0                              |  simul                                                     |  time step when simulation starts (in unit of frequency, see --freq)
  E    |   tend               |      0.0                              |  simul                                                     |  time step when simulation ends (in unit of frequency, see --freq)
  Y    |   eps_abs_integ      |      1e-6                             |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  absolute error for adaptive step-size control 
  Z    |   eps_rel_integ      |      1e-3                             |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  relative error for adaptive step-size control 
  G    |   freeze_forcing     |      -1                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  freeze covariates to their value at specified time
  K    |   like_min           |      1e-17                            |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif         |  if applicable, particles with likelihood smaller than LIKE_MIN are considered lost. Otherwise, lower bound on likelihood.
  U    |   eps_max            |      50                               |  kmcmc, pmcmc                                              |  maximum value allowed for epislon
  S    |   alpha              |      0.02                             |  kmcmc, pmcmc                                              |  smoothing factor of exponential smoothing used to compute smoothed acceptance rate (low values increase degree of smoothing)
  H    |   heat               |      2                                |  mif                                                       |  re-heating accross MIF iterations (scales standard deviatio of proposals)
  L    |   lag                |      0.75                             |  mif                                                       |  lag for fixed-lag smoothing (proportion of the data)
  F    |   freq               |      D                                |  simul                                                     |  print the outputs (and reset incidences to 0 if any) every day (D), week (W), bi-week (B), month (M, taken to be 12.0/365) or year (Y)            |  simplex size used as stopping criteria
  Z    |   size               |      1e-6                             |  ksimplex, simplex                                         |  simplex size used as stopping criteria

### Options without arguments (lower case)
 
 short |     long             |     default                           | algorithm                                                  |                   description
------ | -------------------- | --------------------------------------|------------------------------------------------------------|-------------------------------------------------
  h    |   help               |      -                                |  all                                                       |  print the usage on stdout
  q    |   quiet              |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  no verbosity
  d    |   no_dem_sto         |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off demographic stochasticity
  w    |   no_white_noise     |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off white noises (if any)
  f    |   no_diff            |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off diffusions (if any)
  i    |   interpolation      |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  gsl interpolator for covariates
  t    |   traj               |      -                                |  smc, kalman, mif, simul                                   |  print the trajectories
  r    |   no_filter          |      -                                |  smc                                                       |  do not filter
  c    |   trace              |      -                                |  smc, kalman, ksimplex, simplex                            |  do not write trace output files
  h    |   hat                |      -                                |  smc, kalman                                               |  do not write hat output files
  e    |   pred_res           |      -                                |  smc, kalman                                               |  do not write pred_red output files
  p    |   prior              |      -                                |  smc, kalman, ksimplex, simplex, mif                       |  add log(prior) to the estimated loglikelihood
  s    |   smooth             |      -                                |  kmcmc, pmcmc                                              |  tune epsilon with the value of the acceptance rate obtained with exponential smoothing
  a    |   acc                |      -                                |  kmcmc, pmcmc                                              |  print the acceptance rate 
  z    |   zmq                |      -                                |  pmcmc                                                     |  dispatch particles across machines using a zeromq pipeline
  k    |   chunk              |      -                                |  pmcmc                                                     |  number of particles sent to each machine
  b    |   ic_only            |      -                                |  mif                                                       |  only fit the initial condition using fixed lag smoothing
  l    |   least_squares      |      -                                |  simplex                                                   |  minimize the sum of squared errors instead of maximizing the likelihood
  x    |   transiant          |      -                                |  simul        
  g    |   seed_time          |      -                                |  smc, kmcmc, pmcmc, mif, simul

### outdated options

 short |     long             |     default                           | algorithm                                                  |                   description
------ | -------------------- | --------------------------------------|------------------------------------------------------------|-------------------------------------------------
       |   continue           |      5                                |   simul                                                    |  print the final states in a bifurcation analysis to allow continuation
  b    |   bif                |      -                                |   simul                                                    |  run a bifurcation analysys
  l    |   lyap               |      -                                |   simul                                                    |  compute Lyapunov exponents
  d    |   period_dyn         |      -                                |   simul                                                    |  compute period (dynamical system def)
  u    |   fft                |      -                                |   simul                                                    |  compute period (FFT)
  B    |   block              |      5                                |   simul                                                    |  tuning parameter for max and min detection (has to be an odd number)
  x    |   precision          |      5                                |   simul                                                    |  smallest significant difference to detect variation for min and max detection
  P    |   pipe               |       -                               |   smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  pipe mode (echo theta.json on stdout)
       |   full               |      -                                |   kmcmc, pmcmc                                              |  full update MVN mode



