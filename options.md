 
### Options with arguments (lower case)
 
 short |     long             |     default                           | algorithm                                                  |                   description
------ | -------------------- | --------------------------------------|------------------------------------------------------------|-------------------------------------------------
  d    |   dt                 |      0.0                              |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  integration time step
  i    |   id                 |       0                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  general id (unique integer identifier that will be appended to the output)
  p    |   path               |      '/'                              |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  path where the outputs will be stored
  n    |   n_thread           |      1                                |  smc, pmcmc, mif, simul                                    |  number of threads to be used
  j    |   n_parts            |      1                                |  smc, pmcmc, pmcmc, mif, simul                             |  number of particles
  o    |   n_obs              |      -1                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex              |  number of observations to be fitted (for tempering) *what about mif?*
  a    |   cooling            |      0.999 (kmcmc) / 0.975 (mif)      |  kmcmc, pmcmc, mif                                         |  cooling factor for sampling covariance live tuning
  c    |   cov_switch         |      -1 (kmcmc) / 5 (mif)             |  kmcmc, pmcmc, mif                                         |  select switching iteration from initial covariance to empirical one (mcmc) or to update formula introduced in Ionides et al. 2006 (mif)
  e    |   eps_switch         |      50                               |  kmcmc, pmcmc                                              |  select number of burnin iterations before tuning epsilon
  t    |   n_traj             |      1000                             |  kmcmc, pmcmc                                              |  number of trajectories stored
  m    |   iter               |      10                               |  kmcmc, pmcmc, ksimplex, simplex, mif                      |  number of pmcmc iterations
  b    |   t0                 |      0.0                              |  simul                                                    |  time step when simulation starts (in unit of frequency, see --freq)
  e    |   tend               |      0.0                              |  simul                                                    |  time step when simulation ends (in unit of frequency, see --freq)
  y    |   eps_abs_integ      |      1e-6                             |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  absolute error for adaptive step-size control *does that make sense for other than simplex and ode?*
  z    |   eps_rel_integ      |      1e-3                             |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  relative error for adaptive step-size control *does that make sense for other than simplex and ode?*
  g    |   freeze_forcing     |      -1                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  freeze covariates to their value at specified time
  k    |   like_min           |      1e-17                            |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif         |  particles with likelihood smaller than LIKE_MIN are considered lost *description does not apply to Kalman*
  u    |   eps_max            |      50                               |  kmcmc, pmcmc                                              |  maximum value allowed for epislon
  s    |   alpha              |      0.02                             |  kmcmc, pmcmc                                              |  smoothing factor of exponential smoothing used to compute smoothed acceptance rate (low values increase degree of smoothing)
  h    |   heat               |      2                                |  mif                                                       |  re-heating accross MIF iterations (scales standard deviatio of proposals)
  l    |   lag                |      0.75                             |  mif                                                       |  lag for fixed-lag smoothing (proportion of the data)
  f    |   freq               |      D                                |  simul                                                    |  print the outputs (and reset incidences to 0 if any) every day (D), week (W), bi-week (B), month (M, taken to be 12.0/365) or year (Y)            |  simplex size used as stopping criteria
  z    |   size               |      1e-6                             |  ksimplex, simplex                                         |  simplex size used as stopping criteria

### Options without arguments (upper case)
 
 short |     long             |     default                           | algorithm                                                  |                   description
------ | -------------------- | --------------------------------------|------------------------------------------------------------|-------------------------------------------------
  H    |   help               |      -                                |  all                                                       |  print the usage on stdout
  Q    |   quiet              |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  no verbosity
  D    |   no_dem_sto         |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off demographic stochasticity (if possible) *?*
  W    |   no_white_noise     |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off white noises (if any)
  F    |   no_diff            |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off diffusions (if any)
  I    |   interpolation      |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  gsl interpolator for covariates
  T    |   traj               |      -                                |  smc, kalman, mif, simul                                   |  print the trajectories
  R    |   no_filter          |      -                                |  smc                                                       |  do not filter
  C    |   no_trace           |      -                                |  smc, kalman, ksimplex, simplex                            |  do not write trace output files
  H    |   no_hat             |      -                                |  smc, kalman                                               |  do not write hat output files
  E    |   no_pred_res        |      -                                |  smc, kalman                                               |  do not write pred_red output files
  P    |   prior              |      -                                |  smc, kalman, ksimplex, simplex, mif                       |  add log(prior) to the estimated loglikelihood
  S    |   smooth             |      -                                |  kmcmc, pmcmc                                              |  tune epsilon with the value of the acceptance rate obtained with exponential smoothing
  J    |   transf             |      -                                |  kalman, ksimplex                                          |  add log(JacobianDeterminant(transf)) to the estimated loglikelihood *should we ditch this option??*
  A    |   acc                |      -                                |  kmcmc, pmcmc                                              |  print the acceptance rate *isn't that default mode?*
  Z    |   zmq                |      -                                |  pmcmc                                                     |  dispatch particles across machines using a zeromq pipeline
  K    |   chunk              |      -                                |  pmcmc                                                     |  number of particles sent to each machine
  B    |   ic_only            |      -                                |  mif                                                       |  only fit the initial condition using fixed lag smoothing
  L    |   least_squares      |      -                                |  simplex                                                   |  minimize the sum of squared errors instead of maximizing the likelihood
  X    |   transiant          |      -                                |  simul        

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


### plom-pipe options
 
 short |     long             |     default                           | algorithm                                                  |                   description
------ | -------------------- | --------------------------------------|------------------------------------------------------------|-------------------------------------------------
  q    |   quiet              |       -                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  no verbosity
  P    |   pipe               |       -                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  pipe mode (echo theta.json on stdout)
       |   no_dem_sto         |       -                               |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off demographic stochasticity (if possible) *?*
       |   no_white_noise     |       -                               |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off white noises (if any)
       |   no_diff            |       -                               |  smc, kalman, kmcmc, pmcmc, ksimplex, mif, simul           |  turn off diffusions (if any)
  s    |   DT                 |      0.0                              |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  integration time step
       |   eps_abs            |      1e-6                             |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  absolute error for adaptive step-size control *does that make sense for other than simplex and ode?*
       |   eps_rel            |      1e-3                             |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  relative error for adaptive step-size control *does that make sense for other than simplex and ode?*
  g    |   freeze_forcing     |      -1                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  freeze covariates to their value at specified time
  i    |   id                 |       0                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  general id (unique integer identifier that will be appended to the output)
  p    |   path               |      '/'                              |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  path where the outputs will be stored
  N    |   n_thread           |      1                                |  smc, pmcmc, mif, simul                                    |  number of threads to be used
  l    |   like_min           |      1e-17                            |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif         |  particles with likelihood smaller than LIKE_MIN are considered lost *description does not apply to Kalman*
  J    |                      |      1                                |  smc, pmcmc, pmcmc, mif, simul                             |  number of particles
  o    |   nb_obs             |      -1                               |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex              |  number of observations to be fitted (for tempering) *what about mif?*
  I    |   interpolation      |      -                                |  smc, kalman, kmcmc, pmcmc, ksimplex, simplex, mif, simul  |  gsl interpolator for covariates
  r    |   traj               |      -                                |  smc, kalman, mif, simul                                   |  print the trajectories
  t    |   no_filter          |      -                                |  smc                                                       |  do not filter
  b    |   no_trace           |      -                                |  smc, kalman, ksimplex, simplex                            |  do not write trace output files
  e    |   no_hat             |      -                                |  smc, kalman                                               |  do not write hat output files
  d    |   no_pred_res        |      -                                |  smc, kalman                                               |  do not write pred_red output files
       |   prior              |      -                                |  smc, kalman, ksimplex, simplex, mif                       |  add log(prior) to the estimated loglikelihood
       |   transf             |      -                                |  kalman, ksimplex                                          |  add log(JacobianDeterminant(transf)) to the estimated loglikelihood *should we ditch this option??*
  a    |   cooling            |      0.999 (kmcmc) / 0.975 (mif)      |  kmcmc, pmcmc, mif                                         |  cooling factor for sampling covariance live tuning
  S    |   switch             |      -1 (kmcmc) / 5 (mif)             |  kmcmc, pmcmc, mif                                         |  select switching iteration from initial covariance to empirical one (mcmc) or to update formula introduced in Ionides et al. 2006 (mif)
  E    |   epsilon            |      50                               |  kmcmc, pmcmc                                              |  select number of burnin iterations before tuning epsilon
       |   epsilon_max        |      50                               |  kmcmc, pmcmc                                              |  maximum value allowed for epislon
       |   alpha              |      0.02                             |  kmcmc, pmcmc                                              |  smoothing factor of exponential smoothing used to compute smoothed acceptance rate (low values increase degree of smoothing)
       |   smooth             |      -                                |  kmcmc, pmcmc                                              |  tune epsilon with the value of the acceptance rate obtained with exponential smoothing
  n    |   n_traj             |      1000                             |  kmcmc, pmcmc                                              |  number of trajectories stored
       |   acc                |      -                                |  kmcmc, pmcmc                                              |  print the acceptance rate *isn't that default mode?*
       |   full               |      -                                |  kmcmc, pmcmc                                              |  full update MVN mode
  Z    |   zmq                |      -                                |  pmcmc                                                     |  dispatch particles across machines using a zeromq pipeline
  c    |   chunk              |      -                                |  pmcmc                                                     |  number of particles sent to each machine
  b    |   heat               |      2                                |  mif                                                       |  re-heating accross MIF iterations (scales standard deviatio of proposals)
  L    |   lag                |      0.75                             |  mif                                                       |  lag for fixed-lag smoothing (proportion of the data)
  f    |   ic_only            |      -                                |  mif                                                       |  only fit the initial condition using fixed lag smoothing
  q    |   least_squares      |      -                                |  simplex                                                   |  minimize the sum of squared errors instead of maximizing the likelihood
  M    |   iter               |      10                               |  kmcmc, pmcmc, ksimplex, simplex, mif                      |  number of pmcmc iterations
  S    |   size               |      1e-6                             |  ksimplex, simplex                                         |  simplex size used as stopping criteria
  f    |   freq               |      D                                |  simul                                                    |  print the outputs (and reset incidences to 0 if any) every day (D), week (W), bi-week (B), month (M, taken to be 12.0/365) or year (Y)            |  simplex size used as stopping criteria
  o    |   t0                 |      0.0                              |  simul                                                    |  time step when simulation starts (in unit of frequency, see --freq)
  D    |   tend               |      0.0                              |  simul                                                    |  time step when simulation ends (in unit of frequency, see --freq)
  T    |   transiant          |      -                                |  simul                                                    |  skip a transiant of the specified length ( in number of steps in unit of frequency, see --freq)
  b    |   bif                |      -                                |  simul                                                    |  run a bifurcation analysys
  l    |   lyap               |      -                                |  simul                                                    |  compute Lyapunov exponents
  d    |   period_dyn         |      -                                |  simul                                                    |  compute period (dynamical system def)
  u    |   fft                |      -                                |  simul                                                    |  compute period (FFT)
  B    |   block              |      5                                |  simul                                                    |  tuning parameter for max and min detection (has to be an odd number)
  x    |   precision          |      5                                |  simul                                                    |  smallest significant difference to detect variation for min and max detection
       |   continue           |      5                                |  simul                                                    |  print the final states in a bifurcation analysis to allow continuation
  h    |   help               |      -                                |  all                                                       |  print the usage on stdout
  

