 
 short |     long             |     default     | algorithm                |                   description
------ | -------------------- | ----------------|--------------------------|-------------------------------------------------
  q    |   quiet              |       -         |  smc, kalman, kmcmc, ksimplex  |  no verbosity
  P    |   pipe               |       -         |  smc, kalman, kmcmc, ksimplex      |  pipe mode (echo theta.json on stdout)
       |   no_dem_sto         |       -         |  smc, kalman, kmcmc, ksimplex      |  turn off demographic stochasticity (if possible)*?*
       |   no_white_noise     |       -         |  smc, kalman, kmcmc, ksimplex      |  turn off white noises (if any)
       |   no_diff            |       -         |  smc, kalman, kmcmc, ksimplex      |  turn off diffusions (if any)
  s    |   DT                 |      0.0        |  smc, kalman, kmcmc, ksimplex      |  integration time step
       |   eps_abs            |      1e-6       |  smc, kalman, kmcmc, ksimplex      |  absolute error for adaptive step-size control
       |   eps_rel            |      1e-3       |  smc, kalman, kmcmc, ksimplex      |  relative error for adaptive step-size control
  g    |   freeze_forcing     |      -1         |  smc, kalman, kmcmc, ksimplex      |  freeze covariates to their value at specified time
  i    |   id                 |       0         |  smc, kalman, kmcmc, ksimplex      |  general id (unique integer identifier that will be appended to the output)
  p    |   path               |      '/'        |  smc, kalman, kmcmc, ksimplex      |  path where the outputs will be stored
  N    |   n_thread           |      1          |  smc              |  number of threads to be used
  l    |   LIKE_MIN           |      1e-17      |  smc, kalman, kmcmc, ksimplex      |  particles with likelihood smaller than LIKE_MIN are considered lost *description does not apply to Kalman*
  J    |                      |      1          |  smc              |  number of particles
  o    |   nb_obs             |      -1         |  smc, kalman, kmcmc, ksimplex      |  number of observations to be fitted (for tempering)
  I    |   interpolation      |      -          |  smc, kalman, kmcmc, ksimplex      |  gsl interpolator for covariates
  r    |   traj               |      -          |  smc, kalman      |  print the trajectories
  t    |   no_filter          |      -          |  smc              |  do not filter
  b    |   no_trace           |      -          |  smc, kalman, ksimplex      |  do not write trace output files
  e    |   no_hat             |      -          |  smc, kalman      |  do not write hat output files
  d    |   no_pred_res        |      -          |  smc, kalman      |  do not write pred_red output files
       |   prior              |      -          |  smc, kalman, ksimplex      |  add log(prior) to the estimated loglikelihood
       |   transf             |      -          |  kalman, ksimplex           |  add log(JacobianDeterminant(transf)) to the estimated loglikelihood
       |   full               |      -          |  kmcmc            |  full update MVN mode
  a    |   cooling            |      0.999      |  kmcmc            |  cooling factor for sampling covariance live tuning
  S    |   switch             |      -1         |  kmcmc            |  select switching iteration from initial covariance to empirical one
  E    |   epsilon            |      50         |  kmcmc            |  select number of burnin iterations before tuning epsilon
       |   epsilon_max        |      50         |  kmcmc            |  maximum value allowed for epislon
       |   alpha              |      0.02       |  kmcmc            |  smoothing factor of exponential smoothing used to compute smoothed acceptance rate (low values increase degree of smoothing)
       |   smooth             |      -          |  kmcmc            |  tune epsilon with the value of the acceptance rate obtained with exponential smoothing
  n    |   n_traj             |      1000       |  kmcmc            |  number of trajectories stored
       |   acc                |      -          |  kmcmc            |  print the acceptance rate *isn't that default mode?*
  M    |   iter               |      10         |  kmcmc, ksimplex            |  number of pmcmc iterations
  S    |   size               |      1e-6       |  ksimplex            |  simplex size used as stopping criteria
  h    |   help               |      -          |  all              |  print the usage on stdout
  
