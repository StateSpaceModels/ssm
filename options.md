 
 short |     long             |     default     | algorithm                |                   description
------ | -------------------- | ----------------|--------------------------|-------------------------------------------------
  q    |   quiet              |       -         |  smc, kalman, kmcmc, ksimplex, mif  |  no verbosity
  P    |   pipe               |       -         |  smc, kalman, kmcmc, ksimplex, mif      |  pipe mode (echo theta.json on stdout)
       |   no_dem_sto         |       -         |  smc, kalman, kmcmc, ksimplex, mif      |  turn off demographic stochasticity (if possible)*?*
       |   no_white_noise     |       -         |  smc, kalman, kmcmc, ksimplex, mif      |  turn off white noises (if any)
       |   no_diff            |       -         |  smc, kalman, kmcmc, ksimplex, mif      |  turn off diffusions (if any)
  s    |   DT                 |      0.0        |  smc, kalman, kmcmc, ksimplex, mif      |  integration time step
       |   eps_abs            |      1e-6       |  smc, kalman, kmcmc, ksimplex, mif      |  absolute error for adaptive step-size control *does that make sense for other than simplex and ode?*
       |   eps_rel            |      1e-3       |  smc, kalman, kmcmc, ksimplex, mif      |  relative error for adaptive step-size control *does that make sense for other than simplex and ode?*
  g    |   freeze_forcing     |      -1         |  smc, kalman, kmcmc, ksimplex, mif      |  freeze covariates to their value at specified time
  i    |   id                 |       0         |  smc, kalman, kmcmc, ksimplex, mif      |  general id (unique integer identifier that will be appended to the output)
  p    |   path               |      '/'        |  smc, kalman, kmcmc, ksimplex, mif      |  path where the outputs will be stored
  N    |   n_thread           |      1          |  smc, mif              |  number of threads to be used
  l    |   LIKE_MIN           |      1e-17      |  smc, kalman, kmcmc, ksimplex, mif      |  particles with likelihood smaller than LIKE_MIN are considered lost *description does not apply to Kalman*
  J    |                      |      1          |  smc, mif              |  number of particles
  o    |   nb_obs             |      -1         |  smc, kalman, kmcmc, ksimplex      |  number of observations to be fitted (for tempering) *what about mif?*
  I    |   interpolation      |      -          |  smc, kalman, kmcmc, ksimplex, mif      |  gsl interpolator for covariates
  r    |   traj               |      -          |  smc, kalman, mif      |  print the trajectories
  t    |   no_filter          |      -          |  smc              |  do not filter
  b    |   no_trace           |      -          |  smc, kalman, ksimplex      |  do not write trace output files
  e    |   no_hat             |      -          |  smc, kalman      |  do not write hat output files
  d    |   no_pred_res        |      -          |  smc, kalman      |  do not write pred_red output files
       |   prior              |      -          |  smc, kalman, ksimplex, mif      |  add log(prior) to the estimated loglikelihood
       |   transf             |      -          |  kalman, ksimplex          |  add log(JacobianDeterminant(transf)) to the estimated loglikelihood
  a    |   cooling            |      0.999 (kmcmc) / 0.975 (mif)      |  kmcmc, mif            |  cooling factor for sampling covariance live tuning
  S    |   switch             |      -1 (kmcmc) / 5 (mif)             |  kmcmc, mif            |  select switching iteration from initial covariance to empirical one (mcmc) or to update formula introduced in Ionides et al. 2006 (mif)
  E    |   epsilon            |      50         |  kmcmc            |  select number of burnin iterations before tuning epsilon
       |   epsilon_max        |      50         |  kmcmc            |  maximum value allowed for epislon
       |   alpha              |      0.02       |  kmcmc            |  smoothing factor of exponential smoothing used to compute smoothed acceptance rate (low values increase degree of smoothing)
       |   smooth             |      -          |  kmcmc            |  tune epsilon with the value of the acceptance rate obtained with exponential smoothing
  n    |   n_traj             |      1000       |  kmcmc            |  number of trajectories stored
       |   acc                |      -          |  kmcmc            |  print the acceptance rate *isn't that default mode?*
       |   full               |      -          |  kmcmc            |  full update MVN mode
  b    |   heat               |      2          |  mif            |  re-heating accross MIF iterations (scales standard deviatio of proposals)
  L    |   lag                |      0.75          |  mif            |  lag for fixed-lag smoothing (proportion of the data)
  f    |   ic_only            |      -          |  mif            |  only fit the initial condition using fixed lag smoothing
  M    |   iter               |      10         |  kmcmc, ksimplex, mif            |  number of pmcmc iterations
  S    |   size               |      1e-6       |  ksimplex            |  simplex size used as stopping criteria
  h    |   help               |      -          |  all              |  print the usage on stdout
  
