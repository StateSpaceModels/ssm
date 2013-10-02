 
 short |     long             |     default     | algorithm         |                   description
------ | -------------------- | ----------------|-------------------|-------------------------------------------------
  q    |   quiet              |                 |  smc, kalman      |  no verbosity
  P    |   pipe               |                 |  smc, kalman      |  pipe mode (echo theta.json on stdout)
       |   no_dem_sto         |                 |  smc, kalman      |  turn off demographic stochasticity (if possible)*?*
       |   no_white_noise     |                 |  smc, kalman      |  turn off white noises (if any)
       |   no_diff            |                 |  smc, kalman      |  turn off diffusions (if any)
  s    |   DT                 |                 |  smc, kalman      |  integration time step
       |   eps_abs            |      1e-6       |  smc, kalman      |  absolute error for adaptive step-size control
       |   eps_rel            |      1e-3       |  smc, kalman      |  relative error for adaptive step-size control
  g    |   freeze_forcing     |      -1         |  smc, kalman      |  freeze covariates to their value at specified time
  i    |   id                 |                 |  smc, kalman      |  general id (unique integer identifier that will be appended to the output)
  p    |   path               |                 |  smc, kalman      |  path where the outputs will be stored
  N    |   n_thread           |       1         |  smc              |  number of threads to be used
  l    |   LIKE_MIN           |     1e-17       |  smc, kalman      |  particles with likelihood smaller than LIKE_MIN are considered lost *description does not apply to Kalman*
  J    |                      |                 |  smc              |  number of particles
  o    |   nb_obs             |      -1         |  smc, kalman      |  number of observations to be fitted (for tempering)
  I    |   interpolation      |                 |  smc, kalman      |  gsl interpolator for covariates
  r    |   traj               |                 |  smc, kalman      |  print the trajectories
  t    |   no_filter          |                 |  smc              |  do not filter
  b    |   no_trace           |                 |  smc, kalman      |  do not write trace output files
  e    |   no_hat             |                 |  smc, kalman      |  do not write hat output files
  d    |   no_pred_res        |                 |  smc, kalman      |  do not write pred_red output files
       |   prior              |                 |  smc, kalman      |  add log(prior) to the estimated loglikelihood
       |   transf             |                 |  kalman           |  add log(JacobianDeterminant(transf)) to the estimated loglikelihood
  h    |   help               |                 |  all              |  print the usage on stdout
  
