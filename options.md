
 
 short |     long             |     default     | algorithm |                   description
------ | -------------------- | ----------------|-----------|-------------------------------------------------
  q    |   quiet              |                 |  smc      |  no verbosity
  P    |   pipe               |                 |  smc      |  pipe mode (echo theta.json on stdout)
       |   no_dem_sto         |                 |  smc      |  turn off demographic stochasticity *(if possible)*
       |   no_white_noise     |                 |  smc      |  turn off white noises (if any)
       |   no_diff            |                 |  smc      |  turn off diffusions (if any)
  s    |   DT                 |                 |  smc      |  integration time step
       |   eps_abs            |                 |  smc      |  absolute error for adaptive step-size control
       |   eps_rel            |                 |  smc      |  relative error for adaptive step-size control
  g    |   freeze_forcing     |      -1         |  smc      |  freeze covariates to their value at specified time
  i    |   id                 |                 |  smc      |  general id (unique integer identifier that will be appended to the output)
  p    |   path               |                 |  smc      |  path where the outputs will be stored
  N    |   n_thread           |                 |  smc      |  number of threads to be used
  l    |   LIKE_MIN           |     1e-17       |  smc      |  particles with likelihood smaller than LIKE_MIN are considered lost
  J    |   no_white_noise     |                 |  smc      |  number of particles
  o    |   nb_obs             |      -1         |  smc      |  number of observations to be fitted (for tempering)
  I    |   interpolation      |                 |  smc      |  gsl interpolator for covariates
  r    |   traj               |                 |  smc      |  print the trajectories
  t    |   no_filter          |                 |  smc      |  do not filter
  b    |   no_trace           |                 |  smc      |  do not write trace output files
  e    |   no_hat             |                 |  smc      |  do not write hat output files
  d    |   no_pred_res        |                 |  smc      |  do not write pred_red output files
       |   prior              |                 |  smc      |  add log(prior) to the estimated likelihood
  h    |   help               |                 |  smc      |  print the usage on stdout
  
