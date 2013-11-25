require('pomp')

sir.psr <- function (x, t, params, delta.t, ...) {
    v <- 1/(params["v"]/7)  #time (t) is in week
    beta <- params["r0"]*v
    foi <- beta*x["I"]/(x["S"]+x["I"]+x["R"])
    
    trans <- c(
               reulermultinom(n=1,size=x["S"],rate=c(foi),dt=delta.t), # exits from S
               reulermultinom(n=1,size=x["I"],rate=c(v),dt=delta.t) # exits from I
               )

    ## now connect the compartments
    x[c("S","I","R","cases")]+c(
        -trans[1],
        trans[1]-trans[2],
        trans[2],
        trans[1] # accumulate the S->I
        )
}

sir.ode <- function (x, t, params, delta.t, ...) {
    v <- 1/(params["v"]/7)
    beta <- params["r0"]*v
    foi <- beta*x["I"]/(x["S"]+x["I"]+x["R"])
        
    x[c("S","I","R","cases")]+c(
        -x["S"]*foi,
        x["S"]*foi - v*x["I"],
        v*x["I"],
        x["S"]*foi
        )*delta.t
}


data <- read.csv('node_modules/ssm-tutorial-data/data/data.csv', header=TRUE, na.strings='null')[,2]
#remove NA
times <- (1:length(data))[!is.na(data)]
data <- data[!is.na(data)]


##!!! change step.fun to sir.ode or sir.psr
sir <- pomp(
    data=data.frame(
        time= times,
        reports=data
        ),
    times="time",
    t0=0,
    rprocess=euler.sim(
        step.fun=sir.ode, ##sir.ode or sir.psr
        delta.t=1/14 #time is in weeks
        ),
    parameter.transform=function(params,...){        

        ##TODO inv_logit_ab for r0 to compare with SSM
        params[c('r0', 'v')] <- exp(params[c('r0', 'v')])
        
        params
    },
    parameter.inv.transform=function(params,...){        

        ##TODO logit_ab for r0 to compare with SSM
        params[c('r0', 'v')] <- log(params[c('r0', 'v')])       
       
        params
    },
    rmeasure =  function(x,t,params,...){
        yobs = rnorm(n=1, mean=params['rep']*x['cases'], sd=sqrt((1-params['rep'])*params['rep']*x['cases']))
        if(yobs>0){
            return(yobs)
        } else {
            return(0.0)
        }        
    },
    dmeasure = function(y,x,t,params,log,...){        
        my.mean=params['rep']*x['cases']
        my.sd=sqrt((1-params['rep'])*params['rep']*x['cases'])
        if(y>0.0){
            p = pnorm(y+0.5, mean = my.mean, sd = my.sd, log.p = FALSE)-pnorm(y-0.5, mean = my.mean, sd = my.sd, log.p = FALSE)
        } else{
            p = pnorm(y+0.5, mean = my.mean, sd = my.sd, log.p = FALSE)
        }
        
        if(log){
            return(log(p))
        } else {
            return(p)
        }        
    },
    
    zeronames=c("cases"), # 'cases' is an accumulator variable
    ic.pars=c("S0","I0","R0"), # initial condition parameters
    comp.names=c("S","I","R") # names of the compartments
    )


theta.guess <- c(
    r0=25,
    v=11,
    rep=0.5,
    S.0=70000,
    I.0=10,
    R.0=929990,
    cases.0=0
    )
    
rw.sd <- c(
    r0=0.02,
    v=0.02
    )/sqrt(length(data))

estpars <- c('r0', 'v')

simulate(
    sir,
    params=theta.guess,
    seed=677573454L
    ) -> simul

plot(simul)

pf <- pfilter(sir, params=theta.guess, Np=100, max.fail=1000, tol=1e-17)
logLik(pf)

mif <- mif(
    tol=1e-17,
    sir,
    Nmif=100,
    start=theta.guess,
    transform=TRUE,
    pars=estpars,
    rw.sd=rw.sd,
    Np=1000,
    var.factor=2,
    ic.lag=39,
    method='mif',
    cooling.factor=0.975,
    max.fail=2000,
    verbose=FALSE,
    seed=677573454L
    )

logLik(mif)

res <- conv.rec(mif)

layout(matrix(1:6,3,2))
for (p in estpars){
    plot(exp(res[,p]), type='l', ylab=p)
}

plot(res[,'loglik'], type='l', ylab='loglik')

#interactive plot (good for diagnosis) #plot(mif)

