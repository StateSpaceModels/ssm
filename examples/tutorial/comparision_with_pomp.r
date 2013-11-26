##Special thanks to Aaron King for this template!

require('pomp')

#helper function
invlogitab= function(x, a, b) {
    if (x < 0) {
        return ((b*exp(x)+a)/(1.0+exp(x)));
    } else {
        return ((b+a*exp(-x))/(1.0+exp(-x)) );
    }    
}


model <- "sir_pomp"
pkg <- "pomp"
modelfile <- paste(model,".c",sep="")
solib <- paste(model,.Platform$dynlib.ext,sep="")

cflags <- paste0("PKG_CFLAGS=\"",
                 Sys.getenv("PKG_CFLAGS"),
                 " -I",system.file("include",package="pomp"),"\"")

rv <- system2(
              command=R.home("bin/R"),
              args=c("CMD","SHLIB","-o",solib,modelfile),
              env=cflags
              )


dyn.load(solib) ## load the shared-object library


data <- read.csv('node_modules/ssm-tutorial-data/data/data.csv', header=TRUE, na.strings='null')[,2]
#remove NA
times <- (1:length(data))[!is.na(data)]
data <- data[!is.na(data)]


sir <- pomp(
    data=data.frame(time= times, cases=data),
    times="time",
    t0=0,
    rprocess=euler.sim( step.fun="SIR_stepfn", delta.t=0.25/7 ),
    dmeasure="SIR_dmeasure",
    rmeasure="SIR_rmeasure",
    skeleton.type="vectorfield",
    skeleton="SIR_skelfn",
    parameter.transform="SIR_par_trans",
    parameter.inv.transform="SIR_par_untrans",
    obsnames="cases",
    statenames=c("S","I","R","Inc"),
    paramnames=c("r0","v","rep","S.0","I.0","R.0"), 
    zeronames=c("Inc"),
    comp.names=c("S","I","R"),
    ic.names=c("S.0","I.0","R.0")
    )


theta.guess <- c(
    r0=25,
    v=11,
    rep=0.5,
    S.0=70000,
    I.0=10,
    R.0=929990,
    Inc.0=0
    )
    
rw.sd <- c(
    r0=0.04,
    v=0.02
    )/sqrt(length(data))

simulate(
    sir,
    params=theta.guess,
    seed=677573454L
    ) -> simul

plot(simul)

pf <- pfilter(sir, params=theta.guess, Np=100, max.fail=1000, tol=1e-17)
logLik(pf)


mif <- mif(
    sir,
    tol=1e-17,
    Nmif=100,
    start=theta.guess,
    transform=TRUE,
    pars=c('r0', 'v'),
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


#interactive plot (good for diagnosis)
##plot(mif)


#compare pomp and SSM

layout(matrix(1:6,2,3, byrow=TRUE))
plot(invlogitab(res[,'r0'], 15, 35), type='l', ylab='r0')
plot(exp(res[,'v']), type='l', ylab='v')
plot(res[,'loglik'], type='l', ylab='loglik')


##SSM version

trace <- function(){
    x = read.csv('bin/trace_0.csv')
    plot(x$index, x$r0, type='l')
    plot(x$index, x$pr_v, type='l')
    plot(x$index, x$fitness, type='l')
}

trace()


diag <- function(ind){
    x= read.csv('bin/diag_0.csv')
    layout(matrix(1:6,2,3))

    r0 <- x$r0[x$index==ind]
    for(i in 1:length(r0)){
        r0[i] <- invlogitab(r0[i], 15,35)
    }
       
    plot(as.Date(x$date[x$index==ind]), x$ess[x$index==ind], type='l')
    plot(as.Date(x$date[x$index==ind]), r0, type='l')
    plot(as.Date(x$date[x$index==ind]), x$var_r0[x$index==ind], type='l')
    plot(as.Date(x$date[x$index==ind]), exp(x$pr_v[x$index==ind]), type='l')
    plot(as.Date(x$date[x$index==ind]), x$var_pr_v[x$index==ind], type='l')
}

dyn.unload(solib)
