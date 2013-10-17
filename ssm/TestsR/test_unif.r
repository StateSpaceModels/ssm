
require(coda)



args <- commandArgs()

varname = args[3]
print(varname)

d<-read.csv("trace_0.csv",header=TRUE)

ess = round(effectiveSize(d[,eval(varname)]))

step = floor(length(d[,eval(varname)])/ess)
x = d$r0_paris[seq(1,length(d[,eval(varname)]),step)]
x2 = hist(x,9.5+qunif(seq(0,1,0.1)),plot=FALSE)$counts
t= chisq.test(x2)

print("here")
print(t$p.value) 
if(t$p.value <0.001){
	cat("0",file="outfile.txt")
} else {
	cat("1",file="outfile.txt")
}
