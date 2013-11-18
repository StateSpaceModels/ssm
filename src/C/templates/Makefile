CC=gcc #clang -ferror-limit=2
CFLAGS=-std=gnu99 -O3 -DGSL_RANGE_CHECK_OFF -I $(HOME)/.ssm/include
LDFLAGS=-lssm -lssmtpl -lssm -lm -lgsl -lgslcblas -ljansson -lzmq
LIB=libssmtpl.a
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)

EXEC=simul simplex smc mif pmcmc kalman ksimplex kmcmc worker

all: $(LIB)

install: $(EXEC)
	mv $(EXEC) ../../

$(LIB): $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

simul: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmsimul $(LDFLAGS)

simplex: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmsimplex $(LDFLAGS)

smc: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmsmc $(LDFLAGS) 

mif: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmmif $(LDFLAGS)  

pmcmc: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmpmcmc $(LDFLAGS)  

kalman: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmkalman $(LDFLAGS) 

ksimplex: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmksimplex $(LDFLAGS)

kmcmc: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmkmcmc $(LDFLAGS)

worker: libssmtpl.a
	$(CC) $(CFLAGS) -L. -L $(HOME)/.ssm/lib -o $@ -lssmworker $(LDFLAGS) 

.PHONY: clean

clean:
	rm *.o $(LIB)

uninstall: clean
	rm ../../{simul,simplex,smc,mif,pmcmc,kalman,ksimplex,kmcmc,worker}
