CC=clang -ferror-limit=2
CFLAGS= -std=gnu99 -Wall -O3 -DGSL_RANGE_CHECK_OFF -I kalman -I pmcmc -I simul -I mif -I simplex -I core
LIB=libssm.a libssmsmc.a libssmsimplex.a libssmmif.a libssmpmcmc.a libssmkalman.a libssmksimplex.a libssmkmcmc.a libssmsimul.a #libssmworker.a
ALL_SRC= $(wildcard */*.c)
ALL_SRC_NO_TEMPLATE=$(filter-out templates/input_template.c templates/transform_template.c templates/check_IC_template.c templates/iterator_template.c templates/observed_template.c templates/diff_template.c templates/ode_sde_template.c templates/psr_template.c templates/jac_template.c templates/Ht_template.c templates/Q_template.c templates/step_ekf_template.c, $(ALL_SRC))
SRC=$(filter-out smc/main_smc.c simplex/main_simplex.c mif/main_mif.c worker/main_worker.c pmcmc/main_pmcmc.c kalman/main_kalman.c kalman/main_kmcmc.c kalman/main_ksimplex.c simul/main_simul.c, $(ALL_SRC_NO_TEMPLATE))
INCLUDES=$(wildcard */*.h)
OBJ= $(SRC:.c=.o)
ALL_OBJ= $(ALL_SRC_NO_TEMPLATE:.c=.o)

PREFIX=/usr/local

all: $(LIB)

libssm.a: $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) -o $@ -c $<

libssmsmc.a: smc/main_smc.o
	ar -rcs $@ $^

libssmsimplex.a: simplex/main_simplex.o
	ar -rcs $@ $^

libssmmif.a: mif/main_mif.o
	ar -rcs $@ $^

libssmpmcmc.a: pmcmc/main_pmcmc.o
	ar -rcs $@ $^

libssmkalman.a: kalman/main_kalman.o
	ar -rcs $@ $^

libssmksimplex.a: kalman/main_ksimplex.o
	ar -rcs $@ $^

libssmkmcmc.a: kalman/main_kmcmc.o
	ar -rcs $@ $^

libssmsimul.a: simul/main_simul.o
	ar -rcs $@ $^

#libssmworker.a: worker/main_worker.o
#	ar -rcs $@ $^

.PHONY: clean uninstall install

install:
	test -d $(PREFIX)/lib || mkdir $(PREFIX)/lib; mv $(LIB) $(PREFIX)/lib/; \
	test -d $(PREFIX)/include || mkdir $(PREFIX)/include; cp $(INCLUDES) $(PREFIX)/include/

clean:
	rm $(ALL_OBJ); \
	rm $(LIB)

uninstall:
	rm $(PREFIX)/{lib/libssm.a,lib/libssmsmc.a,lib/libssmsimplex.a,lib/libssmmif.a,lib/libssmpmcmc.a,lib/libssmkalman.a,lib/libssmksimplex.a,lib/libssmkmcmc.a,lib/libssmsimul.a,lib/libssmworker.a,include/ssm.h}
