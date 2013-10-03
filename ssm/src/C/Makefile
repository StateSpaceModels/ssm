CC=clang -ferror-limit=2
CFLAGS= -std=gnu99 -Wall -O3 -DGSL_RANGE_CHECK_OFF -I kalman -I pmcmc -I simulation -I mif -I simplex -I core
LIB=libssm.a libssm_smc.a libssm_simplex.a libssm_mif.a libssm_pmcmc.a libssm_kalman.a libssm_ksimplex.a libssm_kmcmc.a libssm_simul.a libssm_worker.a
ALL_SRC= $(wildcard */*.c)
ALL_SRC_NO_TEMPLATE=$(filter-out templates/input_template.c templates/transform_template.c templates/check_IC_template.c templates/iterator_template.c templates/observed_template.c templates/diff_template.c templates/ode_sde_template.c templates/psr_template.c templates/jac_template.c templates/Ht_template.c templates/Q_template.c templates/step_ekf_template.c, $(ALL_SRC))
SRC=$(filter-out smc/main_smc.c simplex/main_simplex.c mif/main_mif.c worker/main_worker.c pmcmc/main_pmcmc.c kalman/main_kalman.c kalman/main_kmcmc.c kalman/main_ksimplex.c simulation/main_simul.c, $(ALL_SRC_NO_TEMPLATE))
INCLUDES=$(wildcard */*.h)
OBJ= $(SRC:.c=.o)
ALL_OBJ= $(ALL_SRC_NO_TEMPLATE:.c=.o)

PREFIX=/usr/local

all: $(LIB)

libssm.a: $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) -o $@ -c $<

libssm_smc.a: smc/main_smc.o
	ar -rcs $@ $^

libssm_simplex.a: simplex/main_simplex.o
	ar -rcs $@ $^

libssm_mif.a: mif/main_mif.o
	ar -rcs $@ $^

libssm_pmcmc.a: pmcmc/main_pmcmc.o
	ar -rcs $@ $^

libssm_kalman.a: kalman/main_kalman.o
	ar -rcs $@ $^

libssm_ksimplex.a: kalman/main_ksimplex.o
	ar -rcs $@ $^

libssm_kmcmc.a: kalman/main_kmcmc.o
	ar -rcs $@ $^

libssm_simul.a: simulation/main_simul.o
	ar -rcs $@ $^

libssm_worker.a: worker/main_worker.o
	ar -rcs $@ $^

.PHONY: clean uninstall install

install:
	test -d $(PREFIX)/lib || mkdir $(PREFIX)/lib; mv $(LIB) $(PREFIX)/lib/; \
	test -d $(PREFIX)/include || mkdir $(PREFIX)/include; cp $(INCLUDES) $(PREFIX)/include/

clean:
	rm $(ALL_OBJ)

uninstall:
	rm $(PREFIX)/{lib/libssm.a,lib/libssm_smc.a,lib/libssm_simplex.a,lib/libssm_mif.a,lib/libssm_pmcmc.a,lib/libssm_kalman.a,lib/libssm_ksimplex.a,lib/libssm_kmcmc.a,lib/libssm_simul.a,lib/libssm_worker.a,include/ssm.h,include/mcmc_util.h,include/kalman.h,include/simplex.h,include/simulation.h,include/mif.h,include/pmcmc.h}
