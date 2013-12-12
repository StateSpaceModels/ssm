CC=clang

CURRENT_MAKEFILE= $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
TESTS_PATH= $(dir $(abspath $(CURRENT_MAKEFILE)))
CLAR_FIXTURE_PATH= $(TESTS_PATH)resources/

DPKG_ROOT= ../examples/foo/
MODEL_PATH= ../examples/foo/bin

CFLAGS=-g -I.. -I. -I $(HOME)/.ssm/include -Wall -DCLAR_FIXTURE_PATH=\"$(CLAR_FIXTURE_PATH)\"
LDFLAGS=-L$(MODEL_PATH)/C/templates -L $(HOME)/.ssm/lib -lssm -lssmtpl -lssm -lssmtpl -lm -lgsl -lgslcblas -ljansson -lzmq

.PHONY: clean test

# list the objects that go into our test
objects = main.o parameters.o states.o observed.o iterators.o nav.o inputs.o data.o fitness.o calc.o

# build the test executable itself
ssmtest: $(objects) clar.h clar.suite clar.c fixture_data
	$(CC) $(CFLAGS) -o $@ clar.c $(objects) $(LDFLAGS)

# test object files depend on clar macros
$(objects) : clar.h

# build the clar.suite file of test metadata
clar.suite:
	python generate.py .

fixture_data:
	test -d $(CLAR_FIXTURE_PATH) || mkdir -p $(CLAR_FIXTURE_PATH); \
	../bin/ssm install ../examples/foo/package.json -v --src; \
	cp $(MODEL_PATH)/.data.json $(CLAR_FIXTURE_PATH); \
	cp $(DPKG_ROOT)/package.json $(CLAR_FIXTURE_PATH)

# remove all generated files
clean:
	$(RM) -rf *.o clar.suite .clarcache ssmtest ssmtest.dSYM $(CLAR_FIXTURE_PATH)*.json $(MODEL_PATH)

test: ssmtest
	./ssmtest
