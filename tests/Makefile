CC=gcc

CURRENT_MAKEFILE  := $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))
TEST_DIRECTORY    := $(abspath $(CURRENT_MAKEFILE))
CLAR_PATH         := $(dir $(TEST_DIRECTORY))
CLAR_FIXTURE_PATH := $(TEST_DIRECTORY)/resources/

CFLAGS=-g -I.. -I. -Wall -DCLAR_FIXTURE_PATH=\"$(CLAR_FIXTURE_PATH)\"

.PHONY: clean test

# list the objects that go into our test
objects = main.o adding.o

# build the test executable itself
ssmtest: $(objects) clar.h clar.suite $(CLAR_PATH)clar.c
	$(CC) $(CFLAGS) -o $@ "$(CLAR_PATH)clar.c" $(objects)

# test object files depend on clar macros
$(objects) : $(CLAR_PATH)clar.h

# build the clar.suite file of test metadata
clar.suite:
	python "$(CLAR_PATH)generate.py" .

# remove all generated files
clean:
	$(RM) -rf *.o clar.suite .clarcache ssmtest ssmtest.dSYM

test: ssmtest
	./ssmtest
