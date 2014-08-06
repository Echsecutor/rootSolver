# 
# @file Makefile
# @author  Sebastian Schmittner <sebastian@schmittner.pw>
# @version 1.0.2014-06-11
#
# @section Version number format
 #
# The Version number is formatted as "M.S.D" where M is the major
# release branch (backward compatibility to all non-alpha releases of
# the same branch is guaranteed), S is the state of this release (0
# for alpha, 1 for beta, 2 for stable), and D is the date formatted
# as yyyy-mm-dd.)
#
# 
# @section DESCRIPTION
# 
# This Makefile builds several test cases for the multi root solver library. 
# Notice that the actual libraries are header only, nothing needs to be compiled.
# 
# 

export SHELL = /bin/sh
.SUFFIXES:
.SUFFIXES: .cpp .o .out

MAKEFLAGS+= -w -e

export ROOTDIR := $(shell pwd)

export BINDIR = $(ROOTDIR)/bin
export SRCDIR = $(ROOTDIR)/src
export TESTDIR = $(ROOTDIR)/test-out
export TESTCASES = $(ROOTDIR)/test-in


export CC = g++

# notice that "eigen3/..." and mpreal.h need to be included
LOCALINCLUDES=-I$(HOME)/include

DEBUGFLAGS=-D DEBUG=11
#-g
export CFLAGS = -Wall -Werror -ansi -pedantic -std=c++0x -pthread $(LOCALINCLUDES)
#-Wfatal-errors

MPREALLIBS=-lgmpxx -lgmp -lmpfr

export LDFLAGS = $(MPREALLIBS)


export PREFLAGS =
# -D MULTITHREADED=0 -D DEBUG=1

VPATH = $(SRCDIR):$(TESTDIR)

#constructionSite: clean selfconsistencyEquations.out

all: folders testEigen.out testComplex.out testSRS.out testMRS.out testBatch.out testExtraData.out testExtraSolver.out testExtraBatchSolver.out derivativeVerification.out selfconsistencyEquations.out
	@-mv *.dat $(TESTDIR) 2>/dev/null
	@echo 
	@echo "[OK]		All tests Passed! :D"
	@echo 


## production:

SCE-approx: selfconsistencyEquations.cpp
	$(CC) -DUSEEXACTINT=0 -DMULTITHREADED=1 -DDEBUG=3 -O3 $(CFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS)

SCE: selfconsistencyEquations.cpp
	$(CC) -DMULTITHREADED=1 -DDEBUG=3 -O3 $(CFLAGS) -o $(BINDIR)/$@ $^ $(LDFLAGS)



## testcases:

folders:
	@-mkdir $(BINDIR)
	@-mkdir $(TESTDIR)

# Special case: also test the single threaded version.
$(BINDIR)/testBatch: testBatchSolver.cpp
	@echo
	@echo "[...]		Building single threaded $^"
	@echo
	$(CC) -D MULTITHREADED=0 $(CFLAGS) $(DEBUGFLAGS) -o $@-single $^ $(LDFLAGS)
	@echo
	@echo "[...]		Test run for singled threaded batch solver"
	@echo
	$@-single -o $(TESTDIR)/testBatch-single-data > $(TESTDIR)/testBatch-single.out 2>&1
	@echo
	@echo "[OK]		test did not crash ;)"
	@echo
	@echo "[...]		Building multi threaded $^"
	@echo
	$(CC) -D MULTITHREADED=1 $(CFLAGS) $(DEBUGFLAGS) -o $@ $^ $(LDFLAGS)



# The rest should be self explanatory, just build and run all test binaries
# actually comparing results: 2do ;)


$(BINDIR)/%: %.cpp
	@echo
	@echo "[...]		Building $^"
	@echo
	$(CC) $(PREFLAGS) $(CFLAGS) $(DEBUGFLAGS) -o $@ $^ $(LDFLAGS)


selfconsistencyEquations.out: $(BINDIR)/selfconsistencyEquations
	@echo
	@echo "[...]		Test run of $^"
	@echo
	$^ $(TESTCASES)/dummySCE.conf multi-SCE > $(TESTDIR)/$@ 2>&1
	@echo
	@echo "[OK]		$^ multi-threaded test did not crash"
	@echo
	$^ $(TESTCASES)/dummySCE-single.conf single-SCE > $(TESTDIR)/$@-single 2>&1
	@echo
	@echo "[OK]		$^ single-threaded test did not crash"
	@echo "    		$^ you should compare single- with multi-threaded results!"
	@echo



%.out: $(BINDIR)/%
	@echo
	@echo "[...]		Test run of $^"
	@echo
	$^ > $(TESTDIR)/$@ 2>&1
	@echo
	@echo "[OK]		$^ test did not crash ;)"
	@echo


clean:
	@-rm -f $(BINDIR)/* $(TESTDIR)/* *~ $(SRCDIR)/*~ *.o *.dat 2>/dev/null


# actually not used for test cases
%.o : %.cpp
	$(CC) $(PREFLAGS) $(CFLAGS) $(DEBUGFLAGS) -c -o $@ $<

