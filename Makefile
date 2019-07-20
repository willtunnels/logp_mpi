# (c) copyright 1999-2003 by Vrije Universiteit, Amsterdam, The Netherlands.
# For full copyright and restrictions on use see the file COPYRIGHT.

#
# Makefile for LogP/MPI, Version 1.3
#

# Change this to your C compiler
MPICC		= mpicc

# Make this the flags you normally pass to your compilers
CFLAGS		= -O3 -Wall

# Make this the flags you normally pass to your MPI link command (maybe empty)
LDOPTS		=

# Linker flags for the logp library
LDFLAGS		= -L. -llogp

# Define RANLIB as 'echo' if you do not need ranlib
#RANLIB		= echo
RANLIB		= ranlib

############################################################################
# You should not need to change anything below this line
############################################################################

LOGP_SRC	= logp_mpi.c logp_stats.c
LOGP_INC	= logp_mpi.h logp_stats.h
LOGP_OBJ	= $(LOGP_SRC:%.c=%.o)

LOGP_LIB	= liblogp.a

LOGP_TEST_SRC	= logp_test.c
LOGP_TEST_BIN	= logp_test
LOGP_TEST_BIN_LOG = logp_test_mpelog

all: $(LOGP_LIB) $(LOGP_TEST_BIN) # $(LOGP_TEST_BIN_LOG)

clean:
	rm -f *.o

tidy: clean
	rm -f *~
	rm -f $(LOGP_LIB)
	rm -f $(LOGP_TEST_BIN)
	rm -f logp_test.out.*
	rm -f *core

$(LOGP_LIB): $(LOGP_OBJ)
	ar -cr $(LOGP_LIB) $(LOGP_OBJ)
	$(RANLIB) $(LOGP_LIB)

.c.o: $(LOGP_INC)
	$(MPICC) $(CFLAGS) -c $*.c

$(LOGP_TEST_BIN): $(LOGP_TEST_SRC) $(LOGP_LIB) $(LOGP_INC)
	$(MPICC) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDOPTS)

$(LOGP_TEST_BIN_LOG): $(LOGP_TEST_SRC) $(LOGP_LIB) $(LOGP_INC)
	$(MPICC) $(CFLAGS) -mpilog -o $@ $< $(LDFLAGS) $(LDOPTS)
