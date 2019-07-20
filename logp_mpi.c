/*
 * (c) copyright 1999-2003, by the Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include "mpi.h"
#include "logp_mpi.h"
#include "logp_stats.h"

/*
 * logp_mpi.c
 *
 * Interface to the network performance according to the 
 * "parametrized LogP" model, using MPI.
 *
 * The approach described in the paper is still followed, except
 * that we now use confidence interval estimates to determine the
 * stopping criteria.
 */

/* Some search values are made dependent on msg size to be able
 * to speed things up a little.  They're just heuristics.
 */
#define LOGP_SIZE_SMALL(s)	((s) <= 1024)
#define LOGP_SIZE_AVG(s)	(!LOGP_SIZE_SMALL(s) && ((s) <= (64 * 1024)))
#define LOGP_SIZE_BIG(s)	(!LOGP_SIZE_SMALL(s) && !LOGP_SIZE_AVG(s))

/* Extra sanity checks during development: */
/* #define DEBUG_TIMINGS */

#ifdef DEBUG_TIMINGS
#define debug_timings(list) printf list
#else
#define debug_timings(list) do {} while(0)
#endif

/* perform more tests, depending on msg size: */
#define LOGP_ITS_MULT_SMALL	4
#define LOGP_ITS_MULT_AVG	2
#define LOGP_ITS_MULT_MAX	16

#define LOGP_MAX_BUF_EXT	8 /* maxsize extension for g minimum search */

#define LOGP_T_EPS		(1e-8) /* smallest measurable time */

typedef unsigned char byte_type;

#define sec_to_usec(s)	(1E6 * (s))

static void
spin_wait(unsigned long usec)
{
    double now, deadline;

    now = MPI_Wtime();
    deadline = now + ((double)usec) / 1E6;
    while (MPI_Wtime() < deadline)
	;
}

#define TIMING_LOOP 1000

static void
busy(void)
{
    volatile int dummy;
    volatile int i; /* Huh?!: gcc optimization bug/oddity otherwise?! */

    dummy = 42;
    for (i = 0; i < 300; i++) {
	dummy++;
    }
}

static double
measure_clock_ovhd(int diag, double confint)
{
    double start, stop, loop0, loop1;
    double t0[TIMING_LOOP], t1[TIMING_LOOP];
    statistics loop0_stats, loop1_stats;
    int j;

    stats_init(&loop0_stats);
    stats_init(&loop1_stats);

    for (j = 0; j < 5; j++) {
	int i;

	/* first a loop without clock calls */
	start = MPI_Wtime();
	for (i = 0; i < TIMING_LOOP; i++) {
	    busy();
	}
	stop = MPI_Wtime();
	loop0 = (stop - start) / TIMING_LOOP;
	stats_update(&loop0_stats, loop0);

	/* next, the same loop measured using two clock calls */
	for (i = 0; i < TIMING_LOOP; i++) {
	    t0[i] = MPI_Wtime();
	    busy();
	    t1[i] = MPI_Wtime();
	}
	loop1 = 0.0;
	for (i = 0; i < TIMING_LOOP; i++) {
	    loop1 += t1[i] - t0[i];
	}
	loop1 /= TIMING_LOOP;
	stats_update(&loop1_stats, loop1);

	debug_timings(("Busyloop takes %.3f usec,"
		       "measured with timing calls %.3f usec\n",
		       sec_to_usec(loop0), sec_to_usec(loop1)));
	if (diag >= LOGP_DIAG_VERBOSE2) {
	    printf("clock: loop0 %.3f usec, loop1 %.3f usec\n",
		   sec_to_usec(loop0), sec_to_usec(loop1));
	}
    }

    /* The clock overhead is used as a correction on other measurements,
     * so we take the minima as not to overestimate it.
     */
    loop0 = stats_min(&loop0_stats);
    loop1 = stats_min(&loop1_stats);

    if (diag >= LOGP_DIAG_VERBOSE2) {
	printf("clock: loop0 avg %.3f +- %.3f us; loop1 avg %.3f +- %.3f us\n",
	       sec_to_usec(loop0),
	       sec_to_usec(stats_conf_int(&loop0_stats, 1, j, confint)),
	       sec_to_usec(loop1),
	       sec_to_usec(stats_conf_int(&loop1_stats, 1, j, confint)));
    }

    return loop1 - loop0;
}

/* File format: */
#define LOGP_FILE_HEADER    "# LogP network performance data: %s\n"
#define LOGP_FILE_LATENCY_W "# Latency = %9.7f\n"
#define LOGP_FILE_LATENCY_R "# Latency = %lf\n"
#define LOGP_FILE_TABLE	    "# time      bytes os        os_min    os_cnfint or        or_min    or_cnfint g\n"
#define LOGP_FILE_ENTRY_W   "%9ld %7d %9.7f %9.7f %9.7f %9.7f %9.7f %9.7f %9.7f\n"
#define LOGP_FILE_ENTRY_R   "%ld %d %lf %lf %lf %lf %lf %lf %lf\n"
#define LOGP_FILE_ITEMS	    9


int
logp_save(logp_params *logp, char *filename)
{
    FILE *file;
    int i;

    file = fopen(filename, "w");
    if (file == NULL) {
	fprintf(stderr, "LogP: cannot open '%s' for writing\n",	filename);
	return MPI_ERR_INTERN;
    }

    fprintf(file, LOGP_FILE_HEADER,    logp->descr);
    /* L in sec like the others */
    fprintf(file, LOGP_FILE_LATENCY_W, logp->L);
    fprintf(file, LOGP_FILE_TABLE);
    for (i = 0; i < logp->nsizes; i++) {
	logp_size_params *lsp = &logp->params[i];
	fprintf(file, LOGP_FILE_ENTRY_W,
		lsp->time, lsp->size, lsp->os, lsp->os_min, lsp->os_cnfint,
		lsp->or, lsp->or_min, lsp->or_cnfint, lsp->g);    
    }

    fclose(file);
    return MPI_SUCCESS;
}

static void
read_err(char *filename, int line)
{
    fprintf(stderr, "LogP: reading file '%s' failed in line %d\n",
	    filename, line);
}

static void
alloc_err(void)
{
    fprintf(stderr, "LogP: out of memory\n");
}

#define LINESIZE 128

void
logp_free(logp_params *logp)
{
    if (logp == NULL) {
	return;
    }

    if (logp->descr != NULL) free(logp->descr);
    if (logp->params != NULL) free(logp->params);

    free(logp);
}

int
logp_load(logp_params **logpp, char *filename)
{
    FILE *file;
    int i, line, items;
    char linebuffer[LINESIZE];
    char tempdescr[LINESIZE];
    logp_params *logp;

    file = fopen(filename, "r");
    if (file == NULL) {
	fprintf(stderr, "LogP: cannot open '%s' for reading\n",	filename);
	return MPI_ERR_INTERN;
    }

    logp = NULL;

#   define got_read_err(file, filename, line, logp) {	\
	read_err(filename, line);			\
	fclose(file);					\
	logp_free(logp);				\
    }
#   define got_alloc_err(file, logp) {			\
	alloc_err();					\
	fclose(file);					\
	logp_free(logp);				\
    }

    logp = (logp_params *) malloc(sizeof(logp_params));
    if (logp == NULL) {
	got_alloc_err(file, logp);
	return MPI_ERR_INTERN;
    }

    logp->descr  = NULL;
    logp->nsizes = 0;
    logp->L      = 0.0;
    logp->params = NULL;

    line = 0;

    /* ignore first line: */
    line++;
    if (fgets(linebuffer, LINESIZE, file) == NULL) {
	got_read_err(file, filename, line, logp);
	return MPI_ERR_INTERN;
    }
    items = sscanf(linebuffer, LOGP_FILE_HEADER, &tempdescr[0]);
    if (items != 1) {
	got_read_err(file, filename, line, logp);
	return MPI_ERR_INTERN;
    }
    logp->descr = strdup(tempdescr);
    if (logp->descr == NULL) {
	got_alloc_err(file, logp);
	return MPI_ERR_INTERN;
    }

    line++;
    if (fgets(linebuffer, LINESIZE, file) == NULL) {
	got_read_err(file, filename, line, logp);
	return MPI_ERR_INTERN;
    }
    items = sscanf(linebuffer, LOGP_FILE_LATENCY_R, &logp->L);
    if (items != 1) {
	got_read_err(file, filename, line, logp);
	return MPI_ERR_INTERN;
    }
  
    /* ignore third line: */
    line++;
    if (fgets(linebuffer, LINESIZE, file) == NULL) {
	got_read_err(file, filename, line, logp);
	return MPI_ERR_INTERN;
    }

    for (i = 0; ; i++) {
	logp_size_params *lsp;
	int size;

	line++;
	if (fgets(linebuffer, LINESIZE, file) == NULL) {
	    /* EOF: no more lines */
	    break;
	}

	/* extend size_params in logp structure */
	size = i + 1;
	logp->params = (logp_size_params *)
	    realloc(logp->params, size * sizeof(logp_size_params));
	if (logp->params == NULL) {
	    got_alloc_err(file, logp);
	    return MPI_ERR_INTERN;
	}
	logp->nsizes_max = size;
	
	lsp = &logp->params[i];
	items = sscanf(linebuffer, LOGP_FILE_ENTRY_R,
		       &lsp->time, &lsp->size,
		       &lsp->os, &lsp->os_cnfint, &lsp->os_min,
		       &lsp->or, &lsp->or_cnfint, &lsp->or_min, &lsp->g);
	if (items != LOGP_FILE_ITEMS) {
	    got_read_err(file, filename, line, logp);
	    return MPI_ERR_INTERN;
	}
	/* RTT times currently not saved/restored: */
	lsp->RTT = 0.0; 
	lsp->RTT_min = 0.0; 
    }
    logp->nsizes = i;

    fclose(file);

#if 0
    /* Testing: write back to see if was read properly */
    {
	char copyfile[256];
	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	size = snprintf(copyfile, sizeof(copyfile),
			"%s.%d", filename, rank);
	if (size > 0) {
	    logp_save(logp, copyfile);
	}
    }
#endif
    
    *logpp = logp;
    return MPI_SUCCESS;
}

static int
sort_times(const void *v1, const void *v2)
{
    double val1, val2;

    val1 = * (double *) v1;
    val2 = * (double *) v2;

    if (val1 < val2) {
	return -1;
    } else if (val1 > val2) {
	return 1;
    } else {
	return 0;
    }
}


static void
logp_measure_or_os(logp_params *logp,
		   MPI_Comm comm, int rank, int peer,
		   int incr_iters, int max_iters, double confint, double eps,
		   double clock_ovhd,
		   logp_sendmode send_mode, logp_recvmode recv_mode,
		   byte_type *buffer,
		   double *times, double *times_os, int size_i, int diag)
{
    int i, iters, busy;
    double orig_start, start, stop, stop_os;
    double clock_err;
    long delta_theo, delta_theo_min;
    MPI_Status mpi_stat;
    MPI_Request mpi_req;
    logp_size_params *lsp;
    int cur_iters, cur_iters_max;
    int size;
    int ignore_first_runs;
    statistics rtt_stats, os_stats, or_stats;
    double rtt_conf_int, os_conf_int, or_conf_int;

    stats_init(&rtt_stats);
    stats_init(&os_stats);

    /* or/os: delta round trips with varying message sizes
     * we need two rounds here:
     * first: round trips n bytes forward, 0 bytes back
     *        this is also used to determine os
     * then:  same roundtrips 0 bytes forward, n bytes back with a delta
     * Finally we can compute L or g(msgsize) based using RTT's
     * measured and g(0).
     *
     * Assumption: logp->g[0] is already measured.
     */

    lsp = &logp->params[size_i];
    size = lsp->size;

    if (diag >= LOGP_DIAG_DEFAULT) {
	printf("%d: %7d: ", rank, size);
	fflush(stdout);
    }
    if (diag >= LOGP_DIAG_VERBOSE2) {
	printf("\n");
	fflush(stdout);
    }

    lsp->time = time((time_t *) NULL);

    /* first do the round trips without delta, store this in RTT */
    busy = 1;
    iters = 0;
    start = stop = stop_os = 0.0; /* just to avoid compiler warnings */
    rtt_conf_int = os_conf_int = or_conf_int = 0.0; /* ditto */

    /* we use twice the measured time; XXX why? */
    clock_err = 2.0 * clock_ovhd;

    orig_start = MPI_Wtime();

    /* For "small" sizes (<= 1K) increase the amount of tests.
     * This reduces the startup effect per round.
     */
    if (LOGP_SIZE_BIG(size)) {
	cur_iters_max = incr_iters;
    } else {
	cur_iters_max = incr_iters * LOGP_ITS_MULT_MAX;
	if (LOGP_SIZE_SMALL(size)) {
	    max_iters *= LOGP_ITS_MULT_SMALL;
	} else {
	    max_iters *= LOGP_ITS_MULT_AVG;
	}
    }

    /* We don't want delta to be very close to the clock overhead: */
    delta_theo_min = 100.0E6 * clock_ovhd;

    /* The first few runs are usally highest because of cache interactions,
     * especially for small message sizes.
     * Just ignore them for the measurements.
     */
    ignore_first_runs = !LOGP_SIZE_BIG(size);
    if (ignore_first_runs) {
	cur_iters = 2;
    } else {
	cur_iters = incr_iters;
    }

    while (busy) {
	/* tell mirror the amount of tests we're going to do this round */
	MPI_Send(&cur_iters, 1, MPI_INT, peer, 0, comm);

	/* sync */
	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	
	for (i = 0; i < cur_iters; i++) {
	    switch (send_mode) {
	    case LOGP_MPI_SEND:
		start = MPI_Wtime();
		MPI_Send(buffer, size, MPI_BYTE, peer, 0, comm);
		stop_os = MPI_Wtime();
		break;
	    case LOGP_MPI_ISEND:
		start = MPI_Wtime();
		MPI_Isend(buffer, size, MPI_BYTE, peer, 0, comm, &mpi_req);
		stop_os = MPI_Wtime();
		MPI_Wait(&mpi_req, &mpi_stat);
		break;
	    case LOGP_MPI_SSEND:
		start = MPI_Wtime();
		MPI_Ssend(buffer, size, MPI_BYTE, peer, 0, comm);
		stop_os = MPI_Wtime();
		break;
	    }
	    MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	    stop = MPI_Wtime();

	    /* substract the clock_ovhd cost */
	    times   [i] = stop    - start - clock_ovhd;
	    if (times[i] < LOGP_T_EPS) { /* clock_ovhd overestimated? */
		times[i] = LOGP_T_EPS;
	    }
	    times_os[i] = stop_os - start - clock_ovhd;
	    if (times_os[i] < LOGP_T_EPS) { /* clock_ovhd overestimated? */
		times_os[i] = LOGP_T_EPS;
	    }
	}

	if (ignore_first_runs) {
	    /* this was just the first warming-up run */
	    ignore_first_runs = 0;
	    cur_iters = incr_iters;
	    continue;
	}

	iters += cur_iters;

	debug_timings(("%d: os:", rank));
	for (i = 0; i < cur_iters; i++) {
	    stats_update(&rtt_stats, times[i]);
	    stats_update(&os_stats, times_os[i]);

	    debug_timings(("%.3f ", sec_to_usec(times_os[i])));
	}

	/* I think we should be able to use the onesided tail
	 * since the values we are estimating have a lowerbound.
	 */
	rtt_conf_int = stats_conf_int(&rtt_stats, 1, iters, confint);
	os_conf_int  = stats_conf_int(&os_stats,  1, iters, confint);
	  
	if (diag >= LOGP_DIAG_VERBOSE2) {
	    printf("rtt %9.7f +- %9.7f (%3.1f%%); "
		   "os %9.7f +- %9.7f (%3.1f%%)\n",
		   stats_mean(&rtt_stats), rtt_conf_int,
		   100 * rtt_conf_int / stats_mean(&rtt_stats),
		   stats_mean(&os_stats), os_conf_int,
		   100 * os_conf_int / stats_mean(&os_stats));
	}

	/* TODO: take clock_err into account? */
	busy = (iters < max_iters) &&
	    (((os_conf_int / stats_mean(&os_stats)) > 2 * eps) ||
	     ((rtt_conf_int / stats_mean(&rtt_stats)) > 2 * eps));
    }
    /* tell mirror we're ready with round 1 */
    MPI_Send(&busy, 1, MPI_INT, peer, 0, comm);

    lsp->RTT = stats_mean(&rtt_stats);
    lsp->RTT_min = stats_min(&rtt_stats);
    lsp->os  = stats_mean(&os_stats);
    lsp->os_cnfint = os_conf_int;
    lsp->os_min = stats_min(&os_stats);

    if (diag >= LOGP_DIAG_DEFAULT) {
	printf("os %9.7f", lsp->os);
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf(" [%2d] ", iters);
	} else {
	    printf(" ");	    
	}
	debug_timings((", i=%2d;", iters));
	fflush(stdout);
    }
    if (diag >= LOGP_DIAG_VERBOSE2) {
	printf("\n");
	fflush(stdout);
    }

    /* second do the round trips with delta, store this in or */
    busy = 1;
    iters = 0;
    stats_init(&or_stats);


    /* See comments about first runs above */
    ignore_first_runs = !LOGP_SIZE_BIG(size);
    if (ignore_first_runs) {
	cur_iters = 2;
    } else {
	cur_iters = incr_iters;
    }

    while (busy) {
	/* tell mirror the amount of tests we're going to do this round */
	MPI_Send(&cur_iters, 1, MPI_INT, peer, 0, comm);

	/* sync */
	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);

	/* compute a suitable delta for this data size:
	 * add 20% and translate to usec:
	 */
	delta_theo = sec_to_usec(1.2 * lsp->RTT);
	if (delta_theo < delta_theo_min) {
	    delta_theo = delta_theo_min;
	}
	debug_timings((" delta=%ldus;", delta_theo));

	for (i = 0; i < cur_iters; i++) {
	    switch (send_mode) {
	    case LOGP_MPI_SEND:
		MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);
		break;
	    case LOGP_MPI_ISEND:
		MPI_Isend(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_req);
		MPI_Wait(&mpi_req, &mpi_stat);
		break;
	    case LOGP_MPI_SSEND:
		MPI_Ssend(buffer, 0, MPI_BYTE, peer, 0, comm);
		break;
	    }

	    switch (recv_mode) {
	    case LOGP_MPI_RECV:
		spin_wait(delta_theo);
		start = MPI_Wtime();
		MPI_Recv(buffer, size, MPI_BYTE, peer, 0, comm, &mpi_stat);
		stop = MPI_Wtime();
		break;
	    case LOGP_MPI_IRECV:
		MPI_Irecv(buffer, size, MPI_BYTE, peer, 0, comm, &mpi_req);
		spin_wait(delta_theo);
		start = MPI_Wtime();
		MPI_Wait(&mpi_req, &mpi_stat);
		stop = MPI_Wtime();
		break;
	    }

	    /* substract the clock_ovhd cost */
	    times[i] = stop - start - clock_ovhd;
	    if (times[i] < LOGP_T_EPS) {
		times[i] = LOGP_T_EPS;
	    }
	}

	if (ignore_first_runs) {
	    /* this was just the first warming-up run */
	    ignore_first_runs = 0;
	    cur_iters = incr_iters;
	    continue;
	}

	iters += cur_iters;

	debug_timings(("%d: or:", rank));
	for (i = 0; i < cur_iters; i++) {
	    stats_update(&or_stats, times[i]);

	    debug_timings(("%.3f ", sec_to_usec(times[i])));
	}

	or_conf_int  = stats_conf_int(&or_stats,  1, iters, confint);
	  
	if (diag >= LOGP_DIAG_VERBOSE2) {
	    printf("or %9.7f +- %9.7f (%3.1f%%)\n",
		   stats_mean(&or_stats), or_conf_int,
		   100.0 * or_conf_int / stats_mean(&or_stats));
	}

	busy = (iters < max_iters) &&
	    ((or_conf_int / stats_mean(&or_stats)) > 2 * eps);
    }
    /* tell mirror we're ready with round 2 */
    MPI_Send(&busy, 1, MPI_INT, peer, 0, comm);

    lsp->or = stats_mean(&or_stats);
    lsp->or_cnfint = or_conf_int;
    lsp->or_min = stats_min(&or_stats);

    if (diag >= LOGP_DIAG_DEFAULT) {
	printf("or %9.7f", lsp->or);
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf(" [%2d] ", iters);
	} else {
	    printf(" ");
	}
	debug_timings(("(min %.3f us, max %.3f us)",
		       sec_to_usec(stats_min(&or_stats)),
		       sec_to_usec(stats_max(&or_stats))));
	debug_timings((", i=%2d, %.3fs;", iters, MPI_Wtime() - orig_start));
    }

    if (size != 0) {
	/* set g(m) */
	logp_size_params *lsp0 = &logp->params[0];

	if (lsp->RTT > lsp0->RTT) {
	    /* usually */
	    lsp->g = lsp->RTT - lsp0->RTT     + lsp0->g;
	} else {
	    /* Maybe there was a lot of variability while estimating RTT(0);
	     * just take the lowest one measured.
	     */
	    lsp->g = lsp->RTT - lsp0->RTT_min + lsp0->g;
	}
	if (lsp->g <= 0) {
	    printf("cannot determine g(%d): RTT(%d) %9.7f; RTT(0) %9.7f\n",
		   size, size, lsp->RTT, lsp0->RTT);
	    lsp->g = LOGP_T_EPS;
	}
    }

    if (diag >= LOGP_DIAG_DEFAULT) {
	printf(" g %9.7f", lsp->g);
    }

    /* Check os/or/g discrepancies, due to measurement imprecision.
     * For example, direct measurement of g may be smaller than
     * suggested by derived os/or because pipelining messages
     * gives different caching behaviour.  This is hard to avoid..
     */
    if (lsp->os > lsp->g) {
	if (lsp->os > (1.0 + eps) * lsp->g) {
	    if (diag >= LOGP_DIAG_VERBOSE) {
		printf(" os > g!");
	    }
#ifdef FIXUP_G
	    lsp->g = lsp->os;
#endif
	}
    }
    if (lsp->or > lsp->g) {
	if (lsp->or > (1.0 + eps) * lsp->g) {
	    if (diag >= LOGP_DIAG_VERBOSE) {
		printf(" or > g!");
	    }
#ifdef FIXUP_G
	    lsp->g = lsp->or;
#endif
	}
    }

    if (size == 0) {
	/* set L */
	logp_size_params *lsp0;

	lsp0 = &logp->params[0];
	/* L in sec like the others */
	logp->L = (lsp0->RTT / 2.0 - lsp0->g);
	if (logp->L <= 0.0) {
	    logp->L = 0.0;
	    if (diag >= LOGP_DIAG_DEFAULT) {
		printf(" L unknown!");
	    }
	} else {
	    if (diag >= LOGP_DIAG_DEFAULT) {
		printf(" L = %9.2f", logp->L);
	    }
	}
    }

    if (diag >= LOGP_DIAG_DEFAULT) {
	printf("\n");
    }
}


static void
logp_mirror_or_os(MPI_Comm comm, int rank, int peer,
		  byte_type *buffer, int size)
{
    MPI_Status mpi_stat;
    int incr_iters;
    int i;

    /* or/os, first loop: */
    /* printf("%d: or/os with size %d\n", rank, size); */
    for (;;) {
	MPI_Recv(&incr_iters, 1, MPI_INT, peer, 0, comm, &mpi_stat);
	if (incr_iters == 0) {
	    /* done with first loop */
	    break;
	}

	/* sync */
	MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);

	for (i = 0; i < incr_iters; i++) {
	    MPI_Recv(buffer, size, MPI_BYTE, peer, 0, comm, &mpi_stat);
	    MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);
	}
    }

    /* or/os, second loop:  */
    for (;;) {
	MPI_Recv(&incr_iters, 1, MPI_INT, peer, 0, comm, &mpi_stat);
	if (incr_iters == 0) {
	    /* done with second loop */
	    break;
	}

	/* sync */
	MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);

	for (i = 0; i < incr_iters; i++) {
	    MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	    MPI_Send(buffer, size, MPI_BYTE, peer, 0, comm);
	}
    }
}

#define START_RUNS	10

/* Measuring RTT(m) in g(m) quickly: */
#define GM_RTT_COUNT(size) (LOGP_SIZE_SMALL(size) ? 100 : 10)

static void
logp_measure_gm(logp_params *logp,
		int size_i,
		int size,
		MPI_Comm comm, int rank, int peer,
		int incr_iters, int max_iters, double eps,
		double clock_ovhd, logp_sendmode send_mode,
		logp_bufmode buf_mode, byte_type *buffer, int maxsize,
		double *times, double *times_os, int diag)
{
    int i, runs, minruns, maxruns;
    double start, stop, stop_sync;
    double rtt;
    double tot_time, prev_tot_time;
    double avg, prev_avg;
    MPI_Status mpi_stat;
    MPI_Request mpi_req;
    logp_size_params *lsp;
    int firstrun;

    /*
     * Measure g(m) by means of pipelining messages.
     */
    lsp = &logp->params[size_i];

    /* avoid compiler warnings: */
    start = stop = 0.0;
    avg = prev_avg = 0.0;
    prev_tot_time = 0.0;

    if (diag >= LOGP_DIAG_VERBOSE) {
	printf("%d: g(%d)/pipeline\n", rank, size);
    }

    /* Start with a small number of runs.  We will then check if
     * we actually need more runs.
     * Also limit the amount of iterations; especially needed on WAN links.
     */
    runs = START_RUNS;
    minruns = 2 * START_RUNS;
    maxruns = 128 * minruns;

    MPI_Send(&size, 1, MPI_INT, peer, 0, comm);

    /* First measure RTT(m) quickly. */
    /* sync */
    MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	
    start = MPI_Wtime();
    for (i = 0; i < GM_RTT_COUNT(size); i++) {
	MPI_Send(buffer, size, MPI_BYTE, peer, 0, comm);
	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
    }
    stop = MPI_Wtime();
    rtt = (stop - start) / GM_RTT_COUNT(size);

    firstrun = 1;
    for (;;) {
	/* Notify other side about the next number of sends we'll be doing
	 * this round.  By sending 0 we indicate end of this g(m) test.
	 */
	char *bufp, *endp;

	MPI_Send(&runs, 1, MPI_INT, peer, 0, comm);
	if (runs == 0) {
	    break;
	}

	/* sync */
	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	switch (send_mode) {
	case LOGP_MPI_SEND:
	    start = MPI_Wtime();
	    bufp = buffer;
	    endp = buffer + maxsize;
	    for (i = 0; i < runs; i++) {
		MPI_Send(bufp, size, MPI_BYTE, peer, 0, comm);
		if (buf_mode == LOGP_BUF_DYNAMIC) {
		    /* Sending from different place in the buffer may have
		     * noticable influence on performance due to cache.
		     */
		    bufp += 4 * 1024 + (i % 4) * (1024 + 512);
		    if (bufp + size > endp) {
			bufp = buffer;
		    }
		}
	    }
	    stop = MPI_Wtime();
	    break;
	case LOGP_MPI_ISEND:
	    start = MPI_Wtime();
	    for (i = 0; i < runs; i++) {
		MPI_Isend(buffer, size, MPI_BYTE, peer, 0, comm, &mpi_req);
		MPI_Wait(&mpi_req, &mpi_stat);
	    }
	    stop = MPI_Wtime();
	    break;
	case LOGP_MPI_SSEND:
	    start = MPI_Wtime();
	    for (i = 0; i < runs; i++) {
		MPI_Ssend(buffer, size, MPI_BYTE, peer, 0, comm);
	    }
	    stop = MPI_Wtime();
	    break;
	}

	avg = (stop - start) / runs;

	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	stop_sync = MPI_Wtime();
	tot_time = stop_sync - start;

	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("%d: runs = %d: tot_time %.1f"
		   " incr %.1f%% RTT %.1f g_avg %.1f diff %.1f%%\n",
		   rank, runs, sec_to_usec(tot_time),
		   (prev_tot_time != 0.0) ?
		   100 * ((tot_time - prev_tot_time) / prev_tot_time) : 0.0,
		   sec_to_usec(rtt), sec_to_usec(avg),
		   (prev_avg != 0.0) ?
		   100.0 * ((avg - prev_avg) / prev_avg) : 0.0);
	    fflush(stdout);
	}
	    
	if (firstrun) {
	    firstrun = 0;
	    runs *= 2;
	} else {
	    double perc_diff, maxdiff;

	    perc_diff = 100.0 * ((avg - prev_avg) / prev_avg);

	    debug_timings(("%d: g(%d): %d took %.3f sec:"
			   "avg %9.9f, prev %9.9f, diff %.1f%%\n",
			   rank, size, runs, stop - start, avg,
			   prev_avg, perc_diff));

	    /* Be careful about when to stop:
	     * assume saturation is reached when the larger run gets
	     * (almost) the same per byte overhead as the previous one.
	     * Currently assume this if the difference is less than 1%
	     * (1.5% for "big" sizes, or larger even if things get really
	     * out of hand..).
	     */
	    if (runs < maxruns / 4) {
		maxdiff = (LOGP_SIZE_BIG(size) ? 1.5 : 1.0);
	    } else {
		maxdiff = 8.0;
	    }

	    /*
	     * Stop if:
	     *		(max number of runs has been reached)
	     * OR
	     *		(we've done more than the minimum number of runs
	     *    AND   send overhead per byte is within eps of prev measurement
	     *    AND   sync'd total time is not dominated by latency
	     *    AND   sync'd total time is about two times the previous one)
	     */
	    if ((runs >= maxruns) ||
		((runs >= minruns) &&
		 (fabs(perc_diff) < maxdiff) &&
		 (rtt / tot_time < eps / 2) &&
		 (fabs((2*prev_tot_time - rtt/2) / tot_time) < (1.0 + eps))))
	    {
		/* stop */
		runs = 0;
	    } else {
		/* TODO: doubling may be a bit too much if runs > 1000 */
		runs *= 2;
	    }
	}

	prev_avg = avg;
	prev_tot_time = tot_time;
    }

    if (lsp->g != 0.0) {
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("%d: overriding previous g(%d) %9.7f with %9.7f\n",
		   rank, size, lsp->g, avg);
	    fflush(stdout);
	}
    }
    lsp->g = avg;
}

static void
logp_mirror_gm(MPI_Comm comm, int rank, int peer,
	       int incr_iters, logp_sendmode send_mode,
	       logp_bufmode buf_mode, byte_type *buffer, int maxsize)
{
    int i, runs;
    int size;
    MPI_Status mpi_stat;
    double start, stop;

    /* g: 
     * pipeline g(m)
     */

    /* printf("%d: g(m)/pipeline\n", rank); */

    for (;;) {
	char *bufp, *endp;

	MPI_Recv(&size, 1, MPI_INT, peer, 0, comm, &mpi_stat);
    	if (size < 0) {
	    /* no more g(m) measurements */
	    break;
	}

	/* First measure RTT(m) quickly. */

	/* sync */
	MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);
	
	for (i = 0; i < GM_RTT_COUNT(size); i++) {
	    MPI_Recv(buffer, size, MPI_BYTE, peer, 0, comm, &mpi_stat);
	    MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);
	}

	bufp = buffer;
	endp = buffer + maxsize;

	for (;;) {
	    MPI_Recv(&runs, 1, MPI_INT, peer, 0, comm, &mpi_stat);
	    if (runs == 0) {
		break;
	    }

	    MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);

	    start = MPI_Wtime();
	    for (i = 0; i < runs; i++) {
		MPI_Recv(bufp, size, MPI_BYTE, peer, 0, comm, &mpi_stat);
		if (buf_mode == LOGP_BUF_DYNAMIC) {
		    bufp += 4 * 1024 + (i % 4) * (1024 + 512);
		    if (bufp + size > endp) {
			bufp = buffer;
		    }
		}
	    }
	    stop = MPI_Wtime();
	    MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);

	    debug_timings(("%d: avg recv ovhd for size %d < %.2f usec\n",
			   rank, size, sec_to_usec((stop - start) / runs)));
	}
    }
}

static int
logp_discontinuity(int rank, logp_params *logp, int s,
		   double maxdiff, int diag)
{
    double os_exp, os_slope, os_exp_err, os_maxdiff;
    double or_exp, or_slope, or_exp_err, or_maxdiff;
    double  g_exp,  g_slope,  g_exp_err;
    int size_diff, size_diff_1;
    int retval;
    int did_print;

    if (s < 2) {
	/* not enough previous points */
	return 0;
    }

    retval = 0;

    size_diff_1 = logp->params[s - 1].size - logp->params[s - 2].size;
    os_slope = (logp->params[s - 1].os - logp->params[s - 2].os) / size_diff_1;
    or_slope = (logp->params[s - 1].or - logp->params[s - 2].or) / size_diff_1;
    g_slope  = (logp->params[s - 1].g  - logp->params[s - 2].g ) / size_diff_1;

    /* Use linear extrapolation to determine expected value */
    size_diff = logp->params[s].size - logp->params[s - 1].size;
    os_exp = logp->params[s - 1].os + size_diff * os_slope;
    or_exp = logp->params[s - 1].or + size_diff * or_slope;
    g_exp  = logp->params[s - 1].g  + size_diff *  g_slope;

    os_exp_err = or_exp_err = g_exp_err = 0.0;
    if (os_exp > LOGP_T_EPS) {
	os_exp_err = fabs((os_exp - logp->params[s].os) / os_exp);
    }
    if (or_exp > LOGP_T_EPS) {
	or_exp_err = fabs((or_exp - logp->params[s].or) / or_exp);
    }
    if ( g_exp > LOGP_T_EPS) {
	g_exp_err  = fabs(( g_exp - logp->params[s].g ) /  g_exp);
    }

    /* If os/or is very small compared to the size, some non-synchronous
     * mechanism is presumably kicking in, causing the usual overhead to
     * be hidden.  Increase the maxdiff to avoid many useless remeasurements.
     *
     * TODO: find out a reasonable measure "very small"; it's platform
     * dependent of course..  For us,LOGP_T_EPS works fine.
     */
    or_maxdiff = os_maxdiff = maxdiff;
    if (!LOGP_SIZE_SMALL(logp->params[s].size)) {
	if ((os_exp / logp->params[s].size) < LOGP_T_EPS) {
	    os_maxdiff *= 2;
	}
	if ((or_exp / logp->params[s].size) < LOGP_T_EPS) {
	    or_maxdiff *= 2;
	}
    }

    did_print = 0;
    if (os_exp_err > os_maxdiff) {
	if (diag >= LOGP_DIAG_VERBOSE) {
	    if (did_print++ == 0) {
		printf("%d:+%7d: ", rank, logp->params[s].size);
	    }
	    printf("os %9.7f (%.1f%%) ", logp->params[s].os,
		   100 * (logp->params[s].os - os_exp) / os_exp);
	}
	retval = 1;
    }
    if (or_exp_err > or_maxdiff) {
	if (diag >= LOGP_DIAG_VERBOSE) {
	    if (did_print++ == 0) {
		printf("%d:+%7d: ", rank, logp->params[s].size);
	    }
	    printf("or %9.7f (%.1f%%) ", logp->params[s].or,
		   100 * (logp->params[s].or - or_exp) / or_exp);
	}
	retval = 1;
    }
    if (g_exp_err > maxdiff) {
	if (diag >= LOGP_DIAG_VERBOSE) {
	    if (did_print++ == 0) {
		printf("%d:+%7d: ", rank, logp->params[s].size);
	    }
	    printf("g  %9.7f (%.1f%%)", logp->params[s].g,
		   100 * (logp->params[s].g - g_exp) / g_exp);
	}
	retval = 1;
    }
    if (did_print > 0) {
	printf("\n");
    }

    return retval;
}

static double
logp_g_decreasing(int rank, logp_params *logp, int s, int diag)
{
    double g_exp, g_slope;
    int size_diff, size_diff_1;

    if (s < 2) {
	/* not enough previous points */
	return 0.0;
    }

    size_diff_1 = logp->params[s - 1].size - logp->params[s - 2].size;
    g_slope = (logp->params[s - 1].g - logp->params[s - 2].g) / size_diff_1;

    /* Use linear extrapolation to determine expected value */
    size_diff = logp->params[s].size - logp->params[s - 1].size;
    g_exp = logp->params[s - 1].g + size_diff * g_slope;

    if (logp->params[s].g < g_exp) {
	/* g still decreasing */
	double retval;

	retval = (g_exp - logp->params[s].g) / g_exp;
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("%d: size %d: g expected %9.7f actual %9.7f"
		   ", decrease %.1f%%\n",
		   rank, logp->params[s].size, g_exp,
		   logp->params[s].g, retval * 100.0);
	}
	return retval;
    } else {
	return 0.0;
    }
}

static void
logp_insert(logp_params *logp, int s, int size)
{
    int shift;

    /* Shift all parameters starting at s one position up to
     * make room for a new size.
     * Assumption: the caller has checked that it still fits.
     */
    shift = logp->nsizes - s;
    memmove(&logp->params[s + 1], &logp->params[s],
	    shift * sizeof(logp->params[0]));
    
    /* Add the new size */
    logp->params[s].size = size;
    logp->nsizes++;
}

#define DEF_LOG_PROBES	2.0 /* perform log2(maxmsgsize) probes */

int
logp_measure(logp_params **logpp,
	     char *descr,
	     MPI_Comm comm, int measure, int mirror,
	     int logp_g_zero,
	     int incr_iters, int max_iters,
	     double confint, double eps, double max_diff,
	     int minsize, int maxsize, int max_nsizes,
	     logp_sendmode send_mode, logp_recvmode recv_mode,
	     logp_bufmode buf_mode, unsigned int options, int diag)
{
    int retval = MPI_SUCCESS;
    byte_type *buffer;
    int s, rank, peer;
    double start_gm, stop_gm;
    double clock_ovhd;
    double dsize;
    MPI_Status mpi_stat;
    MPI_Request mpi_req;
    double *times, *times_os;
    logp_params *logp;
    int size, maxsize_extended;

    MPI_Comm_rank(comm, &rank);
    if ((rank != measure) && (rank != mirror)) {
	/* sanity check */
	fprintf(stderr, "LogP: my rank %d is not measure %d or mirror %d\n",
		rank, measure, mirror);
	return MPI_ERR_INTERN;
    }
    /* only the processes that measure and mirror continue here */

    /* take possible extension of maximum size into account: */
    maxsize_extended = maxsize * LOGP_MAX_BUF_EXT;
    if ((buffer = (byte_type *) malloc(maxsize_extended)) == NULL) {
	return MPI_ERR_INTERN;
    }
    if ((times = (double *)
	 malloc(incr_iters * LOGP_ITS_MULT_MAX * sizeof(double))) == NULL) {
 	return MPI_ERR_INTERN;
    }
    if ((times_os = (double *)
	 malloc(incr_iters * LOGP_ITS_MULT_MAX * sizeof(double))) == NULL) {
	return MPI_ERR_INTERN;
    }

    if (max_nsizes < 1) {
	return MPI_ERR_INTERN;
    }

    if (rank == measure) {
	/* Allocate logp structure.
	 * To avoid repeatedly having to reallocate all buffers when
	 * adding tests we allocate a reasonable amount in advance.
	 *
	 * TODO: proper error recovery from failed allocs (e.g.,
	 * notify mirror; free resources we did allocate).
	 */
	logp = (logp_params *) malloc(sizeof(logp_params));
	if (logp == NULL) return MPI_ERR_INTERN;
	logp->descr = strdup(descr);
	if (logp->descr == NULL) return MPI_ERR_INTERN;

	logp->params  = (logp_size_params *)
	    malloc(max_nsizes * sizeof(logp_size_params));
	if (logp->params == NULL) {
	    return MPI_ERR_INTERN;
	}
	logp->nsizes_max = max_nsizes;
	logp->nsizes = 0;

	/*
	 * Construct the initial size table.  We'll insert other
	 * values in case we find discontinuities later on.
	 */
	/* always do size 0 */
	logp->params[logp->nsizes++].size = 0;

	/* Then add values from minsize to maxsize.
	 * We prefer doing a log(maxmsgsize) amount of probes.
	 */
	for (dsize = minsize;
	     logp->nsizes < logp->nsizes_max && dsize < maxsize;
	     dsize *= DEF_LOG_PROBES)
	{
	    if (logp->params[logp->nsizes - 1].size == (int) dsize) {
		/* already added this one */
		continue;
	    }
	    logp->params[logp->nsizes++].size = (int) dsize;
	}

	/* if maxsize is not in yet, add it as well */
	if (logp->params[logp->nsizes - 1].size != maxsize) {
	    if (logp->nsizes < logp->nsizes_max) {
		logp->params[logp->nsizes++].size = maxsize;
	    }
	}

	peer = mirror;
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("%d: measuring with peer %d\n", rank, peer);
	}

	/* check the overhead of reading the clock first */
	clock_ovhd = measure_clock_ovhd(diag, confint);
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("Clock overhead is %.3f usec\n", sec_to_usec(clock_ovhd));
	}

	/* a first ping pong for warmup */
	switch (send_mode) {
	case LOGP_MPI_SEND:
	    MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);
	    break;
	case LOGP_MPI_ISEND:
	    MPI_Isend(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_req);
	    break;
	case LOGP_MPI_SSEND:
	    MPI_Ssend(buffer, 0, MPI_BYTE, peer, 0, comm);
	    break;
	}
	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);

	start_gm = MPI_Wtime();

	logp->params[0].g = 0.0;
	logp_measure_gm(logp,
			0 /* only need g(0) */,
			logp_g_zero, /* usually 0, unless problems with that */
			comm, rank, peer,
			incr_iters, max_iters, eps,
			clock_ovhd, send_mode, buf_mode,
			buffer, maxsize_extended,
			times, times_os, diag);
	/* signal end of gm measurements to mirror: */
	size = -1;
	MPI_Send(&size, 1, MPI_INT, peer, 0, comm);

	stop_gm = MPI_Wtime();
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("%d: Measuring g(0) took %.2f sec\n",
		   rank, stop_gm - start_gm);
	}

	for (s = 0; s < logp->nsizes; s++) {
	    /* tell mirror about the size we want to test: */
	    MPI_Send(&logp->params[s].size, 1, MPI_INT, peer, 0, comm);

	    logp_measure_or_os(logp, comm, rank, peer,
			       incr_iters, max_iters, confint, eps,
			       clock_ovhd, send_mode, recv_mode,
			       buffer, times, times_os, s, diag);
	}
    
	/* see if increasing message size might still help */
	while ((logp->nsizes < logp->nsizes_max) &&
	       (logp->params[logp->nsizes - 1].size < maxsize_extended))
	{
	    double decrease;

	    s = logp->nsizes - 1;
	    decrease = logp_g_decreasing(rank, logp, s, diag);
	    if (decrease > eps) {
		logp_insert(logp, s + 1, 2 * logp->params[s].size);

		/* tell mirror about the size we want to test: */
		MPI_Send(&logp->params[s + 1].size, 1, MPI_INT,
			 peer, 0, comm);

		logp_measure_or_os(logp, comm, rank, peer,
				   incr_iters, max_iters, confint, eps,
				   clock_ovhd, send_mode, recv_mode,
				   buffer, times, times_os, s + 1, diag);
	    } else {
		break;
	    }
	}

	/* now go and look for discontinuities: */
	for (s = 2; s < logp->nsizes; s++) {
	    if (logp_discontinuity(rank, logp, s, max_diff, diag)) {
		/* Simply insert a single extra measuring point
		 * halfway between the previous and the current.
		 */
		double minsizediff;
		double sizediff;

		sizediff = (double) (logp->params[s].size - 
				     logp->params[s - 1].size);

		/* We make the minimum size difference between probes
		 * size dependent to avoid doing too many of them.
		 * This could be changed if it is important to find
		 * the exact place where discontinuities occur.
		 */
		minsizediff = eps * (double) logp->params[s].size;
		if (minsizediff < 32) {
		    minsizediff = 32;
		}

		if ((logp->nsizes < logp->nsizes_max) &&
		    (sizediff > minsizediff))
		{
		    logp_insert(logp, s,
			    ((logp->params[s].size +
			      logp->params[s - 1].size) / 2) & ~0x3);

		    /* tell mirror about the size we want to test: */
		    MPI_Send(&logp->params[s].size, 1, MPI_INT,
			     peer, 0, comm);

		    logp_measure_or_os(logp, comm, rank, peer,
				       incr_iters, max_iters, confint, eps,
				       clock_ovhd, send_mode, recv_mode,
				       buffer, times, times_os, s, diag);
		    s--; /* recheck the new point */
		}
	    }
	}

	/* tell mirror we're done */
	size = -1;
	MPI_Send(&size, 1, MPI_INT, peer, 0, comm);

	if (options & LOGP_OPT_GM_DIRECT) {
	    /* in addition, perform direct measurements of g(m),
	     * overriding original estimates.
	     */
	    start_gm = MPI_Wtime();

	    for (s = 0; s < logp->nsizes; s++) {
		logp_measure_gm(logp,
				s,
				logp->params[s].size,
				comm, rank, peer,
				incr_iters, max_iters, eps,
				clock_ovhd, send_mode, buf_mode,
				buffer, maxsize_extended,
				times, times_os, diag);
	    }

	    /* signal end of gm measurements to mirror: */
	    size = -1;
	    MPI_Send(&size, 1, MPI_INT, peer, 0, comm);

	    stop_gm = MPI_Wtime();
	    if (diag >= LOGP_DIAG_VERBOSE) {
		printf("%d: Measuring g(m) took %.2f sec\n",
		       rank, stop_gm - start_gm);
	    }
	}

	/* return the logp result structure to caller */
	*logpp = logp;

    } else if (rank == mirror) {

	/* the process to communicate with: */
	peer = measure;
	if (diag >= LOGP_DIAG_VERBOSE) {
	    printf("%d: mirroring peer %d\n", rank, peer);
	}

	/* a first ping pong for warmup */
	MPI_Recv(buffer, 0, MPI_BYTE, peer, 0, comm, &mpi_stat);
	MPI_Send(buffer, 0, MPI_BYTE, peer, 0, comm);

	logp_mirror_gm(comm, rank, peer, incr_iters, send_mode, buf_mode,
		       buffer, maxsize_extended);

	do {
	    /* wait for size from peer; -1 signals that we're done */
	    MPI_Recv(&size, 1, MPI_INT, peer, 0, comm, &mpi_stat);
	    if (size >= 0) {
		logp_mirror_or_os(comm, rank, peer, buffer, size);
	    }
	} while (size >= 0);

	if (options & LOGP_OPT_GM_DIRECT) {
	    /* in addition, perform direct measurements of g(m) */
	    logp_mirror_gm(comm, rank, peer, incr_iters, send_mode, buf_mode,
			   buffer, maxsize_extended);
	}

	*logpp = NULL; /* we don't have the results at the mirror */
    }

    free(buffer);
    free(times);
    free(times_os);

    return retval;
}

int
logp_init(void)
{
    int retval = 0;

    return retval;
}

int
logp_finalize(void)
{
    int retval = 0;

    return retval;
}

#define logp_param_at_off(logp, size_i, off) \
	(* (double *) ((char *) (& (logp)->params[(size_i)]) + (off)))

static double
extract(logp_params *logp, int offset, int msgsize)
{
    double val_i, val_i_1;
    double slope;
    int i;

    if (msgsize < 0) msgsize = 0;

    /* Just linear search for now */
    for (i = 0; i < logp->nsizes; i++) {
	if (logp->params[i].size >= msgsize) {
	    break;
	}
    }

    if (i >= logp->nsizes) {
	/* msgsize > largest size recorded: extrapolate from last two points */
	int maxsize_i = logp->nsizes - 1;

	val_i   = logp_param_at_off(logp, maxsize_i,     offset);
	val_i_1 = logp_param_at_off(logp, maxsize_i - 1, offset);

	slope = (val_i - val_i_1) /
	    (logp->params[maxsize_i].size - logp->params[maxsize_i - 1].size);

	return val_i + (msgsize - logp->params[maxsize_i].size) * slope;
    } else if (logp->params[i].size == msgsize) {
	/* we found a measured point */
	return logp_param_at_off(logp, i, offset);
    } else {
	/* otherwise: interpolate between size[i - 1] and size[i] */
	val_i   = logp_param_at_off(logp, i,     offset);
	val_i_1 = logp_param_at_off(logp, i - 1, offset);

	slope = (val_i - val_i_1) /
	    (logp->params[i].size - logp->params[i - 1].size);

	return val_i_1 + (msgsize - logp->params[i - 1].size) * slope;
    }
}

time_t logp_time(logp_params *logp, int message_size)
{
    int i;

    /* Just linear search for now */
    for (i = 0; i < logp->nsizes; i++) {
	if (logp->params[i].size >= message_size) {
	    break;
	}
    }

    /* interpolation doesn't make sense here */
    if (i >= logp->nsizes) {
	return -1;
    } else {
	return logp->params[i].time;
    }
}

double
logp_os(logp_params *logp, int message_size)
{
    return extract(logp, offsetof(logp_size_params, os), message_size);
}

double
logp_or(logp_params *logp, int message_size)
{
    return extract(logp, offsetof(logp_size_params, or), message_size);
}

double
logp_g(logp_params *logp, int message_size)
{
    return extract(logp, offsetof(logp_size_params, g), message_size);
}

double
logp_L(logp_params *logp)
{
    /* L in sec like the others */
    return logp->L;
}

