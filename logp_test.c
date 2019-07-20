/*
 * (c) copyright 1999-2003 by Vrije Universiteit, Amsterdam, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include "logp_mpi.h"

static char *prog_name;

static void
usage(void)
{
    fprintf(stderr, "Usage: %s [options]\n", prog_name);
    fprintf(stderr, "options:\n"
    "-flood:         direct measurement of g(m) using flooding\n"
    "-silent:        run silently\n"
    "-verbose:       print details about results of the measurements\n"
    "-static:        send/receive from fixed place in buffer\n"
    "-dynamic:       send/receive from variable place in buffer\n"
    "-min-its <its>: set minimum number of iterations (default %d)\n"
    "-max-its <its>: set maximum number of iterations (default %d)\n"
    "-eps <eps>:     set estimated precision value epsilon (default %.2f)\n"
    "-conf-int <cf>: set confidence interval (default %.2f)\n"
    "-max-diff <df>: set non-linearity difference criterion (default %.2f)\n"
    "-min-size <sz>: set minimum non-zero size to be tested (default %d)\n"
    "-max-size <sz>: set initial maximum size to be tested (default %d)\n"
    "-g-size <sz>:   alternative size for g(0) flooding test (default 0)\n"
    "-Send:          send using MPI_Send\n"
    "-Isend:         send using MPI_Isend\n"
    "-Ssend:         send using MPI_Ssend\n"
    "-Recv:          receive using MPI_Recv\n"
    "-Irecv:         receive using MPI_Irecv\n"
    "-o <file>:      write resulting LogP parameters to file <file>\n",
	    LOGP_DEF_MIN_ITS,
	    LOGP_DEF_MAX_ITS,
	    LOGP_DEF_EPS,
	    LOGP_DEF_CONF_INT,
	    LOGP_DEF_MAX_DIFF,
	    LOGP_DEF_MIN_SIZE,
	    LOGP_DEF_MAX_SIZE
    );
    exit(1);
}

static char *
mpi_recv_name(logp_recvmode recv_mode)
{
    switch (recv_mode) {
    case LOGP_MPI_RECV:	 return "Recv";
    case LOGP_MPI_IRECV: return "Irecv";
    default:             return "UnknownRecvmode";
    }
}

static char *
mpi_send_name(logp_sendmode send_mode)
{
    switch (send_mode) {
    case LOGP_MPI_SEND:	 return "Send";
    case LOGP_MPI_ISEND: return "Isend";
    case LOGP_MPI_SSEND: return "Ssend";
    default:             return "UnknownSendmode";
    }
}

int
main(int argc, char **argv)
{
    int i;
    int err;
    logp_params *logp;
    int rank, nodes;
    int logp_min_its, logp_max_its;
    double logp_conf_int, logp_eps, logp_max_diff;
    int logp_min_size, logp_max_size;
    int logp_g_zero;
    logp_sendmode logp_send_mode;
    logp_recvmode logp_recv_mode;
    logp_bufmode  logp_buf_mode;
    int logp_opt;
    int logp_diag;
    char *logp_outfile;
    char test_name[256];

    prog_name = argv[0];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    if (nodes != 2) {
	/* We might be interested in testing multiple pairs in parallel,
	 * but the names of the output files should be different in that
	 * case, obviously.
	 */
	fprintf(stderr, "%s: should be run on exactly two nodes\n", argv[0]);
	usage();
    }

    /* Set defaults and process options.
     * TODO: check sensibility of parameter values (non-negative etc.)
     */
    logp_g_zero    = 0;
    logp_min_its   = LOGP_DEF_MIN_ITS;
    logp_max_its   = LOGP_DEF_MAX_ITS;
    logp_conf_int  = LOGP_DEF_CONF_INT;
    logp_eps       = LOGP_DEF_EPS;
    logp_max_diff  = LOGP_DEF_MAX_DIFF;
    logp_min_size  = LOGP_DEF_MIN_SIZE;
    logp_max_size  = LOGP_DEF_MAX_SIZE;
    logp_send_mode = LOGP_MPI_SEND;
    logp_recv_mode = LOGP_MPI_RECV;
    logp_buf_mode  = LOGP_BUF_DYNAMIC;
    logp_opt       = LOGP_OPT_NONE;
    logp_diag      = LOGP_DIAG_DEFAULT;
    logp_outfile    = NULL;
    for (i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-flood") == 0) {
	    logp_opt |= LOGP_OPT_GM_DIRECT;
	} else if (strcmp(argv[i], "-silent") == 0) {
	    logp_diag = LOGP_DIAG_SILENT;
	} else if (strcmp(argv[i], "-verbose") == 0) {
	    if (logp_diag < LOGP_DIAG_VERBOSE) {
		logp_diag = LOGP_DIAG_VERBOSE;
	    } else {
		logp_diag++;
	    }
	} else if (strcmp(argv[i], "-min-its") == 0) {
	    if (i + 1 < argc) {
		logp_min_its = atoi(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -min-its: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-max-its") == 0) {
	    if (i + 1 < argc) {
		logp_max_its = atoi(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -max-its: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-eps") == 0) {
	    if (i + 1 < argc) {
		logp_eps = atof(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -eps: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-conf-int") == 0) {
	    if (i + 1 < argc) {
		logp_conf_int = atof(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -conf-int: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-max-diff") == 0) {
	    if (i + 1 < argc) {
		logp_max_diff = atof(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -max-diff: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-g-size") == 0) {
	    if (i + 1 < argc) {
		logp_g_zero = atoi(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -g-size: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-min-size") == 0) {
	    if (i + 1 < argc) {
		logp_min_size = atoi(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -min-size: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-max-size") == 0) {
	    if (i + 1 < argc) {
		logp_max_size = atoi(argv[++i]);
	    } else {
		fprintf(stderr, "%s: -max-size: argument expected\n", argv[0]);
		usage();
	    }
	} else if (strcmp(argv[i], "-static") == 0) {
	    logp_buf_mode = LOGP_BUF_STATIC;
	} else if (strcmp(argv[i], "-dynamic") == 0) {
	    logp_buf_mode = LOGP_BUF_DYNAMIC;
	} else if (strcmp(argv[i], "-Send") == 0) {
	    logp_send_mode = LOGP_MPI_SEND;
	} else if (strcmp(argv[i], "-Isend") == 0) {
	    logp_send_mode = LOGP_MPI_ISEND;
	} else if (strcmp(argv[i], "-Ssend") == 0) {
	    logp_send_mode = LOGP_MPI_SSEND;
	} else if (strcmp(argv[i], "-Recv") == 0) {
	    logp_recv_mode = LOGP_MPI_RECV;
	} else if (strcmp(argv[i], "-Irecv") == 0) {
	    logp_recv_mode = LOGP_MPI_IRECV;
	} else if (strcmp(argv[i], "-o") == 0) {
	    if (i + 1 < argc) {
		logp_outfile = argv[++i];
	    } else {
		fprintf(stderr, "%s: -o: argument expected\n", argv[0]);
		usage();
	    }
	} else {
	    usage();
	}
    }

    sprintf(test_name, "logp_test.%s.%s",
	    mpi_send_name(logp_send_mode), mpi_recv_name(logp_recv_mode));

    err = logp_measure(&logp, test_name,
		       MPI_COMM_WORLD, 0, 1,
		       logp_g_zero,
		       logp_min_its, logp_max_its,
		       logp_conf_int, logp_eps, logp_max_diff,
		       logp_min_size, logp_max_size,
		       LOGP_DEF_DIFF_SIZES,
		       logp_send_mode, logp_recv_mode, logp_buf_mode,
		       logp_opt, logp_diag);

    if (err == MPI_SUCCESS && logp != NULL) {
	if (logp_outfile != NULL) {
	    logp_save(logp, logp_outfile);

#ifdef LOGP_SANITY_CHECK
	    /* A few simple sanity checks of the LogP/MPI API */
	    {
		logp_params *logp_check;
		int size;

		logp_load(&logp_check, logp_outfile);
		logp_save(logp_check, "logp_test.out.check");

		/* and see if retrieval of values works */
		for (size = 1024; size < 1024 * 1024; size *= 5) {
		    printf("%d: os = %9.7f or = %9.7f g = %9.7f\n",
			   size,
			   logp_os(logp_check, size),
			   logp_or(logp_check, size),
			   logp_g (logp_check, size));
		}
	    }
#endif
	}

	logp_free(logp);
    }

    MPI_Finalize();
    return 0;
}
