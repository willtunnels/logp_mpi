#ifndef logp_mpi_h
#define logp_mpi_h

#include <time.h>

/* Note: the next structure should be seen as an opaque type,
 * and not directly be used by clients outside the logp_mpi module.
 */
typedef struct {
    int    size;	/* the sizes (to be) measured */
    time_t time;	/* the time the measurements started */
    double RTT;		/* round trip times */
    double RTT_min;	/* RTT minimum */
    double os;		/* send overheads */
    double os_cnfint;	/* os confidence interval */
    double os_min;	/* os minimum */
    double or;		/* receive overheads */
    double or_cnfint;	/* or confidence interval */
    double or_min;	/* or minimum */
    double g;		/* inverse bandwidths */
} logp_size_params;

typedef struct {
    char   *descr;	/* descriptive name of the measurements/network */
    int     nsizes;	/* current number of sizes */
    int     nsizes_max;	/* max #sizes in following arrays (allocated) */
    double  L;		/* message latency (excluding send/recv overheads) */
    logp_size_params *params;
			/* size dependent parameters, sorted by size */
} logp_params;

typedef enum { LOGP_MPI_SEND, LOGP_MPI_ISEND, LOGP_MPI_SSEND } logp_sendmode;
typedef enum { LOGP_MPI_RECV, LOGP_MPI_IRECV } logp_recvmode;
typedef enum { LOGP_BUF_STATIC, LOGP_BUF_DYNAMIC } logp_bufmode;

/* Reasonable default arguments for the logp_measure() below: */
#define LOGP_DEF_MIN_ITS	3
#define LOGP_DEF_MAX_ITS	18
#define LOGP_DEF_CONF_INT	0.90	/* % confidence interval */
#define LOGP_DEF_EPS		0.03	/* average precision % requested */
#define LOGP_DEF_MAX_DIFF	0.20	/* >20% difference: investigate */
#define LOGP_DEF_MIN_SIZE	1	/* size 0 is always tested */
#define LOGP_DEF_MAX_SIZE	(256 * 1024)
#define LOGP_DEF_DIFF_SIZES	1024
#define LOGP_DEF_SENDMODE	LOGP_MPI_SEND
#define LOGP_DEF_RECVMODE	LOGP_MPI_RECV

#define LOGP_DIAG_SILENT	0
#define LOGP_DIAG_DEFAULT	1
#define LOGP_DIAG_VERBOSE1	2
#define LOGP_DIAG_VERBOSE2	3
#define LOGP_DIAG_VERBOSE	LOGP_DIAG_VERBOSE1

/* Possible values for the "option" argument to logp_measure: */
#define LOGP_OPT_NONE		0x0	/* default "fast" measurements */
#define LOGP_OPT_GM_DIRECT	0x1	/* ditto, but redo g's with flooding */

int logp_init(void);
int logp_finalize(void);
 
int logp_save   (logp_params  *logp, char *filename);
int logp_load   (logp_params **logp, char *filename);
int logp_measure(
    logp_params **logp,
    char *descr,		/* descriptive name of the network	     */
    MPI_Comm comm,		/* MPI communicator used for tests	     */
    int measure,		/* rank of measurement process		     */
    int mirror,			/* rank of mirror process		     */
    int logp_g_zero,		/* msg size to be used for g(0) measurement  */
    int iters,			/* def #iterations,      LOGP_DEF_MIN_ITS    */
    int max_iters,		/* max #iterations,      LOGP_DEF_MAX_ITS    */
    double confint,		/* req confidence interval, LOGP_DEF_CONF_INT*/
    double eps,			/* req error precision,  LOGP_DEF_EPS	     */
    double max_diff,		/* nonlinearity search,  LOGP_DEF_MAX_DIFF   */
    int minsize,		/* min msg size,         LOGP_DEF_MIN_SIZE   */
    int maxsize,		/* initial max msg size, LOGP_DEF_MAX_SIZE   */
    int max_different_sizes,	/* max #sizes,           LOGP_DEF_DIFF_SIZES */
    logp_sendmode send_mode,	/* MPI sendmode,         LOGP_DEF_SENDMODE   */
    logp_recvmode recv_mode,	/* MPI recvmode,         LOGP_DEF_RECVMODE   */
    logp_bufmode buf_mode,	/* how to use Send/receive receive buffers   */
    unsigned int option,	/* additional boolean options		     */
    int diagnostics		/* verbosity level,      LOGP_DIAG_*	     */
    );
void logp_free  (logp_params *logp);

time_t logp_time(logp_params *logp, int message_size);
double logp_os(logp_params *logp, int message_size);
double logp_or(logp_params *logp, int message_size);
double logp_g (logp_params *logp, int message_size);
double logp_L (logp_params *logp);

#endif /* logp_mpi_h */
