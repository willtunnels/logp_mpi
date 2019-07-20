#ifndef LOGP_STATS_H
#define LOGP_STATS_H

typedef struct {
    double num;
    double sum;
    double sum_sq;
    double min;
    double max;
} statistics;

void   stats_init   (statistics *_stats);
void   stats_update (statistics *_stats, double _val);
double stats_mean   (statistics *_stats);
double stats_var    (statistics *_stats);
double stats_min    (statistics *_stats);
double stats_max    (statistics *_stats);
double stats_t      (int _df, int _tail, double _conf);
double stats_conf_int(statistics *_stats, int _tail, int _iters,
		      double _confint);

#endif /* LOGP_STATS_H */
