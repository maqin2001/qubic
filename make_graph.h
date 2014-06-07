#ifndef _MAKE_GRAPH_H
#define _MAKE_GRAPH_H

#include "struct.h"

/* global data */
int col_width;
Edge **edge_list;
Edge *edge_ptr;

/* from cluster */
extern bits16 **profile;
extern int compare_continuous (const void *a, const void *b);

/* global functions */
extern int cluster (FILE *fw, Edge **el, int n);

/* prototypes */
void seed_update (const discrete *s);
static int str_intersect_r (const discrete *s1, const discrete *s2);
/*static int str_negative_intersect_r (const discrete *s1, const discrete *s2);*/
void seed_deduct (const discrete *s);
void make_graph (const char *fn);
static int edge_cmpr(void *a, void *b);
static void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min);
static void fh_dump(struct fibheap *a, Edge **res);
continuous get_spearman (discrete *s1, discrete *s2,int row_1, int row_2, int cnt);		
continuous get_pearson (const discrete *s1, const discrete *s2,int row_1, int row_2, int cnt);		
continuous get_f_socre (continuous a, continuous b, continuous c);

#endif
