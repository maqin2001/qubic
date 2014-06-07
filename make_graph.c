/* Author: Qin Ma <maqin@uga.edu>, Step. 19, 2013
 * Usage: This is part of the bicluster package. Use, redistribute, modify
 *        without limitations.
 *
 * Produces two graphs sequentially, derived from microarray data.
 * 
 * The first graph is generated from the raw data where an edge is defined
 * as two genes having common components from same condition, the score on
 * the edge being the number of the same components. The edges in the first
 * graph, are used as vertices in the second graph, where edges are defined
 * as the common columns between two co-vertex edges in the first graph,
 * with scores defined as the number of same columns for all three genes.
 *
 */

#include "make_graph.h"
/*we can reduce the HEAP_SIZE when the data contain so many genes so that memory is not enough*/
static const int HEAP_SIZE = 20000000;

/**************************************************************************/
/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge in the first graph*/
static int str_intersect_r (const discrete *s1, const discrete *s2)
{
	int common_cnt = 0;
	/* s1 and s2 of equal length, so we check s1 only */
	int i;
	for (i=0;i<cols;i++)
	{
		if (*s1==*s2 && (*s1!=0)) 
			common_cnt++;
		s1++; 
		s2++;
	}
	return common_cnt;
}

/**************************************************************************/
continuous get_pearson (const discrete *s1, const discrete *s2,int row_1, int row_2, int cnt)
{
	continuous ss1[cnt],ss2[cnt];
	continuous score=0, ave1=0, ave2=0, var1=0, var2=0, cntc=cnt;
	int i=0,j=0;
	for (i=0;i<cols;i++)
	{
		if (s1[i]==s2[i] && (s1[i]!=0))
		{
			ss1[j]= arr[row_1][i];	
			ss2[j]= arr[row_2][i];
			j++;
		}
		i++;
	}
	/*get var and ave*/
	for (j=0; j<cnt; j++)
	{
		ave1+=ss1[j];
		ave2+=ss2[j];
	}
	ave1 = ave1/cntc;
	ave2 = ave2/cntc;
	for (j=0; j<cnt; j++)
	{
		var1+= (ss1[j]-ave1)*(ss1[j]-ave1);
		var2+= (ss2[j]-ave2)*(ss2[j]-ave2);
	}
	var1 = sqrt(var1);
        var2 = sqrt(var2);
        for (i=0; i<cnt; i++)
                score += (ss1[i]-ave1)*(ss2[i]-ave2);
        score = fabs(score/(var1*var2));
        return score;
}	
/**************************************************************************/
continuous get_spearman (discrete *s1, discrete *s2,int row_1, int row_2, int cnt)		
{
	discrete ss1[cnt],ss2[cnt];
	continuous ss11[cnt],ss22[cnt],temp1[cnt],temp2[cnt];
	continuous score=0, ave1=0, ave2=0, var1=0, var2=0, cntc = cnt;
	int i=0,j=0;
	for (i=0;i<cols;i++)
	{

		if (symbols[s1[i]]==symbols[s2[i]] && (s1[i]!=0))
		{
			ss11[j]= arr[row_1][i];	
			ss22[j]= arr[row_2][i];
			temp1[j]= arr[row_1][i];	
			temp2[j]= arr[row_2][i];
			j++;
		}
	}
	qsort(temp1, cnt, sizeof *temp1, compare_continuous);	
	for (i=0;i<cnt;i++)
	{	
		ss1[i]=0;
		for(j=0;j<cnt;j++)
		{
			if (ss11[i] == temp1[j])
				ss1[i] = j;
		}
	}
	qsort(temp2, cnt, sizeof *temp2, compare_continuous);	
	for (i=0;i<cnt;i++)
	{	
		ss2[i]=0;
		for(j=0;j<cnt;j++)
			if (ss22[i] == temp2[j])
				ss2[i] = j;
	}
	/*get var and ave*/
	for (j=0; j<cnt; j++)
	{
		ave1+=ss1[j];
		ave2+=ss2[j];
	}
	ave1 = ave1/cntc;
	ave2 = ave2/cntc;
	for (j=0; j<cnt; j++)
	{
		var1+= (ss1[j]-ave1)*(ss1[j]-ave1);
		var2+= (ss2[j]-ave2)*(ss2[j]-ave2);
	}
	var1 = sqrt(var1);
        var2 = sqrt(var2);
        for (i=0; i<cnt; i++)
                score += (ss1[i]-ave1)*(ss2[i]-ave2);
	/*printf ("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ave1, ave2, var1, var2, score);*/
        score = fabs(score/(var1*var2));
        return score;
}
/*************************************************************/
/* String intersection function without string copying, only numbers */
/*caculate the weight of the edge with negative regulation in the first graph*/
/*static int str_negative_intersect_r (const discrete *s1, const discrete *s2)
{
	int common_cnt = 0;
	int i;
	for (i=0;i<cols;i++)
	{
		if ((s1[i]!=0) && (symbols[s1[i]] == -symbols[s2[i]]))
			common_cnt++;
	}
	return common_cnt;
}*/
/**************************************************************************/

/*seed_deduct is not used in the whole program while it may be useful someday*/
void seed_deduct (const discrete *s)
/* remove a row from the profile */
{
	int i;
	discrete ss;
	for (i=0; i<cols; i++)
	{
		ss = s[i];
		profile[i][ss]--;
	}
}

void seed_update (const discrete *s)
{
	int i;
	for (i = 0; i < cols; i++) 
		profile[i][s[i]]++;
}

/*****************************************************************************/

/* Fibonacci heap related subroutines */
static int edge_cmpr(void *a, void *b)
{
	int score_a, score_b;
	score_a = ((Edge *)a)->score;
	score_b = ((Edge *)b)->score;

	if (score_a < score_b) return -1;
	if (score_a == score_b) return 0;
	return 1;
}

/* Maintain a fixed size heap */
static void fh_insert_fixed(struct fibheap *a, Edge *i, Edge **cur_min)
{
	if (a->fh_n < HEAP_SIZE) 
	{
		fh_insert(a, (void *)i);
	}
	else
	{
		if (edge_cmpr(cur_min, i) < 0)
		{
			/* Remove least value and renew */
			fh_extractmin(a);
			fh_insert(a, (void *)i);
			/* Keep a memory of the current min */
			*cur_min = (Edge *)fh_min(a);
		}
	}
}

/*sort the edges in decrease order so that e1 is the largest edge*/
static void fh_dump(struct fibheap *a, Edge **res)
{
	int i;
	int n = a->fh_n;
	for (i=n-1; i>=0; i--)
		res[i] = (Edge *) fh_extractmin(a);
}

/*calculate the F-score*/
continuous get_f_socre (continuous a, continuous b, continuous c)
{
	continuous x, y, z, f, z1, z2;
	x = a;
	y = b;
	z = c;
	f = 0;
	z1 = x/y;
	z2 = x/z;
	f = 2*z1*z2/(z1+z2);
	return f;
}

/**************************************************************************/

void make_graph (const char* fn)
{
	FILE *fw = mustOpen(fn, "w");
	int i, j, cnt, cnt1, cnt2, cnt3;
	int rec_num = 0;
	if (po->COL_WIDTH == 2) 
		po->COL_WIDTH = MAX(cols/20, 2);
	
	/*the virable for spearman calculate*/
	continuous spearman=0; /*int cnt_r=0;*/
	continuous fscore=0; 
	continuous final=0; 

	/* edge_ptr describe edges */
	AllocArray(edge_list, HEAP_SIZE);

	/* Allocating heap structure */
	struct fibheap *heap;
	heap = fh_makeheap();
	fh_setcmp(heap, edge_cmpr);

	/* Generating seed list and push into heap */
	progress("Generating seed list (minimum weight %d)", po->COL_WIDTH);	

	Edge __cur_min = {0, 0, po->COL_WIDTH};
	Edge *_cur_min = &__cur_min;
	Edge **cur_min = & _cur_min;
	/* iterate over all genes to retrieve all edges */
	if (po->IS_TFname)
	{	
		for (i = 0; i < rows; i++)
		{
			cnt = str_intersect_r(arr_c[i], arr_c[TFindex]);
			cnt1 = str_intersect_r(arr_c[i], arr_c[i]);
			cnt2 = str_intersect_r(arr_c[TFindex], arr_c[TFindex]);
			if (po->IS_spearman && cnt>5)
			{
				/*cnt_r = str_negative_intersect_r (arr_c[i], arr_c[j]);		
				cnt = MAX(cnt, cnt_r);*/
				/*get spearman rank corelation*/
				spearman = get_spearman (arr_c[i],arr_c[TFindex],i,TFindex,cnt);
                               	/*printf ("spearman\t%s\t%s\t%d\t%.2f\t%d\t%d\n",genes_n[i],genes_n[j],cnt,spearman,cnt1, cnt2);*/
				fscore = get_f_socre (cnt, cnt1, cnt2);
				cnt3 = cnt;
				cnt = ceil (2*cnt3*spearman*fscore);
				final = 2*cnt3*spearman*fscore;
                       		/*printf ("spearman\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%d\n",i, TFindex, genes_n[i],genes_n[TFindex],cnt3,cnt1, cnt2, fscore,spearman, cnt);*/
                       		printf ("%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",genes_n[i],genes_n[TFindex],cnt3,cnt1, cnt2, fscore,spearman, final);
			}

			if (cnt < (*cur_min)->score) continue;
	
			AllocVar(edge_ptr);
			edge_ptr -> gene_one = i;
			edge_ptr -> gene_two = TFindex;
			edge_ptr -> score = cnt;	
			fh_insert_fixed(heap, edge_ptr, cur_min);
		}
	}
	else{
	for (i = 0; i < rows; i++)
	{
		for (j = i+1; j < rows; j++)
		{
			cnt = str_intersect_r(arr_c[i], arr_c[j]);
			cnt1 = str_intersect_r(arr_c[i], arr_c[i]);
			cnt2 = str_intersect_r(arr_c[j], arr_c[j]);
			if (po->IS_spearman && cnt>5)
			{
				/*cnt_r = str_negative_intersect_r (arr_c[i], arr_c[j]);		
				cnt = MAX(cnt, cnt_r);*/
				/*get spearman rank corelation*/
				spearman = get_spearman (arr_c[i],arr_c[j],i,j,cnt);
                                /*printf ("spearman\t%s\t%s\t%d\t%.2f\t%d\t%d\n",genes_n[i],genes_n[j],cnt,spearman,cnt1, cnt2);*/
				fscore = get_f_socre (cnt, cnt1, cnt2);
				cnt3 = cnt;
				cnt = ceil (2*cnt3*spearman*fscore);
				final = 2*cnt3*spearman*fscore;
                       		printf ("%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n",genes_n[i],genes_n[j],cnt3,cnt1, cnt2, fscore,spearman, final);
			}

			if (cnt < (*cur_min)->score) continue;
		
			AllocVar(edge_ptr);
			edge_ptr -> gene_one = i;
			edge_ptr -> gene_two = j;
			edge_ptr -> score = cnt;	
			fh_insert_fixed(heap, edge_ptr, cur_min);
		}
	}
	}
	rec_num = heap->fh_n;
	if (rec_num == 0)
		errAbort("Not enough overlap between genes");

	/* sort the seeds */
	uglyTime("%d seeds generated", rec_num);
	ReAllocArray(edge_list, rec_num);
	fh_dump(heap, edge_list);

	/* bi-clustering */
        int n_blocks = 0;
	progress("Clustering started");
	n_blocks = cluster(fw, edge_list, rec_num);
	uglyTime("%d clusters are written to %s", n_blocks, fn);

	/* clean up */
	for (i=0; i<rec_num; i++)
		free(edge_list[i]);
	free(edge_list);
}

/***************************************************************************/
