/* 
 * Author:Qin Ma <maqin@csbl.bmb.uga.edu>, Jan. 24, 2010 
 * Usage: This is part of the bicluster package. Use, redistribute, modify without limitations.
 *
 * Include two procedures for file input:
 * o read_continuous() would read a file with this format:
 * ----------------------------------------
 * 	 	cond1	 cond2	 cond3
 * 	gene1    3.14	  -1.2     0.0
 * 	gene2      nd      2.8     4.5
 * ----------------------------------------
 *   values may possibly be any continuous value, e.g. log-ratio of 
 *   lumin intensity for two channels. The procedure then, for each 
 *   row, produces a distribution using method similar to outlier algorithm, 
 *   base on two tails of values (6%),
 *   the middle part, is regarded as insignificant. This would discretize
 *   the continuous value into classes. (If you want divide the data into 
 *   more levels, you can adjust the parameter r and q) See below.
 * 
 * o read_discrete() would read a file with format like:
 * ----------------------------------------
 *  		cond1	 cond2	 cond3
 *  	gene1	    1        1       0
 *  	gene2       1       -1       0
 * ----------------------------------------
 *   the symbols could be any integers (-32768~+32767) and represent distinct classes.
 *   '0', however, will be ignored, and uncounted in the alter algorithms.
 *   since they would represent a no-change class.
 */

#include "read_array.h"

/************************************************************************/
/* Helper variables for tokenizer function */

static char *atom = NULL;
static char delims[] = "\t\r\n";
#define MAXC 100000
/* record the position of each discretized symbol in _symbols_ */
/* an unsigned short can hold all the values between 0 and  USHRT_MAX inclusive. USHRT_MAX must be at least 65535*/
static int bb[USHRT_MAX];

/***********************************************************************/

/* Comparison function for GNU qsort */
int compare_continuous (const void *a, const void *b)
{
    const continuous *da = a;
    const continuous *db = b;
    /*make qsort in the increasing order*/
    return (*da < *db)?-1:(*da != *db);
}
/* emulate gnu gsl quantile function */
/*divide by the number of the data*/
static continuous quantile_from_sorted_data(const continuous sorted_data[], size_t n, double f)
{
	/*floor function returns the largest integral value less than or equal to x*/
	int i = floor((n-1)*f);
	continuous delta = (n-1)*f-i;
	return (1-delta)*sorted_data[i]+delta*sorted_data[i+1];
}
/*divide by the value of the data*/
/*static continuous quantile_from_sorted_data_value(const continuous sorted_data[],size_t n, double f)
{
        return sorted_data[0]+f*(sorted_data[n-1]-sorted_data[0]);
}*/
/***********************************************************************/

static int charset_add(discrete *ar, discrete s)
{
	/*A signed short can hold all the values between SHRT_MIN  and SHRT_MAX inclusive.SHRT_MIN is required to be -32767 or less,SHRT_MAX must be at least 32767*/
	int ps = s + SHRT_MAX;
	if (bb[ps]<0)
	{
		bb[ps] = sigma;
		ar[sigma++] = s;
	}
	return bb[ps];
}

/***********************************************************************/

/* Matrix allocations (continuous and discrete 2d array) */

continuous** alloc2d(int rr, int cc)
{
	continuous** result;
	int i;
	AllocArray(result, rr);
	for (i = 0; i < rr; i++)
		AllocArray(result[i], cc);
	return result;
}

discrete** alloc2c(int rr, int cc)
{
	discrete** result;
	int i;
	AllocArray(result, rr);
	for (i = 0; i < rr; i++)
		AllocArray(result[i], cc);
	return result;
}

/***********************************************************************/

/* Pre-read the datafile, retrieve gene labels and condition labels
 * as well as determine the matrix size
 */
void get_matrix_size (FILE* fp)
{
	/*size_t is the best type to use if you want to represent sizes of objects. 
	 * Using int to represent object sizes is likely to work on most modern systems, but it isn't guaranteed.
	 */
	size_t n = 0;
	char *line;
	/*getline() reads an entire line, storing the address of the buffer containing the text into *line.
	 *the buffer is null-terminated and includes the newline character, if a newline delimiter was found.
	 */
	if (getline(&line, &n, fp)>=0)
	{
		/*strtok function returns a pointer to the next token in str1, where str2 contains the delimiters that determine the token*/
		atom = strtok(line, delims);
		/*delete the first element in atom because the first element corresponding to description column*/
		atom = strtok(NULL, delims);
		while (atom != NULL)
		{
			/*if we do not set atom = strtok(NULL, delims), here is a infinite loop*/
			atom = strtok(NULL, delims);
			cols++;
		}			
	}
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		rows++;
	}
	/*fseed sets the position indicator associated with the stream to a new position defined by adding offset to a reference position specified by origin*/
	fseek(fp, 0, 0);
}

/* Read in the labels on x and y, in microarray terms, genes(rows) and conditions(cols)*/
void read_labels (FILE* fp)
{
	int row = 0;
	int col;
	size_t n = 0;
	char *line;
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		/*currently the first element in atom is the gene name of each row when row>=1, the 0 row corresponding to the line of condition names*/
		if (row >= 1) 
		{
			strcpy(genes_n[row-1], atom);
			/*check if there exist a gene name equals to TFname by -T*/
			if (strcmp (atom, po->TFname)==0)
			{
				TFindex = row-1;
				/*printf ("%d\n",TFindex);*/
			}
		}
		/*delete the first element in atom because the first element corresponding to description column*/
		atom = strtok(NULL, delims);
		col = 0;
		while (atom != NULL)
		{
			if (row == 0) 
				strcpy(conds[col], atom);
			atom = strtok(NULL, delims);
			if (++col == cols) break;
		}
		if (++row == rows+1) break;
	}
	fseek(fp, 0, 0);
}

/*read in the sub-gene list*/
void read_list (FILE* fp)
{
	int i=0, j=0;
	sub_genes_row = 0;
	char line[MAXC];
	while(fgets(line,MAXC,fp)!=NULL)
	{
		atom = strtok(line, delims);
		strcpy(sub_genes[sub_genes_row], atom);	
		sub_genes_row++;
	}

	/*update the sub_list*/
	AllocArray(sublist,rows);
	for (i = 0; i<rows; i++)
		sublist[i] = FALSE;
	for (i=0; i<sub_genes_row; i++)
		for (j=0; j<rows; j++)
			if (strcmp (sub_genes[i], genes_n[j])==0)
				sublist[j] = TRUE;
}

/*read in the f3 value for each gene base on kernel density estimation*/
/*void read_density (FILE* fp)
{
	int k=0;
	char line[MAXC];
	AllocArray(density,rows);
	while(fgets(line,MAXC,fp)!=NULL)
	{
		atom = strtok(line, delims);
		density[k]=atof(atom);	
		k++;
	}
}*/

/* initialize data for discretization */
void init_dis()
{
	int row, col;
	/* store discretized values */
	AllocArray(symbols, USHRT_MAX);
	/* memset sets the first num bytes of the block of memory pointed by ptr to the specified value
	 * memset ( void * ptr, int value, size_t num )*/
	memset(bb, -1, USHRT_MAX*sizeof(*bb));
	/* always add an 'ignore' index so that symbols[0]==0*/
	charset_add(symbols, 0); 
	/*initialize for arr_c*/
	arr_c = alloc2c(rows,cols);
	for (row = 0; row < rows; row++)
        	for (col = 0; col < cols; col++) 
			arr_c[row][col] = 0;
}

void read_discrete (FILE* fp)
{
	int row, col, i;
	init_dis();
	/* import data */
	size_t n = 0;
	char *line;
	row = 1;
	/* Skip first line with condition labels */
	getline(&line, &n, fp);
	/* read the discrete data from the second line */
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		/*skip the first column*/
		atom = strtok(NULL, delims);
		col = 0;
		while (atom != NULL)
		{			
			arr_c[row-1][col] = charset_add(symbols, atoi(atom));
			atom = strtok(NULL, delims);
			if (++col == cols) break;
		}
		if (++row == rows+1) break;
	}
	/* trim the leading spaceholder */
	printf("Discretized data contains %d classes with charset [ ", sigma);
	for (i=0;i<sigma;i++) 
		/*printf("%d ", symbols[i]); printf("]\n");*/
		printf("%d ", i); printf("]\n");
	fseek(fp, 0, 0);
}

void read_continuous (FILE* fp)
{
	int row, col;
	arr = alloc2d(rows,cols);
	/* import data */
	size_t n = 0;
	char *line;
	row = 1;
	/* ignore header line */
	getline(&line, &n, fp);
	while (getline(&line, &n, fp)>=0)
	{
		atom = strtok(line, delims);
		/*skip the first column*/
		atom = strtok(NULL, delims);
		col = 0;
		while (atom != NULL)
		{
			/*we set all the aplha to ignore value 0*/
			/*Checks if parameter atom is either an uppercase or a lowercase alphabetic letter*/
			if (isalpha(*atom)) 
				arr[row-1][col] = 0.0;
			else 
				arr[row-1][col] = atof(atom);
			atom = strtok(NULL, delims);
			if (++col == cols) break;
		}
		if (++row == rows+1) break;
	}
	fseek(fp, 0, 0);
}

/***********************************************************************/

/* Discretize continuous values by revised outlier detection algorithm
 * see details in Algorithm Design section in paper
 */
discrete dis_value(float current, int divided, float *small, int cntl, float *big, int cntu)
{
	int i;
	float d_space = 1.0 / divided;	
	for(i=0; i < divided; i++)
	{		
            if ((cntl > 0) && (current <= quantile_from_sorted_data(small, cntl, d_space * (i+1)))) 
		    return -i-1;
            if ((cntu > 0) && (current >= quantile_from_sorted_data(big, cntu, 1.0 - d_space * (i+1)))) 
		    return i+1;
	}
	return 0;
}

void discretize (const char* stream_nm)
{
	int row, col;
	continuous rowdata[cols];
	float big[cols], small[cols];
	int i,cntu,cntl;
	float f1,f2,f3, upper, lower;
	FILE *fw;
	fw = mustOpen(stream_nm, "w");
	init_dis();
	for (row = 0; row < rows; row++)
	{
		for (col = 0; col < cols; col++) 
			rowdata[col] = arr[row][col];
		qsort(rowdata, cols, sizeof *rowdata, compare_continuous);
		f1 = quantile_from_sorted_data(rowdata,cols,1-po->QUANTILE); 
		f2 = quantile_from_sorted_data(rowdata,cols,po->QUANTILE);
		/*if (po->IS_density)
		{
			f3 = density[row];
		}
		else
		{*/
			f3 = quantile_from_sorted_data(rowdata,cols,0.5);
		/*}*/
		if ((f1-f3)>=(f3-f2))
		{
			upper = 2*f3-f2; 
			lower = f2;
		}
		else
		{
			upper = f1; 
			lower = 2*f3-f1;
		}
		cntu = 0; cntl = 0;
		for (i = 0; i < cols; i++)
		{
			if (rowdata[i] < lower) 
			{ 
				small[cntl] = rowdata[i]; 
				cntl++; 
			}
			if (rowdata[i] > upper) 
			{ 
				big[cntu] = rowdata[i]; 
				cntu++; 
			}
		}		
		for (col = 0; col < cols; col++)
	    		arr_c[row][col] = charset_add(symbols, dis_value(arr[row][col],po->DIVIDED, small, cntl, big, cntu));
		if (abs(cntl-cntu) <= 1)
			fprintf(fw, "%s_unexpressed :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", genes_n[row], lower, upper,cntl, cntu);
		else
			fprintf(fw, "%s :low=%2.5f, up=%2.5f; %d down-regulated,%d up-regulated\n", genes_n[row], lower, upper,cntl, cntu);
    	}
        progress ("Discretization rules are written to %s", stream_nm);
	fclose(fw);
}

/* output the formatted matrix */
void write_imported (const char* stream_nm)
{
	int row, col;
	FILE *fw;
	fw = mustOpen(stream_nm, "w"); 
	fprintf(fw, "o");
	for (col = 0; col < cols; col++)
		fprintf(fw,"\t%s",conds[col]);
	fputc('\n',fw);
	for (row = 0; row < rows; row++)
	{
		fprintf(fw, "%s", genes_n[row]);
		for (col = 0; col < cols; col++)
			fprintf(fw, "\t%d", symbols[arr_c[row][col]]);
		fputc('\n', fw);
	}
	progress("Formatted data are written to %s", stream_nm);
	fclose(fw);
}

/***********************************************************************/
continuous get_KL (discrete *array, discrete *array_background, int a, int b)
{
	int i, j;
       	continuous *num, *num_b;
	continuous IC=0;
	AllocArray (num, sigma);
	AllocArray (num_b, sigma);	
	for (i=0; i<sigma; i++)
	{
		num[i]=0;
		num_b[i]=0;
	}
	for (i=0;i<sigma;i++)
	{
		for (j=0;j<a;j++)
			if (symbols[array[j]] == symbols[i])
				num[i]++;
		for (j=0;j<b;j++)
			if (symbols[array_background[j]] == symbols[i])
				num_b[i]++;
	}	
	for (i=0;i<sigma;i++)
	{
		if (num[i] == 0) continue;
		if (num_b[i] == 0) continue;
		IC += (num[i]/a)*log2((num[i]*b)/(num_b[i]*a));
	}
	return IC;
}
/***********************************************************************/
/*new descretization way based on mixture normal distribution*/
double NormSDist(double x, double a, double b)
{
  /* Cumulative Distribution Function */
  x -= a;
  x /= b;
  if(x > 6)  return 1;
  if(x < -6)  return 0.000001;;
  static const double gamma = 0.231641900, a1 = 0.319381530, a2 = -0.356563782, a3 = 1.781477973, a4 = -1.821255978, a5 = 1.330274429;
  double k = 1.0 / (1 + fabs(x) * gamma);
  double n = k * (a1 + k * (a2 + k * (a3 + k * (a4 + k * a5))));
  a = x;
  a = exp((-1)*a*a/2)*0.39894228040143267793994605993438;
  n = 1 - a * n;
  if(x < 0)  return 1.0 - n; 
  return n;
}

double densityFuction(double x, double a, double d)
{
  /* Probability Density Function */
  x = -1 * (x - a) * (x - a) / (2 * d * d);
  x = exp(x); x *= 0.39894228040143267793994605993438; x /= d;
  return x;
}

void discretize_new (const char* stream_nm)
{
  double results2[3][2], results3[3][3], likelihood[3], table_theta_t1[cols][3], m = 0, d = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, c[3][cols]; 
  int i = 0, j = 0, t[cols], tint = 0, INDEX = 0, num_d = 0, a[cols], b[cols], k = 0, ia = 0, ib = 0, sum = 0, ccc[cols];
  int EM = po->EM; /* parameter with default value being 20 or 150 */
  EM--;
  long long id = 0;
  FILE *fw;  
  fw = mustOpen(stream_nm, "w");  
  init_dis();
  for(id = 0; id < rows; id++)
  {
    for(i = 0; i < cols; i++)
    {
      t[i] = i; /* sort by natural numbers */
      for(j = 0; j < 3; j++)  c[j][i] = 0;
      a[i] = 0;
      b[i] = 0;
      ccc[i] = 0;
    }
    m = 0;
    d = 0;
    for(i = 0; i < cols; i++)
      m += arr[id][i]; /* the summation of one row */ 
    m /= i; /* the mean of one row */
    for(i = 0; i < cols; i++)
    {
      /* the squared of the difference between the sample and the expectation */
	  temp1 = arr[id][i];
      temp1 -= m;
      temp1 *= temp1;
      temp2 += temp1;
    }
    temp2 /= (i-1) ; /* unbiased estimated variance */
    d = sqrt(temp2); /* unbiased estimated standard deviation of one row */  
    likelihood[0] = 0;
    for(i = 0; i < cols; i++)
    {
      temp1 = densityFuction(arr[id][i], m, d);
      temp1 = log(temp1);
      likelihood[0] += temp1;
      for(j = i; j < cols; j++)
      {
        if( arr[id][t[i]] > arr[id][t[j]] )
        {
          tint = t[i];
          t[i] = t[j];
          t[j] = tint;
        }
      }
    }
    likelihood[0] *= -1;
    likelihood[0] += 3;
	/* sort one row */
    results2[0][0] = 0.5; /* default weights */
    results2[0][1] = 0.5; /* default weights */
    tint = cols / 3 - 1; 
    results2[1][0] = arr[id][t[tint]]; 
    tint = cols * 2 / 3 - 1; 
    results2[1][1] = arr[id][t[tint]]; 
	/* Divide-and-Conquer (3) */
    results2[2][0] = d; /* default standard deviation */
    results2[2][1] = d; /* default standard deviation */
    for(INDEX = 0; INDEX < EM; INDEX++)
    {
       if( (results2[2][0]==0) || (results2[2][1]==0) )  break;
	  /* ROUNDs 20 or ROUNDs 150 */
	  for(i = 0; i < cols; i++)
        for(j = 0; j < 2; j++)
        {
          temp = densityFuction(arr[id][i], results2[1][j], results2[2][j]);
          temp *= results2[0][j];
          temp1 = densityFuction(arr[id][i], results2[1][0], results2[2][0]);
          temp2 = densityFuction(arr[id][i], results2[1][1], results2[2][1]);
          temp1 *= results2[0][0];
          temp2 *= results2[0][1];
          temp1 += temp2;
          temp /= temp1;
          table_theta_t1[i][j] = temp;
        }
      for(i = 0; i < 2; i++)
      {
        /* calculate the value of results in the next loop */ 
		temp = 0;
        for(j = 0; j < cols; j++)
        {
          temp += table_theta_t1[j][i];
        }
        temp /= cols;
        results2[0][i] = temp; /* calculate results2[0] */ 
        temp = 0;
        for(j = 0; j < cols; j++)
        {
          temp += ( arr[id][j] * table_theta_t1[j][i] );
        }
        temp /= cols;
        temp /= results2[0][i];
        results2[1][i] = temp; /* calculate results2[1] */ 
        temp = 0;
        for(j = 0; j < cols; j++)
        {
          temp1 = arr[id][j];
          temp1 -= results2[1][i];
          temp1 *= temp1;
          temp1 *= table_theta_t1[j][i];
          temp += temp1;
        }
        temp /= cols;
        temp /= results2[0][i];
        temp = sqrt(temp);
        results2[2][i] = temp; /* calculate results2[2] */ 
      }
    }
    temp = 0;
    for(i = 0; i < cols; i++)
    {
      /* calculate likelihood */ 
	  temp1 = results2[0][0] * densityFuction(arr[id][i], results2[1][0], results2[2][0]);
      temp2 = results2[0][1] * densityFuction(arr[id][i], results2[1][1], results2[2][1]);
      temp1 += temp2;
      temp1 = log(temp1);
      temp += temp1;
    }
    likelihood[1] = -1 * temp + 6;
    num_d = likelihood[0] < likelihood[1] ? 1 : 2;
   /*#################################################################*/
    results3[0][0] = 0.33; /* default weights */
    results3[0][1] = 0.33; /* default weights */
    results3[0][2] = 0.33; /* default weights */
    tint = cols / 4 - 1;
    results3[1][0] = arr[id][t[tint]];
    tint = cols / 2 - 1;
    results3[1][1] = arr[id][t[tint]];
    tint = cols * 3 / 4 - 1;
    results3[1][2] = arr[id][t[tint]];
	/* Divide-and-Conquer (4) */
    results3[2][0] = d; /* default standard deviation */
    results3[2][1] = d; /* default standard deviation */
    results3[2][2] = d; /* default standard deviation */
    for(INDEX = 0; INDEX < EM; INDEX++)
    {
      if( (results3[2][0]==0) || (results3[2][1]==0) || (results3[2][2]==0) )  break;
	  /* ROUNDs 20 or ROUNDs 150 */
	  for(i = 0; i < cols; i++)
        for(j = 0; j < 3; j++)
        {
          temp = densityFuction(arr[id][i], results3[1][j], results3[2][j]);
          temp *= results3[0][j];
          temp1 = densityFuction(arr[id][i], results3[1][0], results3[2][0]);
          temp2 = densityFuction(arr[id][i], results3[1][1], results3[2][1]);
          temp3 = densityFuction(arr[id][i], results3[1][2], results3[2][2]);
          temp1 *= results3[0][0];
          temp2 *= results3[0][1];
          temp3 *= results3[0][2];
          temp1 += temp2;
          temp1 += temp3;
          temp /= temp1;
          table_theta_t1[i][j] = temp;
        }
      for(i = 0; i < 3; i++)
      {
        /* calculate the value of results in the next loop */ 
		temp = 0;
        for(j = 0; j < cols; j++)
        {
          temp += table_theta_t1[j][i];
        }
        temp /= cols;
        results3[0][i] = temp; /* calculate results2[0] */ 
        temp = 0;
        for(j = 0; j < cols; j++)
        {
          temp += ( arr[id][j] * table_theta_t1[j][i] );
        }
        temp /= cols;
        temp /= results3[0][i];
        results3[1][i] = temp; /* calculate results2[1] */ 
        temp = 0;
        for(j = 0; j < cols; j++)
        {
          temp1 = arr[id][j];
          temp1 -= results3[1][i];
          temp1 *= temp1;
          temp1 *= table_theta_t1[j][i];
          temp += temp1;
        }
        temp /= cols;
        temp /= results3[0][i];
        temp = sqrt(temp);
        results3[2][i] = temp; /* calculate results2[2] */ 
      }
    }
    temp = 0;
    for(i = 0; i < cols; i++)
    {
      /* calculate likelihood */ 
	  temp1 = results3[0][0] * densityFuction(arr[id][i], results3[1][0], results3[2][0]);
      temp2 = results3[0][1] * densityFuction(arr[id][i], results3[1][1], results3[2][1]);
      temp3 = results3[0][2] * densityFuction(arr[id][i], results3[1][2], results3[2][2]);
      temp1 += temp2;
      temp1 += temp3;
      temp1 = log(temp1);
      temp += temp1;
    }
    likelihood[2] = -1 * temp + 9;
    num_d = likelihood[num_d-1] < likelihood[2] ? num_d : 3;
    for(i = 0; i < 3; i++)
      for(j = i; j < 3; j++)
      {
        /* sorting the results3 based on the value of results3[2] */ 
		if( results3[1][i] > results3[1][j] )
        {
          temp = results3[1][i];
          results3[1][i] = results3[1][j];
          results3[1][j] = temp;
          temp = results3[0][i];
          results3[0][i] = results3[0][j];
          results3[0][j] = temp;
          temp = results3[2][i];
          results3[2][i] = results3[2][j];
          results3[2][j] = temp;
        }
      }
    temp = 0;
    for(i = 0; i < cols; i++)
    {
      temp1 = arr[id][i] - m;
      temp1 /= d;
      temp1 = temp1 * temp1 * temp1;
      temp += temp1;
    }
  /*#################################################################*/
    if(num_d == 2)
    {
      for(i = 0; i < 2; i++)
      {
        for(j = 0; j < cols; j++)
        {
          c[i][j] = arr[id][j] - results2[1][i];
          c[i][j] *= c[i][j];
          c[i][j] *= -1;
          temp3 = 2 * results2[2][i] * results2[2][i];
          c[i][j] /= temp3;
          c[i][j] = exp(c[i][j]);
          c[i][j] *= results2[0][i];
          c[i][j] /= results2[2][i];
        }
      }
      for(i = 0; i < cols; i++)
      {
        temp3 = c[0][i] > c[1][i] ? c[0][i] / c[1][i] : c[1][i] / c[0][i];
        if(temp3 > 2)  ccc[i] = c[0][i] > c[1][i] ? 1 : 2;
      }
      a[0] = ccc[t[0]];
      ia = 0;
      b[0] = 1;
      ib = 0;
      k = a[0];
      tint = floor(cols/20) + 1;
      for(i = 1; i < cols; i++)
      {
        if(ccc[t[i]] != k)
        {
          ia++;
          a[ia] = ccc[t[i]];        
          k = ccc[t[i]];
          ib++;
          b[ib] = 1;
        }
        else  b[ib]++;
      }
      ia++;
      ib++;
      temp1 = results2[1][0] < results2[1][1] ? results2[1][0] : results2[1][1];
      temp2 = results2[1][0] > results2[1][1] ? results2[1][0] : results2[1][1];
      sum = 0;
      for(i = 0; i < ia; i++)  if(a[i]!=0) sum ++;
      if(sum != 0)
      {
        if(temp < 0)
        {
           if( (a[0] == 0) & (b[0] < tint) ){
             if(b[1] >= tint)
               if(results2[1][a[1]-1] == temp1)  for(i = 0; i < b[0] + b[1]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);}
           else if(a[0] != 0)
             if(results2[1][a[0]-1] == temp1)  for(i = 0; i < b[0]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);
        }
        else
        {
          if( (a[ia-1] == 0) & (b[ia-1] < tint) ){
            if(b[ia-2] >= tint)
              if(results2[1][a[ia-2]-1] == temp2)
                for(i = (cols - b[ia-1] - b[ia-2]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);}
          else if(a[ia-1] != 0)
            if(results2[1][a[ia-1]-1] == temp2)
              for(i = (cols - b[ia-1]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);
        }
      }
    }
  /*#################################################################*/
    else if(num_d == 3)
    {
      for(i = 0; i < 3; i++)
      {
        for(j = 0; j < cols; j++)
        {
          c[i][j] = arr[id][j] - results3[1][i];
          c[i][j] *= c[i][j];
          c[i][j] *= -1;
          temp3 = 2 * results3[2][i] * results3[2][i];
          c[i][j] /= temp3;
          c[i][j] = exp(c[i][j]);
          c[i][j] *= results3[0][i];
          c[i][j] /= results3[2][i];
        }
      }
      for(i = 0; i < cols; i++)
      {
        tint = c[0][i] > c[1][i] ? 1 : 2;
        tint = c[tint-1][i] > c[2][i] ? tint : 3;
        temp3 = c[tint-1][i] / (c[0][i] + c[1][i] + c[2][i] - c[tint-1][i]);
        if(temp3 > 2)  ccc[i] = tint;
      }
      a[0] = ccc[t[0]];
      ia = 0;
      b[0] = 1;
      ib = 0;
      k = a[0];
      tint = floor(cols/20) + 1;
      for(i = 1; i < cols; i++)
      {
        if(ccc[t[i]] != k)
        {
          ia++;
          a[ia] = ccc[t[i]];        
          k = ccc[t[i]];
          ib++;
          b[ib] = 1;        
        }
        else  b[ib]++;
      }
      ia++;
      ib++;
      temp1 = results3[1][0] < results3[1][1] ? results3[1][0] : results3[1][1];
      temp1 = temp1 < results3[1][2] ? temp1 : results3[1][2];
      temp2 = results3[1][0] > results3[1][1] ? results3[1][0] : results3[1][1];
      temp2 = temp2 > results3[1][2] ? temp2 : results3[1][2];
      sum = 0;
      for(i = 0; i < ia; i++)  if(a[i]!=0) sum ++;
      if(sum != 0)
      {
        if( (a[0] == 0) & (b[0] < tint) )
          if(b[1] >= tint)
            if(results3[1][a[1]-1] == temp1)
              for(i = 0; i < b[0] + b[1]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);
        if(a[0] != 0)
          if(results3[1][a[0]-1] == temp1)
            for(i = 0; i < b[0]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);
        if( (a[ia-1] == 0) & (b[ia-1] < tint) )
          if(b[ia-2] >= tint)
            if(results3[1][a[ia-2]-1] == temp2)
              for(i = (cols - b[ia-1] - b[ia-2]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);
        if(a[ia-1] != 0)
          if(results3[1][a[ia-1]-1] == temp2)
            for(i = (cols - b[ia-1]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);
      }
  /*#################################################################*/
    }
    sum = 0;
    for(i = 0; i < cols; i++)
      if(arr_c[id][i] == 0)  sum++;
    if(sum == cols)
    {
      tint = ceil (cols*po->QUANTILE);
      for(i = 0; i < tint; i++)
      {
        arr_c[id][t[i]] = charset_add(symbols,-1);
        arr_c[id][t[cols-i-1]] = charset_add(symbols,1);
      }
	  fprintf(fw, "%s_unexpressed", genes_n[id]);
    }
	else  fprintf(fw, "%s", genes_n[id]);
	fprintf(fw, "\t");
    if(num_d == 2)
      fprintf(fw, "%lf_%lf\t%lf_%lf\n", results2[1][0], results2[2][0], results2[1][1], results2[2][1]);
    else if(num_d == 3)
      fprintf(fw, "%lf_%lf\t%lf_%lf\t%lf_%lf\n", results3[1][0], results3[2][0], results3[1][1], results3[2][1], results3[1][2], results3[2][2]);
    else if(num_d == 1)
      fprintf(fw, "%lf_%lf\n", m, d);
    else
      fprintf(fw, "\n"); 
  }
  progress ("Discretization rules are written to %s", stream_nm); 
  fclose(fw);  
}

void discretize_rpkm (const char* stream_nm)
{
  double results[9][3][9], table_theta_t1[cols][9], m = 0, d = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, temp4 = 0, c[9][cols], cc[9], cut = 0.0009765625, te[9]; 
  int zeroi = 0, i = 0, j = 0, t[cols], tint = 0, INDEX = 0, num_d = 0, a[cols], b[cols], k = 0, ia = 0, ib = 0, sum = 0, ccc[cols], i_cut = 0;
  int UP = 3, DOWN = 2;
  UP++;
  int EM = po->EM; /* parameter with default value being 20 or 150 */
  EM--;
  int sumRPKM = 0;
  double cutRPKM = 0.01;
  long long id = 0;
  FILE *fw;
  fw = mustOpen(stream_nm, "w");
  init_dis();
  for(id = 0; id < rows; id++)
  {
    cut = 0.0009765625;
    i_cut = 0;
    zeroi = 0;
    sum = 0;
    for(i = 0; i < cols; i++)
    {
      t[i] = i; /* sort by natural numbers */
      for(j = 0; j < 9; j++)  c[j][i] = 0;
      a[i] = 0;
      b[i] = 0;
      ccc[i] = 0;
      if(arr[id][i] == 0)  sum++;
	  if(arr[id][i] < cutRPKM)  sumRPKM++;
    }
    m = 0;
    d = 0;
    sumRPKM *= 2;
    if(sumRPKM > cols)
      for(i = 0; i < cols; i++)
        if(arr[id][i] != 0)  arr[id][i] = log10(arr[id][i]);
    for(i = 0; i < cols; i++)
      m += arr[id][i]; /* the summation of one row */ 
    m /= i; /* the mean of one row */
    for(i = 0; i < cols; i++)
    {
      /* the squared of the difference between the sample and the expectation */
	  temp1 = arr[id][i];
      temp1 -= m;
      temp1 *= temp1;
      temp2 += temp1;
    }
    temp2 /= (i-1) ; /* unbiased estimated variance */
    d = sqrt(temp2); /* unbiased estimated standard deviation of one row */  
    for(i = 0; i < cols; i++)
    {
      for(j = i; j < cols; j++)
      {
        if( arr[id][t[i]] > arr[id][t[j]] )
        {
          tint = t[i];
          t[i] = t[j];
          t[j] = tint;
        }
      }
    }
    for(zeroi = 0; zeroi < cols; zeroi++)  if(arr[id][t[zeroi]] == 0)  break;
	/* sort one row */
    /*#################################################################*/
    cut = arr[id][t[0]];
    if(sumRPKM > cols)  cut /= log10(2);
    for(i_cut = 0; i_cut < zeroi; i_cut++)
      if(arr[id][t[i_cut]] >= cut)  break;
    /*#################################################################*/
    for(num_d = DOWN; num_d < UP; num_d++)
    {
      sum *= 2;
      if( sum >= cols )
      {
        num_d = 1;
        break;
      }
      m = 0;
      d = 0;
      for(i = i_cut; i < zeroi; i++)
        m += arr[id][t[i]]; /* the summation of one row */ 
      i -= i_cut;
      m /= i; /* the mean of one row */
      temp2 = 0;
      for(i = i_cut; i < zeroi; i++)
      {
        /* the squared of the difference between the sample and the expectation */
 	    temp1 = arr[id][t[i]];
        temp1 -= m;
        temp1 *= temp1;
        temp2 += temp1;
      }
      i -= i_cut;
      temp2 /= (i-1) ; /* unbiased estimated variance */
      d = sqrt(temp2); /* unbiased estimated standard deviation of one row */  
      for(j = 0; j < 9; j++)
        for(i = 0; i < 9; i++)
        {
          results[j][0][i] = 1; /* default weights */
          results[j][0][i] /= num_d;
          tint = (zeroi - i_cut) * (i + 1) / (num_d + 1) - 1; 
          if(tint >= cols)  tint = zeroi - 1;
          results[j][1][i] = arr[id][t[tint]];
          /* Divide-and-Conquer */
          results[j][2][i] = d; /* default standard deviation */
        }
      for(INDEX = 0; INDEX < EM; INDEX++)
      {
        if(INDEX != 0)
          for(i = 0; i < num_d; i++)
            if( results[num_d][2][i] == 0 )  sum = 999999;
        if(sum == 999999)  break;
    	  /* ROUNDs 20 or ROUNDs 150 */
   	    for(i = i_cut; i < zeroi; i++)
          for(j = 0; j < num_d; j++)
          {
            temp = densityFuction(arr[id][t[i]], results[num_d][1][j], results[num_d][2][j]);
            temp *= results[num_d][0][j];
            temp2 = 0;
            for(tint = 0; tint < num_d; tint++)
            {
              temp1 = densityFuction(arr[id][t[i]], results[num_d][1][tint], results[num_d][2][tint]);
              temp1 *= results[num_d][0][tint];
              temp2 += temp1;
            }
            temp /= temp2;
            table_theta_t1[i][j] = temp;
          }
        temp = 0;
        for(i = 0; i < num_d; i++)
        {
          cc[i] = NormSDist(cut, results[num_d][1][i], results[num_d][2][i]);
          temp += cc[i];
        }
        for(i = 0; i < num_d; i++)  {cc[i] /= temp; cc[i] *= i_cut;}
        temp3 = 0;
        for(i = 0; i < num_d; i++)
        {
          temp = 0;
          for(j = i_cut; j < zeroi; j++)
          {
            temp += table_theta_t1[j][i];
          }
          temp2 = temp;
          temp += cc[i];
          temp /= zeroi;
          te[i] = temp;
          temp3 += temp;
          temp = 0;
          for(j = i_cut; j < zeroi; j++)
          {
            temp += ( arr[id][t[j]] * table_theta_t1[j][i] );
          }
          temp1 = cut - results[num_d][1][i];
          temp1 /= results[num_d][2][i];
          temp1 = densityFuction(temp1, 0, 1);
          temp1 /= NormSDist(temp1, 0, 1);
          temp1 *= results[num_d][2][i];
          temp1 = results[num_d][1][i] - temp1;
          temp1 *= cc[i];
          temp += temp1;
          temp /= temp2;
          results[num_d][1][i] = temp; /* calculate */ 
          temp = cut - results[num_d][1][i];
          temp /= results[num_d][2][i];
          temp1 = densityFuction(temp, 0, 1); 
          temp1 /= NormSDist(temp, 0, 1);
          temp1 *= (cut - results[num_d][1][i]);
          temp1 /= results[num_d][2][i];
          temp1 = 1 - temp1;
          temp1 *= results[num_d][2][i];
          temp1 *= results[num_d][2][i];
          temp1 *= cc[i];
          temp = 0;
          for(j = i_cut; j < zeroi; j++)
          {
            temp1 = arr[id][t[j]];
            temp1 -= results[num_d][1][i];
            temp1 *= temp1;
            temp1 *= table_theta_t1[j][i];
            temp += temp1;
          }
          temp /= temp2;
          temp = sqrt(temp);
          results[num_d][2][i] = temp; /* calculate */ 
        }
        for(i = 0; i < num_d; i++)  results[num_d][0][i] =  te[i] / temp3; /* calculate */ 
      }
    }
    /*#################################################################*/
    if(num_d != 1)
    {
      for(INDEX = DOWN; INDEX < UP; INDEX++)
      {
        temp3 = 0;
        k = INDEX * 2 + 1;
        for(i = 0; i < cols; i++)
        {
          temp2 = 0;
          for(j = 0; j < INDEX; j++)
          {
            temp1 = densityFuction(arr[id][i], results[INDEX][1][j], results[INDEX][2][j]) * results[INDEX][0][j];
            temp2 += temp1;
          }
          temp3 += log10(temp2); 
        }
        k *= log10(cols);
        temp3 *= 2;
        temp3 -= k;
        if(INDEX == DOWN)  temp4 = temp3;
        else
        {
          tint = temp4 > temp3 ? INDEX : tint;
          temp4 = temp4 > temp3 ? temp3 : temp4;
        }
      }
      num_d = tint;
    }
    /*#################################################################*/
    if(num_d == 1)
    {
      for(i = 0; i < cols; i++)  arr_c[id][i] = charset_add(symbols,2);
    }
    else if(num_d == 2)
    {
      for(i = 0; i < 2; i++)
      {
        for(j = 0; j < cols; j++)
        {
          c[i][j] = arr[id][j] - results[num_d][1][i];
          c[i][j] *= c[i][j];
          c[i][j] *= -1;
          temp3 = 2 * results[num_d][2][i] * results[num_d][2][i];
          c[i][j] /= temp3;
          c[i][j] = exp(c[i][j]);
          c[i][j] *= results[num_d][0][i];
          c[i][j] /= results[num_d][2][i];
        }
      }
      for(i = 0; i < cols; i++)
      {
        temp3 = c[0][i] > c[1][i] ? c[0][i] / c[1][i] : c[1][i] / c[0][i];
        if(temp3 > 2)  ccc[i] = c[0][i] > c[1][i] ? 1 : 2;
      }
      a[0] = ccc[t[0]];
      ia = 0;
      b[0] = 1;
      ib = 0;
      k = a[0];
      tint = floor(cols/20) + 1;
      for(i = 1; i < cols; i++)
      {
        if(ccc[t[i]] != k)
        {
          ia++;
          a[ia] = ccc[t[i]];        
          k = ccc[t[i]];
          ib++;
          b[ib] = 1;
        }
        else  b[ib]++;
      }
      ia++;
      ib++;
      temp1 = results[num_d][1][0] < results[num_d][1][1] ? results[num_d][1][0] : results[num_d][1][1];
      temp2 = results[num_d][1][0] > results[num_d][1][1] ? results[num_d][1][0] : results[num_d][1][1];
      sum = 0;
      for(i = 0; i < ia; i++)  if(a[i]!=0) sum ++;
      if(sum != 0)
      {
        if(temp < 0)
        {
           if( (a[0] == 0) & (b[0] < tint) ){
             if(b[1] >= tint)
               if(results[num_d][1][a[1]-1] == temp1)  for(i = 0; i < b[0] + b[1]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);}
           else if(a[0] != 0)
             if(results[num_d][1][a[0]-1] == temp1)  for(i = 0; i < b[0]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);
        }
        else
        {
          if( (a[ia-1] == 0) & (b[ia-1] < tint) ){
            if(b[ia-2] >= tint)
              if(results[num_d][1][a[ia-2]-1] == temp2)
                for(i = (cols - b[ia-1] - b[ia-2]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);}
          else if(a[ia-1] != 0)
            if(results[num_d][1][a[ia-1]-1] == temp2)
              for(i = (cols - b[ia-1]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);
        }
      }
    }
    /*#################################################################*/
    else if(num_d == 3)
    {
      for(i = 0; i < 3; i++)
      {
        for(j = 0; j < cols; j++)
        {
          c[i][j] = arr[id][j] - results[num_d][1][i];
          c[i][j] *= c[i][j];
          c[i][j] *= -1;
          temp3 = 2 * results[num_d][2][i] * results[num_d][2][i];
          c[i][j] /= temp3;
          c[i][j] = exp(c[i][j]);
          c[i][j] *= results[num_d][0][i];
          c[i][j] /= results[num_d][2][i];
        }
      }
      for(i = 0; i < cols; i++)
      {
        tint = c[0][i] > c[1][i] ? 1 : 2;
        tint = c[tint-1][i] > c[2][i] ? tint : 3;
        temp3 = c[tint-1][i] / (c[0][i] + c[1][i] + c[2][i] - c[tint-1][i]);
        if(temp3 > 2)  ccc[i] = tint;
      }
      a[0] = ccc[t[0]];
      ia = 0;
      b[0] = 1;
      ib = 0;
      k = a[0];
      tint = floor(cols/20) + 1;
      for(i = 1; i < cols; i++)
      {
        if(ccc[t[i]] != k)
        {
          ia++;
          a[ia] = ccc[t[i]];        
          k = ccc[t[i]];
          ib++;
          b[ib] = 1;        
        }
        else  b[ib]++;
      }
      ia++;
      ib++;
      temp1 = results[num_d][1][0] < results[num_d][1][1] ? results[num_d][1][0] : results[num_d][1][1];
      temp1 = temp1 < results[num_d][1][2] ? temp1 : results[num_d][1][2];
      temp2 = results[num_d][1][0] > results[num_d][1][1] ? results[num_d][1][0] : results[num_d][1][1];
      temp2 = temp2 > results[num_d][1][2] ? temp2 : results[num_d][1][2];
      sum = 0;
      for(i = 0; i < ia; i++)  if(a[i]!=0) sum ++;
      if(sum != 0)
      {
        if( (a[0] == 0) & (b[0] < tint) )
          if(b[1] >= tint)
            if(results[num_d][1][a[1]-1] == temp1)
              for(i = 0; i < b[0] + b[1]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);
        if(a[0] != 0)
          if(results[num_d][1][a[0]-1] == temp1)
            for(i = 0; i < b[0]; i++)  arr_c[id][t[i]] = charset_add(symbols,-1);
        if( (a[ia-1] == 0) & (b[ia-1] < tint) )
          if(b[ia-2] >= tint)
            if(results[num_d][1][a[ia-2]-1] == temp2)
              for(i = (cols - b[ia-1] - b[ia-2]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);
        if(a[ia-1] != 0)
          if(results[num_d][1][a[ia-1]-1] == temp2)
            for(i = (cols - b[ia-1]); i < cols; i++)  arr_c[id][t[i]] = charset_add(symbols,1);
      }
    /*#################################################################*/
    }
    sum = 0;
    for(i = 0; i < cols; i++)
      if(arr_c[id][i] == 0)  sum++;
    if(sum == cols)
    {
      tint = ceil (cols*po->QUANTILE);
      for(i = 0; i < tint; i++)
      {
        arr_c[id][t[i]] = charset_add(symbols,-1);
        arr_c[id][t[cols-i-1]] = charset_add(symbols,1);
      }
	  fprintf(fw, "%s_unexpressed", genes_n[id]);
    }
	else  fprintf(fw, "%s", genes_n[id]);
	fprintf(fw, "\t");
    if(num_d == 2)
      fprintf(fw, "%lf_%lf\t%lf_%lf\n", results[2][1][0], results[2][2][0], results[2][1][1], results[2][2][1]);
    else if(num_d == 3)
      fprintf(fw, "%lf_%lf\t%lf_%lf\t%lf_%lf\n", results[3][1][0], results[3][2][0], results[3][1][1], results[3][2][1], results[3][1][2], results[3][2][2]);
    else
      fprintf(fw, "\n"); 
  }
  progress ("Discretization rules are written to %s", stream_nm); 
  fclose(fw);  
}


/* discretize the RPKM data into x distributions, x = 2, ..., 9*/
void R2to9(const char* stream_nm)
{
  double results[10][3][10], table_theta_t1[cols][10], m = 0, d = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, c[10][cols], cc[10], sumcut = 0, cut = 0, te[10], BIC[10], kk = 0; 
  int i = 0, j = 0, t[cols], tint = 0, INDEX = 0, k = 0, sum = 0, i_cut = 0;
  int UP = 9, DOWN = 1, num_d = 0;
  UP++;
  int EM = po->EM; /* parameter with default value being 20 or 150 */
  EM--;
  long long id = 0;
  FILE *fw;
  fw = mustOpen(stream_nm, "w");
  init_dis();
  for(id = 0; id < rows; id++)
  {
    for(i = 0; i < cols; i++)
    {
      t[i] = i; /* sort by natural numbers */
    }
    for(i = 0; i < cols; i++)
    {
      for(j = i; j < cols; j++)
      {
        if( arr[id][t[i]] > arr[id][t[j]] )
        {
          tint = t[i];
          t[i] = t[j];
          t[j] = tint;
        }
      }
    }
    cut = 0;
    sum = 0;
    for(i = 0; i < cols; i++)
    {
      sum++;
      if( arr[id][t[i]] != 0)
      {
        cut = log(arr[id][t[i]])/log(2) - 0.00001;
        break;
      }
    }
    for(i = 0; i < cols; i++)
    {
      if( arr[id][t[i]] != 0)  break;
      else  arr[id][t[i]] = pow( 2 ,(cut-2) );
    }
    for(i = 0; i < cols; i++)
      arr[id][t[i]] = log(arr[id][t[i]]) / log(2);
    for(num_d = DOWN; num_d < UP; num_d++)
    {
      for(i = 0; i < 10; i++)
      {
        results[num_d][0][i] = 0;
        results[num_d][1][i] = 0;
        results[num_d][2][i] = 0;
      }
      cc[num_d] = 0;
      BIC[num_d] = 0;
      m = 0;
      d = 0;
      sumcut = 0;
      for(i = 0; i < cols; i++)
      {
        if(arr[id][t[i]] < cut)
        {
          i_cut = i;
          sumcut += arr[id][t[i]];
        }
      }
      for(i = i_cut + 1; i < cols; i++)
        m += arr[id][t[i]]; /* the summation of one row */ 
      m /= (cols - i_cut - 1); /* the mean of one row */
      temp2 = 0;
      for(i = i_cut + 1; i < cols; i++)
      {
        /* the squared of the difference between the sample and the expectation */
        temp1 = arr[id][t[i]];
        temp1 -= m;
        temp1 *= temp1;
        temp2 += temp1;
      }
      temp2 /= (cols - i_cut - 2); /* unbiased estimated variance */
      d = sqrt(temp2); /* unbiased estimated standard deviation of one row */  
      for(j = 0; j < 10; j++)
        for(i = 0; i < 10; i++)
        {
          results[j][0][i] = 1; /* default weights */
          results[j][0][i] /= num_d;
          tint = (cols - i_cut - 1) * (i + 1) / (num_d + 1) - 1; 
          tint += i_cut;
          tint += 1;
          if(tint >= cols)  tint = cols - 1;
          results[j][1][i] = arr[id][t[tint]];
          /* Divide-and-Conquer */
          results[j][2][i] = d; /* default standard deviation */
        }
      for(INDEX = 0; INDEX < EM; INDEX++)
      {
        if( (cols - i_cut - 1) == 0)  sum = 999999;
        if(INDEX != 0)
          for(i = 0; i < num_d; i++)
          {
            if( results[num_d][2][i] == 0 )  sum = 999999;
            else if( results[num_d][0][i] > 999999999 )  sum = 999999;
            else if( results[num_d][1][i] > 999999999 )  sum = 999999;
            else if( results[num_d][2][i] > 999999999 )  sum = 999999;
            else if( results[num_d][0][i] < -999999999 )  sum = 999999;
            else if( results[num_d][1][i] < -999999999 )  sum = 999999;
            else if( results[num_d][2][i] < -999999999 )  sum = 999999;
          }
        if(sum == 999999)  break;
    	  /* ROUNDs 20 or ROUNDs 150 */
   	    for(i = i_cut + 1; i < cols; i++)
          for(j = 0; j < num_d; j++)
          {
            temp = densityFuction(arr[id][t[i]], results[num_d][1][j], results[num_d][2][j]);
            temp *= results[num_d][0][j];
            temp2 = 0;
            for(tint = 0; tint < num_d; tint++)
            {
              temp1 = densityFuction(arr[id][t[i]], results[num_d][1][tint], results[num_d][2][tint]);
              temp1 *= results[num_d][0][tint];
              temp2 += temp1;
            }
            temp /= temp2;
            table_theta_t1[i][j] = temp;
          }
        temp = 0;
        for(i = 0; i < num_d; i++)
        {
          cc[i] = NormSDist(cut, results[num_d][1][i], results[num_d][2][i]);
          temp += cc[i];
        }
        for(i = 0; i < num_d; i++)  {cc[i] /= temp; cc[i] *= (i_cut + 1);}
        temp3 = 0;
        for(i = 0; i < num_d; i++)
        {
          temp = 0;
          for(j = i_cut + 1; j < cols; j++)
          {
            temp += table_theta_t1[j][i];
          }
          temp2 = temp;
          temp += cc[i];
          temp /= cols;
          te[i] = temp;
          temp3 += temp;
          temp = 0;
          for(j = i_cut + 1; j < cols; j++)
          {
            temp += ( arr[id][t[j]] * table_theta_t1[j][i] );
          }
          temp1 = cut - results[num_d][1][i];
          temp1 /= results[num_d][2][i];
          temp1 = densityFuction(temp1, 0, 1);
          temp1 /= NormSDist(temp1, 0, 1);
          temp1 *= results[num_d][2][i];
          temp1 = results[num_d][1][i] - temp1;
          temp1 *= cc[i];
          temp += temp1;
          temp /= temp2;
          results[num_d][1][i] = temp; /* calculate */ 
          temp = cut - results[num_d][1][i];
          temp /= results[num_d][2][i];
          temp1 = densityFuction(temp, 0, 1); 
          temp1 /= NormSDist(temp, 0, 1);
          temp1 *= (cut - results[num_d][1][i]);
          temp1 /= results[num_d][2][i];
          temp1 = 1 - temp1;
          temp1 *= results[num_d][2][i];
          temp1 *= results[num_d][2][i];
          temp1 *= cc[i];
          temp = 0;
          for(j = i_cut + 1; j < cols; j++)
          {
            temp1 = arr[id][t[j]];
            temp1 -= results[num_d][1][i];
            temp1 *= temp1;
            temp1 *= table_theta_t1[j][i];
            temp += temp1;
          }
          temp /= temp2;
          temp = sqrt(temp);
          results[num_d][2][i] = temp; /* calculate */ 
        }
        for(i = 0; i < num_d; i++)  results[num_d][0][i] =  te[i] / temp3; /* calculate */ 
      }
      /*#################################################################*/
      temp3 = 0;
      k = num_d * 2 + 1;
      for(i = i_cut + 1; i < cols; i++)
      {
        temp2 = 0;
        for(j = 0; j < num_d; j++)
        {
          temp1 = densityFuction(arr[id][i], results[num_d][1][j], results[num_d][2][j]) * results[num_d][0][j];
          temp2 += temp1;
        }
        temp3 += log(temp2); 
      }
      for(i = 0; i < i_cut + 1; i++)
      {
        temp2 = 0;
        for(j = 0; j < num_d; j++)
        {
          temp1 = ( 1 - NormSDist(arr[id][i], results[num_d][1][j], results[num_d][2][j]) ) * results[num_d][0][j];
          temp2 += temp1;
        }
        temp3 += log(temp2); 
      }
      kk = k * log(cols);
      temp3 *= 2;
      temp3 -= k;
      BIC[num_d] = temp3;
    }
    /*#################################################################*/
    temp = BIC[1];
    tint = 1;
    for(num_d = DOWN + 1; num_d < UP; num_d++)
      if(temp < BIC[num_d])
      {
        temp = BIC[num_d];
        tint = num_d;
      }
    num_d = tint;
    for(i = 0; i < num_d; i++)
      for(j = 0; j < cols; j++)
      {
        c[i][j] = arr[id][j] - results[num_d][1][i];
        c[i][j] *= c[i][j];
        c[i][j] *= -1;
        temp3 = 2 * results[num_d][2][i] * results[num_d][2][i];
        c[i][j] /= temp3;
        c[i][j] = exp(c[i][j]);
        c[i][j] *= results[num_d][0][i];
        c[i][j] /= results[num_d][2][i];
      }
    for(i = 0; i < cols; i++)
    {
      temp2 = 0;
      temp1 = c[0][i];
      arr_c[id][i] = charset_add(symbols,0);
      for(j = 0; j < num_d; j++)
      {
        temp2 += c[j][i];
        if(temp1 < c[j][i])
        {
          temp1 = c[j][i];
          tint = j + 1;
        }
      }
      temp2 -= temp1;
      if(temp2 == 0)  temp3 = 0;
      else  temp3 = temp1 / temp2;
      if(temp3 > 1.5)  arr_c[id][i] = charset_add(symbols,tint);
    }
    fprintf(fw, "%s", genes_n[id]);
    for(i = 0; i < num_d; i++)
    {
      fprintf(fw, "%lf_%lf\t", results[num_d][1][i], results[num_d][2][i]);
    }
      fprintf(fw, "\n");
  }
  progress ("Discretization rules are written to %s", stream_nm); 
  fclose(fw);
}

/* discretize the microarray data into x distributions, x = 2, ..., 9*/
void N2to9(const char* stream_nm)
{
  double results[10][3][10], table_theta_t1[cols][9], m = 0, d = 0, temp = 0, temp1 = 0, temp2 = 0, temp3 = 0, c[10][cols], cc[10], cut = 0, BIC[10], kk = 0; 
  int i = 0, j = 0, t[cols], tint = 0, INDEX = 0, k = 0, sum = 0;
  int UP = 9, DOWN = 1, num_d = 0;
  UP++;
  int EM = po->EM; /* parameter with default value being 20 or 150 */
  EM--;
  long long id = 0;
  FILE *fw;
  fw = mustOpen(stream_nm, "w");
  init_dis();
  for(id = 0; id < rows; id++)
  {
    for(i = 0; i < cols; i++)
    {
      t[i] = i; /* sort by natural numbers */
    }
    for(i = 0; i < cols; i++)
    {
      for(j = i; j < cols; j++)
      {
        if( arr[id][t[i]] > arr[id][t[j]] )
        {
          tint = t[i];
          t[i] = t[j];
          t[j] = tint;
        }
      }
    }
    cut = 0;
    for(i = 0; i < cols; i++)
      if( arr[id][t[0]] != 0)
        cut = log(arr[id][t[0]])/log(2) - 0.00001;
    for(i = 0; i < cols; i++)
    {
      if( arr[id][t[i]] != 0)  break;
      else  arr[id][t[i]] = pow( 2 ,(cut-2) );
    }
    for(i = 0; i < cols; i++)
      arr[id][t[i]] = log(arr[id][t[i]]) / log(2);
    for(num_d = DOWN; num_d < UP; num_d++)
    {
      for(i = 0; i < 10; i++)
      {
        results[num_d][0][i] = 0;
        results[num_d][1][i] = 0;
        results[num_d][2][i] = 0;
      }
      cc[num_d] = 0;
      BIC[num_d] = 0;
      m = 0;
      d = 0;
      for(i = 0; i < cols; i++)
        m += arr[id][t[i]]; /* the summation of one row */ 
      m /= cols; /* the mean of one row */
      temp2 = 0;
      for(i = 0; i < cols; i++)
      {
        /* the squared of the difference between the sample and the expectation */
        temp1 = arr[id][t[i]];
        temp1 -= m;
        temp1 *= temp1;
        temp2 += temp1;
      }
      temp2 /= (cols - 1); /* unbiased estimated variance */
      d = sqrt(temp2); /* unbiased estimated standard deviation of one row */  
      for(j = 0; j < 10; j++)
        for(i = 0; i < 10; i++)
        {
          results[j][0][i] = 1; /* default weights */
          results[j][0][i] /= num_d;
          tint = cols * (i + 1) / (num_d + 1) - 1; 
          if(tint >= cols)  tint = cols - 1;
          results[j][1][i] = arr[id][t[tint]];
          /* Divide-and-Conquer */
          results[j][2][i] = d; /* default standard deviation */
        }
      for(INDEX = 0; INDEX < EM; INDEX++)
      {
        if(INDEX != 0)
          for(i = 0; i < num_d; i++)
            if( results[num_d][2][i] == 0 )  sum = 999999;
        if(sum == 999999)  break;
    	  /* ROUNDs 20 or ROUNDs 150 */
   	    for(i = 0; i < cols; i++)
          for(j = 0; j < num_d; j++)
          {
            temp = densityFuction(arr[id][t[i]], results[num_d][1][j], results[num_d][2][j]);
            temp *= results[num_d][0][j];
            temp2 = 0;
            for(tint = 0; tint < num_d; tint++)
            {
              temp1 = densityFuction(arr[id][t[i]], results[num_d][1][tint], results[num_d][2][tint]);
              temp1 *= results[num_d][0][tint];
              temp2 += temp1;
            }
            temp /= temp2;
            table_theta_t1[i][j] = temp;
          }
        for(i = 0; i < num_d; i++)
        {
          temp = 0;
          for(j = 0; j < cols; j++)
          {
            temp += table_theta_t1[j][i];
          }
          temp /= cols;
          results[num_d][0][i] = temp; /* calculate */ 
          temp = 0;
          for(j = 0; j < cols; j++)
          {
            temp += ( arr[id][t[j]] * table_theta_t1[j][i] );
          }
          temp /= cols;
          temp /= results[num_d][0][i];
          results[num_d][1][i] = temp; /* calculate */ 
          temp = 0;
          for(j = 0; j < cols; j++)
          {
            temp1 = arr[id][t[j]];
            temp1 -= results[num_d][1][i];
            temp1 *= temp1;
            temp1 *= table_theta_t1[j][i];
            temp += temp1;
          }
          temp /= cols;
          temp /= results[num_d][0][i];
          temp = sqrt(temp);
          results[num_d][2][i] = temp; /* calculate */ 
        }
      }
      /*#################################################################*/
      temp3 = 0;
      k = num_d * 2 + 1;
      for(i = 0; i < cols; i++)
      {
        temp2 = 0;
        for(j = 0; j < num_d; j++)
        {
          temp1 = densityFuction(arr[id][i], results[num_d][1][j], results[num_d][2][j]) * results[num_d][0][j];
          temp2 += temp1;
        }
        temp3 += log(temp2); 
      }
      kk = k * log(cols);
      temp3 *= 2;
      temp3 -= k;
      BIC[num_d] = temp3;
    }
    /*#################################################################*/
    temp = BIC[1];
    tint = 1;
    for(num_d = DOWN + 1; num_d < UP; num_d++)
      if(temp < BIC[num_d])
      {
        temp = BIC[num_d];
        tint = num_d;
      }
    num_d = tint;
    for(i = 0; i < num_d; i++)
      for(j = 0; j < cols; j++)
      {
        c[i][j] = arr[id][j] - results[num_d][1][i];
        c[i][j] *= c[i][j];
        c[i][j] *= -1;
        temp3 = 2 * results[num_d][2][i] * results[num_d][2][i];
        c[i][j] /= temp3;
        c[i][j] = exp(c[i][j]);
        c[i][j] *= results[num_d][0][i];
        c[i][j] /= results[num_d][2][i];
      }
    for(i = 0; i < cols; i++)
    {
      temp2 = 0;
      temp1 = c[0][i];
      arr_c[id][i] = charset_add(symbols,0);
      for(j = 0; j < num_d; j++)
      {
        temp2 += c[j][i];
        if(temp1 < c[j][i])
        {
          temp1 = c[j][i];
          tint = j + 1;
        }
      }
      temp2 -= temp1;
      if(temp2 == 0)  temp3 = 0;
      else  temp3 = temp1 / temp2;
      if(temp3 > 1.5)  arr_c[id][i] = charset_add(symbols,tint);
    }
    fprintf(fw, "%s", genes_n[id]);
    for(i = 0; i < num_d; i++)
    {
      fprintf(fw, "%lf_%lf\t", results[num_d][1][i], results[num_d][2][i]);
    }
      fprintf(fw,"\n");
  }
  progress ("Discretization rules are written to %s", stream_nm); 
  fclose(fw);
}
