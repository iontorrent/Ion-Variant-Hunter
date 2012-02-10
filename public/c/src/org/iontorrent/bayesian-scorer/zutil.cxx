/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "zutil.h"

long long filesize(const char *filename)
{
    FILE *fp = ckopen(filename, "rb");
    fseeko(fp, 0, SEEK_END);
    off_t s = ftello(fp);
    fclose(fp);
    return s;
} 

int MBofFile(const char *filename)
{
    long long s = filesize(filename);
    int fs = (int) (s  >> 20);
    return fs;
}

char *reverse (const char *from, long long len, char *to, char *_compl)
{
	const char *s;
	char *t;

	for (s = from+len-1, t = to; s >= from; --s, ++t)
		if ((*t = _compl[*s]) == 0)
			*t = 'N';
	*t = '\0';
	return to;
}  

char *reverse (const char *from, long long len, char *to)
{
        const char *s;
        char *t;

        for (s = from+len-1, t = to; s >= from; --s, ++t)
                *t = *s;
        *t = '\0';
        return to;
}

char *reverse_seq(char *seq, long long len, char *_compl)
{
    char *start, *end;
    start = seq; 
    end = start+len-1;
    while (start < end) {
	char t = _compl[*start];
	*start = _compl[*end];
	*end = t;
	start++;
	end--;
    }
    if (start == end) {
	*start = _compl[*start];
    }
    return seq;
}

char *reverse_seq(char *seq, long long len)
{
    char *start, *end;
    start = seq;
    end = start+len-1;
    while (start < end) {
        char t = *start;
        *start = *end;
        *end = t;
        start++;
        end--;
    }
    return seq;
}

char *find_string(char *s, const char *dl)
{
    int i = strlen(dl);
    char *t;

    for (t = s; *t; t++) {
	if (*t == *dl && strncmp(t, dl, i) == 0) return t;
    }
    return NULL;
}

char *ckalloc(long long s)
{
    char *a = new char[s];
    if (a == NULL) {
	printf("out of memory\n");
	exit(1);
    }
    return a;
}

FILE *ckopen(const char *file, const char *format)
{
    FILE *fp = fopen(file, format);
    if (fp == NULL) {
	fprintf(stderr, "Can not open %s %s\n", file, format);
	exit(1);
    }
    return fp;
}

char *strsave(const char *a)
{
    char *b;
    b = ckalloc(strlen(a)+1);
    strcpy(b, a);
    return b;
}

char *strsave(const char *s, long long l)
{
    char *a = new char[l+1];
    strncpy(a, s, l);
    a[l] = '\0';
    return a;
}


char *process_name(char *name) 
{
    char *b = strchr(name, '/');
    if (b) {
	*b = '.';
	b++;
    }
    int len = strcspn(name, " .:\n");
/*
    printf("%s %d\n", name, len);
*/
    name[len] = '\0';
/*
    printf("%s\n", name);
*/
/*
    char *a = strchr(name, '.');
    char *d = strchr(name,':');
    if (a || d) {
        if (a) *a = 0;
	if (d) *d = 0;
    } else {
	int i = strlen(name);
	if (name[i-1] == '\n') {
	    name[i-1] = '\0';
	}
    }
*/
    return b;
}

int  chomp(char *s)
{
    int i = strlen(s);
    if (s[i-1] == '\n') { s[i-1] = '\0'; i--;}
    return i;
}
/*
void fatal(const char *s)
{
    fprintf(stderr, s);
    exit(1);
}
*/

void fatal(const char *format,...) { 
  va_list argList;                   
  va_start(argList,format);
    char buf[1048];
    vsnprintf(buf, sizeof(buf), format, argList);
    fprintf(stderr, "%s", buf);
    va_end(argList);
    exit(1);
}


int getassay(FILE *fp, char *pm1, char *pm2, double &tm1, double &tm2, int &pos1, int &pos2, double &pset)
{
    char line[1000];
    int tr, t2;
    float s;
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, "ForSeq", 6) == 0) {
	    sscanf(line, "ForSeq    %s", pm1);
	} else if (strncmp(line, "RevSeq", 6) == 0) {
	    sscanf(line, "RevSeq    %s", pm2);
	} else if (strncmp(line, "ForTm", 5)==0) {
	    sscanf(line, "ForTm     %f", &s);
	    tm1 = s;
	} else if (strncmp(line, "RevTm", 5)==0) {
	    sscanf(line, "RevTm     %f", &s);
	    tm2 = s;
	    return 1;
	} else if (strncmp(line, "ForDims", 7)==0) {
	    sscanf(line, "ForDims   %d %d", &tr, &t2);
	    pos1 = t2;
	} else if (strncmp(line, "RevDims", 7)==0) {
	    sscanf(line, "RevDims   %d %d", &tr, &t2);
	    pos2 = t2;
	} else if (strncmp(line, "Pset", 4)==0) {
	    sscanf(line, "%*s %f", &s);
	    pset = s;
	}
    }
    return 0;
}

int getassay(FILE *fp, char *pm1, char *pm2, double &tm1, double &tm2, int &pos1, int &pos2)
{
    char line[1000];
    int tr, t2;
    float s;
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, "ForSeq", 6) == 0) {
	    sscanf(line, "ForSeq    %s", pm1);
	} else if (strncmp(line, "RevSeq", 6) == 0) {
	    sscanf(line, "RevSeq    %s", pm2);
	} else if (strncmp(line, "ForTm", 5)==0) {
	    sscanf(line, "ForTm     %f", &s);
	    tm1 = s;
	} else if (strncmp(line, "RevTm", 5)==0) {
	    sscanf(line, "RevTm     %f", &s);
	    tm2 = s;
	    return 1;
	} else if (strncmp(line, "ForDims", 7)==0) {
	    sscanf(line, "ForDims   %d %d", &tr, &t2);
	    pos1 = t2;
	} else if (strncmp(line, "RevDims", 7)==0) {
	    sscanf(line, "RevDims   %d %d", &tr, &t2);
	    pos2 = t2;
	}
    }
    return 0;
}

int getassay(FILE *fp, char *pm1, char *pm2, char *p)
{
    char line[1000];
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, "ForSeq", 6) == 0) {
	    sscanf(line, "ForSeq    %s", pm1);
	} else if (strncmp(line, "RevSeq", 6) == 0) {
	    sscanf(line, "RevSeq    %s", pm2);
	    return 1;
	} else if (strncmp(line, "ProbeSeq", 8) == 0) {
	    sscanf(line, "ProbeSeq  %s", p);
	}
    }
    return 0;
}

void init_scode (char *_scode)
{
		int i;
		for (i=0; i< 128; ++i)
			_scode[i] = AMBIG;

		_scode['A'] = 0;
		_scode['C'] = 1;
		_scode['G'] = 2;
		_scode['T'] = 3;
		_scode['a'] = 0;
		_scode['c'] = 1;
		_scode['g'] = 2;
		_scode['t'] = 3;
}

void init_compl (char *_compl)
{
	int i;
	for (i = 0; i < 128;i++) _compl[i] = 'X';
    _compl['A'] = 'T';
    _compl['C'] = 'G';
    _compl['G'] = 'C';
    _compl['T'] = 'A';
    _compl['a'] = 'T';
    _compl['c'] = 'G';
    _compl['g'] = 'C';
    _compl['t'] = 'A';
    
    /* ambiguity codes */
    _compl['B'] = 'V';	/* B = C, G or T */
    _compl['D'] = 'H';	/* D = A, G or T */
    _compl['H'] = 'D';	/* H = A, C or T */
    _compl['K'] = 'M';	/* K = G or T */
    _compl['M'] = 'K';	/* M = A or C */
    _compl['N'] = 'N';	/* N = A, C, G or T */
    _compl['R'] = 'Y';	/* R = A or G (purines) */
    _compl['S'] = 'S';	/* S = C or G */
    _compl['V'] = 'B';	/* V = A, C or G */
    _compl['W'] = 'W';	/* W = A or T */
    _compl['X'] = 'X';	/* X = A, C, G or T */
    _compl['Y'] = 'R';	/* Y = C or T (pyrimidines) */

    _compl['b'] = 'v';  /* B = C, G or T */
    _compl['d'] = 'h';  /* D = A, G or T */
    _compl['h'] = 'd';  /* H = A, C or T */
    _compl['k'] = 'm';  /* K = G or T */
    _compl['m'] = 'k';  /* M = A or C */
    _compl['n'] = 'n';  /* N = A, C, G or T */
    _compl['r'] = 'y';  /* R = A or G (purines) */
    _compl['s'] = 's';  /* S = C or G */
    _compl['v'] = 'b';  /* V = A, C or G */
    _compl['w'] = 'w';  /* W = A or T */
    _compl['y'] = 'r';  /* Y = C or T (pyrimidines) */

    _compl['.'] = '.';
    
}

// a is the array to store the matrix, cur and cur2 are the indices of the last
// matrix stored, and c is the direction of left or right, in the regular
// probe case sym is 1 and c does not matter.
// order is the mapping of ACGT to 0-3.
int read_one_matrix(int *a, int &cur, int &cur2, char c, char mismatch_neg, int n, int sym, char *order, char *line, FILE *fp)
{
    int i;
    int sign = 1;
     int gap  = 0;
     int ATorAA = 0;
    if (mismatch_neg) sign = -1;
    for (i = 0; i < 16; i++) a[i] = 1;
    a[0]=a[5]=a[10]=a[15]=0;
    int find = 0;
    if (c == 'l' || c == 'r') {
	ATorAA = 1;
	c = toupper(c);
    }
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, "Default",  7) == 0) {
	    int d, e;
	    sscanf(line+7, "%d %d", &d, &e);
	    int i;
	    for (i = 0; i < 16; i++) {
		a[i] = sign*e;
	    }
	    if (ATorAA) {
		a[3]=a[6]= a[9] = a[12] = sign*d;
	    } else {
		a[0]=a[5]=a[10]=a[15]=sign*d;
	    }
	} else if (strncmp(line, "Gap", 3) == 0) {
	    sscanf(line+3, "%d", &gap); 
	} else if (strncmp(line, "End", 3)== 0) {
	    find =1;
	    int j = n;
	    int b = cur;
	    if (b == 0 && j > 0) fatal("The first matrix is missing\n");
	    if (!sym && c == 'R') b = cur2;
	    a -= 16*(j-b);
	    while (b < j) {
		int i;
		for (i = 0; i < 16; i++) {
		    a[i] = a[i-16];
		}
		a+= 16;
		b++;
	    }
	    if (!sym && c == 'R') { cur2 = n+1;}
	    else {cur = n+1;}
	    break;
	} else {
	    char tmp[100];
	    int d;
	    sscanf(line, "%s %d", tmp, &d);
	    char c = tmp[0];
	    char c1 = tmp[1];
	    //printf("%c %c\n", c, c1);
	    if (order[c] ==AMBIG || order[c1] == AMBIG) {
		fprintf(stderr, "Matrix Ingnore line \"%s\"", line);
	    } else {
		if (ATorAA) {
		    a[order[c]*4+order[c1]] = sign*d;
		} else {
		    a[order[c]*4+3-order[c1]]=sign*d;
		}
	    }
	}
    }
    if (!find) {fatal("A matrix is not ended");}
    return gap;
}

// just read in one matrix and stored it in a, line is a temparory
// memory.

int read_m(int *a, char *order, char *line, FILE *fp)
{
    int cur = 0;
    int n = 0;
    return read_one_matrix(a, cur, cur, 'L', 0, n, 1, order, line, fp);
}


int find_string(FILE *fp, char *str, char *line)
{  
    int i = strlen(str);
    while (fgets(line, 1000, fp)) {
	if (strncmp(line, str, i)==0) {
	    return 1; 
	}
    }
    return 0;
}

void zero_diagnal(int *d, int *e, int num)
{
    int each = 128;
    int i, j, k;
    for (i = 0; i < num; i++, d+=each*each, e+=each*each) {
	for (j = 0; j < each; j++) {
	    int same = d[each*j+j];
	    int same2 = e[each*j+j];
	    for (k = 0; k < each; k++) {
		d[each*j+k] -= same;
		e[each*j+k] -= same2;
	    } 
	}
    }
}

int the_same_letter(char a, char b)
{
    return (toupper(a) == toupper(b));
}

static int band2_galign(char *seq, int **psi, int length, int threshold)
{
    int a1, a2, a3,a4, a5, j;
    char *s;
    a1 = a2 = a3 = a4=a5= 0;
    for (j = 1; j <= length; j++) {
	int *pp = psi[j-1];
	int gap2 = pp['X']*2+1;
	int gap1;
	if (j == length) gap1 = gap2;
	else {
	    gap1 = psi[j]['X']*2+1;
	    if (gap1 > gap2) gap1 = gap2;
	}
	s = seq+j-1;
	int e;
	a1 += pp[*s++];
	e = a2+gap2;
	if (e < a1) a1 = e;
	else e = a1;
	a2 += pp[*s++];
	e+=gap1;
	int f = a3+gap2;
	if (e>f) e = f;
	if (e>a2) e = a2;
	else a2 = e;
	e+=gap1;
	a3 += pp[*s++];
	f = a4+gap2;
	if (e>f) e = f;
	if (e>a3) e = a3;
	else a3 = e;
	e+=gap1;
	a4 += pp[*s++];
	f = a5+gap2;
	if (e>f) e = f;
	if (e>a4) e = a4;
	else a4 = e;
	e += gap1;
	a5 += pp[*s];
	if (a5 > e) a5 = e;
	if (a3 > threshold && a4 > threshold && a2 > threshold && a5 > threshold && a1 > threshold) 
	    return -1;
    }
    if (a1 > a2) a1 = a2;
    if (a1 > a3) a1 = a3;
    if (a1 > a4) a1 = a4;
    if (a1 > a5) return a5;
    return a1;
}

static int band1_galign(char *seq, int **psi, int length, int threshold)
{
    int a1, a2, a3, j;
    char *s;
    a1 = a2 = a3 = 0;
    int gap_next = psi[0]['X']*2+1;
    for (j = 1; j < length; j++) {
	int *pp = psi[j-1];
	int gap2 = gap_next;
	int gap1;
	gap_next = psi[j]['X']*2+1;
	if (gap_next > gap2) gap1 = gap2;
	else gap1 = gap_next;
	s = seq+j-1;
	int e;
	a1 += pp[*s++];
	e = a2+gap2;
	if (e < a1) a1 = e;
	else e = a1;
	a2 += pp[*s++];
	e+=gap1;
	int f = a3+gap2;
	if (e>f) e = f;
	if (e>a2) e = a2;
	else a2 = e;
	e += gap1;
	a3 += pp[*s];
	if (a3 > e) a3 = e;
	if (a2 > threshold && a3 > threshold && a1 > threshold) return -1;
    }
    int *pp = psi[j-1];
    s = seq+j-1;
    a1 += pp[*s++];
    if (a2+gap_next < a1) a1 = a2+gap_next;
    a2+= pp[*s++];
    if (a3+gap_next < a2) a2 = a3+gap_next;
    a3 += pp[*s];
    if (a1 > a2) a1 = a2;
    if (a1 > a3) a1 = a3;
    if (a1 > threshold) return -1;
    return a1;
}

// seq is of length length+band
int banded_gapalign(char *seq, int **psi, int length, int threshold, int band, int align, int gpenalty, char *primer_seq)
{
    if (align == 0) {
	if (band == 3) return band1_galign(seq, psi, length, threshold);
	if (band == 5) return band2_galign(seq, psi, length, threshold);
    } 
    int i;
    char *s;
    double aa[(length+1)*(band+1)];
    for (i = 0; i < band; i++) aa[i] = (double) i/ 1000.0;
    aa[i] = 10000.0;
    int j;
    double min_row=0, *last, *cur;
    last = cur = aa;
    for (j = 1; j <= length; j++) {
	int *pp = psi[j-1];
	int gap2 = pp['X']*2+1;
	int gap1;
	if (j == length) gap1 = gap2;
	else {
	    gap1 = psi[j]['X']*2+1;
	    if (gap1 > gap2) gap1 = gap2;
	}
	s = seq+j-1;
	double e = 10000;
	min_row = e;
	if (align) {
	    cur = last+band+1;
	    last[band] = 10000;
	}
	for (i = 0; i < band; i++, s++) {
	    e += gap1;
	    double f = last[i+1]+gap2;
	    if (e > f) e = f;
	    double b = last[i]+pp[*s];
	    if (b < e) e = b;
	    cur[i] = e;
	    //fprintf(stderr, "%f ", e);
	    if (e < min_row) min_row = e;
	}
	//fprintf(stderr, "\n");
	if (((int)min_row) > threshold) return -1;
	last = cur;
    }
    for (i = 0; i < band; i++) {
	if (min_row ==  last[i]) {
	    j = i;
	}
    }
    double score = min_row;
    if (align) {
	i = j;
	j = length;
	char line1[1000], line2[1000], linem[1000];
	int p = 999;
	last = cur-band-1;
	line1[p] = line2[p] = linem[p] = 0;
	while (cur > aa) {
	    //fprintf(stderr, "%d %d\n", i,j);
	    p--;
	    linem[p] = ' ';
	    int *pp = psi[j-1];
	    int gap2 = pp['X']*2+1;
	    int gap1;
	    if (j == length) gap1 = gap2;
	    else {
		gap1 = psi[j]['X']*2+1;
		if (gap1 > gap2) gap1 = gap2;
	    }
	    if ((int)cur[i] == (int) last[i+1]+gap2) {
		line1[p] = '-';
		cur = last; last -= band+1; j--; 
		line2[p] = primer_seq[j]; i++;
	    } else if ((int) cur[i] == (int) cur[i-1]+gap1) {
		line1[p] = seq[j-1+i];
		line2[p] = '-';
		i--;
	    } else {
		line1[p] = seq[j-1+i];
		cur = last; last -= band+1; j--; 
		line2[p] = primer_seq[j];
		if (the_same_letter(line1[p], line2[p])) 
		    linem[p] = '|';
	    }
	}
	printf("%s\n%s\n%s\n", line1+p, linem+p, line2+p);
    }
    if (((int) score) <= threshold)
	return (int) score;
    else return -1;

}

//match >0 mismatch <0 gap>0 align=0/1
int local_align(const char *seq1, const char *seq2, int len1, int len2, int ma, int mi, int ga, int align, int display_cut, int upper_triangle)
{
    int i, j;
    int ss = 0;
    int e, d, u, best1=-1, best2=-1;
    int s[len1];
    int tb[len2][len1];
    for (i = 0; i < len1; i++) s[i] = 0;
    for (i = 0; i < len2; i++) {
	int code = 1;
	e = d=0;
	if (upper_triangle == 1) {
	    j = i+1;
	    d = s[i];
	} else j = 0;
	for (; j < len1; j++) {
	    code =1;
	    u = s[j];
	    if (e < u) {
		e = u;
		code = -1;
	    }
	    e -= ga;
	    if (the_same_letter(seq2[i], seq1[j])) d+=ma; else d+=mi;
	    if (e < d) {e = d; code = 0;}
	    if (e < 0) {
		e = 0; code =-2;
	    }
	    d = u;
	    s[j] = e;
	    tb[i][j] = code;
	    if (e > ss) {
		ss = e;
		best1 = i;
		best2 = j;
	    }
	}
    }
    
    if (align && ss > display_cut) {
	printf("score %d\n", ss);
	char line1[len1+len2];
	char line2[len1+len2];
	char linem[len1+len2];
	int p = len1+len2-1;
	line1[p] = line2[p] = linem[p] = 0;
	i = best1; j = best2;
	while (i >= 0 && j >=0) {
	    p--;
	    linem[p] = ' ';
	    if (tb[i][j] == 0) {
		line1[p] = seq2[i];
		line2[p] = seq1[j];
		if (the_same_letter(line1[p],line2[p])) linem[p] = '|';
		j--; i--;
	    } else if (tb[i][j] == -1) {
		line1[p] = seq2[i];
		line2[p] = '-';
		i--;
	    } else if (tb[i][j] == 1) {
		line1[p] = '-';
		line2[p] = seq1[j];
		j--;
	    } else { // -2 start of alignment, so display end here
		p++;
		break;
	    }
	}
	printf("%s\n%s\n%s\n", line1+p, linem+p, line2+p);
    }
    //printf("%d \n", ss);
    return ss;
} 

unsigned dwordmap(unsigned a, int wsize)
{
    //return ((a >>6) << 4 | (a & 15)); //1^(w-2)01^2
    //return ((a >>8) <<6 |(a & 63));  // 1^(w-3)01^3
    // return ((a >>4) << 2 | (a &3)); // 1^(w-1)01
    return ((a >> 10) << 8 | ( a & 255));
}

// delete the base n of sequence s. 0-based.
void delete_one_base(const char *s, char *d, int len, int n)
{
    if (n >= len) {
	strcpy(d, s);
	return;
    }
    strncpy(d, s, n);
    strcpy(d+n, s+n+1);
}


void printn()
{
    int outsize = printf("\n") ;
    check_printf;
}

void print_line(const char *line)
{
    int outsize = printf("%s", line);
    check_printf;
}

int compare_id(const char *id, const  char *id1)
{
    int a, b, c, d=0, e=0;
    int x,y , z, u=0, v=0;
    sscanf(id, ">%d_%d_%d_%d_%d", &a, &b, &c, &d, &e);
    sscanf(id1, ">%d_%d_%d_%d_%d", &x, &y, &z, &u, &v);

//    sscanf(id, ">%d_%d_%d", &a, &b, &c);
//    sscanf(id1, ">%d_%d_%d", &x, &y, &z);

    if (a < x) return -1;
    if (a > x) return 1;
    if (b < y) return -1;
    if (b > y) return 1;
    if (c > z) return 1;
    if (c < z) return -1;
    if (d > u) return 1;
     if (d < u) return -1;
    if (e > v) return 1;
     if (e < v) return -1;
    return 0;
}

int panel_number(char *id)
{
    int p;
    sscanf(id, ">%d", &p);
    return p;
}

int  skip2end(FILE *fp, char *s, int ll, FILE *fp1)
{
    do {
	if (fp1) fprintf(fp1, "%s", s);
	int i = strlen(s);
	if (s[i-1] == '\n') return 1;
    } while (fgets(s, ll, fp));
    return 0;
}

void mini_cut_cal(int &mc, int a, int b)
{ if (a < b) a = b; if (a < mc) mc = a;} // mc = min(sc, max(a, b)) 

int count_symbol(char *s, char a) {
    if (s == NULL || *s == 0) return 0; 
    int x = 1;
    while ((s = strchr(s, a))) {
	x++;
	s++;
    }
     return x;
}

void output_progress(char *progress_file, double a, const char *name)
{
    if (progress_file) {
        FILE *fp = ckopen(progress_file, "w");
        fprintf(fp, "%d,%s\n",(int)  (a*100.0), name);
        fclose(fp);
    } else {
        fprintf(stderr, "%d,%s\n",(int)  (a*100.0), name);
    }
}
void output_progress(char *progress_file, double a)
{
    output_progress(progress_file, a, "mapreads");
}
int lenofread(char *a)
{
    return lenofread(a, 1);
}

int lenofread(char *a, int colorspace)
{
    int i = strcspn(a, "\r\n")-colorspace;
    return i;
}
 
