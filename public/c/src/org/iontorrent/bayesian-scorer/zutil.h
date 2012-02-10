/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

#ifndef __zutil_h__
#define __zutil_h__ 1

#include <stdarg.h>

char *reverse (const char *from, long long len, char *to, char *co);
char *reverse (const char *from, long long len, char *to);
char *reverse_seq(char *seq, long long len, char *comp);
char *reverse_seq(char *seq, long long len);
char *find_string(char *seq, const char *limit);
char *strsave(const char *s, long long l);
char *strsave(const char *s);
char *ckalloc(long long s);
FILE *ckopen(const char *, const char *);
char *process_name(char *name);
int chomp(char *s);
//void fatal(const char *s);
void fatal(const char *format,...);
int getassay(FILE *ftaqout, char *pm1, char *pm2, double &tm1, double &tm2, int &pos1, int &pos2);
int getassay(FILE *ftaqout, char *pm1, char *pm2, double &tm1, double &tm2, int &pos1, int &pos2, double &pset);
int getassay(FILE *ftaqout, char *pm1, char *pm2, char *p);
void init_scode (char *_scode);
void init_compl (char *_compl);
int read_one_matrix(int *a, int &cur, int &cur2, char c, char mismatch_neg, int n, int sym, char *order, char *line, FILE *fp);
int read_m(int *a, char *order, char *line, FILE *fp);
int find_string(FILE *fp, char *str, char *line);
void zero_diagnal(int *d, int *e, int num);
int banded_gapalign(char *seq, int **seq2, int len2, int threshold, int band,int align, int gpenalty, char *primer_seq);
int local_align(const char *seq1, const char *seq2, int len1, int len2, int ma, int mi, int ga, int align, int dcut, int ut);
int the_same_letter(char a, char b);
unsigned dwordmap(unsigned a, int w);
void delete_one_base(const char *s, char *d, int len, int n);
void printn();
void print_line(const char *c);
int compare_id(const char *id1, const char *id2);
int panel_number(char *id);
int MBofFile(const char *f);
long long filesize(const char *f);
int skip2end(FILE *fp, char *s, int ll, FILE *fp1);
void mini_cut_cal(int &mc, int a, int b);
void output_progress(char *f, double a);
void output_progress(char *f, double a, const char *n);
int count_symbol(char *, char);
int lenofread(char *);
int lenofread(char *, int);


#define AMBIG 100
#define CHECK_MEM(a) {if ((a)==NULL) {fatal("out of memory");}}
#define CHECK_DEL_ARRAY(a) {if ((a)) { delete [] a; a = NULL;}}
#define CHECK_DEL(a) {if ((a)) { delete a; a = NULL;}}
#define check_printf { if (outsize <0) {perror("Error Output Results"); exit(-1);}}

#endif /* !READSMAP_ZUTIL_H */
