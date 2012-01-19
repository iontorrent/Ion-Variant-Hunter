
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "fasta-io.h"
#include "zutil.h"

#include "posterior_flow.h"
#define EC_CEILING 500.0
static FILE *outp = stdout;
static FILE *logfile = stderr;
static int num_error = 0;
static int cur_pos;
static int debug = 0;

double poisson (double x, double Lam)
{
	if (Lam > EC_CEILING) {
	    x *= EC_CEILING/Lam;
	    Lam = EC_CEILING;
	}
        if (x < 0.0)
        {
                return 0.0;

        }
 	if (x <= Lam) return 0;

        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        /*
         * float iFac = 1.0;
         * float LamPow = 1.0;
         */
        double                   ret = 0.0;
        double                  log_iFac = 0.0;
        double                  log_LamPow = 0.0;

        unsigned int    i, iMax = (unsigned int) x;
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

        for (i = 0; i <= iMax; i++)
        {

                /*
                 * if ((LamPow / iFac) < UNDER_LIMIT) { break;
                 * }
                 */
                //if ((log_LamPow - log_iFac) <= -20)
                //{
                  //      break;
                //}

                /*
                 * ret += LamPow/iFac;
                 */
//                ret += exp(log_LamPow - log_iFac);

                /*
                 * iFac *= (float) (i + 1);
                 * LamPow *= Lam;
                 */
                log_iFac += log((float) (i + 1));
                log_LamPow += log(Lam);

        }

        //return ret * (float) exp(-Lam);
        return -(log_LamPow-log_iFac-Lam)/10.0;
}


double flow_poster_prob_calc::calc_post_prob(int num, unsigned char *flow_order, int *flow_sig, int dir)
{
    int i, j;
    unsigned char *temp_r;
    int *temp_s;
    if (dir==-1) {
	temp_r = ref; ref = ref_rev;
	temp_s = ref_r; ref_r = ref_r_rev;
    }
    if (do_split) {
	temp[0] = 0;
	for (i = 1; i < num; i++) temp[i] = temp[i-1] + prob_ins(ref[0], ref_r[0], flow_order[i-1], flow_sig[i-1]);
	for (i = 0; i < rlen; i++) {
	    double p = temp[0];
	    temp[0] += prob_del(ref_r[i]);
	    for (j = 0; j < 4; j++){  split[j] = 0.0; c_sign[j] = 0;}
	    for (j = 0; j < num; j++) {
		double m = p; 
		p += prob_match(ref[i],ref_r[i], flow_order[j], flow_sig[j]);
		double x = prob_ins(ref[i+1], ref_r[i+1], flow_order[j], flow_sig[j]);
		double q = temp[j]+x;
		if (q > p) q = p;
		p = temp[j+1];
		double q1 = p+prob_del(ref_r[i]);
		if (q1 < q) q = q1;
		if (ref[i] == flow_order[j]) {
		    q1 = split[flow_order[j]]+prob_match(ref[i], ref_r[i], flow_order[j], flow_sig[j]+c_sign[flow_order[j]]);
		    if (q1 < q) q = q1;
		}
		split[0] += x;
		split[1] += x;
		split[2] += x;
		split[3] += x;
		split[flow_order[j]] = m; 
		c_sign[flow_order[j]] = flow_sig[j];
		temp[j+1] = q;
	    }
	}

	if (dir == -1) {
	    ref = temp_r; ref_r = temp_s;
	}
	return temp[num];
    }
    temp[0] = 0;
    for (i = 1; i <= rlen; i++) temp[i] = temp[i-1]+prob_del(ref_r[i-1]);
    if (debug) {
	for (i = 0; i < rlen; i++) printf("%d:%d ", ref[i],ref_r[i]);
    	printf("read ");
    	for (i = 0; i < num; i++) printf("%d:%d ", flow_order[i], flow_sig[i]);
    }
    double bscore = 1e20;
    int end_p = 0;
    ref[rlen] = 5;
    for (i = 0; i < num; i++) {
	double p = temp[0];
	if (i >=LF) temp[0] += prob_ins(ref[0], ref_r[0], flow_order[i], flow_sig[i]); //do local temp[0] is always zero.
	for (j = 0; j < rlen; j++) {
	    p += prob_match(ref[j],ref_r[j], flow_order[i], flow_sig[i]);
	    double q = temp[j]+prob_del(ref_r[j]);
	    if (q > p) q = p;
	    p = temp[j+1];
	    double x = prob_ins(ref[j+1], ref_r[j+1], flow_order[i], flow_sig[i]);
	    double q1 = p+x;
	    temp[j+1] = (q > q1)? q1:q;
	}
	if (i >= num-LF && temp[rlen] < bscore) { end_p = i;bscore = temp[rlen]; }//local
    }
        if (dir == -1) {
            ref = temp_r; ref_r = temp_s;
        }
    //return temp[rlen]; local 
    if (debug) printf("score %f end p %d\n", (float) bscore, end_p);
    return bscore; // local
}


double flow_poster_prob_calc::prob_match(int ref_len, int flow_len)
{
    double x = abs(ref_len*100-flow_len)+MINI_PEN;
    double y = length_factor(ref_len);

    //if (x <= 60 && x >=40) {
//	return 3.0;
  //  }
    if (x < 50) {
    	x = -log10(1.0- pow(x/64.5, 3.0))*10.0;
    	y = y/ length_factor(1);
    }
    /*    
    if (flow_len > 500) {
        return 5+x*y*.7;
    }
    */

    return x*y;
}

double flow_poster_prob_calc::prob_match(unsigned char refb, int ref_len, unsigned char flow_base, int flow_sign)
{
    if (refb == flow_base) return prob_match(ref_len, flow_sign);
    return prob_ins(refb, ref_len, flow_base, flow_sign)+prob_del(ref_len);
}

double flow_poster_prob_calc::prob_ins(unsigned char refbase, int len, unsigned char insbase, int flow_len)
{
    int x;
    if (refbase == insbase) {
	x = len*100-flow_len;
	if (x < flow_len) x = flow_len;
    } else x = flow_len;
    if (x < 50) 
        return -log10(1.0- pow(x/64.5, 3.0))*10.0;
    double y = length_factor(1);
    return x*y;
}

double flow_poster_prob_calc::best_hyp(flow_list *alter_hyp, flow_list *read_list, int cov, float vh)
{
    int num_hyp = alter_hyp->num_list();
    int best = -1;
    double bscore = 1e20;
    double sbscore = 1e20;
    int i, num_r = read_list->num_list();
    unsigned char *s;
    int *f;
    int len;
    if (num_hyp == 0) return 0.0;
    for (i = 0; i < num_hyp; i++) {
	int j;
	double prior;
	alter_hyp->get_list(i, s, f, len, prior);
	setRef(len, s, f);
	double score = 0.0;
    	unsigned char *s1;
    	int *f1;
    	int len1;
	if (debug) alter_hyp->output_hyp(i);
	for (j = 0; j < num_r; j++) {
	    double p1;
	    read_list->get_list(j, s1, f1, len1, p1);
	    double x = calc_post_prob(len1, s1, f1, (p1 >0.0)? 1:-1);
	    score += x;
	}
	score -= log10(prior)*10.0;

	if (score < bscore) {
	    sbscore = bscore;
	    bscore = score;
	    best = i;
	} else if (score < sbscore) {
	    sbscore  = score;
	}
    }
    if (debug) printf("best =%f  second=%f\n", bscore, sbscore);
    double ex =-(sbscore-bscore)/num_r/10.0;
    ex = pow(10.0, ex);
    if (num_r > cov) ex = ex*(cov+num_r);
    else ex = ex*cov;
    double xxx = poisson(num_r, ex);
    if (num_r > cov) xxx = xxx*cov/num_r;
    if (vh > 75.0 && (best == 0 || xxx < 1.0)) {
	best = 1;
	xxx = vh*2.0;
    }
    //xxx = -log10(xxx)*10.0;
    //printf("PPP=log odd p-value=%f %f %f %f %f\n", (float) (sbscore-bscore), (float) sbscore, (float) bscore, (float) xxx, (float) ex);
    if (best != 0) {
	//printf("PPP=Refe="); alter_hyp->output_hyp(0); 
    	//printf("PPP=Call="); alter_hyp->output_hyp(best);
	xxx = xxx*0.63;  // normalization
	fprintf(outp, "Bayesian_Score=%f\n", (float) xxx);
    } else {
	//printf("PPP1=Refe="); alter_hyp->output_hyp(0);
	fprintf(outp, "Bayesian_Score=0\n");
    }
}


double flow_poster_prob_calc::prob_del(int len)
{
    return DEL_PENALTY;
}


// Below are just some testing code, not intended for final tools.

static unsigned char code(char a)
{
    a = toupper(a);
    if (a == 'A') return 0;
    if (a == 'C') return 1;
    if (a == 'G') return 2;
    if (a == 'T') return 3;
    if (a == 'N') return 0;
    fatal("%c not ACGT\n", a);
}

static int getflow(unsigned char *f, int *fs, const char *s, int left, int right, char *seq, int co)
{
    int i;
    const char *ss = s+left;
    //printf("%s %d %d %d\n", seq, left, right, co);
    i = co-1;
    char a = *ss;
    char ac[] = "ACGT";
    int count = 0;
    while (i >= 0) {
	if (toupper(*ss) != toupper(a)) {
	    fs[i] = count;
	    count = 1;
	    f[i] = code(a);
	    i--;
	    a = *ss;
	} else count++;
	ss--;
    }
    ss = seq;
    i = co-1; a = ac[f[i]]; count = fs[i];
    while (*ss && *ss != '.') {
	if (toupper(*ss) != toupper(a)) {
            fs[i] = count;
            count = 1;
            f[i] = code(a);
            i++;
            a = *ss;
        } else count++;
	ss++;
    }
    ss = s+right;
    if (toupper(*ss) == toupper(a)) co--;
    while (co >= 0) {
        if (toupper(*ss) != toupper(a)) {
            fs[i] = count;
            count = 1;
            f[i] = code(a);
            i++;
	    co--;
            a = *ss;
        } else count++;
        ss++;
    }
    return i;
}

static unsigned char total_flow_order[5000];

static int count_of_flow(char *s)
{
    int i = 1;
    s++;
    while (*s) {
	if (*s != *(s-1)) i++;
	s++;
    }
    return i;
}

static int set_hyp(char *line, flow_list *hy, const char *s)  // 2 hypothesis first.
{
    int pos;
    char ref[1000], targ[1000];
    unsigned char f[100];
    int fs[100];
    //printf("PPP=%s", line);
    chomp(line);
    //fprintf(outp, "%s", line);
    sscanf(line, "%*s %d %*s %s %s", &pos, ref, targ);
    cur_pos = pos;
    pos--; // to o-offset
    int i = getflow(f, fs, s, pos-1, pos+strlen(ref), ref, CONTEXT);
    hy->add_list(f, fs, i, 1.0);
    i = getflow(f, fs, s, pos-1, pos+strlen(ref), targ, CONTEXT);
    int rv = 0;
    if (strlen(ref) == strlen(targ)) {
	hy->add_list(f, fs, i, 0.005);
    } else {
	hy->add_list(f, fs, i, 0.001);
	if (strlen(ref) == 2 && strlen(targ) ==1) {  // for deletion, try alternative SNP cases.
	    int y = strlen(ref)-strlen(targ);
	    fs[CONTEXT-1]+=y;
	    hy->add_list(f, fs, i, 0.005);
	    fs[CONTEXT-1] -= y;
	    fs[i-CONTEXT] += y;
	    hy->add_list(f, fs, i, 0.005);
	}  
    }
    int b= count_of_flow(targ);
    rv = b-1;
    return rv;
}

static void expand_flow(unsigned char *s, int ss, char *fl)
{
    int i;
    char *f = fl;
    for (i = 0 ; i < ss; i++, f++) {
	if (*f == 0) f = fl;
	s[i] = code(*f);
    }
}

static int add_read(FILE *fp, char *line, flow_list *rs, int adj) 
{
    char line2[10000];
    while (fgets(line2, sizeof line2, fp)) {
	if (strncmp(line2, "FLOW", 4)==0) break;
    } 
    int fs[10000], i = 0;
    char *s = line2;
    while (s = strchr(s, ',')) {
	s++;
	fs[i] = atoi(s);
    	i++;
    }
    fs[7] -= 100;
    if (line[0] == 'I') return 1;
    int x, y, dir;
    int c, fi, pos;
    char d, dd;
    sscanf(line, "%*s %*s %*s %d %c %*s %*s %*s %c:%d:%*c:%*d:%d", &pos, &d, &dd, &x, &fi);
    if (d == '-') dir = -1; else dir = 1;
    if (dir == -1) {
	char *ss = line, *s;
 	while (s = strchr(ss,'\t')) {
	    ss = s+1;
	}
	sscanf(ss, "%c:%d:%*c:%*d:%d", &dd, &x, &fi);
    }

    unsigned char d1 = code(dd);
    if (debug) printf("%c  %d %d %d %d %d\t", d, x, total_flow_order[x], fs[x], d1, fi);
    int lr = pos-cur_pos;
    if (dir == -1) lr = -lr;
    
    if (x !=7 && fi != fs[x] || d1 != total_flow_order[x]) {
	int t;
   	for (t =1; t < 6; t++) {
	    if (fs[x+t] == fi && total_flow_order[x+t] == d1) {
		x = x+t; break;
	    }
	    if (fs[x-t] == fi && total_flow_order[x-t] == d1) {
                x = x-t; break;
            }
	}
	if (t == 6) {/*fprintf(logfile, "WARNING, cannot find the correct flow intensities from %sWARNING. %s", line, line2);*/ num_error++; return 1;}
	if (debug) printf("%c  %d %d %d %d %d\t", d, x, total_flow_order[x], fs[x], d1, fi);
    }

    y = x+1;
    x--;
    if (dir == 1) c = CONTEXT+1+adj/2; else c = CONTEXT+1+adj;
    c += lr;
    while (x >= 7 && c > 0) {
	if (fs[x] > 60) c--;
	x--;
    }
    x++;
    if (dir == 1) c = CONTEXT+1+adj; else c = CONTEXT+1+adj/2;
    c -= lr;
    while (y < i && c > 0) {
	if (fs[y] > 60) c--;
	y++;
    }
    if (debug) {
	printf("%d %d %d %s", x, y, i, line);
    	for (i = x; i < y; i++) printf("%d:%d ", total_flow_order[i], fs[i]);
    	printf("\n");
    }
    rs->add_list(total_flow_order+x, fs+x, y-x, (double) dir);
    return 1;
}

static int get_cov(char *c)
{
    char *s = strstr(c, "Num-spanning-reads=");
    if (!s)  return 0;
    return atoi(s+19);
}

static void get_allele_freq(char *s, float &a, int &b)
{
    float x;
    int n;
    a = 0; b = 0;
    s = strstr(s, "Num-spanning-reads");
    if (!s) return;
    s = strchr(s, '=');
    n = atoi(s+1);
    s = strstr(s, "Num-variant-reads");
    if (!s) return;
    s = strchr(s, '=');
    b = atoi(s+1);
    a = (float) b / (float) n;
}

static int need_combine(char *s, char *t)
{
    char ref1[100], targ1[100];
    int pos1;
    char id1[100];
    sscanf(s, "%s %d %*s %s %s", id1, &pos1, ref1, targ1);
    char ref2[100], targ2[100];
    int pos2;
    char id2[100];
    sscanf(t, "%s %d %*s %s %s", id2, &pos2, ref2, targ2);
    if (strcmp(id1, id2) != 0) return 0;
    if (pos1 != pos2-1) return 0;
    if (strlen(targ1) != 1 || strlen(targ2) != 1) return 0;
    if (strlen(ref1) != strlen(ref2)+1) return 0;
    if (strcmp(ref1+1, ref2) != 0) return 0;
    float p1, p2;
    int n1, n2;
    get_allele_freq(s, p1, n1);
    get_allele_freq(t, p2, n2);
    if (p1 < p2) {
	if (p1 > 0.2) return 0;
	if (n2 < n1*3/2) return 0;
	return 1;
    } else {
	if (p2 > 0.2) return 0;
	if (n1 < n2*3/2) return 0;
	return 2;
    }
}

static void split(char *s, char *&s1, char *&s2, float &cov, int &ns_ref, int &nv)
{
    char *t = strstr(s, "Variant-freq");
    if (!t) fatal("No Variant-freq field in vcf line.\n");
    *t = 0;
    t = strchr(t+1, '=');
    cov = atof(t+1);
    t = strchr(t, ';');
    s1 = t;
    t = strstr(t, "Num-spanning-ref-reads");
    if (t) {
    	t = strchr(t+1, '=');
    	ns_ref = atoi(t+1);
    } else {
	ns_ref = 0;
	t = s1;
    }
    t = strstr(t, "Num-variant-reads");
    if (!t) fatal("No Num-variant-reads field in vcf line.\n");
    *t = 0;
    t = strchr(t+1, '=');
    nv = atoi(t+1);
    t = strchr(t, ';');
    s2 = t;
}

static void combine_str(char *s, char *t)
{
    char *s1, *s2;
    float cov1, cov2;
    int ns_ref1, ns_ref2, nv1, nv2;
    split(s, s1, s2, cov1, ns_ref1, nv1);
    split(t, s1, s2, cov2, ns_ref2, nv2);
    if (ns_ref1 < nv2) ns_ref1=nv2;
    sprintf(s, "%sVariant-freq=%f%sNum-spanning-ref-reads=%d;Num-variant-reads=%d%s", t, cov1+cov2, s1, ns_ref1-nv2, nv1+nv2, s2);
}

static float vhscore(char *line)
{
    char *s = strstr(line, "Score");
    if (!s) return -1.0;
    s = strchr(s, '=');
    float x = atof(s+1);
    return x;
}

main(int argc, char *argv[])
{
    flow_poster_prob_calc aa;
    unsigned char ref[10]={0,1,2,1,0,3,1,2,1,0};
    int flow[10] = {1,3,3,1,1,1,2,2,1,1};
    aa.setRef(10, ref, flow); 
    unsigned char f[10] = {0,1,2,1,0,3,1,2,1,0};
    int fm[10] = {100,300,300,100,100,100,200,200,100,100};
    double x = aa.calc_post_prob(10, f, fm);
    
    unsigned char f1[21] = {0,   1,   2,3,0,  1,2, 3,  0,1,2,  3,0, 1,  2, 3,0, 1, 2,3, 0};
    int          fm1[21] = {100,300,300,3,3,100,4,11,100,1,1,100,1,200,200,1,1,100,1,1,100};
    x = aa.calc_post_prob(21, f1, fm1);
    unsigned char f2[25] = {0,   1, 2, 3, 0,   1,  2, 3,0,  1,2, 3,  0,1,2,  3,0, 1,  2, 3,0, 1, 2,3, 0};
    int          fm2[25] = {100,300,0, 0, 100, 0, 300,3,3,100,4,11,100,1,1,100,1,200,200,1,1,100,1,1,100};
    x = aa.calc_post_prob(25, f2, fm2);
    unsigned char f3[21] = {0,   1,   2,3,0,  1,  2, 3,  0,1,2,  3,0, 1,  2, 3,0, 1, 2,3, 0};
    int          fm3[21] = {100,300,300,3,100,100,4,11,100,1,1,100,1,200,200,1,1,100,1,1,100};
    x = aa.calc_post_prob(21, f3, fm3);
    unsigned char ref1[11]={0,1,0,2,1,0,3,1,2,1,0};
    int flow1[11] = {1,3,1,3,1,1,1,2,2,1,1};
    aa.setRef(11, ref1, flow1);
    x = aa.calc_post_prob(25, f2, fm2);
    unsigned char ref2[11]={0,1,2,0,1,0,3,1,2,1,0};
    int flow2[11] = {1,3,3,1,1,1,1,2,2,1,1};

    aa.setRef(11, ref2, flow2);
    x = aa.calc_post_prob(21, f3, fm3);

    /*
    unsigned char ref4[4]={0,2,0,3};
    int flow4[4] = {1,1,1,1};
    unsigned char ref5[5]={0,2,1,0,3};
    int flow5[5] = {1,1,1,1,1};
    unsigned char ref6[5]={0,2,3,0,3};
    int flow6[5] = {1,1,1,1,1};
    unsigned char f5[12] = {0,1,2,3,0,1,2,3,0,1,2,3};
    int fm5[12] = {100,0,100,0,0,100,0,0,100,0,0,100};
    unsigned char f6[8] = {0,1,2,3,0,1,2,3};
    int fm6[8]={100,0,100,100,100,0,0,100};
    aa.setRef(4, ref4, flow4);
    x = aa.calc_post_prob(12,f5,fm5);
    printf("example1 base %f\n", (float) x);
    x = aa.calc_post_prob(8, f6,fm6);
    printf("example2 base %f\n", (float) x);
    aa.setRef(5, ref5, flow5);
    x = aa.calc_post_prob(12,f5,fm5);
    printf("example1 alter %f\n", (float) x);
    x = aa.calc_post_prob(8, f6,fm6);
    printf("ref5 against f6 %f\n", (float) x);
    aa.setRef(5, ref6, flow6);
    x = aa.calc_post_prob(8, f6,fm6);
    printf("example2 alter %f\n", (float) x);
    */
    flow_list *hy = new flow_list(), *rs = new flow_list();
    /*
    hy->add_list(ref4, flow4, 4, 1.0);
    hy->add_list(ref5, flow5, 5, 0.001);
    rs->add_list(f5, fm5, 12, 1.0);
    rs->add_list(f6,fm6, 8, 1.0);
    aa.best_hyp(hy, rs);  
    */

    if (argc < 4) {
	fprintf(stderr, "%s ref_file devfile flow_order_string output_vcffile [input_vcffile] [-log logfile]", argv[0]);
	exit(1);
    }
    FastaFile *fafile = new FastaFile(SEQTYPE_NT);
    if (!fafile->Open(argv[1],"r"))
        return 0;
    FastaSeq faseq;
    if (!fafile->Read(faseq))
    {
	fatal("cannot open \n");
    }
    const char *s = faseq.Sequence();
    FILE *fp = ckopen(argv[2], "r");
    char line[100000], last_line[100000];
    last_line[0] = 0;
    expand_flow(total_flow_order, 5000, argv[3]);
    outp = ckopen(argv[4], "w");
    int print_header = 1;
    int out_header = 1;
    if (argc > 5) {
	int i;
	for (i = 5; i < argc; i++) {
	    if (argv[i][0] == '-') {
		if (strcmp(argv[i]+1, "log") == 0) {
		    logfile = ckopen(argv[i+1], "w");
		} else if (strcmp(argv[i]+1, "debug")==0|| strcmp(argv[i]+1, "d")) {
		    debug = atoi(argv[i+1]);
		} else {
		    fprintf(stderr, "wrong option %s\n", argv[i]);
		    exit(1);
		}
		i++;
	    } else {
        	FILE *input_dev = ckopen(argv[i], "r");
		print_header = 0;
		while (fgets(line, sizeof line, input_dev)) {
	    	    if (line[0] == '#') {
			if (out_header && line[1] != '#') {
		   	    fprintf(outp, "##FILTER=<ID=Bayesian_Score,Description=\"Using Bayesian model to re-evaluate the quality of the variant prediction.\">\n");
		    	    out_header = 0;
		    	}  
                    	fprintf(outp, line);
	    	    } else break;
		}
	    }
	}
    }
    int diff, coverage;
    int need_check = 0;
    float vh;
    while (fgets(line, sizeof line, fp)) {
	if (line[0] == '#') { 
	    if (print_header) fprintf(outp, line); 
	    continue;
	}
	if (line[0] == '\n') break;
	if (out_header) {
	    out_header = 0;
	    fprintf(outp, "##FILTER=<ID=Bayesian_Score,Number=1,Type=float,Description=\"Using Bayesian model to re-evaluate the quality of the variant prediction.\">\n");
	}
	char sid[100];
	sscanf(line, "%s", sid);
	//fprintf(stderr, "seq %s\n", sid);
	while (strcmp(sid, faseq.Label())!=0) {
	   if (!fafile->Read(faseq)) {
        	fatal("cannot open \n");
    	    }
	    //fprintf(stderr, "from reference %s\n", faseq.Label());
	    s = faseq.Sequence();
	}
	if (need_check && (x = need_combine(last_line, line))) {
	    if (x == 1) {
		hy->reset();
		diff = set_hyp(line, hy, s);
		coverage = get_cov(line);
		combine_str(last_line, line);
		chomp(last_line);
	    }
	    need_check = 0;
	} else {
	    if (last_line[0] != 0) {
		fprintf(outp, last_line);
		vh = vhscore(last_line);
		aa.best_hyp(hy, rs, coverage, vh);
	    }
            hy->reset();
            rs->reset();
	    strcpy(last_line, line);
	    chomp(last_line);        
	    diff = set_hyp(line, hy, s);
            coverage = get_cov(line);
	    need_check = 1;
	}

  	//printf("coverage=%d\n", coverage);
	fgets(line, sizeof line, fp); // definition line skip
	while (fgets(line, sizeof line, fp)) {
	    if (line[0] == '\n') break;
	    if (line[0] == 'I') continue;
	    if (!add_read(fp, line, rs, diff)) exit(1);
	}
	//aa.best_hyp(hy, rs, coverage);
    }
    fprintf(outp, last_line);
    vh = vhscore(last_line);
    aa.best_hyp(hy, rs, coverage, vh);
    if (num_error > 0) {
	fprintf(stderr, "Detecting %d reads with wrong flow information that can not be resolved by rescorer. Please check log file.\n", num_error);
    }
}


