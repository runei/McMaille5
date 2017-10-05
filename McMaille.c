/* ../../McMaille/McMaille.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "f2c.h"

/* Common Block Declarations */

union {
	struct {
	integer nhkl0, lhkl, ndat;
	real dmin__, slabda2;
	integer ihh[30000]	/* was [3][10000] */;
	real al[9]	/* was [3][3] */, pi, cri[10000], difp[10000], difm[
		10000], th2[10000], fobs[10000], sum_f__;
	integer nind;
	real w2[10000];
	integer nmx, ndat10;
	} _1;
	struct {
	integer nhkl0, lhkl, ndat;
	real dmin__, slabda2;
	integer ihh[30000]	/* was [3][10000] */;
	real al[9]	/* was [3][3] */, pi, cri[10000], difp[10000], difm[
		10000], th2[10000], fobs[10000], sum_f__;
	integer nind;
	real w2[10000];
	integer nhkl, ndat10;
	} _2;
} cal_;

#define cal_1 (cal_._1)
#define cal_2 (cal_._2)

struct {
	integer ind[1000000]	/* was [100][10000] */;
} cal2_;

#define cal2_1 cal2_

struct troc_1_ {
	integer iwr, irid;
};

#define troc_1 (*(struct troc_1_ *) &troc_)

struct {
	real qq[2000]	/* was [200][10] */, bb[10], b[10];
	integer h__[200], k[200], l[200], npaf;
	real afi[10];
	integer nr, indic;
	real pds[200];
	integer npaf2;
} truc_;

#define truc_1 truc_

/* Initialized data */

struct {
	integer e_1;
	integer fill_2[1];
	} troc_ = { 20 };


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__21 = 21;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__4 = 4;
static integer c__28 = 28;
static integer c__5 = 5;
static integer c__35 = 35;
static doublereal c_b1413 = 3e4;
static integer c__24 = 24;
static integer c__25 = 25;
static integer c__12 = 12;
static integer c__0 = 0;

int MAIN__(void)
{
	/* Format strings */
	static char fmt_1[] = "(\002  Entry file (no extension) ??\002,$)";
	static char fmt_2[] = "(a20)";
	static char fmt_5[] = "(a80)";
	static char fmt_3000[] = "(/1x,\002====================================="
		"========\r        ====================\002/1x,\002McMaille versi"
		"on \002,a,\002 by A. Le Bail - 2006 -\r         alb@cristal.or"
		"g\002/1x,\002=============================================\r    "
		"    ====================\002//1x,\002Using generic filename :"
		" \002,a)";
	static char fmt_8[] = "(20a4)";
	static char fmt_9660[] = "(\002! Wavelength, zeropoint, Ngrid\002)";
	static char fmt_9661[] = "(f9.6,2x,f7.4,\002 0\002)";
	static char fmt_9662[] = "(\002! Codes for symmetry\002)";
	static char fmt_9663[] = "(\002! W, Nind\002)";
	static char fmt_9664[] = "(f5.3,i3)";
	static char fmt_9[] = "(\002 Wavelength : \002,f9.5,\002 Zeropoint : "
		"\002,f8.4/)";
	static char fmt_10[] = "(\002 Width of the columnar profile shape, W  ="
		" \002,f9.4)";
	static char fmt_331[] = "(\002 Max non-indexed reflections, NIND  = \002"
		",i4)";
	static char fmt_9668[] = "(\002!Pmin, Pmax, Vmin, Vmax, Rmin, Rmax, Rmax"
		"ref\002)";
	static char fmt_9666[] = "(\002! Ntests, Nruns\002)";
	static char fmt_9667[] = "(\002!  2-theta   Intensity\002)";
	static char fmt_22[] = "(2x,f10.3,f9.4,2f10.3)";
	static char fmt_180[] = "(\002  Results in cubic, run, tests :\002,i3,f1"
		"2.0)";
	static char fmt_1115[] = "(14x,f5.3,f8.4,f9.1,i3)";
	static char fmt_115[] = "(f5.3,f12.0,f8.4,f9.1,2i3)";
	static char fmt_280[] = "(\002  Results in hexagonal, run, tests :\002,i"
		"3,f12.0)";
	static char fmt_281[] = "(\002  Results in rhombohedral, run, tests :"
		"\002,i3,f12.0)";
	static char fmt_1215[] = "(14x,f5.3,2f8.4,f9.1,i3)";
	static char fmt_215[] = "(f5.3,f12.0,2f8.4,f9.1,2i3)";
	static char fmt_380[] = "(\002  Results in tetragonal, run, tests :\002,"
		"i3,f12.0)";
	static char fmt_480[] = "(\002  Results in orthorhombic, run, tests :"
		"\002,i3,f12.0)";
	static char fmt_1415[] = "(14x,f5.3,3f8.4,f9.1,i3)";
	static char fmt_415[] = "(f5.3,f12.0,3f8.4,f9.1,2i3)";
	static char fmt_580[] = "(\002  Results in monoclinic, run, tests :\002,"
		"i3,f12.0)";
	static char fmt_1515[] = "(14x,f5.3,3f8.4,f7.2,f9.1,i3)";
	static char fmt_515[] = "(f5.3,f12.0,3f8.4,f7.2,f9.1,2i3)";
	static char fmt_680[] = "(\002  Results in triclinic, run, tests :\002,i"
		"3,f12.0)";
	static char fmt_1615[] = "(14x,f5.3,3f8.4,3f7.2,f9.1,i3)";
	static char fmt_615[] = "(f5.3,f12.0,3f8.4,3f7.2,f9.1,2i3)";
	static char fmt_7000[] = "(\002   M(20) = \002,f9.2)";
	static char fmt_7001[] = "(\002   F(20) = \002,f9.2,\002 (\002,f8.4"
		",\002,\002,i4,\002)\002)";
	static char fmt_19993[] = "(\002IN  F.o.M.    Volume   V/V1      a      "
		"  b        c       alpha   beta    gamma   Bravais lattice\002)";
	static char fmt_2002[] = "(i2,f8.2,f11.3,f6.2,39x,3f9.4,3f8.3)";
	static char fmt_20031[] = "(i2,f8.2,f11.3,f6.2,3x,3f9.4,3f8.3,4x,a1,2x,a"
		"11)";
	static char fmt_2069[] = "(i2,f8.2,f11.3,f6.2,1x,a1,1x,a12,1x,a7,1x,a8,i"
		"7,3f9.4,3f8.3,6f10.4)";
	static char fmt_2001[] = "(f8.3,f11.3,f6.2,2i4,3f9.4,3f8.3)";
	static char fmt_19994[] = "(\002IN  F(20)     Volume   V/V1      a      "
		"  b        c       alpha   beta    gamma   Bravais lattice\002)";
	static char fmt_19995[] = "(\002IN  M(20)     Volume   V/V1      a      "
		"  b        c       alpha   beta    gamma   Bravais lattice\002)";
	static char fmt_1999[] = "(\002    Rp     Vol     Vol/V1 Ind Nsol    a  "
		"      b     \r        c      alpha  beta  gamma\002)";
	static char fmt_19992[] = "(\002    Rp2    Vol     Vol/V1 Ind Nsol    a "
		"       b     \r        c      alpha  beta  gamma\002)";
	static char fmt_19996[] = "(\002IN  F.o.M.    Volume             a      "
		"  b        c       alpha   beta    gamma   Bravais lattice\002)";
	static char fmt_20032[] = "(i2,f8.2,f11.3,9x,3f9.4,3f8.3,4x,a1,2x,a11)";
	static char fmt_3955[] = "(//\002  Type any character and a RETURN to co"
		"ntinue : \002,$)";

	/* System generated locals */
	address a__1[2];
	integer i__1, i__2[2], i__3, i__4, i__5;
	real r__1;
	cllist cl__1;
	inlist ioin__1;

	/* Builtin functions */
	integer s_wsle(cilist *), e_wsle(void), do_lio(integer *, integer *, char 
		*, ftnlen), s_cmp(char *, char *, ftnlen, ftnlen);
	int s_copy(char *, char *, ftnlen, ftnlen);
	integer s_wsfe(cilist *), e_wsfe(void), s_rsfe(cilist *), do_fio(integer *
		, char *, ftnlen), e_rsfe(void), i_len(char *, ftnlen);
	int s_cat(char *, char **, integer *, integer *, ftnlen);
	integer f_inqu(inlist *), f_clos(cllist *), s_wsli(icilist *), e_wsli(
		void), s_rsle(cilist *), e_rsle(void);
	int s_stop(char *, ftnlen);
	double asin(doublereal), sin(doublereal), d_mod(doublereal *, doublereal *
		), log(doublereal), sqrt(doublereal), atan(doublereal), exp(
		doublereal);

	/* Local variables */
	static time_t time_end__;
	extern int cpu_time__(real *), esp_init__(integer *);
	static doublereal ntimelim[6];
	static logical pressedk;
	static char indxprog[12];
	static integer interest;
	static real a, b, c__, d__[10000];
	static integer i__, j, k, l;
	static real w, x;
	static integer i1, i2, i3, n1, n2;
	static real v1, w1[10000], v2, v3, x2, x3, totaltime;
	extern int open_read1__(integer *, char *, ftnlen);
	static real bb[8], dd;
	static integer ib[10000];
	static char bl[1];
	static real am, hh[3], ap;
	static integer na;
	static real bp, bm;
	static integer nb, im[10], kk, ll[10000], km[10000], ln, nm;
	static real ql[6], cy, hw[10000], qo[10000];
	static integer nr;
	static real rp[10000], vm;
	static integer ip;
	static real cp, vp, cm;
	static integer nc;
	static real vr;
	static integer jj;
	static real po;
	static integer jh, ii;
	static time_t time_begin__;
	static integer id1[100000], id2[100000], lc0, km2[10000], km3[100000], ll2[100000], ip2;
	static real th3[10000], xf1, hw4[10000], bbb[10000], ff20[10000], afi[8];
	static integer rp2[10000];
	static integer lca, nda;
	static real cel[60000]	/* was [6][10000] */, fm20[10000];
	static integer igc, nha[16000], ifi[10000], nhb[16000], jhh[30000]	/* 
		was [3][10000] */;
	static char nam[80];
	static real pma[6], vgc[10000];
	static integer igm, ihr;
	static doublereal igt;
	static real pmi[6];
	static integer imn[10], lll[10000];
	static real rmi;
	extern int open_write1__(integer *, char *, ftnlen);
	static real del;
	static integer ibr;
	static logical qex;
	static integer lpr;
	static real ddt, ddq, pos[16000], bet, ang, alp, gam;
	static integer jjj;
	static real sss, dmi;
	static integer icn, icz, nvk;
	static real yoy;
	static integer igc2;
	static real fcal[10000], diff, deld;
	static char fend[1], file[80];
	static integer ncel, jifi;
	static real bpar[6];
	static integer iref, ihkl[30000]	/* was [3][10000] */;
	extern int datn_(char *, char *, ftnlen, ftnlen);
	static real sang;
	static integer nhkl, ipen, isee;
	static char more[11];
	static real pmin;
	static integer nmax;
	static real dump[10000], pmax, rmax, spar, xfom[10000], vmin;
	static integer nsol[10000];
	static real yobs[16000], vmax, rmin, zero, tmax;
	extern int brav_(integer *, integer *, integer *);
	static integer nrun;
	static real vmon, text[20];
	static integer nout;
	static real betp, betm;
	static integer ifin, imem;
	static real diff2;
	extern int sort_(integer *, real *, integer *);
	static integer nsys[6];
	static real hmin, flog, slog, step;
	static integer npts;
	static real delt;
	static integer ifin1, ifin2, ndat2, ifin3;
	static real dmax1, dmax2, dmax3, omeg;
	static integer kpos;
	static real amda1, amda2;
	static integer npat1;
	static real rmax0[6], rmax2;
	static integer nrun2;
	extern int sort2_(integer *, integer *, integer *), 
		sort3_(integer *, integer *, integer *);
	static integer ncalc;
	extern int dcell_(real *, real *, real *);
	static integer icode, indic, ifile;
	static real delta[3];
	static integer iseed;
	static real ycalc[16000];
	static integer lfile;
	extern int progressview_(char *, ftnlen);
	static real deltc;
	extern doublereal randi_(integer *);
	static integer nglob;
	static real adelt;
	static integer ngrid;
	static real bdelt, theta[10000];
	static integer llhkl, iverb, irefs[10000];
	static real pndat;
	extern int killk_(logical *);
	static real rglob, cdelt, vdelt, tkill;
	static integer ibest;
	static real costh, tanth, astep;
	static char tempo[80];
	static real xnhkl, procs, vtric, sinth;
	static integer istep;
	static real omegt;
	static integer kiter;
	static real sum_y__, ttmax;
	static integer ivers;
	static real thmax, thmin, excrg, vorth;
	static integer imemt, nruns, isyst[7];
	static real cncalc[10000], slabda;
	static integer nruns2, nblack;
	static real deltab, deltad, escape;
	extern int celref_(integer *, real *, real *, integer *, 
		real *, integer *, real *, real *);
	static real celold[6];
	static integer iiseed, ipedig;
	extern int mcmnam_(integer *, char *, ftnlen);
	static char buffer[79];
	static real celpre[6], somega[10000], deltap, deltct[3];
	static integer ncells;
	static doublereal ntried;
	static integer nexcrg;
	static real timlim;
	extern int supcel_(integer *, integer *, real *, integer 
		*, real *, integer *);
	static integer notric, iprocs;
	static real astart, pstart[3];
	extern int celref2_(integer *, real *, real *, integer *,
		 real *, integer *, real *, real *), calcul1_(real *, real *), 
		calcul2_(real *, integer *, real *, integer *, integer *);
	static char select1[80];
	extern int filedel_(integer *, char *, ftnlen);
	static real betdelt;
	static doublereal ntriedb, ncycles;
	static char datenow[7];
	static real rmaxref;
	static doublereal ntriedt;
	static real rpsmall, pstartb[6];
	static integer nsolmax;
	static char timenow[8];
	static real astartt[3];

	/* Fortran I/O blocks */
	static cilist io___2 = { 0, 6, 0, 0, 0 };
	static cilist io___3 = { 0, 6, 0, 0, 0 };
	static cilist io___4 = { 0, 6, 0, 0, 0 };
	static cilist io___13 = { 0, 6, 0, fmt_1, 0 };
	static cilist io___14 = { 0, 5, 0, fmt_2, 0 };
	static cilist io___15 = { 0, 6, 0, 0, 0 };
	static cilist io___20 = { 0, 6, 0, 0, 0 };
	static cilist io___21 = { 0, 19, 1, fmt_5, 0 };
	static cilist io___23 = { 0, 21, 0, fmt_5, 0 };
	static icilist io___26 = { 0, buffer, 0, 0, 79, 1 };
	static icilist io___27 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___28 = { 0, 20, 0, fmt_3000, 0 };
	static cilist io___29 = { 0, 20, 0, 0, 0 };
	static cilist io___30 = { 0, 20, 0, 0, 0 };
	static cilist io___31 = { 0, 20, 0, 0, 0 };
	static cilist io___35 = { 1, 0, 0, fmt_8, 0 };
	static cilist io___38 = { 0, 20, 0, fmt_8, 0 };
	static cilist io___39 = { 1, 0, 0, 0, 0 };
	static cilist io___48 = { 0, 28, 0, fmt_8, 0 };
	static cilist io___49 = { 0, 28, 0, fmt_9660, 0 };
	static cilist io___50 = { 0, 28, 0, fmt_9661, 0 };
	static cilist io___51 = { 0, 28, 0, fmt_9662, 0 };
	static cilist io___52 = { 0, 28, 0, 0, 0 };
	static cilist io___54 = { 1, 0, 0, 0, 0 };
	static cilist io___56 = { 0, 28, 0, fmt_9663, 0 };
	static cilist io___57 = { 0, 28, 0, fmt_9664, 0 };
	static cilist io___58 = { 1, 0, 0, 0, 0 };
	static cilist io___59 = { 0, 20, 0, 0, 0 };
	static cilist io___60 = { 0, 20, 0, fmt_9, 0 };
	static cilist io___61 = { 0, 20, 0, 0, 0 };
	static cilist io___62 = { 0, 20, 0, 0, 0 };
	static cilist io___63 = { 0, 20, 0, 0, 0 };
	static cilist io___64 = { 0, 20, 0, 0, 0 };
	static cilist io___65 = { 0, 20, 0, 0, 0 };
	static cilist io___66 = { 0, 20, 0, 0, 0 };
	static cilist io___67 = { 0, 20, 0, 0, 0 };
	static cilist io___68 = { 0, 20, 0, 0, 0 };
	static cilist io___69 = { 0, 20, 0, 0, 0 };
	static cilist io___70 = { 0, 20, 0, fmt_10, 0 };
	static cilist io___71 = { 0, 20, 0, 0, 0 };
	static cilist io___72 = { 0, 20, 0, fmt_331, 0 };
	static cilist io___80 = { 0, 28, 0, fmt_9668, 0 };
	static cilist io___81 = { 0, 28, 0, 0, 0 };
	static cilist io___82 = { 1, 0, 0, 0, 0 };
	static cilist io___83 = { 1, 0, 0, 0, 0 };
	static cilist io___86 = { 0, 20, 0, 0, 0 };
	static cilist io___88 = { 0, 20, 0, 0, 0 };
	static cilist io___89 = { 0, 20, 0, 0, 0 };
	static cilist io___90 = { 0, 20, 0, 0, 0 };
	static cilist io___91 = { 0, 20, 0, 0, 0 };
	static cilist io___92 = { 0, 20, 0, 0, 0 };
	static cilist io___93 = { 0, 20, 0, 0, 0 };
	static cilist io___95 = { 0, 20, 0, 0, 0 };
	static cilist io___96 = { 0, 20, 0, 0, 0 };
	static cilist io___99 = { 0, 20, 0, 0, 0 };
	static cilist io___100 = { 0, 20, 0, 0, 0 };
	static cilist io___101 = { 0, 28, 0, fmt_9666, 0 };
	static cilist io___102 = { 0, 28, 0, 0, 0 };
	static cilist io___103 = { 0, 28, 0, fmt_9667, 0 };
	static cilist io___104 = { 0, 20, 0, 0, 0 };
	static cilist io___105 = { 0, 20, 0, 0, 0 };
	static cilist io___106 = { 1, 0, 0, 0, 0 };
	static cilist io___107 = { 0, 20, 0, 0, 0 };
	static cilist io___108 = { 0, 20, 0, 0, 0 };
	static cilist io___109 = { 1, 0, 0, 0, 0 };
	static cilist io___112 = { 0, 20, 0, 0, 0 };
	static cilist io___113 = { 0, 20, 0, 0, 0 };
	static cilist io___115 = { 0, 20, 0, 0, 0 };
	static cilist io___116 = { 0, 20, 0, 0, 0 };
	static cilist io___117 = { 0, 20, 0, 0, 0 };
	static cilist io___118 = { 0, 20, 0, 0, 0 };
	static cilist io___119 = { 0, 20, 0, 0, 0 };
	static cilist io___120 = { 0, 20, 0, 0, 0 };
	static cilist io___121 = { 0, 20, 0, 0, 0 };
	static cilist io___122 = { 0, 6, 0, 0, 0 };
	static cilist io___123 = { 0, 6, 0, 0, 0 };
	static cilist io___124 = { 0, 6, 0, 0, 0 };
	static cilist io___125 = { 0, 6, 0, 0, 0 };
	static cilist io___126 = { 0, 6, 0, 0, 0 };
	static cilist io___127 = { 0, 6, 0, 0, 0 };
	static cilist io___128 = { 0, 6, 0, 0, 0 };
	static cilist io___129 = { 0, 0, 1, 0, 0 };
	static cilist io___130 = { 0, 0, 1, 0, 0 };
	static cilist io___132 = { 0, 20, 0, 0, 0 };
	static cilist io___133 = { 0, 20, 0, 0, 0 };
	static cilist io___134 = { 0, 20, 0, 0, 0 };
	static cilist io___135 = { 0, 20, 0, 0, 0 };
	static cilist io___136 = { 0, 6, 0, 0, 0 };
	static cilist io___137 = { 0, 6, 0, 0, 0 };
	static cilist io___138 = { 0, 6, 0, 0, 0 };
	static cilist io___140 = { 0, 28, 0, 0, 0 };
	static cilist io___146 = { 0, 20, 0, 0, 0 };
	static cilist io___147 = { 0, 20, 0, 0, 0 };
	static cilist io___149 = { 0, 20, 0, fmt_22, 0 };
	static cilist io___150 = { 0, 20, 0, 0, 0 };
	static cilist io___162 = { 0, 20, 0, 0, 0 };
	static cilist io___163 = { 0, 20, 0, 0, 0 };
	static cilist io___164 = { 0, 20, 0, 0, 0 };
	static cilist io___165 = { 0, 20, 0, 0, 0 };
	static cilist io___166 = { 0, 20, 0, 0, 0 };
	static cilist io___167 = { 0, 6, 0, 0, 0 };
	static cilist io___168 = { 0, 6, 0, 0, 0 };
	static cilist io___169 = { 0, 6, 0, 0, 0 };
	static cilist io___170 = { 0, 6, 0, 0, 0 };
	static cilist io___171 = { 0, 6, 0, 0, 0 };
	static cilist io___172 = { 0, 20, 0, 0, 0 };
	static cilist io___173 = { 0, 20, 0, 0, 0 };
	static cilist io___174 = { 0, 20, 0, 0, 0 };
	static cilist io___175 = { 0, 20, 0, 0, 0 };
	static cilist io___176 = { 0, 20, 0, 0, 0 };
	static cilist io___177 = { 0, 6, 0, 0, 0 };
	static cilist io___178 = { 0, 6, 0, 0, 0 };
	static cilist io___179 = { 0, 6, 0, 0, 0 };
	static cilist io___180 = { 0, 6, 0, 0, 0 };
	static cilist io___181 = { 0, 6, 0, 0, 0 };
	static cilist io___182 = { 0, 20, 0, 0, 0 };
	static cilist io___183 = { 0, 20, 0, 0, 0 };
	static cilist io___193 = { 0, 20, 0, 0, 0 };
	static cilist io___194 = { 0, 20, 0, 0, 0 };
	static cilist io___195 = { 0, 20, 0, 0, 0 };
	static cilist io___196 = { 0, 20, 0, 0, 0 };
	static cilist io___197 = { 0, 20, 0, 0, 0 };
	static cilist io___198 = { 0, 20, 0, 0, 0 };
	static cilist io___199 = { 0, 20, 0, 0, 0 };
	static cilist io___201 = { 0, 6, 0, 0, 0 };
	static cilist io___202 = { 0, 6, 0, 0, 0 };
	static cilist io___203 = { 0, 6, 0, 0, 0 };
	static cilist io___208 = { 0, 20, 0, 0, 0 };
	static cilist io___209 = { 0, 20, 0, 0, 0 };
	static cilist io___210 = { 0, 20, 0, 0, 0 };
	static cilist io___211 = { 0, 20, 0, 0, 0 };
	static cilist io___212 = { 0, 20, 0, fmt_180, 0 };
	static cilist io___214 = { 0, 20, 0, 0, 0 };
	static cilist io___215 = { 0, 20, 0, 0, 0 };
	static cilist io___216 = { 0, 20, 0, 0, 0 };
	static cilist io___217 = { 0, 35, 0, 0, 0 };
	static cilist io___218 = { 0, 35, 0, 0, 0 };
	static cilist io___244 = { 0, 20, 0, 0, 0 };
	static cilist io___245 = { 0, 6, 0, 0, 0 };
	static cilist io___246 = { 0, 20, 0, 0, 0 };
	static cilist io___247 = { 0, 6, 0, 0, 0 };
	static cilist io___263 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___264 = { 0, 20, 0, 0, 0 };
	static cilist io___265 = { 0, 20, 0, 0, 0 };
	static cilist io___266 = { 0, 20, 0, 0, 0 };
	static cilist io___267 = { 0, 20, 0, 0, 0 };
	static cilist io___268 = { 0, 20, 0, 0, 0 };
	static cilist io___269 = { 0, 6, 0, 0, 0 };
	static cilist io___270 = { 0, 6, 0, 0, 0 };
	static cilist io___271 = { 0, 6, 0, 0, 0 };
	static cilist io___281 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___282 = { 0, 20, 0, fmt_115, 0 };
	static cilist io___283 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___284 = { 0, 20, 0, fmt_115, 0 };
	static cilist io___285 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___287 = { 0, 20, 0, 0, 0 };
	static cilist io___288 = { 0, 20, 0, 0, 0 };
	static cilist io___289 = { 0, 20, 0, 0, 0 };
	static cilist io___290 = { 0, 6, 0, 0, 0 };
	static cilist io___291 = { 0, 6, 0, 0, 0 };
	static cilist io___292 = { 0, 20, 0, 0, 0 };
	static cilist io___293 = { 0, 20, 0, 0, 0 };
	static cilist io___294 = { 0, 20, 0, 0, 0 };
	static cilist io___295 = { 0, 20, 0, 0, 0 };
	static cilist io___296 = { 0, 20, 0, fmt_280, 0 };
	static cilist io___297 = { 0, 20, 0, fmt_281, 0 };
	static cilist io___298 = { 0, 20, 0, 0, 0 };
	static cilist io___299 = { 0, 20, 0, 0, 0 };
	static cilist io___300 = { 0, 20, 0, 0, 0 };
	static cilist io___301 = { 0, 35, 0, 0, 0 };
	static cilist io___302 = { 0, 35, 0, 0, 0 };
	static cilist io___308 = { 0, 20, 0, 0, 0 };
	static cilist io___309 = { 0, 6, 0, 0, 0 };
	static cilist io___310 = { 0, 20, 0, 0, 0 };
	static cilist io___311 = { 0, 6, 0, 0, 0 };
	static cilist io___312 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___313 = { 0, 20, 0, 0, 0 };
	static cilist io___314 = { 0, 20, 0, 0, 0 };
	static cilist io___315 = { 0, 20, 0, 0, 0 };
	static cilist io___316 = { 0, 20, 0, 0, 0 };
	static cilist io___317 = { 0, 6, 0, 0, 0 };
	static cilist io___318 = { 0, 6, 0, 0, 0 };
	static cilist io___322 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___323 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___324 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___325 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___326 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___327 = { 0, 20, 0, 0, 0 };
	static cilist io___328 = { 0, 20, 0, 0, 0 };
	static cilist io___329 = { 0, 20, 0, 0, 0 };
	static cilist io___330 = { 0, 20, 0, 0, 0 };
	static cilist io___331 = { 0, 6, 0, 0, 0 };
	static cilist io___332 = { 0, 20, 0, 0, 0 };
	static cilist io___333 = { 0, 20, 0, 0, 0 };
	static cilist io___334 = { 0, 20, 0, 0, 0 };
	static cilist io___335 = { 0, 20, 0, 0, 0 };
	static cilist io___336 = { 0, 20, 0, fmt_380, 0 };
	static cilist io___337 = { 0, 20, 0, 0, 0 };
	static cilist io___338 = { 0, 20, 0, 0, 0 };
	static cilist io___339 = { 0, 20, 0, 0, 0 };
	static cilist io___340 = { 0, 35, 0, 0, 0 };
	static cilist io___341 = { 0, 35, 0, 0, 0 };
	static cilist io___342 = { 0, 20, 0, 0, 0 };
	static cilist io___343 = { 0, 6, 0, 0, 0 };
	static cilist io___344 = { 0, 20, 0, 0, 0 };
	static cilist io___345 = { 0, 6, 0, 0, 0 };
	static cilist io___346 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___347 = { 0, 20, 0, 0, 0 };
	static cilist io___348 = { 0, 20, 0, 0, 0 };
	static cilist io___349 = { 0, 20, 0, 0, 0 };
	static cilist io___350 = { 0, 20, 0, 0, 0 };
	static cilist io___351 = { 0, 6, 0, 0, 0 };
	static cilist io___352 = { 0, 6, 0, 0, 0 };
	static cilist io___353 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___354 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___355 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___356 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___357 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___358 = { 0, 20, 0, 0, 0 };
	static cilist io___359 = { 0, 20, 0, 0, 0 };
	static cilist io___360 = { 0, 20, 0, 0, 0 };
	static cilist io___361 = { 0, 20, 0, 0, 0 };
	static cilist io___362 = { 0, 6, 0, 0, 0 };
	static cilist io___364 = { 0, 20, 0, 0, 0 };
	static cilist io___365 = { 0, 20, 0, 0, 0 };
	static cilist io___366 = { 0, 20, 0, 0, 0 };
	static cilist io___367 = { 0, 20, 0, 0, 0 };
	static cilist io___368 = { 0, 20, 0, fmt_480, 0 };
	static cilist io___369 = { 0, 20, 0, 0, 0 };
	static cilist io___370 = { 0, 20, 0, 0, 0 };
	static cilist io___371 = { 0, 20, 0, 0, 0 };
	static cilist io___372 = { 0, 35, 0, 0, 0 };
	static cilist io___373 = { 0, 35, 0, 0, 0 };
	static cilist io___376 = { 0, 20, 0, 0, 0 };
	static cilist io___377 = { 0, 6, 0, 0, 0 };
	static cilist io___378 = { 0, 20, 0, 0, 0 };
	static cilist io___379 = { 0, 6, 0, 0, 0 };
	static cilist io___380 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___381 = { 0, 20, 0, 0, 0 };
	static cilist io___382 = { 0, 20, 0, 0, 0 };
	static cilist io___383 = { 0, 20, 0, 0, 0 };
	static cilist io___384 = { 0, 20, 0, 0, 0 };
	static cilist io___385 = { 0, 6, 0, 0, 0 };
	static cilist io___386 = { 0, 6, 0, 0, 0 };
	static cilist io___396 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___397 = { 0, 20, 0, fmt_415, 0 };
	static cilist io___398 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___399 = { 0, 20, 0, fmt_415, 0 };
	static cilist io___400 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___401 = { 0, 20, 0, 0, 0 };
	static cilist io___402 = { 0, 20, 0, 0, 0 };
	static cilist io___403 = { 0, 20, 0, 0, 0 };
	static cilist io___404 = { 0, 20, 0, 0, 0 };
	static cilist io___405 = { 0, 20, 0, 0, 0 };
	static cilist io___406 = { 0, 6, 0, 0, 0 };
	static cilist io___408 = { 0, 20, 0, 0, 0 };
	static cilist io___409 = { 0, 20, 0, 0, 0 };
	static cilist io___410 = { 0, 20, 0, 0, 0 };
	static cilist io___411 = { 0, 20, 0, 0, 0 };
	static cilist io___412 = { 0, 20, 0, fmt_580, 0 };
	static cilist io___413 = { 0, 20, 0, 0, 0 };
	static cilist io___414 = { 0, 20, 0, 0, 0 };
	static cilist io___415 = { 0, 20, 0, 0, 0 };
	static cilist io___416 = { 0, 35, 0, 0, 0 };
	static cilist io___417 = { 0, 35, 0, 0, 0 };
	static cilist io___420 = { 0, 20, 0, 0, 0 };
	static cilist io___421 = { 0, 6, 0, 0, 0 };
	static cilist io___422 = { 0, 20, 0, 0, 0 };
	static cilist io___423 = { 0, 6, 0, 0, 0 };
	static cilist io___424 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___425 = { 0, 20, 0, 0, 0 };
	static cilist io___426 = { 0, 20, 0, 0, 0 };
	static cilist io___427 = { 0, 20, 0, 0, 0 };
	static cilist io___428 = { 0, 20, 0, 0, 0 };
	static cilist io___429 = { 0, 6, 0, 0, 0 };
	static cilist io___430 = { 0, 6, 0, 0, 0 };
	static cilist io___434 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___435 = { 0, 20, 0, fmt_515, 0 };
	static cilist io___436 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___437 = { 0, 20, 0, fmt_515, 0 };
	static cilist io___438 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___439 = { 0, 20, 0, 0, 0 };
	static cilist io___440 = { 0, 20, 0, 0, 0 };
	static cilist io___441 = { 0, 20, 0, 0, 0 };
	static cilist io___442 = { 0, 20, 0, 0, 0 };
	static cilist io___443 = { 0, 20, 0, 0, 0 };
	static cilist io___444 = { 0, 20, 0, 0, 0 };
	static cilist io___445 = { 0, 6, 0, 0, 0 };
	static cilist io___447 = { 0, 20, 0, 0, 0 };
	static cilist io___448 = { 0, 20, 0, 0, 0 };
	static cilist io___449 = { 0, 20, 0, 0, 0 };
	static cilist io___450 = { 0, 20, 0, 0, 0 };
	static cilist io___451 = { 0, 20, 0, fmt_680, 0 };
	static cilist io___452 = { 0, 20, 0, 0, 0 };
	static cilist io___453 = { 0, 20, 0, 0, 0 };
	static cilist io___454 = { 0, 20, 0, 0, 0 };
	static cilist io___455 = { 0, 35, 0, 0, 0 };
	static cilist io___456 = { 0, 35, 0, 0, 0 };
	static cilist io___461 = { 0, 20, 0, 0, 0 };
	static cilist io___462 = { 0, 6, 0, 0, 0 };
	static cilist io___463 = { 0, 20, 0, 0, 0 };
	static cilist io___464 = { 0, 6, 0, 0, 0 };
	static cilist io___465 = { 0, 6, 0, fmt_1615, 0 };
	static cilist io___466 = { 0, 20, 0, 0, 0 };
	static cilist io___467 = { 0, 20, 0, 0, 0 };
	static cilist io___468 = { 0, 20, 0, 0, 0 };
	static cilist io___469 = { 0, 20, 0, 0, 0 };
	static cilist io___470 = { 0, 6, 0, 0, 0 };
	static cilist io___471 = { 0, 6, 0, 0, 0 };
	static cilist io___472 = { 0, 6, 0, fmt_1615, 0 };
	static cilist io___473 = { 0, 20, 0, fmt_615, 0 };
	static cilist io___474 = { 0, 6, 0, fmt_1615, 0 };
	static cilist io___475 = { 0, 20, 0, fmt_615, 0 };
	static cilist io___476 = { 0, 6, 0, fmt_1615, 0 };
	static cilist io___477 = { 0, 20, 0, 0, 0 };
	static cilist io___478 = { 0, 20, 0, 0, 0 };
	static cilist io___479 = { 0, 20, 0, 0, 0 };
	static cilist io___480 = { 0, 20, 0, 0, 0 };
	static cilist io___481 = { 0, 20, 0, 0, 0 };
	static cilist io___482 = { 0, 20, 0, 0, 0 };
	static cilist io___483 = { 0, 20, 0, 0, 0 };
	static cilist io___484 = { 0, 20, 0, 0, 0 };
	static cilist io___485 = { 0, 6, 0, 0, 0 };
	static cilist io___486 = { 0, 6, 0, 0, 0 };
	static cilist io___487 = { 0, 6, 0, 0, 0 };
	static cilist io___488 = { 0, 20, 0, 0, 0 };
	static cilist io___489 = { 0, 20, 0, 0, 0 };
	static cilist io___490 = { 0, 20, 0, 0, 0 };
	static cilist io___491 = { 0, 20, 0, 0, 0 };
	static cilist io___492 = { 0, 20, 0, 0, 0 };
	static cilist io___493 = { 0, 20, 0, 0, 0 };
	static cilist io___494 = { 0, 20, 0, 0, 0 };
	static cilist io___495 = { 0, 20, 0, 0, 0 };
	static cilist io___496 = { 0, 35, 0, 0, 0 };
	static cilist io___497 = { 0, 35, 0, 0, 0 };
	static cilist io___498 = { 0, 20, 0, 0, 0 };
	static cilist io___499 = { 0, 6, 0, 0, 0 };
	static cilist io___500 = { 0, 20, 0, 0, 0 };
	static cilist io___501 = { 0, 6, 0, 0, 0 };
	static cilist io___502 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___503 = { 0, 20, 0, 0, 0 };
	static cilist io___504 = { 0, 20, 0, 0, 0 };
	static cilist io___505 = { 0, 20, 0, 0, 0 };
	static cilist io___506 = { 0, 20, 0, 0, 0 };
	static cilist io___507 = { 0, 6, 0, 0, 0 };
	static cilist io___508 = { 0, 6, 0, 0, 0 };
	static cilist io___509 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___510 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___511 = { 0, 20, 0, 0, 0 };
	static cilist io___512 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___513 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___514 = { 0, 6, 0, 0, 0 };
	static cilist io___515 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___516 = { 0, 20, 0, fmt_115, 0 };
	static cilist io___517 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___518 = { 0, 20, 0, fmt_115, 0 };
	static cilist io___519 = { 0, 6, 0, fmt_1115, 0 };
	static cilist io___521 = { 0, 20, 0, 0, 0 };
	static cilist io___522 = { 0, 20, 0, 0, 0 };
	static cilist io___523 = { 0, 20, 0, 0, 0 };
	static cilist io___524 = { 0, 6, 0, 0, 0 };
	static cilist io___525 = { 0, 6, 0, 0, 0 };
	static cilist io___526 = { 0, 20, 0, 0, 0 };
	static cilist io___527 = { 0, 20, 0, 0, 0 };
	static cilist io___528 = { 0, 20, 0, 0, 0 };
	static cilist io___529 = { 0, 20, 0, 0, 0 };
	static cilist io___530 = { 0, 20, 0, 0, 0 };
	static cilist io___531 = { 0, 20, 0, 0, 0 };
	static cilist io___532 = { 0, 20, 0, 0, 0 };
	static cilist io___533 = { 0, 20, 0, 0, 0 };
	static cilist io___534 = { 0, 20, 0, 0, 0 };
	static cilist io___535 = { 0, 35, 0, 0, 0 };
	static cilist io___536 = { 0, 35, 0, 0, 0 };
	static cilist io___538 = { 0, 20, 0, 0, 0 };
	static cilist io___539 = { 0, 6, 0, 0, 0 };
	static cilist io___540 = { 0, 20, 0, 0, 0 };
	static cilist io___541 = { 0, 6, 0, 0, 0 };
	static cilist io___542 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___543 = { 0, 20, 0, 0, 0 };
	static cilist io___544 = { 0, 20, 0, 0, 0 };
	static cilist io___545 = { 0, 20, 0, 0, 0 };
	static cilist io___546 = { 0, 20, 0, 0, 0 };
	static cilist io___547 = { 0, 6, 0, 0, 0 };
	static cilist io___548 = { 0, 6, 0, 0, 0 };
	static cilist io___549 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___550 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___551 = { 0, 20, 0, 0, 0 };
	static cilist io___552 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___553 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___554 = { 0, 6, 0, 0, 0 };
	static cilist io___555 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___556 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___557 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___558 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___559 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___560 = { 0, 20, 0, 0, 0 };
	static cilist io___561 = { 0, 20, 0, 0, 0 };
	static cilist io___562 = { 0, 20, 0, 0, 0 };
	static cilist io___563 = { 0, 20, 0, 0, 0 };
	static cilist io___564 = { 0, 6, 0, 0, 0 };
	static cilist io___565 = { 0, 20, 0, 0, 0 };
	static cilist io___566 = { 0, 20, 0, 0, 0 };
	static cilist io___567 = { 0, 20, 0, 0, 0 };
	static cilist io___568 = { 0, 20, 0, 0, 0 };
	static cilist io___569 = { 0, 20, 0, 0, 0 };
	static cilist io___570 = { 0, 20, 0, 0, 0 };
	static cilist io___571 = { 0, 20, 0, 0, 0 };
	static cilist io___572 = { 0, 20, 0, 0, 0 };
	static cilist io___573 = { 0, 35, 0, 0, 0 };
	static cilist io___574 = { 0, 35, 0, 0, 0 };
	static cilist io___575 = { 0, 20, 0, 0, 0 };
	static cilist io___576 = { 0, 6, 0, 0, 0 };
	static cilist io___577 = { 0, 20, 0, 0, 0 };
	static cilist io___578 = { 0, 6, 0, 0, 0 };
	static cilist io___579 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___580 = { 0, 20, 0, 0, 0 };
	static cilist io___581 = { 0, 20, 0, 0, 0 };
	static cilist io___582 = { 0, 20, 0, 0, 0 };
	static cilist io___583 = { 0, 20, 0, 0, 0 };
	static cilist io___584 = { 0, 6, 0, 0, 0 };
	static cilist io___585 = { 0, 6, 0, 0, 0 };
	static cilist io___586 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___587 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___588 = { 0, 20, 0, 0, 0 };
	static cilist io___589 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___590 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___591 = { 0, 6, 0, 0, 0 };
	static cilist io___592 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___593 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___594 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___595 = { 0, 20, 0, fmt_215, 0 };
	static cilist io___596 = { 0, 6, 0, fmt_1215, 0 };
	static cilist io___597 = { 0, 20, 0, 0, 0 };
	static cilist io___598 = { 0, 20, 0, 0, 0 };
	static cilist io___599 = { 0, 20, 0, 0, 0 };
	static cilist io___600 = { 0, 20, 0, 0, 0 };
	static cilist io___601 = { 0, 6, 0, 0, 0 };
	static cilist io___602 = { 0, 20, 0, 0, 0 };
	static cilist io___603 = { 0, 20, 0, 0, 0 };
	static cilist io___604 = { 0, 20, 0, 0, 0 };
	static cilist io___605 = { 0, 20, 0, 0, 0 };
	static cilist io___606 = { 0, 20, 0, 0, 0 };
	static cilist io___607 = { 0, 20, 0, 0, 0 };
	static cilist io___608 = { 0, 20, 0, 0, 0 };
	static cilist io___609 = { 0, 20, 0, 0, 0 };
	static cilist io___610 = { 0, 35, 0, 0, 0 };
	static cilist io___611 = { 0, 35, 0, 0, 0 };
	static cilist io___613 = { 0, 6, 0, 0, 0 };
	static cilist io___614 = { 0, 20, 0, 0, 0 };
	static cilist io___615 = { 0, 6, 0, 0, 0 };
	static cilist io___616 = { 0, 20, 0, 0, 0 };
	static cilist io___617 = { 0, 6, 0, 0, 0 };
	static cilist io___618 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___619 = { 0, 20, 0, 0, 0 };
	static cilist io___620 = { 0, 20, 0, 0, 0 };
	static cilist io___621 = { 0, 20, 0, 0, 0 };
	static cilist io___622 = { 0, 20, 0, 0, 0 };
	static cilist io___623 = { 0, 6, 0, 0, 0 };
	static cilist io___624 = { 0, 6, 0, 0, 0 };
	static cilist io___625 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___626 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___627 = { 0, 20, 0, 0, 0 };
	static cilist io___628 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___629 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___630 = { 0, 6, 0, 0, 0 };
	static cilist io___631 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___632 = { 0, 20, 0, fmt_415, 0 };
	static cilist io___633 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___634 = { 0, 20, 0, fmt_415, 0 };
	static cilist io___635 = { 0, 6, 0, fmt_1415, 0 };
	static cilist io___636 = { 0, 20, 0, 0, 0 };
	static cilist io___637 = { 0, 20, 0, 0, 0 };
	static cilist io___638 = { 0, 20, 0, 0, 0 };
	static cilist io___639 = { 0, 20, 0, 0, 0 };
	static cilist io___640 = { 0, 20, 0, 0, 0 };
	static cilist io___641 = { 0, 6, 0, 0, 0 };
	static cilist io___642 = { 0, 20, 0, 0, 0 };
	static cilist io___643 = { 0, 20, 0, 0, 0 };
	static cilist io___644 = { 0, 20, 0, 0, 0 };
	static cilist io___645 = { 0, 20, 0, 0, 0 };
	static cilist io___646 = { 0, 20, 0, 0, 0 };
	static cilist io___647 = { 0, 20, 0, 0, 0 };
	static cilist io___648 = { 0, 20, 0, 0, 0 };
	static cilist io___649 = { 0, 20, 0, 0, 0 };
	static cilist io___650 = { 0, 35, 0, 0, 0 };
	static cilist io___651 = { 0, 35, 0, 0, 0 };
	static cilist io___654 = { 0, 6, 0, 0, 0 };
	static cilist io___655 = { 0, 20, 0, 0, 0 };
	static cilist io___656 = { 0, 6, 0, 0, 0 };
	static cilist io___657 = { 0, 20, 0, 0, 0 };
	static cilist io___658 = { 0, 6, 0, 0, 0 };
	static cilist io___659 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___660 = { 0, 20, 0, 0, 0 };
	static cilist io___661 = { 0, 20, 0, 0, 0 };
	static cilist io___662 = { 0, 20, 0, 0, 0 };
	static cilist io___663 = { 0, 20, 0, 0, 0 };
	static cilist io___664 = { 0, 6, 0, 0, 0 };
	static cilist io___665 = { 0, 6, 0, 0, 0 };
	static cilist io___666 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___667 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___668 = { 0, 20, 0, 0, 0 };
	static cilist io___669 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___670 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___671 = { 0, 6, 0, 0, 0 };
	static cilist io___672 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___673 = { 0, 20, 0, fmt_515, 0 };
	static cilist io___674 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___675 = { 0, 20, 0, fmt_515, 0 };
	static cilist io___676 = { 0, 6, 0, fmt_1515, 0 };
	static cilist io___677 = { 0, 20, 0, 0, 0 };
	static cilist io___678 = { 0, 20, 0, 0, 0 };
	static cilist io___679 = { 0, 20, 0, 0, 0 };
	static cilist io___680 = { 0, 20, 0, 0, 0 };
	static cilist io___681 = { 0, 20, 0, 0, 0 };
	static cilist io___682 = { 0, 20, 0, 0, 0 };
	static cilist io___684 = { 0, 20, 0, 0, 0 };
	static cilist io___685 = { 0, 20, 0, 0, 0 };
	static cilist io___686 = { 0, 20, 0, 0, 0 };
	static cilist io___687 = { 0, 20, 0, 0, 0 };
	static cilist io___688 = { 0, 20, 0, 0, 0 };
	static cilist io___689 = { 0, 20, 0, 0, 0 };
	static cilist io___690 = { 0, 20, 0, 0, 0 };
	static cilist io___691 = { 0, 20, 0, 0, 0 };
	static cilist io___692 = { 0, 20, 0, 0, 0 };
	static cilist io___693 = { 0, 20, 0, fmt_19993, 0 };
	static cilist io___694 = { 0, 20, 0, 0, 0 };
	static cilist io___708 = { 0, 24, 0, fmt_2002, 0 };
	static cilist io___709 = { 0, 20, 0, fmt_20031, 0 };
	static cilist io___710 = { 0, 25, 0, fmt_2069, 0 };
	static cilist io___711 = { 0, 20, 0, 0, 0 };
	static cilist io___712 = { 0, 20, 0, 0, 0 };
	static cilist io___713 = { 0, 20, 0, 0, 0 };
	static cilist io___714 = { 0, 20, 0, 0, 0 };
	static cilist io___715 = { 0, 20, 0, 0, 0 };
	static cilist io___716 = { 0, 20, 0, 0, 0 };
	static cilist io___717 = { 0, 20, 0, 0, 0 };
	static cilist io___718 = { 0, 20, 0, 0, 0 };
	static cilist io___719 = { 0, 20, 0, fmt_19994, 0 };
	static cilist io___720 = { 0, 20, 0, 0, 0 };
	static cilist io___721 = { 0, 20, 0, fmt_20031, 0 };
	static cilist io___722 = { 0, 20, 0, 0, 0 };
	static cilist io___723 = { 0, 20, 0, 0, 0 };
	static cilist io___724 = { 0, 20, 0, 0, 0 };
	static cilist io___725 = { 0, 20, 0, 0, 0 };
	static cilist io___726 = { 0, 20, 0, 0, 0 };
	static cilist io___727 = { 0, 20, 0, 0, 0 };
	static cilist io___728 = { 0, 20, 0, 0, 0 };
	static cilist io___729 = { 0, 20, 0, 0, 0 };
	static cilist io___730 = { 0, 20, 0, fmt_19995, 0 };
	static cilist io___731 = { 0, 20, 0, 0, 0 };
	static cilist io___732 = { 0, 20, 0, fmt_20031, 0 };
	static cilist io___733 = { 0, 20, 0, 0, 0 };
	static cilist io___734 = { 0, 20, 0, 0, 0 };
	static cilist io___735 = { 0, 20, 0, 0, 0 };
	static cilist io___736 = { 0, 20, 0, 0, 0 };
	static cilist io___737 = { 0, 20, 0, 0, 0 };
	static cilist io___738 = { 0, 20, 0, 0, 0 };
	static cilist io___739 = { 0, 20, 0, 0, 0 };
	static cilist io___740 = { 0, 20, 0, 0, 0 };
	static cilist io___741 = { 0, 20, 0, 0, 0 };
	static cilist io___742 = { 0, 20, 0, fmt_1999, 0 };
	static cilist io___745 = { 0, 20, 0, 0, 0 };
	static cilist io___746 = { 0, 20, 0, 0, 0 };
	static cilist io___747 = { 0, 20, 0, 0, 0 };
	static cilist io___748 = { 0, 20, 0, 0, 0 };
	static cilist io___749 = { 0, 20, 0, 0, 0 };
	static cilist io___750 = { 0, 20, 0, 0, 0 };
	static cilist io___751 = { 0, 20, 0, 0, 0 };
	static cilist io___752 = { 0, 20, 0, 0, 0 };
	static cilist io___753 = { 0, 20, 0, 0, 0 };
	static cilist io___756 = { 0, 20, 0, fmt_2001, 0 };
	static cilist io___757 = { 0, 24, 0, fmt_2002, 0 };
	static cilist io___759 = { 0, 20, 0, 0, 0 };
	static cilist io___760 = { 0, 20, 0, 0, 0 };
	static cilist io___761 = { 0, 20, 0, 0, 0 };
	static cilist io___762 = { 0, 20, 0, 0, 0 };
	static cilist io___763 = { 0, 20, 0, 0, 0 };
	static cilist io___764 = { 0, 35, 0, 0, 0 };
	static cilist io___765 = { 0, 35, 0, 0, 0 };
	static cilist io___766 = { 0, 35, 0, 0, 0 };
	static cilist io___767 = { 0, 35, 0, 0, 0 };
	static cilist io___768 = { 0, 35, 0, 0, 0 };
	static cilist io___769 = { 0, 35, 0, 0, 0 };
	static cilist io___770 = { 0, 35, 0, 0, 0 };
	static cilist io___771 = { 0, 35, 0, 0, 0 };
	static cilist io___772 = { 0, 35, 0, 0, 0 };
	static cilist io___773 = { 0, 35, 0, 0, 0 };
	static cilist io___774 = { 0, 35, 0, 0, 0 };
	static cilist io___775 = { 0, 35, 0, 0, 0 };
	static cilist io___776 = { 0, 35, 0, 0, 0 };
	static cilist io___777 = { 0, 35, 0, 0, 0 };
	static cilist io___778 = { 0, 20, 0, 0, 0 };
	static cilist io___779 = { 0, 20, 0, 0, 0 };
	static cilist io___780 = { 0, 20, 0, 0, 0 };
	static cilist io___781 = { 0, 20, 0, 0, 0 };
	static cilist io___783 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___784 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___785 = { 0, 20, 0, 0, 0 };
	static cilist io___786 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___787 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___788 = { 0, 6, 0, 0, 0 };
	static cilist io___828 = { 0, 20, 0, 0, 0 };
	static cilist io___829 = { 0, 20, 0, 0, 0 };
	static cilist io___830 = { 0, 20, 0, 0, 0 };
	static cilist io___832 = { 0, 12, 0, "(20A4)", 0 };
	static cilist io___833 = { 0, 12, 0, "(A,F9.4)", 0 };
	static cilist io___842 = { 0, 12, 0, "(2F8.3,F8.5,I8)", 0 };
	static cilist io___847 = { 0, 12, 0, "(I3,I7,5f10.5)", 0 };
	static cilist io___848 = { 0, 12, 0, "(8I5,8i3)", 0 };
	static cilist io___849 = { 0, 12, 0, "(8(F7.0,1X))", 0 };
	static cilist io___850 = { 0, 12, 0, "(8(F7.0,1X))", 0 };
	static cilist io___851 = { 0, 12, 0, "(8I10)", 0 };
	static cilist io___852 = { 0, 12, 0, "(10F8.3)", 0 };
	static cilist io___855 = { 0, 12, 0, "(I8)", 0 };
	static cilist io___856 = { 0, 12, 0, "(2F8.2)", 0 };
	static cilist io___857 = { 0, 20, 0, 0, 0 };
	static cilist io___858 = { 0, 20, 0, 0, 0 };
	static cilist io___859 = { 0, 20, 0, 0, 0 };
	static cilist io___860 = { 0, 20, 0, 0, 0 };
	static cilist io___861 = { 0, 20, 0, 0, 0 };
	static cilist io___862 = { 0, 20, 0, 0, 0 };
	static cilist io___863 = { 0, 20, 0, fmt_1999, 0 };
	static cilist io___864 = { 0, 20, 0, fmt_2001, 0 };
	static cilist io___865 = { 0, 20, 0, 0, 0 };
	static cilist io___866 = { 0, 20, 0, 0, 0 };
	static cilist io___867 = { 0, 20, 0, 0, 0 };
	static cilist io___868 = { 0, 20, 0, 0, 0 };
	static cilist io___869 = { 0, 20, 0, 0, 0 };
	static cilist io___870 = { 0, 35, 0, 0, 0 };
	static cilist io___871 = { 0, 35, 0, 0, 0 };
	static cilist io___872 = { 0, 35, 0, 0, 0 };
	static cilist io___873 = { 0, 35, 0, 0, 0 };
	static cilist io___874 = { 0, 35, 0, 0, 0 };
	static cilist io___875 = { 0, 35, 0, 0, 0 };
	static cilist io___876 = { 0, 35, 0, 0, 0 };
	static cilist io___877 = { 0, 35, 0, 0, 0 };
	static cilist io___878 = { 0, 35, 0, 0, 0 };
	static cilist io___879 = { 0, 35, 0, 0, 0 };
	static cilist io___880 = { 0, 35, 0, 0, 0 };
	static cilist io___881 = { 0, 35, 0, 0, 0 };
	static cilist io___882 = { 0, 35, 0, 0, 0 };
	static cilist io___883 = { 0, 35, 0, 0, 0 };
	static cilist io___884 = { 0, 20, 0, 0, 0 };
	static cilist io___885 = { 0, 20, 0, 0, 0 };
	static cilist io___886 = { 0, 20, 0, 0, 0 };
	static cilist io___887 = { 0, 20, 0, 0, 0 };
	static cilist io___888 = { 0, 20, 0, fmt_7000, 0 };
	static cilist io___889 = { 0, 20, 0, fmt_7001, 0 };
	static cilist io___890 = { 0, 20, 0, 0, 0 };
	static cilist io___891 = { 0, 6, 0, fmt_7000, 0 };
	static cilist io___892 = { 0, 6, 0, fmt_7001, 0 };
	static cilist io___893 = { 0, 6, 0, 0, 0 };
	static cilist io___894 = { 0, 20, 0, 0, 0 };
	static cilist io___895 = { 0, 20, 0, 0, 0 };
	static cilist io___896 = { 0, 20, 0, 0, 0 };
	static cilist io___897 = { 0, 20, 0, 0, 0 };
	static cilist io___898 = { 0, 20, 0, 0, 0 };
	static cilist io___899 = { 0, 20, 0, 0, 0 };
	static cilist io___900 = { 0, 20, 0, fmt_1999, 0 };
	static cilist io___901 = { 0, 20, 0, fmt_2001, 0 };
	static cilist io___902 = { 0, 20, 0, 0, 0 };
	static cilist io___903 = { 0, 20, 0, 0, 0 };
	static cilist io___904 = { 0, 20, 0, 0, 0 };
	static cilist io___905 = { 0, 20, 0, 0, 0 };
	static cilist io___906 = { 0, 20, 0, 0, 0 };
	static cilist io___907 = { 0, 20, 0, 0, 0 };
	static cilist io___908 = { 0, 20, 0, 0, 0 };
	static cilist io___909 = { 0, 20, 0, fmt_19992, 0 };
	static cilist io___910 = { 0, 20, 0, fmt_2001, 0 };
	static cilist io___911 = { 0, 20, 0, 0, 0 };
	static cilist io___912 = { 0, 20, 0, 0, 0 };
	static cilist io___913 = { 0, 20, 0, 0, 0 };
	static cilist io___914 = { 0, 20, 0, 0, 0 };
	static cilist io___915 = { 0, 20, 0, 0, 0 };
	static cilist io___916 = { 0, 20, 0, 0, 0 };
	static cilist io___917 = { 0, 20, 0, 0, 0 };
	static cilist io___918 = { 0, 20, 0, 0, 0 };
	static cilist io___919 = { 0, 20, 0, 0, 0 };
	static cilist io___920 = { 0, 20, 0, 0, 0 };
	static cilist io___927 = { 0, 20, 0, fmt_19992, 0 };
	static cilist io___928 = { 0, 20, 0, fmt_2001, 0 };
	static cilist io___929 = { 0, 24, 0, fmt_2002, 0 };
	static cilist io___930 = { 0, 20, 0, fmt_2001, 0 };
	static cilist io___931 = { 0, 24, 0, fmt_2002, 0 };
	static cilist io___932 = { 0, 20, 0, 0, 0 };
	static cilist io___933 = { 0, 20, 0, 0, 0 };
	static cilist io___934 = { 0, 20, 0, 0, 0 };
	static cilist io___935 = { 0, 20, 0, 0, 0 };
	static cilist io___936 = { 0, 20, 0, 0, 0 };
	static cilist io___937 = { 0, 20, 0, 0, 0 };
	static cilist io___938 = { 0, 20, 0, 0, 0 };
	static cilist io___939 = { 0, 20, 0, 0, 0 };
	static cilist io___940 = { 0, 20, 0, 0, 0 };
	static cilist io___941 = { 0, 20, 0, 0, 0 };
	static cilist io___942 = { 0, 20, 0, 0, 0 };
	static cilist io___943 = { 0, 20, 0, 0, 0 };
	static cilist io___944 = { 0, 20, 0, 0, 0 };
	static cilist io___945 = { 0, 20, 0, 0, 0 };
	static cilist io___946 = { 0, 20, 0, 0, 0 };
	static cilist io___947 = { 0, 20, 0, 0, 0 };
	static cilist io___948 = { 0, 20, 0, 0, 0 };
	static cilist io___949 = { 0, 20, 0, 0, 0 };
	static cilist io___950 = { 0, 20, 0, 0, 0 };
	static cilist io___951 = { 0, 20, 0, 0, 0 };
	static cilist io___952 = { 0, 20, 0, 0, 0 };
	static cilist io___953 = { 0, 20, 0, 0, 0 };
	static cilist io___954 = { 0, 20, 0, 0, 0 };
	static cilist io___955 = { 0, 20, 0, 0, 0 };
	static cilist io___956 = { 0, 20, 0, 0, 0 };
	static cilist io___957 = { 0, 20, 0, 0, 0 };
	static cilist io___958 = { 0, 20, 0, 0, 0 };
	static cilist io___959 = { 0, 20, 0, 0, 0 };
	static cilist io___960 = { 0, 20, 0, 0, 0 };
	static cilist io___961 = { 0, 20, 0, fmt_19996, 0 };
	static cilist io___963 = { 0, 20, 0, fmt_20032, 0 };
	static cilist io___964 = { 0, 20, 0, 0, 0 };
	static cilist io___967 = { 0, 20, 0, 0, 0 };
	static cilist io___968 = { 0, 20, 0, 0, 0 };
	static cilist io___969 = { 0, 20, 0, 0, 0 };
	static cilist io___970 = { 0, 20, 0, 0, 0 };
	static cilist io___971 = { 0, 20, 0, fmt_20032, 0 };
	static cilist io___972 = { 0, 20, 0, 0, 0 };
	static cilist io___973 = { 0, 20, 0, 0, 0 };
	static cilist io___974 = { 0, 20, 0, 0, 0 };
	static cilist io___975 = { 0, 20, 0, 0, 0 };
	static cilist io___976 = { 0, 20, 0, 0, 0 };
	static cilist io___977 = { 0, 20, 0, 0, 0 };
	static cilist io___978 = { 0, 20, 0, 0, 0 };
	static cilist io___979 = { 0, 20, 0, 0, 0 };
	static cilist io___980 = { 0, 20, 0, 0, 0 };
	static cilist io___981 = { 0, 20, 0, 0, 0 };
	static icilist io___982 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___983 = { 0, 20, 0, 0, 0 };
	static icilist io___984 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___985 = { 0, 20, 0, 0, 0 };
	static icilist io___986 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___987 = { 0, 20, 0, 0, 0 };
	static icilist io___988 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___989 = { 0, 20, 0, 0, 0 };
	static icilist io___990 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___991 = { 0, 20, 0, 0, 0 };
	static icilist io___992 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___993 = { 0, 20, 0, 0, 0 };
	static icilist io___994 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___995 = { 0, 20, 0, 0, 0 };
	static icilist io___996 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___997 = { 0, 20, 0, 0, 0 };
	static icilist io___998 = { 0, buffer, 0, 0, 79, 1 };
	static cilist io___1001 = { 0, 20, 0, 0, 0 };
	static cilist io___1002 = { 0, 6, 0, fmt_3955, 0 };
	static cilist io___1003 = { 0, 5, 0, 0, 0 };	



/*     Version 4.00 parallelized with OpenMP for multi-core processors */
/*       but this version is slightly modified for monoprocessors */

/*     MAILLE in french = CELL in english */
/*     Mc for Monte Carlo */
/*     Pronounce : MacMy */

/* *********************************************************************** */

/*     A Monte Carlo and grid search code for indexing powder patterns */

/*     For more details see the documentation at */
/*                   http://www.cristal.org/McMaille/ */
/*             or    http://sdpd.univ-lemans.fr/McMaille/ */

/*              by A. Le Bail - September 2002 for version 0.9 */
/*                              as well as for versions 1.0, 2.0 and 3.0 */
/*                              October 2006 for version 4.00 */
/*                        alb@cristal.org */
/*                        http://www.cristal.org/ */

/*                        Résidence Cristal - Appt 213 */
/*                        2, rue de Gasperi */
/*                        72100 Le Mans */
/*                        FRANCE */

/*   Versions 0.9 : cubic only */
/*            1.0 : hexagonal/trigonal/rhombohedral, tetragonal, */
/*                  orthorhombic added, plus .ckm and .prf files */
/*            2.0 : monoclinic and triclinic added in MC */
/*                  but not in grid search (too long) */
/*            3.0 : columnar peak shapes instead of Gaussian */
/*                  in versions 0.9-2.0 */
/*                  no Le Bail fit contrarily to versions 0.9-2.0 */
/*                  only fit by percentage of inclusion of the */
/*                  calculated column into the observed one */
/*            3.02: black box mode */
/*            3.03: improved Monte Carlo */
/*            3.04: two-phases mode */
/*            4.00: automatisation improved : more chances to identify */
/*                   the correct cell in "black box" mode */
/*                  Identification of the Bravais lattice */
/*                  Parallelization by using OpenMP directives */
/*                  improving the speed with multicore processors */
/*                  speed x1.7 to 1.8 with "dual core" or "core duo" */
/*                  speed x3.6 expected with the quad core in 2007 */
/*                  speed x79 expected with the 80-core in 2012...;-) */



/* *********************************************************************** */

/*    Copyright (C) 2002-2006 Armel Le Bail */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., */
/* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. */


/* ********************************************************************** */




/*      USE OMP_LIB */




/* $OMP THREADPRIVATE(/cal/,/cal2/) */

/*      CALL KMP_GET_STACKSIZE(ISIZE) */
/*      print *,ISIZE */
/*      ISIZE=200000000 */
/*      CALL KMP_SET_STACKSIZE(ISIZE) */
/*      CALL KMP_GET_STACKSIZE(ISIZE) */
/*      print *,ISIZE */

/*  Search for the number of processors available */

	iprocs = 1;
/*      IPROCS=OMP_GET_NUM_PROCS() */
	s_wsle(&io___2);
	e_wsle();
	s_wsle(&io___3);
	do_lio(&c__9, &c__1, "Number of used processors :", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&iprocs, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___4);
	e_wsle();
	procs = (real) iprocs;
	if (procs < 2.f) {
	procs = 1.f;
	}

	cal_1.pi = 114.59156f;
/* 360./3.1415926 */
	n1 = 1;
	n2 = 2;
	x = 0.f;

/* Open files using name from command line and standard extensions. */
/* The OPEN statements may need to be changed for some computers. */
/* Subroutine SXNM gets the generic filename from the command line. */
/* If nothing is found, then the user is prompted for the filename. */
	mcmnam_(&ln, nam, (ftnlen)80);
/* 	PRINT *,LN,NAM */
	i__1 = ln - 4;
	if (s_cmp(nam + i__1, ".exe", ln - i__1, (ftnlen)4) == 0) {
	goto L334;
	}
	s_copy(file, nam, (ftnlen)80, ln);
	printf("jshsosflkmklsf\n");
	lfile = ln;
	goto L335;

L334:
	s_wsfe(&io___13);
	e_wsfe();
	s_rsfe(&io___14);
	do_fio(&c__1, file, (ftnlen)80);
	e_rsfe();
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, file, (ftnlen)80);
	e_wsle();
	lfile = i_len(file, (ftnlen)80);
	while(*(unsigned char *)&file[lfile - 1] == ' ') {
	--lfile;
	}

L335:
	lc0 = 19;
	lca = 21;
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".inp";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	// snprintf(tempo, 80, "%s%s", file, ".inp");
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
	goto L3;
	}
	filedel_(&c__21, tempo, (ftnlen)80);
L3:
	open_write1__(&c__21, tempo, (ftnlen)80);
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".dat";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (qex) {
	goto L333;
	}
	s_wsle(&io___20);
	do_lio(&c__9, &c__1, "That file does not exist, try again...", (ftnlen)38)
		;
	e_wsle();
	goto L334;
L333:
	open_read1__(&c__19, tempo, (ftnlen)80);
L4:
	i__1 = s_rsfe(&io___21);
	if (i__1 != 0) {
	goto L6;
	}
	i__1 = do_fio(&c__1, select1, (ftnlen)80);
	if (i__1 != 0) {
	goto L6;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	goto L6;
	}
	if (*(unsigned char *)select1 != '!') {
	s_wsfe(&io___23);
	do_fio(&c__1, select1, (ftnlen)80);
	e_wsfe();
	}
	goto L4;
L6:
	cl__1.cerr = 0;
	cl__1.cunit = lc0;
	cl__1.csta = 0;
	f_clos(&cl__1);
	cl__1.cerr = 0;
	cl__1.cunit = lca;
	cl__1.csta = 0;
	f_clos(&cl__1);
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".inp";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	open_read1__(&c__21, tempo, (ftnlen)80);
	lpr = 20;
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".imp";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl  = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
	goto L7;
	}
	filedel_(&c__20, tempo, (ftnlen)80);
L7:
	open_write1__(&c__20, tempo, (ftnlen)80);
	s_wsli(&io___26);
	do_lio(&c__9, &c__1, "McMaille version ", (ftnlen)17);
	do_lio(&c__9, &c__1, "4.00    ", (ftnlen)8);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	s_wsli(&io___27);
	do_lio(&c__9, &c__1, "Data file : ", (ftnlen)12);
	do_lio(&c__9, &c__1, file, lfile);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	s_wsfe(&io___28);
	do_fio(&c__1, "4.00    ", (ftnlen)8);
	do_fio(&c__1, file, lfile);
	e_wsfe();
	s_wsle(&io___29);
	e_wsle();
	s_wsle(&io___30);
	do_lio(&c__9, &c__1, "  Number of Processors :", (ftnlen)24);
	do_lio(&c__3, &c__1, (char *)&iprocs, (ftnlen)sizeof(integer));
	e_wsle();
	s_wsle(&io___31);
	e_wsle();
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	time(&time_begin__);

/* .....READ PROBLEM IDENTIFICATION */

	io___35.ciunit = lca;
	i__1 = s_rsfe(&io___35);
	if (i__1 != 0) {
	goto L3500;
	}
	for (i__ = 1; i__ <= 20; ++i__) {
	i__1 = do_fio(&c__1, (char *)&text[i__ - 1], (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3500;
	}
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	goto L3500;
	}
	s_wsfe(&io___38);
	for (i__ = 1; i__ <= 20; ++i__) {
	do_fio(&c__1, (char *)&text[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();

/* .....READ wavelength and type of calculation */
/*     SLABDA = wavelength */
/*     ZERO   = zeropoint to be added at thr beginning */
/*     NGRID  if = 1  : grid cell generation */
/*            if = 0  : Monte Carlo cell generation */
/*            if = 2  : both */
/*            if = 3  : black box Monte Carlo only... */
/*            if = -3 : black box Monte Carlo only, without triclinic... */
/*            if = 4  : black box Monte Carlo + grid search... */

	io___39.ciunit = lca;
	i__1 = s_rsle(&io___39);
	if (i__1 != 0) {
	goto L3501;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&slabda, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	goto L3501;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&zero, (ftnlen)sizeof(real));
	if (i__1 != 0) {
	goto L3501;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&ngrid, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	goto L3501;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
	goto L3501;
	}
	iverb = 0;
	if (slabda < 0.f) {
	iverb = 1;
	slabda = -slabda;
	}
	bb[1] = slabda;
	bb[0] = 0.f;
/* zeropoint after correction... = 0. */
	afi[1] = 0.f;
/* code for wavelength refinement */
	afi[0] = 1.f;
/* code for zeropoint refinement */
	notric = 0;
	if (ngrid == -3) {
	notric = 1;
	ngrid = 3;
	}

/* .....READ codes for search in crystalline systems */

/*    NSYS(n) */
/*     n */
/*     1   Cubic */
/*     2   Hexagonal/trigonal/Rhombohedral */
/*     3   Tetragonal */
/*     4   Orthorhombic */
/*     5   Monoclinic */
/*     6   Triclinic */

/*     if NSYS(n)=0 : no search */
/*     if NSYS(n)=1 : search */
/*     if NSYS(2)=2 : search in rhombohedral */

	nblack = 0;
	if (ngrid == 4) {
	nblack = 1;
	}
	if (ngrid == 4) {
	ngrid = 3;
	}
	if (ngrid == 3) {
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "-new.dat";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L9325;
	}
	filedel_(&c__28, tempo, (ftnlen)80);
L9325:
	open_write1__(&c__28, tempo, (ftnlen)80);
	s_wsfe(&io___48);
	for (i__ = 1; i__ <= 20; ++i__) {
		do_fio(&c__1, (char *)&text[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___49);
	e_wsfe();
	s_wsfe(&io___50);
	do_fio(&c__1, (char *)&slabda, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&zero, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___51);
	e_wsfe();
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, "1 0 0 0 0 0", (ftnlen)11);
	e_wsle();
	for (i__ = 1; i__ <= 6; ++i__) {
/* L9300: */
		nsys[i__ - 1] = 1;
	}
	if (notric == 1) {
		nsys[5] = 0;
	}
	} else {
	io___54.ciunit = lca;
	i__1 = s_rsle(&io___54);
	if (i__1 != 0) {
		goto L3502;
	}
	for (i__ = 1; i__ <= 6; ++i__) {
		i__1 = do_lio(&c__3, &c__1, (char *)&nsys[i__ - 1], (ftnlen)
			sizeof(integer));
		if (i__1 != 0) {
		goto L3502;
		}
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L3502;
	}
	}

/* .....Read the tolerated error on 2-theta values W */
/*          which is also the column width of the */
/*          columnar profile shape */
/*          and how many non-indexed reflections NIND */

/*     If W is given as negative, then individual W */
/*     values will be read later (triplets : 2-theta, I, Width) */
/*          moreover, Width will be multiplied by -W */
/*          (use W = -1 for no change...) */

	if (ngrid == 3) {
	cal_1.nind = 3;
	w = slabda * .3f / 1.54056f;
	s_wsfe(&io___56);
	e_wsfe();
	s_wsfe(&io___57);
	do_fio(&c__1, (char *)&w, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cal_1.nind, (ftnlen)sizeof(integer));
	e_wsfe();
	} else {
	io___58.ciunit = lca;
	i__1 = s_rsle(&io___58);
	if (i__1 != 0) {
		goto L3503;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&w, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3503;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&cal_1.nind, (ftnlen)sizeof(
		integer));
	if (i__1 != 0) {
		goto L3503;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L3503;
	}
	}

/* .....Some printing */

	s_wsle(&io___59);
	e_wsle();
	s_wsfe(&io___60);
	do_fio(&c__1, (char *)&slabda, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&zero, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsle(&io___61);
	e_wsle();
	s_wsle(&io___62);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	if (ngrid == 0) {
	s_wsle(&io___63);
	do_lio(&c__9, &c__1, "    Monte Carlo cell generation", (ftnlen)31);
	e_wsle();
	}
	if (ngrid == 1) {
	s_wsle(&io___64);
	do_lio(&c__9, &c__1, "        Grid cell generation", (ftnlen)28);
	e_wsle();
	}
	if (ngrid == 2) {
	s_wsle(&io___65);
	do_lio(&c__9, &c__1, "Both searches - Monte Carlo AND grid", (ftnlen)
		36);
	e_wsle();
	}
	if (ngrid == 3) {
	s_wsle(&io___66);
	do_lio(&c__9, &c__1, " -- Black box Monte Carlo search --", (ftnlen)
		35);
	e_wsle();
	}
	if (nblack == 1) {
	s_wsle(&io___67);
	do_lio(&c__9, &c__1, " Black box Monte Carlo + grid search ", (ftnlen)
		37);
	e_wsle();
	}
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___69);
	e_wsle();
	s_wsfe(&io___70);
	do_fio(&c__1, (char *)&w, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsle(&io___71);
	e_wsle();
	s_wsfe(&io___72);
	do_fio(&c__1, (char *)&cal_1.nind, (ftnlen)sizeof(integer));
	e_wsfe();

/* .....READ Min/Max parameters, volume, Rp */
/* ...  RMI : If Rp < RMI, stop, should be the good cell */
/* ...  RMAX : Keep a refined cell if Rp < Rmax */
/* ...  RMAXREF : if Rp < RMAXREF, refine that cell by Monte Carlo */


	if (ngrid == 3) {
	pmin = 2.f;
	pmax = 30.f;
	vmin = 8.f;
	vmax = 2.7e4f;
	rmi = .02f;
	rmax = .15f;
	rmaxref = .4f;
	s_wsfe(&io___80);
	e_wsfe();
	s_wsle(&io___81);
	do_lio(&c__9, &c__1, " 2. 50. 8. 125000. 0.05 0.15 0.50", (ftnlen)33);
	e_wsle();
	} else {
	io___82.ciunit = lca;
	i__1 = s_rsle(&io___82);
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&pmin, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&pmax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&vmin, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&rmi, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&rmax, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&rmaxref, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3504;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L3504;
	}
	}
	if (pmin < 0.f) {
	io___83.ciunit = lca;
	i__1 = s_rsle(&io___83);
	if (i__1 != 0) {
		goto L3504;
	}
	for (i__ = 1; i__ <= 6; ++i__) {
		i__1 = do_lio(&c__4, &c__1, (char *)&pmi[i__ - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		goto L3504;
		}
		i__1 = do_lio(&c__4, &c__1, (char *)&pma[i__ - 1], (ftnlen)sizeof(
			real));
		if (i__1 != 0) {
		goto L3504;
		}
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L3504;
	}
	}
	s_wsle(&io___86);
	e_wsle();
	if (pmin > 0.f) {
/* 	WRITE(20,*)' Min/Max a,b,c, V ',PMIN,PMAX,VMIN,VMAX */
	for (i__ = 1; i__ <= 3; ++i__) {
		j = i__ + 3;
		pmi[i__ - 1] = pmin;
		pma[i__ - 1] = pmax;
		pmi[j - 1] = 60.f;
		pma[j - 1] = 120.f;
/* L30: */
	}
	} else {
	s_wsle(&io___88);
	do_lio(&c__9, &c__1, " Min/Max a cell parameter ", (ftnlen)26);
	do_lio(&c__4, &c__1, (char *)&pmi[0], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___89);
	do_lio(&c__9, &c__1, " Min/Max b cell parameter ", (ftnlen)26);
	do_lio(&c__4, &c__1, (char *)&pmi[1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[1], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___90);
	do_lio(&c__9, &c__1, " Min/Max c cell parameter ", (ftnlen)26);
	do_lio(&c__4, &c__1, (char *)&pmi[2], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___91);
	do_lio(&c__9, &c__1, " Min/Max alpha cell parameter ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&pmi[3], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[3], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___92);
	do_lio(&c__9, &c__1, " Min/Max beta  cell parameter ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&pmi[4], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[4], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___93);
	do_lio(&c__9, &c__1, " Min/Max gamma cell parameter ", (ftnlen)30);
	do_lio(&c__4, &c__1, (char *)&pmi[5], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[5], (ftnlen)sizeof(real));
	e_wsle();
	}

/*  NR is test for automatic Rmax decrease */

	nr = 0;
	if (rmax < 0.f) {
	nr = 1;
	rmax = -rmax;
	}

/* 	WRITE(20,*) */
/* 	WRITE(20,*)' Min/Max volumes ',VMIN,VMAX */
	s_wsle(&io___95);
	e_wsle();
	s_wsle(&io___96);
	do_lio(&c__9, &c__1, " Min/Max Rp, Rmaxref ", (ftnlen)21);
	do_lio(&c__4, &c__1, (char *)&rmi, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&rmax, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&rmaxref, (ftnlen)sizeof(real));
	e_wsle();

/* .....According to NGRID, read either grid steps */
/*                              or Monte Carlo parameters */
/*                              or both */
	if (ngrid == 3) {
	spar = .02f;
	sang = .05f;
/*      WRITE(28,9665) */
/* 9665  FORMAT('! Spar, Sang') */
/*      WRITE(28,*)'0.02 0.2' */
	s_wsle(&io___99);
	e_wsle();
	s_wsle(&io___100);
	do_lio(&c__9, &c__1, " Steps on (a,b,c) and angles ", (ftnlen)29);
	do_lio(&c__4, &c__1, (char *)&spar, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&sang, (ftnlen)sizeof(real));
	e_wsle();
	s_wsfe(&io___101);
	e_wsfe();
	s_wsle(&io___102);
	do_lio(&c__9, &c__1, "-100 20", (ftnlen)7);
	e_wsle();
	s_wsfe(&io___103);
	e_wsfe();
	} else {
	if (ngrid == 0) {
		goto L12;
	}
	if (ngrid == 1) {
		goto L11;
	}
	if (ngrid == 2) {
		goto L11;
	}
	s_wsle(&io___104);
	e_wsle();
	s_wsle(&io___105);
	do_lio(&c__9, &c__1, " UNKNOWN NGRID PARAMETER : STOP", (ftnlen)31);
	e_wsle();
	s_stop("", (ftnlen)0);
L11:

/* .....READ grid steps */
/*     SPAR = step on cell parameters */
/*     SANG = step on angles */

	io___106.ciunit = lca;
	i__1 = s_rsle(&io___106);
	if (i__1 != 0) {
		goto L3505;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&spar, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3505;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&sang, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3505;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L3505;
	}
	s_wsle(&io___107);
	e_wsle();
	s_wsle(&io___108);
	do_lio(&c__9, &c__1, " Steps on (a,b,c) and angles ", (ftnlen)29);
	do_lio(&c__4, &c__1, (char *)&spar, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&sang, (ftnlen)sizeof(real));
	e_wsle();
	if (ngrid == 2) {
		goto L12;
	}
	goto L13;
L12:

/*  Continue up to NTIMELIM tests */
/*  Save parameters if Rp < Rmax */
/*  Make NRUNS times those NTIMELIM tests */
/*  If NTIMELIM is negative, |NTIMELIM| will apply to cubic */
/*             and |NTIMELIM|*50. for tetragonal, hexagonal, */
/*                 etc */

	io___109.ciunit = lca;
	i__1 = s_rsle(&io___109);
	if (i__1 != 0) {
		goto L3506;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&timlim, (ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L3506;
	}
	i__1 = do_lio(&c__3, &c__1, (char *)&nruns, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
		goto L3506;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L3506;
	}
	if (timlim < 0.f) {
		timlim = -timlim;
		s_wsle(&io___112);
		e_wsle();
		s_wsle(&io___113);
		do_lio(&c__9, &c__1, " N of runs ", (ftnlen)11);
		do_lio(&c__3, &c__1, (char *)&nruns, (ftnlen)sizeof(integer));
		e_wsle();
		ntimelim[0] = timlim;
		s_wsle(&io___115);
		do_lio(&c__9, &c__1, " N of MC events in cubic        ", (ftnlen)
			32);
		do_lio(&c__5, &c__1, (char *)&ntimelim[0], (ftnlen)sizeof(
			doublereal));
		e_wsle();
		ntimelim[1] = ntimelim[0] * 20.f;
		s_wsle(&io___116);
		do_lio(&c__9, &c__1, " N of MC events in tetra/hexa   ", (ftnlen)
			32);
		do_lio(&c__5, &c__1, (char *)&ntimelim[1], (ftnlen)sizeof(
			doublereal));
		e_wsle();
		ntimelim[2] = ntimelim[1];
		ntimelim[3] = ntimelim[2] * 20.f;
		s_wsle(&io___117);
		do_lio(&c__9, &c__1, " N of MC events in orthorhombic ", (ftnlen)
			32);
		do_lio(&c__5, &c__1, (char *)&ntimelim[3], (ftnlen)sizeof(
			doublereal));
		e_wsle();
		ntimelim[4] = ntimelim[3] * 20.f;
		s_wsle(&io___118);
		do_lio(&c__9, &c__1, " N of MC events in monoclinic   ", (ftnlen)
			32);
		do_lio(&c__5, &c__1, (char *)&ntimelim[4], (ftnlen)sizeof(
			doublereal));
		e_wsle();
		ntimelim[5] = ntimelim[4] * 20.f;
		s_wsle(&io___119);
		do_lio(&c__9, &c__1, " N of MC events in triclinic    ", (ftnlen)
			32);
		do_lio(&c__5, &c__1, (char *)&ntimelim[5], (ftnlen)sizeof(
			doublereal));
		e_wsle();
	} else {
		for (i__ = 1; i__ <= 6; ++i__) {
		ntimelim[i__ - 1] = timlim;
/* L8221: */
		}
		s_wsle(&io___120);
		e_wsle();
		s_wsle(&io___121);
		do_lio(&c__9, &c__1, " N of MC events, N of runs ", (ftnlen)27);
		do_lio(&c__4, &c__1, (char *)&timlim, (ftnlen)sizeof(real));
		do_lio(&c__3, &c__1, (char *)&nruns, (ftnlen)sizeof(integer));
		e_wsle();
	}

L13:
	;
	}

/* ... Make a WARNING */

	if (ngrid == 3) {
	s_wsle(&io___122);
	e_wsle();
	s_wsle(&io___123);
	do_lio(&c__9, &c__1, "      This is the black box mode, can be long."
		"..", (ftnlen)48);
	e_wsle();
	if (notric == 1) {
		s_wsle(&io___124);
		e_wsle();
		s_wsle(&io___125);
		do_lio(&c__9, &c__1, "               No triclinic search.", (
			ftnlen)35);
		e_wsle();
	}
	}
	s_wsle(&io___126);
	e_wsle();
	s_wsle(&io___127);
	do_lio(&c__9, &c__1, "  To cancel and save, type K (capital letter) anyt"
		"ime", (ftnlen)53);
	e_wsle();
	s_wsle(&io___128);
	e_wsle();

/* .....And now, read couples of 2-theta and I values */

	cal_1.sum_f__ = 0.f;
	cal_1.ndat = 1;
L14:
	if (w > 0.f) {
	io___129.ciunit = lca;
	i__1 = s_rsle(&io___129);
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&cal_1.th2[cal_1.ndat - 1], (
		ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&cal_1.fobs[cal_1.ndat - 1], (
		ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L15;
	}
	} else {
	io___130.ciunit = lca;
	i__1 = s_rsle(&io___130);
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&cal_1.th2[cal_1.ndat - 1], (
		ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&cal_1.fobs[cal_1.ndat - 1], (
		ftnlen)sizeof(real));
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = do_lio(&c__4, &c__1, (char *)&w1[cal_1.ndat - 1], (ftnlen)
		sizeof(real));
	if (i__1 != 0) {
		goto L15;
	}
	i__1 = e_rsle();
	if (i__1 != 0) {
		goto L15;
	}
	}
	if (cal_1.th2[cal_1.ndat - 1] >= 180.f) {
	goto L3508;
	}
	++cal_1.ndat;
	if (cal_1.ndat > 100) {
	s_wsle(&io___132);
	do_lio(&c__9, &c__1, "Max data = 100 !", (ftnlen)16);
	e_wsle();
	}
	if (cal_1.ndat > 100) {
	s_stop("", (ftnlen)0);
	}
	goto L14;
L15:
	--cal_1.ndat;

/* .... Verify if these are d values or 2-theta */

	if (cal_1.th2[1] < cal_1.th2[0]) {
	s_wsle(&io___133);
	e_wsle();
	s_wsle(&io___134);
	do_lio(&c__9, &c__1, "    WARNING : DATA were given as d(A) values", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___135);
	e_wsle();
	s_wsle(&io___136);
	e_wsle();
	s_wsle(&io___137);
	do_lio(&c__9, &c__1, "    WARNING : DATA were given as d(A) values", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___138);
	e_wsle();
	i__1 = cal_1.ndat;
	for (nda = 1; nda <= i__1; ++nda) {
/* L1412: */
		cal_1.th2[nda - 1] = asin(slabda / (cal_1.th2[nda - 1] * 2.f)) * 
			cal_1.pi;
	}
	}

	i__1 = cal_1.ndat;
	for (nda = 1; nda <= i__1; ++nda) {
	if (w > 0.f) {
		w1[nda - 1] = w;
		if (ngrid == 3) {
		s_wsle(&io___140);
		do_lio(&c__4, &c__1, (char *)&cal_1.th2[nda - 1], (ftnlen)
			sizeof(real));
		do_lio(&c__4, &c__1, (char *)&cal_1.fobs[nda - 1], (ftnlen)
			sizeof(real));
		e_wsle();
		}
	} else {
		w1[nda - 1] *= -w;
	}
	if (cal_1.th2[cal_1.ndat - 1] >= 180.f) {
		goto L3508;
	}

/*     Addition of the Zeropoint */

	cal_1.th2[nda - 1] += zero;
	d__[nda - 1] = slabda / (sin(cal_1.th2[nda - 1] / cal_1.pi) * 2.f);
/* Computing 2nd power */
	r__1 = d__[nda - 1];
	qo[nda - 1] = 1 / (r__1 * r__1);
/* L1413: */
	}
	if (ngrid == 3 && cal_1.ndat >= 20) {
	cal_1.ndat = 20;
	}
	if (ngrid == 3) {
	cl__1.cerr = 0;
	cl__1.cunit = 28;
	cl__1.csta = 0;
	f_clos(&cl__1);
	}
	nhkl = cal_1.ndat;

/* ... NDAT10 is the max limit for the number of calculated */
/*           peak positions = 10 times the number of */
/*           observed peak positions */

	cal_1.ndat10 = cal_1.ndat * 10;
	ndat2 = cal_1.ndat << 1;

	i__1 = nhkl;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L9400: */
	cal_1.sum_f__ += cal_1.fobs[i__ - 1];
	}
	nmax = cal_1.ndat - cal_1.nind;

/* .....END OF DATA READING */

	cl__1.cerr = 0;
	cl__1.cunit = lca;
	cl__1.csta = "DELETE";
	f_clos(&cl__1);

/*  Output of some data */

	s_wsle(&io___146);
	e_wsle();
	s_wsle(&io___147);
	do_lio(&c__9, &c__1, "    2-THETA     d(A)    Iobs       W", (ftnlen)36);
	e_wsle();
	i__1 = nhkl;
	for (nm = 1; nm <= i__1; ++nm) {
	s_wsfe(&io___149);
	do_fio(&c__1, (char *)&cal_1.th2[nm - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&d__[nm - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&cal_1.fobs[nm - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&w1[nm - 1], (ftnlen)sizeof(real));
	e_wsfe();
/* L21: */
	}
	s_wsle(&io___150);
	e_wsle();

/* ...  Various starting values initialized */

/* Computing 2nd power */
	r__1 = slabda;
	cal_1.slabda2 = r__1 * r__1 / 4.f;
	if (pmin > 0.f) {
	for (i__ = 1; i__ <= 3; ++i__) {
		deltct[i__ - 1] = 30.f;
		astartt[i__ - 1] = 60.f;
		delta[i__ - 1] = (pmax - pmin) / 2.f;
/* L23: */
		pstart[i__ - 1] = pmin;
	}
	deltc = 15.f;
	astart = 90.f;
	} else {
	for (i__ = 1; i__ <= 3; ++i__) {
		j = i__ + 3;
		deltct[i__ - 1] = (pma[j - 1] - pmi[j - 1]) / 2.f;
		astartt[i__ - 1] = pmi[j - 1];
		delta[i__ - 1] = (pma[i__ - 1] - pmi[i__ - 1]) / 2.f;
/* L24: */
		pstart[i__ - 1] = pmi[i__ - 1];
	}
	deltc = (pma[4] - pmi[4]) / 2.f;
	astart = pmi[4];
	}
	for (i__ = 1; i__ <= 6; ++i__) {
/* L1925: */
	rmax0[i__ - 1] = rmax;
	}
	i__1 = cal_1.ndat;
	for (j = 1; j <= i__1; ++j) {
	cal_1.w2[j - 1] = w1[j - 1] / 2.f;
	cal_1.difp[j - 1] = cal_1.th2[j - 1] + cal_1.w2[j - 1];
	cal_1.difm[j - 1] = cal_1.th2[j - 1] - cal_1.w2[j - 1];
/* L25: */
	}
/*   DELTAB = zone explored in angstrom around a good cell */
/*   DMIN acts as a lower d limit for keeping reflections */
	deltab = .02f;
/*      DMIN=SLABDA/(2.*SIN(TH2(NHKL)/PI))-DELTAB */
	cal_1.dmin__ = d__[nhkl - 1] - deltab;
/*  Dmax values will help to determine max cell parameters */
/*    in black box mode */
	dmax1 = d__[0] + deltab * 2.f;
	dmax2 = d__[1] + deltab * 2.f;
	dmax3 = d__[2] + deltab * 2.f;

/*  Warning on the wavelength... */

	s_wsle(&io___162);
	e_wsle();
	s_wsle(&io___163);
	do_lio(&c__9, &c__1, " dmax = ", (ftnlen)8);
	do_lio(&c__4, &c__1, (char *)&d__[0], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___164);
	do_lio(&c__9, &c__1, " Be sure that your choice of max cell parameters", (
		ftnlen)48);
	e_wsle();
	s_wsle(&io___165);
	do_lio(&c__9, &c__1, "  ensures the exploration of all possibilities.", (
		ftnlen)47);
	e_wsle();
	s_wsle(&io___166);
	e_wsle();
	s_wsle(&io___167);
	e_wsle();
	s_wsle(&io___168);
	do_lio(&c__9, &c__1, " dmax = ", (ftnlen)8);
	do_lio(&c__4, &c__1, (char *)&d__[0], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___169);
	do_lio(&c__9, &c__1, " Be sure that your choice of max cell parameters", (
		ftnlen)48);
	e_wsle();
	s_wsle(&io___170);
	do_lio(&c__9, &c__1, "  ensures the exploration of all possibilities.", (
		ftnlen)47);
	e_wsle();
	s_wsle(&io___171);
	e_wsle();
	if (dmax1 > 30.f) {
	s_wsle(&io___172);
	e_wsle();
	s_wsle(&io___173);
	do_lio(&c__9, &c__1, " WARNING : dmax > 30 A.", (ftnlen)23);
	e_wsle();
	s_wsle(&io___174);
	do_lio(&c__9, &c__1, " Divide the wavelength by 10 and try again...", 
		(ftnlen)45);
	e_wsle();
	s_wsle(&io___175);
	do_lio(&c__9, &c__1, " and then, multiply the cell parameters by 10.",
		 (ftnlen)46);
	e_wsle();
	s_wsle(&io___176);
	e_wsle();
	s_wsle(&io___177);
	e_wsle();
	s_wsle(&io___178);
	do_lio(&c__9, &c__1, " WARNING : dmax > 30 A.", (ftnlen)23);
	e_wsle();
	s_wsle(&io___179);
	do_lio(&c__9, &c__1, " Divide the wavelength by 10 and try again...", 
		(ftnlen)45);
	e_wsle();
	s_wsle(&io___180);
	do_lio(&c__9, &c__1, " and then, multiply the cell parameters by 10.",
		 (ftnlen)46);
	e_wsle();
	s_wsle(&io___181);
	e_wsle();
	}

	s_wsle(&io___182);
	e_wsle();
	s_wsle(&io___183);
	do_lio(&c__9, &c__1, "   Dmin = ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&cal_1.dmin__, (ftnlen)sizeof(real));
	e_wsle();
/* Computing 2nd power */
	r__1 = cal_1.dmin__;
	cal_1.dmin__ = 1 / (r__1 * r__1);
/*   DELTAD = zone explored in 1/100 of degrees around a good cell */
	deltad = .2f;
/*  IGC = number of retained cells, IGM = max IGC */
/*  IREF = code for refining the best cell if it had */
/*         not Rp < Rmin */
/*  IGT = total number of cells satisfying to Rp < Rmax */
/*        including multiple identical cells */
	igc = 0;
	igt = 0.f;
	iref = 0;
	igm = 10000;
	ihr = 0;
	nruns2 = 1;
	rpsmall = 1.f;

/*  ESCAPE : a value of 0.15 means that in 15% of the tests, */
/*            a parameter change may be accepted even if that */
/*            change does not lead to any Rp or number of */
/*            indexed reflections improvement */

	escape = .15f;

	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	s_wsle(&io___193);
	e_wsle();
	s_wsle(&io___194);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___195);
	do_lio(&c__9, &c__1, "                RESULTS - RESULTS - RESULTS - RESU"
		"LTS", (ftnlen)53);
	e_wsle();
	s_wsle(&io___196);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___197);
	do_lio(&c__9, &c__1, "                EXPLORED CELL PARAMETERS AND VOLUM"
		"ES:", (ftnlen)53);
	e_wsle();
	s_wsle(&io___198);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___199);
	e_wsle();
/* ------------------------------------------------------------------------- */
/*     Initialisation */

	esp_init__(&iseed);

/* ------------------------------------------------------------------------- */


/* .....AND NOW : Generate the cell, either by Monte Carlo */
/*               or by grid search */

	if (ngrid == 1) {
	goto L700;
	}

/* ...  Cell generation by Monte Carlo */

	s_wsle(&io___201);
	e_wsle();
	s_wsle(&io___202);
	do_lio(&c__9, &c__1, "Monte Carlo search :", (ftnlen)20);
	e_wsle();

	if (nsys[0] == 0) {
	goto L200;
	}

/*    Cubic case */

	s_wsle(&io___203);
	do_lio(&c__9, &c__1, "Cubic:        Rp     a        V     Nind", (ftnlen)
		40);
	e_wsle();

	ifile = 1;
	ncycles = 200.f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 90.f;

	if (ngrid == 3) {
	nruns = 1;
	pmin = dmax1 * .9f;
	pmax = dmax1 * 3.1f;
	vmin = pmin * pmin * pmin;
	vmax = pmax * pmax * pmax;
	ntimelim[0] = (vmax - vmin) * .5f;
	if (ntimelim[0] > 1e4f) {
		ntimelim[0] = 1e4f;
	}
	delta[0] = (pmax - pmin) / 2.f;
	pstart[0] = pmin;
	}

	s_wsle(&io___208);
	e_wsle();
	s_wsle(&io___209);
	do_lio(&c__9, &c__1, "Cubic Monte Carlo search :", (ftnlen)26);
	e_wsle();
	s_wsle(&io___210);
	do_lio(&c__9, &c__1, " Max a, V ", (ftnlen)10);
	do_lio(&c__4, &c__1, (char *)&pmax, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___211);
	e_wsle();

	if (iverb == 1) {
	s_wsfe(&io___212);
	do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ntimelim[0], (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsle(&io___214);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___215);
	do_lio(&c__9, &c__1, " Rp  Trial number   a        V    Nind Icod", (
		ftnlen)43);
	e_wsle();
	s_wsle(&io___216);
	e_wsle();
	}

/*     READ hkl Miller indices in cub.hkl */

	s_copy(tempo, "cub.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___217);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 6;
	if (cal_1.nhkl0 > 400) {
	cal_1.nhkl0 = 400;
	}
	i__1 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L101: */
	s_rsle(&io___218);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

	i__1 = nruns;
	for (nrun = 1; nrun <= i__1; ++nrun) {
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	rmax = rmaxref;
	rmin = rmax;

/* ...  here starts the loop */

	interest = 0;
	tmax = ntimelim[0] / procs;
	ttmax = ntimelim[0] * 10.f;
	ncells = (integer) ntimelim[0];
	iiseed = 0;
	ntried = 0.f;
	ntriedt = 0.f;
	nout = 0;

/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC, */
/* $OMP& RMAX2,A,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,DDT,DDQ, */
/* $OMP& DIFF,DIFF2) */
/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
/* $OMP& celpre,rmin,rmax,bb,afi) */
/* $OMP DO */

	i__3 = ncells;
	for (ncel = 1; ncel <= i__3; ++ncel) {
		if (nout >= 1) {
		goto L196;
		}
		if (interest >= 1) {
		goto L196;
		}
		++iiseed;
		if (iiseed == 1) {
		iseed = ((iseed - ncel * nrun) / 2 << 1) + 1;
		}

L102:

		ntriedb = 0.f;
		celpre[0] = pstart[0] + delta[0] * 2.f * randi_(&iseed);
		ntried += 1.f;
		goto L104;
L103:
		del = deltab * (1.f - ntriedb / cy);
		celpre[0] = pstartb[0] + del * (randi_(&iseed) - .5f) * 2.f;
		ntriedb += 1.f;
L104:
		celpre[1] = celpre[0];
		celpre[2] = celpre[0];
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L105: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}
		dcell_(celpre, cal_1.al, &v1);
		if (ntried > tmax) {
		++nout;
		goto L196;
		}
		if (ntriedb != 0.f) {
		goto L106;
		}
		if (v1 > vmax || v1 < vmin) {
		ntried += -1.f;
		ntriedt += 1.f;
		if (ntriedt > ttmax) {
			++nout;
			goto L196;
		}
		goto L102;
		}

L106:
		calcul1_(&diff, &diff2);
		if (cal_1.nmx > cal_1.ndat10) {
		ntried += -1;
		goto L102;
		}
		if (ntriedb != 0.f) {
		goto L114;
		}

/* ... Here are the 2 criteria for selecting a cell to be */
/* ...      "refined" by Monte Carlo (NCYCLES events) */
/* ... Rp value satisfying (<0.5) ??? or enough hkl explained ??? */

		if (cal_1.lhkl >= nmax) {
		rmax = diff;
		icode = 2;
		if (diff <= rmaxref) {
			icode = 1;
		}
		} else {
		icode = 1;
		}
		if (diff > rmax) {
		goto L117;
		}
		if (cal_1.lhkl < nmax) {
		goto L117;
		}
L114:
		if (diff <= rmax) {
		llhkl = cal_1.lhkl;
		rmax = diff;
		rmax2 = diff2;
		a = celpre[0];
		v2 = v1;
		if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			v3 = v1;
		}

/* ... "Refine" that cell (by Monte Carlo too...) */

		pstartb[0] = celpre[0];
		}
		if (ntriedb <= ncycles) {
		goto L103;
		}
		ntriedb = 0.f;
		if (rmax >= rmax0[0]) {
		goto L117;
		}
		if (rmax2 >= .15f) {
		goto L117;
		}
		ipen = cal_1.ndat - llhkl;
		if (ipen > cal_1.nind) {
		goto L117;
		}

/* $OMP CRITICAL(STORE1) */

		++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

		igt += 1.f;
		if (nr == 1) {
		if (igt > 50.f) {
			if (ntried / igt < 100.f) {
			if (rmax0[0] > .2f) {
				rmax0[0] -= rmax0[0] * .05f;
				s_wsle(&io___244);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now "
					"Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[0], (ftnlen)
					sizeof(real));
				e_wsle();
				s_wsle(&io___245);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now "
					"Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[0], (ftnlen)
					sizeof(real));
				e_wsle();
			}
			}
		}
		}

		if (igc > 10000) {
		s_wsle(&io___246);
		do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (
			ftnlen)36);
		e_wsle();
		s_wsle(&io___247);
		do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (
			ftnlen)36);
		e_wsle();
		--igc;
		++interest;
		}
		cel[igc * 6 - 6] = a;
		cel[igc * 6 - 5] = a;
		cel[igc * 6 - 4] = a;
		cel[igc * 6 - 3] = 90.f;
		cel[igc * 6 - 2] = 90.f;
		cel[igc * 6 - 1] = 90.f;

/* $OMP END CRITICAL(STORE1) */

/* ... Check for supercell */

		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = a;
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L140: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}
		dcell_(celpre, cal_1.al, &v1);

/* $OMP CRITICAL(STORE2) */

		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		km[igc - 1] = llhkl;
		km2[igc - 1] = cal_1.lhkl;
		ifi[igc - 1] = ifile;
		nsol[igc - 1] = 1;
		vgc[igc - 1] = v1;
		rp[igc - 1] = rmax;
		rp2[igc - 1] = diff;
		if (rp[igc - 1] < rpsmall) {
		rpsmall = rp[igc - 1];
		isee = 1;
		} else {
		isee = 0;
		}
		supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__3);
		brav_(&cal_1.lhkl, ihkl, &ibr);
		ib[igc - 1] = ibr;
		a = cel[igc * 6 - 6];
		cel[igc * 6 - 5] = a;
		cel[igc * 6 - 4] = a;
		v2 = vgc[igc - 1];

/* $OMP END CRITICAL(STORE2) */

/* ... Check for interesting result */

/*      IF(INTEREST.GE.1)GO TO 196 */
		indic = 1;
		bb[2] = a;
		bb[3] = a;
		bb[4] = a;
		bb[5] = 90.f;
		bb[6] = 90.f;
		bb[7] = 90.f;
		afi[2] = 1.f;
		afi[3] = 1.f;
		afi[4] = 1.f;
		afi[5] = 0.f;
		afi[6] = 0.f;
		afi[7] = 0.f;
		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = a;
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L110: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}

/* $OMP CRITICAL(FOUND) */

		if (rp[igc - 1] < rmi) {
		++interest;
		s_wsfe(&io___263);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___264);
		e_wsle();
		s_wsle(&io___265);
		do_lio(&c__9, &c__1, "======================================"
			"==============\r ===========================", (
			ftnlen)81);
		e_wsle();
		s_wsle(&io___266);
		e_wsle();
		s_wsle(&io___267);
		do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT "
			": Rp < Rmin !", (ftnlen)51);
		e_wsle();
		s_wsle(&io___268);
		e_wsle();
		s_wsle(&io___269);
		e_wsle();
		s_wsle(&io___270);
		do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT "
			": Rp < Rmin !", (ftnlen)51);
		e_wsle();
		s_wsle(&io___271);
		e_wsle();

/* ... Refine that cell */

		dcell_(celpre, cal_1.al, &v1);
		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
		cncalc[igc - 1] = (real) ncalc;
		if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] * 
				2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
		}
/*      WRITE(20,7000)FM20(IGC) */
/*      WRITE(20,7001)FF20(IGC),DDT,NCALC */
/* 	WRITE(20,*) */
/*      PRINT 7000,FM20(IGC) */
/*      PRINT 7001,FF20(IGC),DDT,NCALC */
/* 	PRINT * */
		iref = 1;
		goto L197;
		} else {

/*  Anyway, calculate the M20 and F20 values */

		dcell_(celpre, cal_1.al, &v1);
		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
		cncalc[igc - 1] = (real) ncalc;
		if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] * 
				2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
		}
		}

/* Test if cell already found */

		if (igc > 1) {
		i__4 = igc - 1;
		for (i__ = 1; i__ <= i__4; ++i__) {
			if (ifi[i__ - 1] != ifile) {
			goto L118;
			}
			vdelt = vgc[igc - 1] / 300.f;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
			goto L118;
			}
			++nsol[i__ - 1];
			if (rp[igc - 1] < rp[i__ - 1]) {
			if (isee == 1) {
				s_wsfe(&io___281);
				do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real))
					;
				do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(
					integer));
				e_wsfe();
			}
			km[i__ - 1] = km[igc - 1];
			vgc[i__ - 1] = vgc[igc - 1];
			rp[i__ - 1] = rp[igc - 1];
			cel[i__ * 6 - 6] = cel[igc * 6 - 6];
			cel[i__ * 6 - 5] = cel[igc * 6 - 5];
			cel[i__ * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			if (nsol[i__ - 1] > 5) {
			ntried = tmax + 1.f;
			++nout;
			}
			goto L119;
L118:
			;
		}
		if (iverb == 1) {
			s_wsfe(&io___282);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		if (isee == 1) {
			s_wsfe(&io___283);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
L119:
		;
		} else {
		if (iverb == 1) {
			s_wsfe(&io___284);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		if (isee == 1) {
			s_wsfe(&io___285);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		}
L197:

/* $OMP END CRITICAL(FOUND) */




/* First criterium reinitialized to Rmaxref */

L117:
		rmax = rmaxref;

/* ... Stop if max limit of Monte Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

/*   END ON MC tests */

L196:
		;
	}

/* $OMP END DO NOWAIT */
/* $OMP END PARALLEL */

	if (interest >= 1) {
		goto L5000;
	}
	killk_(&pressedk);
	if (rmin == rmax) {
		goto L198;
	}
	if (iverb == 1) {
		s_wsle(&io___287);
		e_wsle();
		s_wsle(&io___288);
		do_lio(&c__9, &c__1, "Best result : a=", (ftnlen)16);
		do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, "V=", (ftnlen)2);
		do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " Rp=", (ftnlen)4);
		do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___289);
		e_wsle();
		datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L198:
	if (pressedk) {
		goto L5000;
	}

/*  END ON NRUNS */

/* L199: */
	}

L200:
	if (nsys[1] == 0) {
	goto L300;
	}

/*    Hexagonal case */


	if (ngrid == 3) {
	ihr = 1;
	}
L290:
	if (ihr == 2) {
	nsys[1] = 2;
	}
	if (nsys[1] == 1) {
	s_wsle(&io___290);
	do_lio(&c__9, &c__1, "Hexagonal:    Rp     a       c        V     Ni"
		"nd", (ftnlen)48);
	e_wsle();
	rpsmall = 1.f;
	} else {
	s_wsle(&io___291);
	do_lio(&c__9, &c__1, "Rhombohedral: Rp     a       c        V     Ni"
		"nd", (ftnlen)48);
	e_wsle();
	rpsmall = 1.f;
	}

	ifile = 2;
	ncycles = 500.f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 120.f;

	if (ngrid == 3) {
	nruns = 10;
	pmin = 2.f;
	pmax = 30.f;
	if (nsys[1] == 2) {
		pmax = 60.f;
	}
	pma[2] = dmax1 * 3.1f;
	if (nsys[1] == 2) {
		pma[2] = dmax1 * 6.1f;
	}
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	pma[0] = dmax1 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	pma[1] = pma[0];
	vmin = 8.f;
	for (i__ = 1; i__ <= 3; ++i__) {
		pmi[i__ - 1] = pmin;
		delta[i__ - 1] = (pma[i__ - 1] - pmi[i__ - 1]) / 2.f;
/* L223: */
		pstart[i__ - 1] = pmi[i__ - 1];
	}
	vmax = pma[0] * pma[1] * pma[2];
	if (vmax > 4e3f) {
		vmax = 4e3f;
	}
	ntimelim[1] = vmax * 5.f;
	}

	s_wsle(&io___292);
	e_wsle();
	s_wsle(&io___293);
	do_lio(&c__9, &c__1, "Hexagonal/Trigonal/Rhomboedral Monte Carlo search :"
		, (ftnlen)51);
	e_wsle();
	s_wsle(&io___294);
	do_lio(&c__9, &c__1, " Max(a,c), V ", (ftnlen)13);
	do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___295);
	e_wsle();

	if (iverb == 1) {
	if (nsys[1] == 1) {
		s_wsfe(&io___296);
		do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ntimelim[1], (ftnlen)sizeof(doublereal));
		e_wsfe();
	}
	if (nsys[1] == 2) {
		s_wsfe(&io___297);
		do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ntimelim[1], (ftnlen)sizeof(doublereal));
		e_wsfe();
	}
	s_wsle(&io___298);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___299);
	do_lio(&c__9, &c__1, " Rp  Trial number    a      c        V  Nind I"
		"cod", (ftnlen)49);
	e_wsle();
	s_wsle(&io___300);
	e_wsle();
	}

/*     READ hkl Miller indices in hex.hkl */

	if (nsys[1] == 2) {
	goto L260;
	}
	s_copy(tempo, "hex.hkl", (ftnlen)80, (ftnlen)7);
	goto L261;
L260:
	s_copy(tempo, "rho.hkl", (ftnlen)80, (ftnlen)7);
	ifile = 7;
L261:
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___301);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	if (nsys[1] == 2) {
	goto L262;
	}
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	goto L263;
L262:
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 600) {
	cal_1.nhkl0 = 600;
	}
L263:
	i__1 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L201: */
	s_rsle(&io___302);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

	i__1 = nruns;
	for (nrun = 1; nrun <= i__1; ++nrun) {
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	rmax = rmaxref;
	rmin = rmax;

/* ...  here starts the loop */

	interest = 0;
	tmax = ntimelim[1] / procs;
	ttmax = ntimelim[1] * 10.f;
	ncells = (integer) ntimelim[1];
	iiseed = 0;
	ntried = 0.;
	ntriedt = 0.;
	nout = 0;
	celpre[0] = pstart[0] + delta[0] * 2.f * randi_(&iseed);
	celpre[1] = celpre[0];
	celpre[2] = pstart[2] + delta[2] * 2.f * randi_(&iseed);
	celold[0] = celpre[0];
	celold[2] = celpre[2];
	rglob = 1.f;
	nglob = 0;

/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC, */
/* $OMP& RMAX2,A,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP, */
/* $OMP& DIFF,DIFF2,DDT,DDQ) */
/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
/* $OMP DO */

	i__3 = ncells;
	for (ncel = 1; ncel <= i__3; ++ncel) {
		if (nout >= 1) {
		goto L296;
		}
		if (interest >= 1) {
		goto L296;
		}
		++iiseed;
		if (iiseed == 1) {
		iseed = ((iseed - ncel * nrun) / 2 << 1) + 1;
		}

L202:

/*     Which parameter to vary ? a or c ? */

		ntriedb = 0.f;
		ip = 3;
		if (randi_(&iseed) > .5f) {
		ip = 1;
		}
		celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.f * randi_(&
			iseed);
		ntried += 1.f;
		goto L204;
L203:
		del = deltab * (1.f - ntriedb / cy);
		i__ = 3;
		if (randi_(&iseed) > .5f) {
		i__ = 1;
		}
		celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - .5f) 
			* 2.f;
		ntriedb += 1.f;
L204:
		celpre[1] = celpre[0];
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L205: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}
		dcell_(celpre, cal_1.al, &v1);
		if (ntried > tmax) {
		++nout;
		goto L296;
		}
		if (ntriedb != 0.f) {
		goto L206;
		}
		if (v1 > vmax || v1 < vmin) {
		ntried += -1.f;
		ntriedt += 1.f;
		if (ntriedt > ttmax) {
			++nout;
			goto L296;
		}
		goto L202;
		}

L206:
		calcul1_(&diff, &diff2);
		if (cal_1.nmx > cal_1.ndat10) {
		ntried += -1;
		goto L202;
		}
		if (ntriedb != 0.f) {
		goto L214;
		}

/* ... Rp value satisfying ??? */

		if (diff < rglob || cal_1.lhkl > nglob) {
		rglob = diff;
		nglob = cal_1.lhkl;
		celold[ip - 1] = celpre[ip - 1];
		}
		if (cal_1.lhkl >= nmax) {
		rmax = diff;
		icode = 2;
		if (diff <= rmaxref) {
			icode = 1;
		}
		} else {
		icode = 1;
		}
		if (diff > rmax) {
		goto L217;
		}
		if (cal_1.lhkl < nmax) {
		goto L217;
		}
L214:
		if (diff <= rmax) {
		llhkl = cal_1.lhkl;
		rmax = diff;
		rmax2 = diff2;
		a = celpre[0];
		c__ = celpre[2];
		v2 = v1;
		if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[2] = c__;
			v3 = v1;
		}

/* ... "Refine" that cell (by Monte Carlo too...) */

		pstartb[0] = celpre[0];
		pstartb[2] = celpre[2];
		}
		if (ntriedb <= ncycles) {
		goto L203;
		}
		rglob = .5f;
		nglob = ndat2;
		if (ip == 1) {
		celold[ip - 1] = a;
		}
		if (ip == 3) {
		celold[ip - 1] = c__;
		}
		ntriedb = 0.f;
		if (rmax >= rmax0[1]) {
		goto L217;
		}
		if (rmax2 >= .15f) {
		goto L217;
		}
		ipen = cal_1.ndat - llhkl;
		if (ipen > cal_1.nind) {
		goto L217;
		}

/* $OMP CRITICAL (STORE1) */

		++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

		igt += 1.f;
		if (nr == 1) {
		if (igt > 50.f) {
			if (ntried / igt < 1e3f) {
			if (rmax0[1] > .2f) {
				rmax0[1] -= rmax0[1] * .05f;
				s_wsle(&io___308);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now "
					"Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[1], (ftnlen)
					sizeof(real));
				e_wsle();
				s_wsle(&io___309);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now "
					"Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[1], (ftnlen)
					sizeof(real));
				e_wsle();
			}
			}
		}
		}

		if (igc > 10000) {
		s_wsle(&io___310);
		do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (
			ftnlen)36);
		e_wsle();
		s_wsle(&io___311);
		do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (
			ftnlen)36);
		e_wsle();
		--igc;
		++interest;
/*      GO TO 5000 */
		}
		cel[igc * 6 - 6] = a;
		cel[igc * 6 - 5] = a;
		cel[igc * 6 - 4] = c__;
		cel[igc * 6 - 3] = 90.f;
		cel[igc * 6 - 2] = 90.f;
		cel[igc * 6 - 1] = 120.f;

/* $OMP END CRITICAL(STORE1) */

/* ... Check for supercell */

		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = c__;
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L240: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}
		dcell_(celpre, cal_1.al, &v1);

/* $OMP CRITICAL(STORE2) */

		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		km[igc - 1] = llhkl;
		km2[igc - 1] = cal_1.lhkl;
		ifi[igc - 1] = ifile;
		nsol[igc - 1] = 1;
		vgc[igc - 1] = v1;
		rp[igc - 1] = rmax;
		rp2[igc - 1] = diff;
		if (rp[igc - 1] < rpsmall) {
		rpsmall = rp[igc - 1];
		isee = 1;
		} else {
		isee = 0;
		}
		supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__2);
		brav_(&cal_1.lhkl, ihkl, &ibr);
		ib[igc - 1] = ibr;
		a = cel[igc * 6 - 6];
		cel[igc * 6 - 5] = a;
		c__ = cel[igc * 6 - 4];
		v2 = vgc[igc - 1];

/* $OMP END CRITICAL(STORE2) */

/* ... Check for interesting result */

/*      IF(INTEREST.GE.1)GO TO 296 */
		indic = 2;
		bb[2] = a;
		bb[3] = a;
		bb[4] = c__;
		bb[5] = 90.f;
		bb[6] = 90.f;
		bb[7] = 120.f;
		afi[2] = 1.f;
		afi[3] = 1.f;
		afi[4] = 1.f;
		afi[5] = 0.f;
		afi[6] = 0.f;
		afi[7] = 0.f;
		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = c__;
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L210: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}

/* $OMP CRITICAL(FOUND) */

		if (rp[igc - 1] < rmi) {
		++interest;
		s_wsfe(&io___312);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___313);
		e_wsle();
		s_wsle(&io___314);
		do_lio(&c__9, &c__1, "======================================"
			"==============\r ===========================", (
			ftnlen)81);
		e_wsle();
		s_wsle(&io___315);
		e_wsle();
		s_wsle(&io___316);
		do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT "
			": Rp < Rmin!", (ftnlen)50);
		e_wsle();
		s_wsle(&io___317);
		e_wsle();
		s_wsle(&io___318);
		do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT "
			": Rp < Rmin!", (ftnlen)50);
		e_wsle();

/* ... Refine that cell */

		dcell_(celpre, cal_1.al, &v1);
		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
		cncalc[igc - 1] = (real) ncalc;
		if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] * 
				2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
		}
/*      WRITE(20,7000)FM20(IGC) */
/*      WRITE(20,7001)FF20(IGC),DDT,NCALC */
/* 	WRITE(20,*) */
/*      PRINT 7000,FM20(IGC) */
/*      PRINT 7001,FF20(IGC),DDT,NCALC */
/* 	PRINT * */
		iref = 1;
		goto L297;
		} else {

/*  Anyway, calculate the M20 and F20 values */

		dcell_(celpre, cal_1.al, &v1);
		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
		cncalc[igc - 1] = (real) ncalc;
		if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] * 
				2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
		}
		}

/* Test if cell already found */

		if (igc > 1) {
		i__4 = igc - 1;
		for (i__ = 1; i__ <= i__4; ++i__) {
			if (ifi[i__ - 1] != ifile) {
			goto L218;
			}
			vdelt = vgc[igc - 1] / 300.f;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
			goto L218;
			}
			adelt = cel[igc * 6 - 6] / 500.f;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
			goto L218;
			}
			++nsol[i__ - 1];
			if (rp[igc - 1] < rp[i__ - 1]) {
			if (isee == 1) {
				s_wsfe(&io___322);
				do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real))
					;
				do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(
					integer));
				e_wsfe();
			}
			km[i__ - 1] = km[igc - 1];
			vgc[i__ - 1] = vgc[igc - 1];
			rp[i__ - 1] = rp[igc - 1];
			cel[i__ * 6 - 6] = cel[igc * 6 - 6];
			cel[i__ * 6 - 5] = cel[igc * 6 - 5];
			cel[i__ * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			if (nsol[i__ - 1] > 5) {
			ntried = tmax + 1.f;
			++nout;
			}
			goto L219;
L218:
			;
		}
		if (iverb == 1) {
			s_wsfe(&io___323);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		if (isee == 1) {
			s_wsfe(&io___324);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
L219:
		;
		} else {
		if (iverb == 1) {
			s_wsfe(&io___325);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		if (isee == 1) {
			s_wsfe(&io___326);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		}
L297:

/* $OMP END CRITICAL(FOUND) */

L217:
		rmax = rmaxref;
		if (randi_(&iseed) > escape) {
		celpre[ip - 1] = celold[ip - 1];
		}

/*  END ON MC TESTS */

L296:
		;
	}

/* $OMP END DO NOWAIT */
/* $OMP END PARALLEL */


/* ... Stop if max limit of Monte Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (interest >= 1) {
		goto L5000;
	}
	killk_(&pressedk);
	if (rmin == rmax) {
		goto L298;
	}
	if (iverb == 1) {
		s_wsle(&io___327);
		e_wsle();
		s_wsle(&io___328);
		do_lio(&c__9, &c__1, "Best result : a = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___329);
		do_lio(&c__9, &c__1, "Best result : c = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, "V = ", (ftnlen)4);
		do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___330);
		e_wsle();
		datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L298:
	if (pressedk) {
		goto L5000;
	}

/*  END ON NRUNS */

/* L299: */
	}


	++ihr;
	if (ihr == 2) {
	goto L290;
	}
L300:
	if (nsys[2] == 0) {
	goto L400;
	}

/*    Tetragonal case */


	rpsmall = 1.f;
	s_wsle(&io___331);
	do_lio(&c__9, &c__1, "Tetragonal:   Rp     a       c        V     Nind", (
		ftnlen)48);
	e_wsle();

	ifile = 3;
	ncycles = 500.f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 90.f;

	if (ngrid == 3) {
	nruns = 10;
	pmin = 2.f;
	pmax = 30.f;
	pma[2] = dmax1 * 4.1f;
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	pma[0] = dmax1 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	pma[1] = pma[0];
	vmin = 8.f;
	for (i__ = 1; i__ <= 3; ++i__) {
		pmi[i__ - 1] = pmin;
		delta[i__ - 1] = (pma[i__ - 1] - pmi[i__ - 1]) / 2.f;
/* L323: */
		pstart[i__ - 1] = pmi[i__ - 1];
	}
	vmax = pma[0] * pma[1] * pma[2];
	if (vmax > 4e3f) {
		vmax = 4e3f;
	}
	ntimelim[2] = vmax * 5.f;
	}

	s_wsle(&io___332);
	e_wsle();
	s_wsle(&io___333);
	do_lio(&c__9, &c__1, "Tetragonal Monte Carlo search :", (ftnlen)31);
	e_wsle();
	s_wsle(&io___334);
	do_lio(&c__9, &c__1, " Max(a,c), V ", (ftnlen)13);
	do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___335);
	e_wsle();

	if (iverb == 1) {
	s_wsfe(&io___336);
	do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ntimelim[2], (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsle(&io___337);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___338);
	do_lio(&c__9, &c__1, " Rp  Trial number    a      c        V  Nind I"
		"cod", (ftnlen)49);
	e_wsle();
	s_wsle(&io___339);
	e_wsle();
	}

/*     READ hkl Miller indices in tet.hkl */

	s_copy(tempo, "tet.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___340);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	i__1 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L301: */
	s_rsle(&io___341);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

	i__1 = nruns;
	for (nrun = 1; nrun <= i__1; ++nrun) {
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	rmax = rmaxref;
	rmin = rmax;

/* ...  here starts the loop */

	interest = 0;
	tmax = ntimelim[2] / procs;
	ttmax = ntimelim[2] * 10.f;
	ncells = (integer) ntimelim[2];
	iiseed = 0;
	ntried = 0.f;
	ntriedt = 0.f;
	nout = 0;
	celpre[0] = pstart[0] + delta[0] * 2.f * randi_(&iseed);
	celpre[1] = celpre[0];
	celpre[2] = pstart[2] + delta[2] * 2.f * randi_(&iseed);
	celold[0] = celpre[0];
	celold[2] = celpre[2];
	rglob = 1.f;
	nglob = 0;

/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC, */
/* $OMP& RMAX2,A,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP, */
/* $OMP& DIFF,DIFF2,DDT,DDQ) */
/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
/* $OMP DO */

	i__3 = ncells;
	for (ncel = 1; ncel <= i__3; ++ncel) {
		if (nout >= 1) {
		goto L396;
		}
		if (interest >= 1) {
		goto L396;
		}
		++iiseed;
		if (iiseed == 1) {
		iseed = ((iseed - ncel * nrun) / 2 << 1) + 1;
		}

L302:

/*     Which parameter to vary ? a or c ? */

		ntriedb = 0.f;
		ip = 3;
		if (randi_(&iseed) > .5f) {
		ip = 1;
		}
		celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.f * randi_(&
			iseed);
		ntried += 1.f;
		goto L304;
L303:
		del = deltab * (1.f - ntriedb / cy);
		i__ = 3;
		if (randi_(&iseed) > .5f) {
		i__ = 1;
		}
		celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - .5f) 
			* 2.f;
		ntriedb += 1.f;
L304:
		celpre[1] = celpre[0];
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L305: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}
		dcell_(celpre, cal_1.al, &v1);
		if (ntried > tmax) {
		++nout;
		goto L396;
		}
		if (ntriedb != 0.f) {
		goto L306;
		}
		if (v1 > vmax || v1 < vmin) {
		ntried += -1.f;
		ntriedt += 1.f;
		if (ntriedt > ttmax) {
			++nout;
			goto L396;
		}
		goto L302;
		}

L306:
		calcul1_(&diff, &diff2);
		if (cal_1.nmx > cal_1.ndat10) {
		ntried += -1;
		goto L302;
		}
		if (ntriedb != 0.f) {
		goto L314;
		}

/* ... Rp value satisfying ??? */

		if (diff < rglob || cal_1.lhkl > nglob) {
		rglob = diff;
		nglob = cal_1.lhkl;
		celold[ip - 1] = celpre[ip - 1];
		}
		if (cal_1.lhkl >= nmax) {
		rmax = diff;
		icode = 2;
		if (diff <= rmaxref) {
			icode = 1;
		}
		} else {
		icode = 1;
		}
		if (diff > rmax) {
		goto L317;
		}
		if (cal_1.lhkl < nmax) {
		goto L317;
		}
L314:
		if (diff <= rmax) {
		llhkl = cal_1.lhkl;
		rmax = diff;
		rmax2 = diff2;
		a = celpre[0];
		c__ = celpre[2];
		v2 = v1;
		if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[2] = c__;
			v3 = v1;
		}

/* ... "Refine" that cell (by Monte Carlo too...) */

		pstartb[0] = celpre[0];
		pstartb[2] = celpre[2];
		}
		if (ntriedb <= ncycles) {
		goto L303;
		}
		rglob = .5f;
		nglob = ndat2;
		if (ip == 1) {
		celold[ip - 1] = a;
		}
		if (ip == 3) {
		celold[ip - 1] = c__;
		}
		ntriedb = 0.f;
		if (rmax >= rmax0[2]) {
		goto L317;
		}
		if (rmax2 >= .15f) {
		goto L317;
		}
		ipen = cal_1.ndat - llhkl;
		if (ipen > cal_1.nind) {
		goto L317;
		}

/* $OMP CRITICAL(STORE1) */

		++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

		igt += 1.f;
		if (nr == 1) {
		if (igt > 50.f) {
			if (ntried / igt < 1e3f) {
			if (rmax0[2] > .2f) {
				rmax0[2] -= rmax0[2] * .05f;
				s_wsle(&io___342);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now "
					"Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[2], (ftnlen)
					sizeof(real));
				e_wsle();
				s_wsle(&io___343);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now "
					"Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[2], (ftnlen)
					sizeof(real));
				e_wsle();
			}
			}
		}
		}

		if (igc > 10000) {
		s_wsle(&io___344);
		do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (
			ftnlen)36);
		e_wsle();
		s_wsle(&io___345);
		do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (
			ftnlen)36);
		e_wsle();
		--igc;
		++interest;
/*      GO TO 5000 */
		}
		cel[igc * 6 - 6] = a;
		cel[igc * 6 - 5] = a;
		cel[igc * 6 - 4] = c__;
		cel[igc * 6 - 3] = 90.f;
		cel[igc * 6 - 2] = 90.f;
		cel[igc * 6 - 1] = 90.f;

/* $OMP END CRITICAL(STORE1) */

/* ... Check for supercell */

		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = c__;
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L340: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}
		dcell_(celpre, cal_1.al, &v1);

/* $OMP CRITICAL(STORE2) */

		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		km[igc - 1] = llhkl;
		km2[igc - 1] = cal_1.lhkl;
		ifi[igc - 1] = ifile;
		nsol[igc - 1] = 1;
		vgc[igc - 1] = v1;
		rp[igc - 1] = rmax;
		rp2[igc - 1] = diff;
		if (rp[igc - 1] < rpsmall) {
		rpsmall = rp[igc - 1];
		isee = 1;
		} else {
		isee = 0;
		}
		supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__4);
		brav_(&cal_1.lhkl, ihkl, &ibr);
		ib[igc - 1] = ibr;
		a = cel[igc * 6 - 6];
		cel[igc * 6 - 5] = a;
		c__ = cel[igc * 6 - 4];
		v2 = vgc[igc - 1];

/* $OMP END CRITICAL(STORE2) */

/* ... Check for interesting result */

/*      IF(INTEREST.GE.1)GO TO 396 */
		indic = 2;
		bb[2] = a;
		bb[3] = a;
		bb[4] = c__;
		bb[5] = 90.f;
		bb[6] = 90.f;
		bb[7] = 90.f;
		afi[2] = 1.f;
		afi[3] = 1.f;
		afi[4] = 1.f;
		afi[5] = 0.f;
		afi[6] = 0.f;
		afi[7] = 0.f;
		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = c__;
		for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L310: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
		}

/* $OMP CRITICAL(FOUND) */

		if (rp[igc - 1] < rmi) {
		++interest;
		s_wsfe(&io___346);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___347);
		e_wsle();
		s_wsle(&io___348);
		do_lio(&c__9, &c__1, "======================================"
			"==============\r ===========================", (
			ftnlen)81);
		e_wsle();
		s_wsle(&io___349);
		e_wsle();
		s_wsle(&io___350);
		do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT "
			": Rp < Rmin!", (ftnlen)50);
		e_wsle();
		s_wsle(&io___351);
		e_wsle();
		s_wsle(&io___352);
		do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT "
			": Rp < Rmin!", (ftnlen)50);
		e_wsle();

/* ... Refine that cell */

		dcell_(celpre, cal_1.al, &v1);
		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
		cncalc[igc - 1] = (real) ncalc;
		if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] * 
				2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
		}
/*      WRITE(20,7000)FM20(IGC) */
/*      WRITE(20,7001)FF20(IGC),DDT,NCALC */
/* 	WRITE(20,*) */
/*      PRINT 7000,FM20(IGC) */
/*      PRINT 7001,FF20(IGC),DDT,NCALC */
/* 	PRINT * */
		iref = 1;
		goto L397;
		} else {

/*  Anyway, calculate the M20 and F20 values */

		dcell_(celpre, cal_1.al, &v1);
		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
		cncalc[igc - 1] = (real) ncalc;
		if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] * 
				2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
		}
		}

/* Test if cell already found */

		if (igc > 1) {
		i__4 = igc - 1;
		for (i__ = 1; i__ <= i__4; ++i__) {
			if (ifi[i__ - 1] != ifile) {
			goto L318;
			}
			vdelt = vgc[igc - 1] / 300.f;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
			goto L318;
			}
			adelt = cel[igc * 6 - 6] / 500.f;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
			goto L318;
			}
			++nsol[i__ - 1];
			if (rp[igc - 1] < rp[i__ - 1]) {
			if (isee == 1) {
				s_wsfe(&io___353);
				do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real))
					;
				do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
				do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(
					integer));
				e_wsfe();
			}
			km[i__ - 1] = km[igc - 1];
			vgc[i__ - 1] = vgc[igc - 1];
			rp[i__ - 1] = rp[igc - 1];
			cel[i__ * 6 - 6] = cel[igc * 6 - 6];
			cel[i__ * 6 - 5] = cel[igc * 6 - 5];
			cel[i__ * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			if (nsol[i__ - 1] > 5) {
			ntried = tmax + 1.f;
			++nout;
			}
			goto L319;
L318:
			;
		}
		if (iverb == 1) {
			s_wsfe(&io___354);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		if (isee == 1) {
			s_wsfe(&io___355);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
L319:
		;
		} else {
		if (iverb == 1) {
			s_wsfe(&io___356);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal))
				;
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		if (isee == 1) {
			s_wsfe(&io___357);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		}
L397:

/* $OMP END CRITICAL(FOUND) */


L317:
		rmax = rmaxref;
		if (randi_(&iseed) > escape) {
		celpre[ip - 1] = celold[ip - 1];
		}

/*   END ON 4000 tests */

L396:
		;
	}

/* $OMP END DO NOWAIT */
/* $OMP END PARALLEL */


/* ... Stop if max limit of Monte Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (interest >= 1) {
		goto L5000;
	}
	killk_(&pressedk);
	if (rmin == rmax) {
		goto L398;
	}
	if (iverb == 1) {
		s_wsle(&io___358);
		e_wsle();
		s_wsle(&io___359);
		do_lio(&c__9, &c__1, "Best result : a = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___360);
		do_lio(&c__9, &c__1, "Best result : c = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, "V = ", (ftnlen)4);
		do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___361);
		e_wsle();
		datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L398:
	if (pressedk) {
		goto L5000;
	}

/*  END ON NRUNS */

/* L399: */
	}


L400:
	if (nsys[3] == 0) {
	goto L500;
	}

/*    Orthorhombic case */


	rpsmall = 1.f;
	s_wsle(&io___362);
	do_lio(&c__9, &c__1, "Orthorhombic: Rp     a       b       c        V   "
		"  Nind", (ftnlen)56);
	e_wsle();
	ifile = 4;
	ncycles = 1e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 90.f;

	if (ngrid == 3) {
	nruns2 = 6;
	nruns = 10;
	pmin = 2.f;
	pmax = 20.f;
	pma[0] = dmax1 * 2.1f;
	pma[1] = dmax2 * 2.1f;
	pma[2] = dmax3 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	if (pma[1] > pmax) {
		pma[1] = pmax;
	}
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	vorth = pma[0] * pma[1] * pma[2];
	if (vorth > 3e3f) {
		vorth = 3e3f;
	}
	vmax = vorth;
	for (i__ = 1; i__ <= 3; ++i__) {
		pmi[i__ - 1] = pmin;
		delta[i__ - 1] = (pma[i__ - 1] - pmi[i__ - 1]) / 2.f;
/* L423: */
		pstart[i__ - 1] = pmi[i__ - 1];
	}
	}

	s_wsle(&io___364);
	e_wsle();
	s_wsle(&io___365);
	do_lio(&c__9, &c__1, "Orthorhombic Monte Carlo search :", (ftnlen)33);
	e_wsle();
	s_wsle(&io___366);
	do_lio(&c__9, &c__1, " Max(a,b,c), V ", (ftnlen)15);
	do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___367);
	e_wsle();

	if (iverb == 1) {
	s_wsfe(&io___368);
	do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ntimelim[3], (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsle(&io___369);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___370);
	do_lio(&c__9, &c__1, " Rp  Trial number    a   b   c       V  Nind i"
		"cod", (ftnlen)49);
	e_wsle();
	s_wsle(&io___371);
	e_wsle();
	}

/*     READ hkl Miller indices in ort.hkl */

	s_copy(tempo, "ort.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___372);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__1 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L401: */
	s_rsle(&io___373);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

	i__1 = nruns2;
	for (nrun2 = 1; nrun2 <= i__1; ++nrun2) {

	if (ngrid == 3) {
		if (nrun2 == 1) {
		vmin = 8.f;
		vmax = 500.f;
		if (vorth < 500.f) {
			vmax = vorth;
		}
		ntimelim[3] = (vmax - vmin) * 20.f;
		}
		if (nrun2 == 2) {
		if (vorth < 500.f) {
			goto L500;
		}
		vmin = 500.f;
		vmax = 1e3f;
		if (vorth < 1e3f) {
			vmax = vorth;
		}
		ntimelim[3] = (vmax - vmin) * 20.f;
		}
		if (nrun2 == 3) {
		if (vorth < 1e3f) {
			goto L500;
		}
		vmin = 1e3f;
		vmax = 1500.f;
		if (vorth < 1500.f) {
			vmax = vorth;
		}
		ntimelim[3] = (vmax - vmin) * 20.f;
		}
		if (nrun2 == 4) {
		if (vorth < 1500.f) {
			goto L500;
		}
		vmin = 1500.f;
		vmax = 2e3f;
		if (vorth < 2e3f) {
			vmax = vorth;
		}
		ntimelim[3] = (vmax - vmin) * 20.f;
		}
		if (nrun2 == 5) {
		if (vorth < 2e3f) {
			goto L500;
		}
		vmin = 2e3f;
		vmax = 2500.f;
		if (vorth < 2500.f) {
			vmax = vorth;
		}
		ntimelim[3] = (vmax - vmin) * 20.f;
		}
		if (nrun2 == 6) {
		if (vorth < 2500.f) {
			goto L500;
		}
		vmin = 2500.f;
		vmax = 3e3f;
		if (vorth < 3e3f) {
			vmax = vorth;
		}
		ntimelim[3] = (vmax - vmin) * 20.f;
		}
	}

	i__3 = nruns;
	for (nrun = 1; nrun <= i__3; ++nrun) {
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
		rmax = rmaxref;
		rmin = rmax;

/* ...  here starts the loop */

		interest = 0;
		tmax = ntimelim[3] / procs;
		ttmax = ntimelim[3] * 10.f;
		ncells = (integer) ntimelim[3];
		iiseed = 0;
		ntried = 0.f;
		ntriedt = 0.f;
		nout = 0;
		celpre[0] = pstart[0] + delta[0] * 2.f * randi_(&iseed);
		celpre[1] = pstart[1] + delta[1] * 2.f * randi_(&iseed);
		celpre[2] = pstart[2] + delta[2] * 2.f * randi_(&iseed);
		celold[0] = celpre[0];
		celold[1] = celpre[1];
		celold[2] = celpre[2];
		rglob = 1.f;
		nglob = 0;

/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC, */
/* $OMP& RMAX2,A,B,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X, */
/* $OMP& DIFF,DIFF2,DDT,DDQ) */
/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
/* $OMP DO */

		i__4 = ncells;
		for (ncel = 1; ncel <= i__4; ++ncel) {
		if (nout >= 1) {
			goto L496;
		}
		if (interest >= 1) {
			goto L496;
		}
		++iiseed;
		if (iiseed == 1) {
			iseed = ((iseed - ncel * nrun) / 2 << 1) + 1;
		}

L402:

/*     Which parameter to vary ? a or b or c ? */

		ntriedb = 0.f;
		x = randi_(&iseed);
		if (x >= 0.f && x < .33333f) {
			ip = 1;
		}
		if (x >= .33333f && x < .66666f) {
			ip = 2;
		}
		if (x >= .66666f && x <= 1.f) {
			ip = 3;
		}
		celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2 * randi_(&
			iseed);
		ntried += 1.f;
		goto L404;
L403:
		del = deltab * (1.f - ntriedb / cy);
		x = randi_(&iseed);
		if (x >= 0.f && x < .33333f) {
			i__ = 1;
		}
		if (x >= .33333f && x < .66666f) {
			i__ = 2;
		}
		if (x >= .66666f && x < 1.f) {
			i__ = 3;
		}
		celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - 
			.5f) * 2.f;
		ntriedb += 1.f;
L404:
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L405: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}
		dcell_(celpre, cal_1.al, &v1);
		if (ntried > tmax) {
			++nout;
			goto L496;
		}
		if (ntriedb != 0.f) {
			goto L406;
		}
		if (v1 > vmax || v1 < vmin) {
			ntried += -1.f;
			ntriedt += 1.f;
			if (ntriedt > ttmax) {
			++nout;
			goto L496;
			}
			goto L402;
		}

L406:
		calcul1_(&diff, &diff2);
		if (cal_1.nmx > cal_1.ndat10) {
			ntried += -1;
			goto L402;
		}
		if (ntriedb != 0.f) {
			goto L414;
		}

/* ... Rp value satisfying ??? */

		if (diff < rglob || cal_1.lhkl > nglob) {
			rglob = diff;
			nglob = cal_1.lhkl;
			celold[ip - 1] = celpre[ip - 1];
		}
		if (cal_1.lhkl >= nmax) {
			rmax = diff;
			icode = 2;
			if (diff <= rmaxref) {
			icode = 1;
			}
		} else {
			icode = 1;
		}
		if (diff > rmax) {
			goto L417;
		}
		if (cal_1.lhkl < nmax) {
			goto L417;
		}
L414:
		if (diff <= rmax) {
			llhkl = cal_1.lhkl;
			rmax = diff;
			rmax2 = diff2;
			a = celpre[0];
			b = celpre[1];
			c__ = celpre[2];
			v2 = v1;
			if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[1] = b;
			bpar[2] = c__;
			v3 = v1;
			}

/* ... "Refine" that cell (by Monte Carlo too...) */

			pstartb[0] = celpre[0];
			pstartb[1] = celpre[1];
			pstartb[2] = celpre[2];
		}
		if (ntriedb <= ncycles) {
			goto L403;
		}
		rglob = .5f;
		nglob = ndat2;
		if (ip == 1) {
			celold[ip - 1] = a;
		}
		if (ip == 2) {
			celold[ip - 1] = b;
		}
		if (ip == 3) {
			celold[ip - 1] = c__;
		}
		ntriedb = 0.f;
		if (rmax >= rmax0[3]) {
			goto L417;
		}
		if (rmax2 >= .15f) {
			goto L417;
		}
		ipen = cal_1.ndat - llhkl;
		if (ipen > cal_1.nind) {
			goto L417;
		}

/* $OMP CRITICAL(STORE1) */

		++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

		igt += 1.f;
		if (nr == 1) {
			if (igt > 50.f) {
			if (ntried / igt < 1e4f) {
				if (rmax0[3] > .2f) {
				rmax0[3] -= rmax0[3] * .05f;
				s_wsle(&io___376);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, "
					"now Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[3], (
					ftnlen)sizeof(real));
				e_wsle();
				s_wsle(&io___377);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, "
					"now Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[3], (
					ftnlen)sizeof(real));
				e_wsle();
				}
			}
			}
		}

		if (igc > 10000) {
			s_wsle(&io___378);
			do_lio(&c__9, &c__1, "   More than 10000 good cells = ST"
				"OP", (ftnlen)36);
			e_wsle();
			s_wsle(&io___379);
			do_lio(&c__9, &c__1, "   More than 10000 good cells = ST"
				"OP", (ftnlen)36);
			e_wsle();
			--igc;
			++interest;
		}
		cel[igc * 6 - 6] = a;
		cel[igc * 6 - 5] = b;
		cel[igc * 6 - 4] = c__;
		cel[igc * 6 - 3] = 90.f;
		cel[igc * 6 - 2] = 90.f;
		cel[igc * 6 - 1] = 90.f;

/* $OMP END CRITICAL(STORE1) */

/* ... Check for supercell */

		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c__;
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L440: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}
		dcell_(celpre, cal_1.al, &v1);

/* $OMP CRITICAL(STORE2) */

		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		km[igc - 1] = llhkl;
		km2[igc - 1] = cal_1.lhkl;
		ifi[igc - 1] = ifile;
		nsol[igc - 1] = 1;
		vgc[igc - 1] = v1;
		rp[igc - 1] = rmax;
		rp2[igc - 1] = diff;
		if (rp[igc - 1] < rpsmall) {
			rpsmall = rp[igc - 1];
			isee = 1;
		} else {
			isee = 0;
		}
		supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__1);
		brav_(&cal_1.lhkl, ihkl, &ibr);
		ib[igc - 1] = ibr;
		a = cel[igc * 6 - 6];
		b = cel[igc * 6 - 5];
		c__ = cel[igc * 6 - 4];
		v2 = vgc[igc - 1];

/* $OMP END CRITICAL(STORE2) */

/* ... Check for interesting result */

/*      IF(INTEREST.GE.1)GO TO 496 */
		indic = 0;
		bb[2] = a;
		bb[3] = b;
		bb[4] = c__;
		bb[5] = 90.f;
		bb[6] = 90.f;
		bb[7] = 90.f;
		afi[2] = 1.f;
		afi[3] = 1.f;
		afi[4] = 1.f;
		afi[5] = 0.f;
		afi[6] = 0.f;
		afi[7] = 0.f;
		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c__;
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L410: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}

/* $OMP CRITICAL(FOUND) */

		if (rp[igc - 1] < rmi) {
			++interest;
			s_wsfe(&io___380);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			s_wsle(&io___381);
			e_wsle();
			s_wsle(&io___382);
			do_lio(&c__9, &c__1, "=================================="
				"==================\r ===========================",
				 (ftnlen)81);
			e_wsle();
			s_wsle(&io___383);
			e_wsle();
			s_wsle(&io___384);
			do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RES"
				"ULT : Rp < Rmin!", (ftnlen)50);
			e_wsle();
			s_wsle(&io___385);
			e_wsle();
			s_wsle(&io___386);
			do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RES"
				"ULT : Rp < Rmin!", (ftnlen)50);
			e_wsle();

/* ... Refine that cell */

			dcell_(celpre, cal_1.al, &v1);
			calcul2_(&diff, ihkl, th3, &ncalc, &igc);
			celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &
				ddq);
			cncalc[igc - 1] = (real) ncalc;
			if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq)
				;
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
			} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] 
				* 2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
			}
/*      WRITE(20,7000)FM20(IGC) */
/*      WRITE(20,7001)FF20(IGC),DDT,NCALC */
/* 	WRITE(20,*) */
/*      PRINT 7000,FM20(IGC) */
/*      PRINT 7001,FF20(IGC),DDT,NCALC */
/* 	PRINT * */
			iref = 1;
			goto L497;
		} else {

/*  Anyway, calculate the M20 and F20 values */

			dcell_(celpre, cal_1.al, &v1);
			calcul2_(&diff, ihkl, th3, &ncalc, &igc);
			celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &
				ddq);
			cncalc[igc - 1] = (real) ncalc;
			if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq)
				;
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
			} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] 
				* 2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
			}
		}

/* Test if cell already found */

		if (igc > 1) {
			i__5 = igc - 1;
			for (i__ = 1; i__ <= i__5; ++i__) {
			if (ifi[i__ - 1] != ifile) {
				goto L418;
			}
			vdelt = vgc[igc - 1] / 300.f;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
				goto L418;
			}
			adelt = cel[igc * 6 - 6] / 500.f;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			na = 0;
			if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
				na = 1;
			}
			bdelt = cel[igc * 6 - 5] / 500.f;
			bp = cel[igc * 6 - 5] + bdelt;
			bm = cel[igc * 6 - 5] - bdelt;
			nb = 0;
			if (cel[i__ * 6 - 6] > bp || cel[i__ * 6 - 6] < bm) {
				nb = 1;
			}
			cdelt = cel[igc * 6 - 4] / 500.f;
			cp = cel[igc * 6 - 4] + cdelt;
			cm = cel[igc * 6 - 4] - cdelt;
			nc = 0;
			if (cel[i__ * 6 - 6] > cp || cel[i__ * 6 - 6] < cm) {
				nc = 1;
			}
			if (na == 1 && nb == 1 && nc == 1) {
				goto L418;
			}
			++nsol[i__ - 1];
			if (rp[igc - 1] < rp[i__ - 1]) {
				if (isee == 1) {
				s_wsfe(&io___396);
				do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real)
					);
				do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real)
					);
				do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(
					integer));
				e_wsfe();
				}
				km[i__ - 1] = km[igc - 1];
				vgc[i__ - 1] = vgc[igc - 1];
				rp[i__ - 1] = rp[igc - 1];
				cel[i__ * 6 - 6] = cel[igc * 6 - 6];
				cel[i__ * 6 - 5] = cel[igc * 6 - 5];
				cel[i__ * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			if (nsol[i__ - 1] > 5) {
				ntried = tmax + 1.f;
				++nout;
			}
			goto L419;
L418:
			;
			}
			if (iverb == 1) {
			s_wsfe(&io___397);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer))
				;
			e_wsfe();
			}
			if (isee == 1) {
			s_wsfe(&io___398);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			}
L419:
			;
		} else {
			if (iverb == 1) {
			s_wsfe(&io___399);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer))
				;
			e_wsfe();
			}
			if (isee == 1) {
			s_wsfe(&io___400);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			}
		}
L497:

/* $OMP END CRITICAL(FOUND) */

L417:
		rmax = rmaxref;
		if (randi_(&iseed) > escape) {
			celpre[ip - 1] = celold[ip - 1];
		}

/*  END ON MC tests */

L496:
		;
		}

/* $OMP END DO NOWAIT */
/* $OMP END PARALLEL */

/* ... Stop if max limit of Monte Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

		if (interest >= 1) {
		goto L5000;
		}
		killk_(&pressedk);
		if (rmin == rmax) {
		goto L498;
		}
		if (iverb == 1) {
		s_wsle(&io___401);
		e_wsle();
		s_wsle(&io___402);
		do_lio(&c__9, &c__1, "Best result : a = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___403);
		do_lio(&c__9, &c__1, "Best result : b = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[1], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___404);
		do_lio(&c__9, &c__1, "Best result : c = ", (ftnlen)18);
		do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, "V = ", (ftnlen)4);
		do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___405);
		e_wsle();
		datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
		}
		rmin = rmax;
L498:
		if (pressedk) {
		goto L5000;
		}

/*  END ON NRUNS */

/* L499: */
	}
	}

L500:
	if (nsys[4] == 0) {
	goto L600;
	}

/*    Monoclinic case */


	rpsmall = 1.f;
	s_wsle(&io___406);
	do_lio(&c__9, &c__1, "Monoclinic:   Rp     a       b       c       bet  "
		"   V     Nind", (ftnlen)63);
	e_wsle();
	ifile = 5;
	ncycles = 2e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[5] = 90.f;

	if (ngrid == 3) {
	nruns = 20;
	nruns2 = 6;
	pmin = 2.f;
	pmax = 20.f;
	pma[0] = dmax1 * 2.1f;
	pma[1] = dmax1 * 2.1f;
	pma[2] = dmax2 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	if (pma[1] > pmax) {
		pma[1] = pmax;
	}
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	vmon = pma[0] * pma[1] * pma[2];
	if (vmon > 3e3f) {
		vmon = 3e3f;
	}
	for (i__ = 1; i__ <= 3; ++i__) {
		pmi[i__ - 1] = pmin;
		delta[i__ - 1] = (pma[i__ - 1] - pmi[i__ - 1]) / 2.f;
/* L523: */
		pstart[i__ - 1] = pmi[i__ - 1];
	}
	}

	s_wsle(&io___408);
	e_wsle();
	s_wsle(&io___409);
	do_lio(&c__9, &c__1, "Monoclinic Monte Carlo search :", (ftnlen)31);
	e_wsle();
	s_wsle(&io___410);
	do_lio(&c__9, &c__1, " Max(a,b,c), V ", (ftnlen)15);
	do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___411);
	e_wsle();

	if (iverb == 1) {
	s_wsfe(&io___412);
	do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ntimelim[4], (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsle(&io___413);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___414);
	do_lio(&c__9, &c__1, " Rp  Trial number    a   b   c    bet   V  Nin"
		"d Icod", (ftnlen)52);
	e_wsle();
	s_wsle(&io___415);
	e_wsle();
	}

/*     READ hkl Miller indices in mon.hkl */

	s_copy(tempo, "mon.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___416);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L501: */
	s_rsle(&io___417);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

	i__3 = nruns2;
	for (nrun2 = 1; nrun2 <= i__3; ++nrun2) {

	if (ngrid == 3) {
		if (nrun2 == 1) {
		vmin = 8.f;
		vmax = 500.f;
		if (vmon < 500.f) {
			vmax = vmon;
		}
		ntimelim[4] = (vmax - vmin) * 200.f;
		}
		if (nrun2 == 2) {
		if (vmon < 500.f) {
			goto L600;
		}
		vmin = 500.f;
		vmax = 1e3f;
		if (vmon < 1e3f) {
			vmax = vmon;
		}
		ntimelim[4] = (vmax - vmin) * 200.f;
		}
		if (nrun2 == 3) {
		if (vmon < 1e3f) {
			goto L600;
		}
		vmin = 1e3f;
		vmax = 1500.f;
		if (vmon < 1500.f) {
			vmax = vmon;
		}
		ntimelim[4] = (vmax - vmin) * 200.f;
		}
		if (nrun2 == 4) {
		if (vmon < 1500.f) {
			goto L600;
		}
		vmin = 1500.f;
		vmax = 2e3f;
		if (vmon < 2e3f) {
			vmax = vmon;
		}
		ntimelim[4] = (vmax - vmin) * 200.f;
		}
		if (nrun2 == 5) {
		if (vmon < 2e3f) {
			goto L600;
		}
		vmin = 2e3f;
		vmax = 2500.f;
		if (vmon < 2500.f) {
			vmax = vmon;
		}
		ntimelim[4] = (vmax - vmin) * 200.f;
		}
		if (nrun2 == 6) {
		if (vmon < 2500.f) {
			goto L600;
		}
		vmin = 2500.f;
		vmax = 3e3f;
		if (vmon < 3e3f) {
			vmax = vmon;
		}
		ntimelim[4] = (vmax - vmin) * 200.f;
		}
	}

	i__1 = nruns;
	for (nrun = 1; nrun <= i__1; ++nrun) {
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
		rmax = rmaxref;
		rmin = rmax;

/* ...  here starts the loop */

		interest = 0;
		tmax = ntimelim[4] / procs;
		ttmax = ntimelim[4] * 10.f;
		ncells = (integer) ntimelim[4];
		iiseed = 0;
		ntried = 0.f;
		ntriedt = 0.f;
		nout = 0;
		celpre[0] = pstart[0] + delta[0] * 2.f * randi_(&iseed);
		celpre[1] = pstart[1] + delta[1] * 2.f * randi_(&iseed);
		celpre[2] = pstart[2] + delta[2] * 2.f * randi_(&iseed);
		celpre[4] = astart + deltc * 2.f * randi_(&iseed);
		celold[0] = celpre[0];
		celold[1] = celpre[1];
		celold[2] = celpre[2];
		celold[4] = celpre[4];
		rglob = 1.f;
		nglob = 0;
/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,DELD,V1,ICODE,LLHKL,IHKL,TH3, */
/* $OMP& RMAX2,A,B,C,BET,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X, */
/* $OMP& DIFF,DIFF2,NCALC,DDT,DDQ) */
/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
/* $OMP& celpre,celold,rglob,nold,rmin,rmax,bb,afi) */
/* $OMP DO */

		i__4 = ncells;
		for (ncel = 1; ncel <= i__4; ++ncel) {
		if (nout >= 1) {
			goto L596;
		}
		if (interest >= 1) {
			goto L596;
		}
		++iiseed;
		if (iiseed == 1) {
			iseed = ((iseed - ncel * nrun) / 2 << 1) + 1;
		}

L502:

/*     Which parameter to vary ? a or b or c or beta ? */

		ntriedb = 0.f;
		x = randi_(&iseed);
		if (x >= 0.f && x < .25f) {
			ip = 1;
		}
		if (x >= .25f && x < .5f) {
			ip = 2;
		}
		if (x >= .5f && x < .75f) {
			ip = 3;
		}
		if (x >= .75f && x <= 1.f) {
			ip = 5;
		}
		if (ip != 5) {
			celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.f * 
				randi_(&iseed);
		} else {
			celpre[ip - 1] = astart + deltc * 2.f * randi_(&iseed);
		}
		ntried += 1.f;
		goto L504;
L503:
		del = deltab * (1.f - ntriedb / cy);
		deld = deltad * (1.f - ntriedb / cy);
		x = randi_(&iseed);
		if (x >= 0.f && x < .25f) {
			i__ = 1;
		}
		if (x >= .25f && x < .5f) {
			i__ = 2;
		}
		if (x >= .5f && x < .75f) {
			i__ = 3;
		}
		if (x >= .75f && x <= 1.f) {
			i__ = 5;
		}
		if (i__ != 5) {
			celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed)
				 - .5f) * 2.f;
		} else {
			celpre[i__ - 1] = pstartb[i__ - 1] + deld * (randi_(&
				iseed) - .5f) * 2.f;
		}
		ntriedb += 1.f;
L504:
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L505: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}
		dcell_(celpre, cal_1.al, &v1);
		if (ntried > tmax) {
			++nout;
			goto L596;
		}
		if (ntriedb != 0.f) {
			goto L506;
		}
		if (v1 > vmax || v1 < vmin) {
			ntried += -1.f;
			ntriedt += 1.f;
			if (ntriedt > ttmax) {
			++nout;
			goto L596;
			}
			goto L502;
		}

L506:
		calcul1_(&diff, &diff2);
		if (cal_1.nmx > cal_1.ndat10) {
			ntried += -1;
			goto L502;
		}
		if (ntriedb != 0.f) {
			goto L514;
		}

/* ... Rp value satisfying ??? */

		if (diff < rglob || cal_1.lhkl > nglob) {
			rglob = diff;
			nglob = cal_1.lhkl;
			celold[ip - 1] = celpre[ip - 1];
		}
		if (cal_1.lhkl >= nmax) {
			rmax = diff;
			icode = 2;
			if (diff <= rmaxref) {
			icode = 1;
			}
		} else {
			icode = 1;
		}
		if (diff > rmax) {
			goto L517;
		}
		if (cal_1.lhkl < nmax) {
			goto L517;
		}
L514:
		if (diff <= rmax) {
			llhkl = cal_1.lhkl;
			rmax = diff;
			rmax2 = diff2;
			a = celpre[0];
			b = celpre[1];
			c__ = celpre[2];
			bet = celpre[4];
			v2 = v1;
			if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[1] = b;
			bpar[2] = c__;
			bpar[4] = bet;
			v3 = v1;
			}

/* ... "Refine" that cell (by Monte Carlo too...) */

			pstartb[0] = celpre[0];
			pstartb[1] = celpre[1];
			pstartb[2] = celpre[2];
			pstartb[4] = celpre[4];
		}
		if (ntriedb <= ncycles) {
			goto L503;
		}
		rglob = .5f;
		nglob = ndat2;
		if (ip == 1) {
			celold[ip - 1] = a;
		}
		if (ip == 2) {
			celold[ip - 1] = b;
		}
		if (ip == 3) {
			celold[ip - 1] = c__;
		}
		if (ip == 5) {
			celold[ip - 1] = bet;
		}
		ntriedb = 0.f;
		if (rmax >= rmax0[4]) {
			goto L517;
		}
		if (rmax2 >= .15f) {
			goto L517;
		}
		ipen = cal_1.ndat - llhkl;
		if (ipen > cal_1.nind) {
			goto L517;
		}

/* $OMP CRITICAL(STORE1) */

		++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

		igt += 1.f;
		if (nr == 1) {
			if (igt > 50.f) {
			if (ntried / igt < 1e5f) {
				if (rmax0[4] > .2f) {
				rmax0[4] -= rmax0[4] * .05f;
				s_wsle(&io___420);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, "
					"now Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[4], (
					ftnlen)sizeof(real));
				e_wsle();
				s_wsle(&io___421);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, "
					"now Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[4], (
					ftnlen)sizeof(real));
				e_wsle();
				}
			}
			}
		}

		if (igc > 10000) {
			s_wsle(&io___422);
			do_lio(&c__9, &c__1, "   More than 10000 good cells = ST"
				"OP", (ftnlen)36);
			e_wsle();
			s_wsle(&io___423);
			do_lio(&c__9, &c__1, "   More than 10000 good cells = ST"
				"OP", (ftnlen)36);
			e_wsle();
			--igc;
			++interest;
/*      GO TO 5000 */
		}
		cel[igc * 6 - 6] = a;
		cel[igc * 6 - 5] = b;
		cel[igc * 6 - 4] = c__;
		cel[igc * 6 - 3] = 90.f;
		cel[igc * 6 - 2] = bet;
		cel[igc * 6 - 1] = 90.f;

/* $OMP END CRITICAL(STORE1) */

/* ... Check for supercell */

		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c__;
		celpre[4] = bet;
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L540: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}
		dcell_(celpre, cal_1.al, &v1);

/* $OMP CRITICAL(STORE2) */

		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		km[igc - 1] = llhkl;
		km2[igc - 1] = cal_1.lhkl;
		ifi[igc - 1] = ifile;
		nsol[igc - 1] = 1;
		vgc[igc - 1] = v1;
		rp[igc - 1] = rmax;
		rp2[igc - 1] = diff;
		if (rp[igc - 1] < rpsmall) {
			rpsmall = rp[igc - 1];
			isee = 1;
		} else {
			isee = 0;
		}
		supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__1);
		brav_(&cal_1.lhkl, ihkl, &ibr);
		ib[igc - 1] = ibr;
		a = cel[igc * 6 - 6];
		b = cel[igc * 6 - 5];
		c__ = cel[igc * 6 - 4];
		v2 = vgc[igc - 1];

/* $OMP END CRITICAL(STORE2) */

/* ... Check for interesting result */

/*      IF(INTEREST.GE.1)GO TO 596 */
		indic = 0;
		bb[2] = a;
		bb[3] = b;
		bb[4] = c__;
		bb[5] = 90.f;
		bb[6] = bet;
		bb[7] = 90.f;
		afi[2] = 1.f;
		afi[3] = 1.f;
		afi[4] = 1.f;
		afi[5] = 0.f;
		afi[6] = 1.f;
		afi[7] = 0.f;
		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c__;
		celpre[4] = bet;
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L510: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}

/* $OMP CRITICAL(FOUND) */

		if (rp[igc - 1] < rmi) {
			++interest;
			s_wsfe(&io___424);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			s_wsle(&io___425);
			e_wsle();
			s_wsle(&io___426);
			do_lio(&c__9, &c__1, "=================================="
				"==================\r ===========================",
				 (ftnlen)81);
			e_wsle();
			s_wsle(&io___427);
			e_wsle();
			s_wsle(&io___428);
			do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RES"
				"ULT : Rp < Rmin!", (ftnlen)50);
			e_wsle();
			s_wsle(&io___429);
			e_wsle();
			s_wsle(&io___430);
			do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RES"
				"ULT : Rp < Rmin!", (ftnlen)50);
			e_wsle();

/* ... Refine that cell */

			dcell_(celpre, cal_1.al, &v1);
			calcul2_(&diff, ihkl, th3, &ncalc, &igc);
			celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &
				ddq);
			cncalc[igc - 1] = (real) ncalc;
			if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq)
				;
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
			} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] 
				* 2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
			}
/*      WRITE(20,7000)FM20(IGC) */
/*      WRITE(20,7001)FF20(IGC),DDT,NCALC */
/* 	WRITE(20,*) */
/*      PRINT 7000,FM20(IGC) */
/*      PRINT 7001,FF20(IGC),DDT,NCALC */
/* 	PRINT * */
			iref = 1;
			goto L597;
		} else {

/*  Anyway, calculate the M20 and F20 values */

			dcell_(celpre, cal_1.al, &v1);
			calcul2_(&diff, ihkl, th3, &ncalc, &igc);
			celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &
				ddq);
			cncalc[igc - 1] = (real) ncalc;
			if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq)
				;
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
			} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] 
				* 2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
			}
		}

/* Test if cell already found */


		if (igc > 1) {
			i__5 = igc - 1;
			for (i__ = 1; i__ <= i__5; ++i__) {
			if (ifi[i__ - 1] != ifile) {
				goto L518;
			}
			vdelt = vgc[igc - 1] / 300.f;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
				goto L518;
			}
			bdelt = cel[igc * 6 - 5] / 500.f;
			bp = cel[igc * 6 - 5] + bdelt;
			bm = cel[igc * 6 - 5] - bdelt;
			if (cel[i__ * 6 - 5] > bp || cel[i__ * 6 - 5] < bm) {
				goto L518;
			}
			betdelt = cel[igc * 6 - 2] / 500.f;
			betp = cel[igc * 6 - 2] + betdelt;
			betm = cel[igc * 6 - 2] - betdelt;
			if (cel[i__ * 6 - 2] > betp || cel[i__ * 6 - 2] < 
				betm) {
				goto L518;
			}
			adelt = cel[igc * 6 - 6] / 500.f;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			na = 0;
			if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
				na = 1;
			}
			cdelt = cel[igc * 6 - 4] / 500.f;
			cp = cel[igc * 6 - 4] + cdelt;
			cm = cel[igc * 6 - 4] - cdelt;
			nc = 0;
			if (cel[i__ * 6 - 6] > cp || cel[i__ * 6 - 6] < cm) {
				nc = 1;
			}
			if (na == 1 && nc == 1) {
				goto L518;
			}
			++nsol[i__ - 1];
			if (rp[igc - 1] < rp[i__ - 1]) {
				if (isee == 1) {
				s_wsfe(&io___434);
				do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real)
					);
				do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real)
					);
				do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(
					integer));
				e_wsfe();
				}
				km[i__ - 1] = km[igc - 1];
				vgc[i__ - 1] = vgc[igc - 1];
				rp[i__ - 1] = rp[igc - 1];
				cel[i__ * 6 - 6] = cel[igc * 6 - 6];
				cel[i__ * 6 - 5] = cel[igc * 6 - 5];
				cel[i__ * 6 - 4] = cel[igc * 6 - 4];
				cel[i__ * 6 - 2] = cel[igc * 6 - 2];
			}
			--igc;
			if (nsol[i__ - 1] > 5) {
				ntried = tmax + 1.f;
				++nout;
			}
			goto L519;
L518:
			;
			}
			if (iverb == 1) {
			s_wsfe(&io___435);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer))
				;
			e_wsfe();
			}
			if (isee == 1) {
			s_wsfe(&io___436);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			}
L519:
			;
		} else {
			if (iverb == 1) {
			s_wsfe(&io___437);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer))
				;
			e_wsfe();
			}
			if (isee == 1) {
			s_wsfe(&io___438);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			}
		}
L597:

/* $OMP END CRITICAL(FOUND) */

L517:
		rmax = rmaxref;
		if (randi_(&iseed) > escape) {
			celpre[ip - 1] = celold[ip - 1];
		}

/*  END ON MC tests */

L596:
		;
		}

/* $OMP END DO NOWAIT */
/* $OMP END PARALLEL */

/* ... Stop if max limit of Monte Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

		if (interest >= 1) {
		goto L5000;
		}
		killk_(&pressedk);
		if (rmin == rmax) {
		goto L598;
		}
		if (iverb == 1) {
		s_wsle(&io___439);
		e_wsle();
		s_wsle(&io___440);
		do_lio(&c__9, &c__1, "Best result : a =    ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___441);
		do_lio(&c__9, &c__1, "Best result : b =    ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[1], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___442);
		do_lio(&c__9, &c__1, "Best result : c =    ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___443);
		do_lio(&c__9, &c__1, "Best result : beta = ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[4], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, "  V = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___444);
		e_wsle();
		datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
		}
		rmin = rmax;
L598:
		if (pressedk) {
		goto L5000;
		}

/*  END ON NRUNS */

/* L599: */
	}
	}

L600:
	if (nsys[5] == 0) {
	goto L700;
	}

/*    Triclinic case */


	rpsmall = 1.f;
	s_wsle(&io___445);
	do_lio(&c__9, &c__1, "Triclinic:    Rp\r     a       b       c       alp"
		"    bet    gam     V     Nind", (ftnlen)78);
	e_wsle();
	ifile = 6;
	ncycles = 5e3f;
	cy = ncycles * 1.1f;

	if (ngrid == 3) {
	nruns = 20;
	nruns2 = 8;
	pmin = 2.f;
	pmax = 20.f;
	pma[0] = dmax1 * 1.5f;
	pma[1] = dmax2 * 1.5f;
	pma[2] = dmax3 * 1.5f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	if (pma[1] > pmax) {
		pma[1] = pmax;
	}
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	vtric = pma[0] * pma[1] * pma[2];
	if (vtric > 2e3f) {
		vtric = 2e3f;
	}
	vmax = vtric;
	for (i__ = 1; i__ <= 3; ++i__) {
		pmi[i__ - 1] = pmin;
		delta[i__ - 1] = (pma[i__ - 1] - pmi[i__ - 1]) / 2.f;
/* L623: */
		pstart[i__ - 1] = pmi[i__ - 1];
	}
	}

	s_wsle(&io___447);
	e_wsle();
	s_wsle(&io___448);
	do_lio(&c__9, &c__1, "Triclinic Monte Carlo search :", (ftnlen)30);
	e_wsle();
	s_wsle(&io___449);
	do_lio(&c__9, &c__1, " Max(a,b,c), V ", (ftnlen)15);
	do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[1], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&vtric, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___450);
	e_wsle();

	if (iverb == 1) {
	s_wsfe(&io___451);
	do_fio(&c__1, (char *)&nrun, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ntimelim[5], (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsle(&io___452);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___453);
	do_lio(&c__9, &c__1, " Rp Trial number  a  b  c  alp bet gam  V  Nin"
		"d icod", (ftnlen)52);
	e_wsle();
	s_wsle(&io___454);
	e_wsle();
	}

/*     READ hkl Miller indices in tri.hkl */

	s_copy(tempo, "tri.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___455);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__1 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L601: */
	s_rsle(&io___456);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

	i__1 = nruns2;
	for (nrun2 = 1; nrun2 <= i__1; ++nrun2) {

	if (ngrid == 3) {
		if (nrun2 == 1) {
		vmin = 8.f;
		vmax = 250.f;
		if (vtric < 250.f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 2) {
		if (vtric < 250.f) {
			goto L700;
		}
		vmin = 250.f;
		vmax = 500.f;
		if (vtric < 500.f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 3) {
		if (vtric < 500.f) {
			goto L700;
		}
		vmin = 500.f;
		vmax = 750.f;
		if (vtric < 750.f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 4) {
		if (vtric < 750.f) {
			goto L700;
		}
		vmin = 750.f;
		vmax = 1e3f;
		if (vtric < 1e3f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 5) {
		if (vtric < 1e3f) {
			goto L700;
		}
		vmin = 1e3f;
		vmax = 1250.f;
		if (vtric < 1250.f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 6) {
		if (vtric < 1250.f) {
			goto L700;
		}
		vmin = 1250.f;
		vmax = 1500.f;
		if (vtric < 1500.f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 7) {
		if (vtric < 1500.f) {
			goto L700;
		}
		vmin = 1500.f;
		vmax = 1750.f;
		if (vtric < 1750.f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
		if (nrun2 == 8) {
		if (vtric < 1750.f) {
			goto L700;
		}
		vmin = 1750.f;
		vmax = 2e3f;
		if (vtric < 2e3f) {
			vmax = vtric;
		}
		ntimelim[5] = (vmax - vmin) * 4e3f;
		}
	}

	i__3 = nruns;
	for (nrun = 1; nrun <= i__3; ++nrun) {
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
		rmax = rmaxref;
		rmin = rmax;

/* ...  here starts the loop */

		interest = 0;
		tmax = ntimelim[5] / procs;
		ttmax = ntimelim[5] * 10.f;
		ncells = (integer) ntimelim[5];
		iiseed = 0;
		ntried = 0.f;
		ntriedt = 0.f;
		nout = 0;

/* ...  here starts the loop */

		celpre[0] = pstart[0] + delta[0] * 2.f * randi_(&iseed);
		celpre[1] = pstart[1] + delta[1] * 2.f * randi_(&iseed);
		celpre[2] = pstart[2] + delta[2] * 2.f * randi_(&iseed);
		celpre[3] = astartt[0] + deltct[0] * 2.f * randi_(&iseed);
		celpre[4] = astartt[1] + deltct[1] * 2.f * randi_(&iseed);
		celpre[5] = astartt[2] + deltct[2] * 2.f * randi_(&iseed);
		celold[0] = celpre[0];
		celold[1] = celpre[1];
		celold[2] = celpre[2];
		celold[3] = celpre[3];
		celold[4] = celpre[4];
		celold[5] = celpre[5];
		rglob = 1.f;
		nglob = 0;

/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,DELD,V1,ICODE,LLHKL,IHKL,TH3, */
/* $OMP& RMAX2,A,B,C,ALP,BET,GAM,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X, */
/* $OMP& DIFF,DIFF2,IP2,ANG,NCALC,DDT,DDQ) */
/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
/* $OMP DO */

		i__4 = ncells;
		for (ncel = 1; ncel <= i__4; ++ncel) {
		if (nout >= 1) {
			goto L696;
		}
		if (interest >= 1) {
			goto L696;
		}
		++iiseed;
		if (iiseed == 1) {
			iseed = ((iseed - ncel * nrun) / 2 << 1) + 1;
		}

L602:

/*     Which parameter to vary ? a or b or c or alpha or beta or gamma ? */

		ntriedb = 0.f;
		x = randi_(&iseed);
		if (x >= 0.f && x < .16666f) {
			ip = 1;
		}
		if (x >= .16666f && x < .33333f) {
			ip = 2;
		}
		if (x >= .33333f && x < .5f) {
			ip = 3;
		}
		if (x >= .5f && x < .66666f) {
			ip = 4;
		}
		if (x >= .66666f && x < .83333f) {
			ip = 5;
		}
		if (x >= .83333f && x <= 1.f) {
			ip = 6;
		}
		if (ip != 4 && ip != 5 && ip != 6) {
			celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.f * 
				randi_(&iseed);
		} else {
			celpre[ip - 1] = astartt[ip - 4] + deltct[ip - 4] * 2.f * 
				randi_(&iseed);
		}
		ang = celpre[3] + celpre[4] + celpre[5];
		if (ang >= 360.f && ang <= 180.f) {
			goto L696;
		}
		ntried += 1.f;
		goto L604;
L603:
		del = deltab * (1.f - ntriedb / cy);
		deld = deltad * (1.f - ntriedb / cy);
		x = randi_(&iseed);
		if (x >= 0.f && x < .16666f) {
			ip2 = 1;
		}
		if (x >= .16666f && x < .33333f) {
			ip2 = 2;
		}
		if (x >= .33333f && x < .5f) {
			ip2 = 3;
		}
		if (x >= .5f && x < .66666f) {
			ip2 = 4;
		}
		if (x >= .66666f && x < .83333f) {
			ip2 = 5;
		}
		if (x >= .83333f && x <= 1.f) {
			ip2 = 6;
		}
		if (ip2 != 4 && ip2 != 5 && ip2 != 6) {
			celpre[ip2 - 1] = pstartb[ip2 - 1] + del * (randi_(&iseed)
				 - .5f) * 2.f;
		} else {
			celpre[ip2 - 1] = pstartb[ip2 - 1] + deld * (randi_(&
				iseed) - .5f) * 2.f;
		}
		ang = celpre[3] + celpre[4] + celpre[5];
		if (ang >= 360.f && ang <= 180.f) {
			goto L603;
		}
		ntriedb += 1.f;
L604:
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L605: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}
		dcell_(celpre, cal_1.al, &v1);
		if (ntried > tmax) {
			++nout;
			goto L696;
		}
		if (ntriedb != 0.f) {
			goto L606;
		}
		if (v1 > vmax || v1 < vmin) {
			ntried += -1.f;
			ntriedt += 1.f;
			if (ntriedt > ttmax) {
			++nout;
			goto L696;
			}
			goto L602;
		}

L606:
		calcul1_(&diff, &diff2);
		if (cal_1.nmx > cal_1.ndat10) {
			ntried += -1;
			goto L602;
		}
		if (ntriedb != 0.f) {
			goto L614;
		}

/* ... Rp value satisfying ??? */

		if (diff < rglob || cal_1.lhkl > nglob) {
			rglob = diff;
			nglob = cal_1.lhkl;
			celold[ip - 1] = celpre[ip - 1];
		}
		if (cal_1.lhkl >= nmax) {
			rmax = diff;
			icode = 2;
			if (diff <= rmaxref) {
			icode = 1;
			}
		} else {
			icode = 1;
		}
		if (diff > rmax) {
			goto L617;
		}
		if (cal_1.lhkl < nmax) {
			goto L617;
		}
L614:
		if (diff <= rmax) {
			llhkl = cal_1.lhkl;
			rmax = diff;
			rmax2 = diff2;
			a = celpre[0];
			b = celpre[1];
			c__ = celpre[2];
			alp = celpre[3];
			bet = celpre[4];
			gam = celpre[5];
			v2 = v1;
			if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[1] = b;
			bpar[2] = c__;
			bpar[3] = alp;
			bpar[4] = bet;
			bpar[5] = gam;
			v3 = v1;
			}

/* ... "Refine" that cell (by Monte Carlo too...) */

			pstartb[0] = celpre[0];
			pstartb[1] = celpre[1];
			pstartb[2] = celpre[2];
			pstartb[3] = celpre[3];
			pstartb[4] = celpre[4];
			pstartb[5] = celpre[5];
		}
		if (ntriedb <= ncycles) {
			goto L603;
		}
		rglob = .5f;
		nglob = ndat2;
		if (ip == 1) {
			celold[ip - 1] = a;
		}
		if (ip == 2) {
			celold[ip - 1] = b;
		}
		if (ip == 3) {
			celold[ip - 1] = c__;
		}
		if (ip == 4) {
			celold[ip - 1] = alp;
		}
		if (ip == 5) {
			celold[ip - 1] = bet;
		}
		if (ip == 6) {
			celold[ip - 1] = gam;
		}
		ntriedb = 0.f;
		if (rmax >= rmax0[5]) {
			goto L617;
		}
		if (rmax2 >= .15f) {
			goto L617;
		}
		ipen = cal_1.ndat - llhkl;
		if (ipen > cal_1.nind) {
			goto L617;
		}

/* $OMP CRITICAL(STORE1) */

		++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

		igt += 1.f;
		if (nr == 1) {
			if (igt > 50.f) {
			if (ntried / igt < 1e5f) {
				if (rmax0[5] > .2f) {
				rmax0[5] -= rmax0[5] * .05f;
				s_wsle(&io___461);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, "
					"now Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[5], (
					ftnlen)sizeof(real));
				e_wsle();
				s_wsle(&io___462);
				do_lio(&c__9, &c__1, "  Rmax reduced by 5%, "
					"now Rmax = ", (ftnlen)33);
				do_lio(&c__4, &c__1, (char *)&rmax0[5], (
					ftnlen)sizeof(real));
				e_wsle();
				}
			}
			}
		}

		if (igc > 10000) {
			s_wsle(&io___463);
			do_lio(&c__9, &c__1, "   More than 10000 good cells = ST"
				"OP", (ftnlen)36);
			e_wsle();
			s_wsle(&io___464);
			do_lio(&c__9, &c__1, "   More than 10000 good cells = ST"
				"OP", (ftnlen)36);
			e_wsle();
			--igc;
			++interest;
/*      GO TO 5000 */
		}
		cel[igc * 6 - 6] = a;
		cel[igc * 6 - 5] = b;
		cel[igc * 6 - 4] = c__;
		cel[igc * 6 - 3] = alp;
		cel[igc * 6 - 2] = bet;
		cel[igc * 6 - 1] = gam;

/* $OMP END CRITICAL(STORE1) */

/* ... Check for supercell */

		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c__;
		celpre[3] = alp;
		celpre[4] = bet;
		celpre[5] = gam;
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L640: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}
		dcell_(celpre, cal_1.al, &v1);

/* $OMP CRITICAL(STORE2) */

		calcul2_(&diff, ihkl, th3, &ncalc, &igc);
		km[igc - 1] = llhkl;
		km2[igc - 1] = cal_1.lhkl;
		ifi[igc - 1] = ifile;
		nsol[igc - 1] = 1;
		vgc[igc - 1] = v1;
		rp[igc - 1] = rmax;
		rp2[igc - 1] = diff;
		if (rp[igc - 1] < rpsmall) {
			rpsmall = rp[igc - 1];
			isee = 1;
		} else {
			isee = 0;
		}
		supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__1);
		brav_(&cal_1.lhkl, ihkl, &ibr);
		ib[igc - 1] = ibr;
		a = cel[igc * 6 - 6];
		b = cel[igc * 6 - 5];
		c__ = cel[igc * 6 - 4];
		v2 = vgc[igc - 1];

/* $OMP END CRITICAL(STORE2) */

/* ... Check for interesting result */

/*      IF(INTEREST.GE.1)GO TO 696 */
		indic = 0;
		bb[2] = a;
		bb[3] = b;
		bb[4] = c__;
		bb[5] = alp;
		bb[6] = bet;
		bb[7] = gam;
		afi[2] = 1.f;
		afi[3] = 1.f;
		afi[4] = 1.f;
		afi[5] = 1.f;
		afi[6] = 1.f;
		afi[7] = 1.f;
		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c__;
		celpre[3] = alp;
		celpre[4] = bet;
		celpre[5] = gam;
		for (i__ = 1; i__ <= 3; ++i__) {
			for (j = 1; j <= 3; ++j) {
/* L610: */
			cal_1.al[i__ + j * 3 - 4] = 0.f;
			}
		}

/* $OMP CRITICAL(FOUND) */

		if (rp[igc - 1] < rmi) {
			++interest;
			s_wsfe(&io___465);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&alp, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&gam, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			s_wsle(&io___466);
			e_wsle();
			s_wsle(&io___467);
			do_lio(&c__9, &c__1, "=================================="
				"==================\r ===========================",
				 (ftnlen)81);
			e_wsle();
			s_wsle(&io___468);
			e_wsle();
			s_wsle(&io___469);
			do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RES"
				"ULT : Rp < Rmin!", (ftnlen)50);
			e_wsle();
			s_wsle(&io___470);
			e_wsle();
			s_wsle(&io___471);
			do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RES"
				"ULT : Rp < Rmin!", (ftnlen)50);
			e_wsle();

/* ... Refine that cell */

			dcell_(celpre, cal_1.al, &v1);
			calcul2_(&diff, ihkl, th3, &ncalc, &igc);
			celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &
				ddq);
			cncalc[igc - 1] = (real) ncalc;
			if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq)
				;
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
			} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] 
				* 2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
			}
/*      WRITE(20,7000)FM20(IGC) */
/*      WRITE(20,7001)FF20(IGC),DDT,NCALC */
/* 	WRITE(20,*) */
/*      PRINT 7000,FM20(IGC) */
/*      PRINT 7001,FF20(IGC),DDT,NCALC */
/* 	PRINT * */
			iref = 1;
			goto L697;
		} else {

/*  Anyway, calculate the M20 and F20 values */

			dcell_(celpre, cal_1.al, &v1);
			calcul2_(&diff, ihkl, th3, &ncalc, &igc);
			celref2_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &
				ddq);
			cncalc[igc - 1] = (real) ncalc;
			if (cal_1.ndat >= 20) {
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq)
				;
			ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
			} else {
			pndat = (real) cal_1.ndat;
			fm20[igc - 1] = qo[cal_1.ndat - 1] / (cncalc[igc - 1] 
				* 2.f * ddq);
			ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
			}
		}

/* Test if cell already found */


		if (igc > 1) {
			i__5 = igc - 1;
			for (i__ = 1; i__ <= i__5; ++i__) {
			if (ifi[i__ - 1] != ifile) {
				goto L618;
			}
			vdelt = vgc[igc - 1] / 300.f;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
				goto L618;
			}
			adelt = cel[igc * 6 - 6] / 500.f;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			na = 0;
			if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
				na = 1;
			}
			bdelt = cel[igc * 6 - 5] / 500.f;
			bp = cel[igc * 6 - 5] + bdelt;
			bm = cel[igc * 6 - 5] - bdelt;
			nb = 0;
			if (cel[i__ * 6 - 6] > bp || cel[i__ * 6 - 6] < bm) {
				nb = 1;
			}
			cdelt = cel[igc * 6 - 4] / 500.f;
			cp = cel[igc * 6 - 4] + cdelt;
			cm = cel[igc * 6 - 4] - cdelt;
			nc = 0;
			if (cel[i__ * 6 - 6] > cp || cel[i__ * 6 - 6] < cm) {
				nc = 1;
			}
			if (na == 1 && nb == 1 && nc == 1) {
				goto L618;
			}
			na = 0;
			if (cel[i__ * 6 - 5] > ap || cel[i__ * 6 - 5] < am) {
				na = 1;
			}
			nb = 0;
			if (cel[i__ * 6 - 5] > bp || cel[i__ * 6 - 5] < bm) {
				nb = 1;
			}
			nc = 0;
			if (cel[i__ * 6 - 5] > cp || cel[i__ * 6 - 5] < cm) {
				nc = 1;
			}
			if (na == 1 && nb == 1 && nc == 1) {
				goto L618;
			}
			++nsol[i__ - 1];
			if (rp[igc - 1] < rp[i__ - 1]) {
				if (isee == 1) {
				s_wsfe(&io___472);
				do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real)
					);
				do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real)
					);
				do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&alp, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&gam, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(
					real));
				do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(
					integer));
				e_wsfe();
				}
/*     1WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN */
				km[i__ - 1] = km[igc - 1];
				vgc[i__ - 1] = vgc[igc - 1];
				rp[i__ - 1] = rp[igc - 1];
				cel[i__ * 6 - 6] = cel[igc * 6 - 6];
				cel[i__ * 6 - 5] = cel[igc * 6 - 5];
				cel[i__ * 6 - 4] = cel[igc * 6 - 4];
				cel[i__ * 6 - 3] = cel[igc * 6 - 3];
				cel[i__ * 6 - 2] = cel[igc * 6 - 2];
				cel[i__ * 6 - 1] = cel[igc * 6 - 1];
			}
			--igc;
			if (nsol[i__ - 1] > 5) {
				ntried = tmax + 1.f;
				++nout;
			}
			goto L619;
L618:
			;
			}
			if (iverb == 1) {
			s_wsfe(&io___473);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&alp, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&gam, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer))
				;
			e_wsfe();
			}
			if (isee == 1) {
			s_wsfe(&io___474);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&alp, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&gam, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			}
L619:
			;
		} else {
			if (iverb == 1) {
			s_wsfe(&io___475);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(
				doublereal));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&alp, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&gam, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer))
				;
			e_wsfe();
			}
			if (isee == 1) {
			s_wsfe(&io___476);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&alp, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&gam, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
			}
		}
L697:

/* $OMP END CRITICAL(FOUND) */

L617:
		rmax = rmaxref;
		if (randi_(&iseed) > escape) {
			celpre[ip - 1] = celold[ip - 1];
		}

/*  END ON MC tests */

L696:
		;
		}

/* $OMP END DO NOWAIT */
/* $OMP END PARALLEL */

/* ... Stop if max limit of Monte Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

		if (interest >= 1) {
		goto L5000;
		}
		killk_(&pressedk);
/* L616: */
		if (rmin == rmax) {
		goto L698;
		}
		if (iverb == 1) {
		s_wsle(&io___477);
		e_wsle();
		s_wsle(&io___478);
		do_lio(&c__9, &c__1, "Best result : a =    ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___479);
		do_lio(&c__9, &c__1, "Best result : b =    ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[1], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___480);
		do_lio(&c__9, &c__1, "Best result : c =    ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___481);
		do_lio(&c__9, &c__1, "Best result : alph = ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[3], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___482);
		do_lio(&c__9, &c__1, "Best result : beta = ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[4], (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___483);
		do_lio(&c__9, &c__1, "Best result : gamm = ", (ftnlen)21);
		do_lio(&c__4, &c__1, (char *)&bpar[5], (ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, "  V = ", (ftnlen)6);
		do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
		e_wsle();
		s_wsle(&io___484);
		e_wsle();
		datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
		}
		rmin = rmax;
L698:
		if (pressedk) {
		goto L5000;
		}

/*  END ON NRUNS */

/* L699: */
	}
	}


/*   SECOND PART OF McMAILLE : GRID SEARCH... */

L700:
	if (ngrid == 0) {
	goto L5000;
	}
	if (nblack == 0 && ngrid == 3) {
	goto L5000;
	}
	s_wsle(&io___485);
	e_wsle();
	s_wsle(&io___486);
	do_lio(&c__9, &c__1, "Grid search :", (ftnlen)13);
	e_wsle();

/* ...  Cell generation by systematic grid */


	if (nsys[0] == 0) {
	goto L1200;
	}

/*    Cubic case */

	rpsmall = 1.f;
	s_wsle(&io___487);
	do_lio(&c__9, &c__1, "Cubic:        Rp     a       V     Nind", (ftnlen)
		39);
	e_wsle();
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 1;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.f;
	ncycles = 200.f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 90.f;
	if (iverb == 1) {
	s_wsle(&io___488);
	e_wsle();
	s_wsle(&io___489);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___490);
	do_lio(&c__9, &c__1, "Grid search :", (ftnlen)13);
	e_wsle();
	s_wsle(&io___491);
	do_lio(&c__9, &c__1, "   Results in cubic :", (ftnlen)21);
	e_wsle();
	s_wsle(&io___492);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	}

	if (ngrid == 3) {
	pmin = 2.f;
	pmax = dmax1 * 3.1f;
	pmi[0] = pmin;
	pma[0] = pmax;
	vmin = 8.f;
	vmax = pmax * pmax * pmax;
	if (iverb == 1) {
		s_wsle(&io___493);
		do_lio(&c__9, &c__1, " Max a, V ", (ftnlen)10);
		do_lio(&c__4, &c__1, (char *)&pmax, (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
		e_wsle();
	}
	}

	if (iverb == 1) {
	s_wsle(&io___494);
	do_lio(&c__9, &c__1, " Rp  Trial number    a         V  Nind Icod", (
		ftnlen)43);
	e_wsle();
	s_wsle(&io___495);
	e_wsle();
	}

/*     READ hkl Miller indices in cub.hkl */

	s_copy(tempo, "cub.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___496);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 6;
	if (cal_1.nhkl0 > 400) {
	cal_1.nhkl0 = 400;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1101: */
	s_rsle(&io___497);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
L1102:
	ntriedb = 0.f;
	celpre[0] += spar;
	ntried += 1.f;
	goto L1104;
L1103:
	del = deltab * (1.f - ntriedb / cy);
	celpre[0] = pstartb[0] + del * (randi_(&iseed) - .5f) * 2.f;
	ntriedb += 1.f;
L1104:
	celpre[1] = celpre[0];
	celpre[2] = celpre[0];
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1105: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.f) {
	goto L1116;
	}
	if (ntriedb != 0.f) {
	goto L1106;
	}
	if (v1 > vmax || v1 < vmin) {
	ntried += -1.f;
	goto L1102;
	}

L1106:
	calcul1_(&diff, &diff2);
	if (cal_1.nmx > cal_1.ndat10) {
	ntried += -1;
	goto L1102;
	}
	if (ntriedb != 0.f) {
	goto L1114;
	}

/* ... Rp value satisfying ??? */

	if (cal_1.lhkl >= nmax) {
	rmax = diff;
	icode = 2;
	if (diff <= rmaxref) {
		icode = 1;
	}
	} else {
	icode = 1;
	}
	celold[0] = celpre[0];
	if (diff > rmax) {
	goto L1117;
	}
	if (cal_1.lhkl < nmax) {
	goto L1117;
	}
L1114:
	if (diff <= rmax) {
	llhkl = cal_1.lhkl;
	rmax = diff;
	rmax2 = diff2;
	a = celpre[0];
	v2 = v1;
	if (diff < rmin) {
		rmin = diff;
		bpar[0] = a;
		v3 = v1;
	}

/* ... "Refine" that cell (by Monte Carlo too...) */

	pstartb[0] = celpre[0];
	}
	if (ntriedb <= ncycles) {
	goto L1103;
	}
	ntriedb = 0.f;
	if (rmax >= rmax0[0]) {
	goto L1117;
	}
	if (rmax2 >= .15f) {
	goto L1117;
	}
	ipen = cal_1.ndat - llhkl;
	if (ipen > cal_1.nind) {
	goto L1117;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.f;
	if (nr == 1) {
	if (igt > 50.f) {
		if (ntried / igt < 100.f) {
		if (rmax0[0] > .1f) {
			rmax0[0] -= rmax0[0] * .05f;
			s_wsle(&io___498);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[0], (ftnlen)sizeof(
				real));
			e_wsle();
			s_wsle(&io___499);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[0], (ftnlen)sizeof(
				real));
			e_wsle();
		}
		}
	}
	}

	if (igc > 10000) {
	s_wsle(&io___500);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___501);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	--igc;
	goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = a;
	cel[igc * 6 - 3] = 90.f;
	cel[igc * 6 - 2] = 90.f;
	cel[igc * 6 - 1] = 90.f;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = a;
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1140: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	km[igc - 1] = llhkl;
	km2[igc - 1] = cal_1.lhkl;
	ifi[igc - 1] = ifile;
	nsol[igc - 1] = 1;
	vgc[igc - 1] = v1;
	rp[igc - 1] = rmax;
	rp2[igc - 1] = diff;
	if (rp[igc - 1] < rpsmall) {
	rpsmall = rp[igc - 1];
	isee = 1;
	} else {
	isee = 0;
	}
	supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__3);
	brav_(&cal_1.lhkl, ihkl, &ibr);
	ib[igc - 1] = ibr;
	a = cel[igc * 6 - 6];
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = a;
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
	s_wsfe(&io___502);
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___503);
	e_wsle();
	s_wsle(&io___504);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___505);
	e_wsle();
	s_wsle(&io___506);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();
	s_wsle(&io___507);
	e_wsle();
	s_wsle(&io___508);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();

/* ... Refine that cell */

	indic = 1;
	bb[2] = a;
	bb[3] = a;
	bb[4] = a;
	bb[5] = 90.f;
	bb[6] = 90.f;
	bb[7] = 90.f;
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 0.f;
	afi[6] = 0.f;
	afi[7] = 0.f;
	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = a;
	for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L1110: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
		cncalc[igc - 1] = (real) ncalc;
		fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
		ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		s_wsfe(&io___509);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___510);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___511);
		e_wsle();
		s_wsfe(&io___512);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___513);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___514);
		e_wsle();
	}
	iref = 1;
	goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
	i__3 = igc - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (ifi[i__ - 1] != ifile) {
		goto L1118;
		}
		vdelt = vgc[igc - 1] / 300.f;
		vp = vgc[igc - 1] + vdelt;
		vm = vgc[igc - 1] - vdelt;
		if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
		goto L1118;
		}
		++nsol[i__ - 1];
		if (rp[igc - 1] < rp[i__ - 1]) {
		if (isee == 1) {
			s_wsfe(&io___515);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		km[i__ - 1] = km[igc - 1];
		vgc[i__ - 1] = vgc[igc - 1];
		rp[i__ - 1] = rp[igc - 1];
		cel[i__ * 6 - 6] = cel[igc * 6 - 6];
		cel[i__ * 6 - 5] = cel[igc * 6 - 5];
		cel[i__ * 6 - 4] = cel[igc * 6 - 4];
		}
		--igc;
		goto L1119;
L1118:
		;
	}
	if (iverb == 1) {
		s_wsfe(&io___516);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___517);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
L1119:
	;
	} else {
	if (iverb == 1) {
		s_wsfe(&io___518);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___519);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	}



L1117:
	rmax = rmaxref;
	celpre[0] = celold[0];

/* ... Stop if max limit of grid tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (celpre[0] > pma[0]) {
	goto L1116;
	}
	tkill = d_mod(&ntried, &c_b1413);
	if (tkill >= 0.f) {
	killk_(&pressedk);
	if (pressedk) {
		goto L1116;
	}
	}
	goto L1102;
L1116:
	if (rmin == rmax) {
	goto L1198;
	}
	if (iverb == 1) {
	s_wsle(&io___521);
	e_wsle();
	s_wsle(&io___522);
	do_lio(&c__9, &c__1, "Best result : a=", (ftnlen)16);
	do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, "V=", (ftnlen)2);
	do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, " Rp= ", (ftnlen)5);
	do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___523);
	e_wsle();
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L1198:
	if (pressedk) {
	goto L5000;
	}


L1200:
	if (nsys[1] == 0) {
	goto L1300;
	}


/*    Hexagonal case */


	ihr = 0;
	if (ngrid == 3) {
	ihr = 1;
	}
L1290:
	if (ihr == 2) {
	nsys[1] = 2;
	}
	if (iverb == 1) {
	if (nsys[1] == 1) {
		s_wsle(&io___524);
		do_lio(&c__9, &c__1, "Hexagonal:    Rp     a      c       V     "
			"Nind", (ftnlen)46);
		e_wsle();
	}
	if (nsys[1] == 2) {
		s_wsle(&io___525);
		do_lio(&c__9, &c__1, "Rhombohedral: Rp     a      c       V     "
			"Nind", (ftnlen)46);
		e_wsle();
	}
	}
	rpsmall = 1.f;
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 2;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.f;
	ncycles = 500.f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 120.f;
	if (iverb == 1) {
	s_wsle(&io___526);
	e_wsle();
	s_wsle(&io___527);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___528);
	do_lio(&c__9, &c__1, "Grid search :", (ftnlen)13);
	e_wsle();
	if (nsys[1] == 1) {
		s_wsle(&io___529);
		do_lio(&c__9, &c__1, "   Results in hexagonal", (ftnlen)23);
		e_wsle();
	}
	if (nsys[1] == 2) {
		s_wsle(&io___530);
		do_lio(&c__9, &c__1, "   Results in rhombohedral", (ftnlen)26);
		e_wsle();
	}
	s_wsle(&io___531);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	}

	if (ngrid == 3) {
	pmin = 2.f;
	pmax = 30.f;
	pmi[0] = pmin;
	pmi[2] = pmin;
	pma[2] = dmax1 * 6.1f;
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	pma[0] = dmax1 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	pma[1] = pma[0];
	vmin = 8.f;
	vmax = pma[0] * pma[1] * pma[2];
	if (vmax > 4e3f) {
		vmax = 4e3f;
	}
	if (iverb == 1) {
		s_wsle(&io___532);
		do_lio(&c__9, &c__1, " Max(a,c), V ", (ftnlen)13);
		do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
		e_wsle();
	}
	}

	if (iverb == 1) {
	s_wsle(&io___533);
	do_lio(&c__9, &c__1, " Rp  Trial number    a      c        V  Nind I"
		"cod", (ftnlen)49);
	e_wsle();
	s_wsle(&io___534);
	e_wsle();
	}

/*     READ hkl Miller indices in hex.hkl */

	if (nsys[1] == 2) {
	goto L1260;
	}
	s_copy(tempo, "hex.hkl", (ftnlen)80, (ftnlen)7);
	goto L1261;
L1260:
	s_copy(tempo, "rho.hkl", (ftnlen)80, (ftnlen)7);
	ifile = 7;
L1261:
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___535);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	if (nsys[1] == 2) {
	goto L1262;
	}
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	goto L1263;
L1262:
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 600) {
	cal_1.nhkl0 = 600;
	}
L1263:
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1201: */
	s_rsle(&io___536);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = celpre[0];
	celpre[2] = pmi[2] - spar;
	ifin = 1;

L1202:

/*     Which parameter to vary ? a or c ? */

	ntriedb = 0.f;
	if (ifin == 1) {
	celpre[0] += spar;
	if (celpre[0] > pma[0]) {
		goto L1216;
	}
	ifin = 0;
	ntried += 1.f;
	}
	celpre[2] += spar;
	if (celpre[2] > pma[2]) {
	celpre[2] = pmi[2] - spar;
	ifin = 1;
	goto L1202;
	}
	ntried += 1.f;
	goto L1204;
L1203:
	del = deltab * (1.f - ntriedb / cy);
	i__ = 3;
	if (randi_(&iseed) > .5f) {
	i__ = 1;
	}
	celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - .5f) * 2.f;
	ntriedb += 1.f;
L1204:
	celpre[1] = celpre[0];
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1205: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.f) {
	goto L1216;
	}
	if (ntriedb != 0.f) {
	goto L1206;
	}
	if (v1 > vmax || v1 < vmin) {
	ntried += -1.f;
	goto L1202;
	}

L1206:
	calcul1_(&diff, &diff2);
	if (cal_1.nmx > cal_1.ndat10) {
	ntried += -1;
	goto L1202;
	}
	if (ntriedb != 0.f) {
	goto L1214;
	}

/* ... Rp value satisfying ??? */

	if (cal_1.lhkl >= nmax) {
	rmax = diff;
	icode = 2;
	if (diff <= rmaxref) {
		icode = 1;
	}
	} else {
	icode = 1;
	}
	celold[0] = celpre[0];
	celold[2] = celpre[2];
	if (diff > rmax) {
	goto L1217;
	}
	if (cal_1.lhkl < nmax) {
	goto L1217;
	}
L1214:
	if (diff <= rmax) {
	llhkl = cal_1.lhkl;
	rmax = diff;
	rmax2 = diff2;
	a = celpre[0];
	c__ = celpre[2];
	v2 = v1;
	if (diff < rmin) {
		rmin = diff;
		bpar[0] = a;
		bpar[2] = c__;
		v3 = v1;
	}

/* ... "Refine" that cell (by Monte Carlo too...) */

	pstartb[0] = celpre[0];
	pstartb[2] = celpre[2];
	}
	if (ntriedb <= ncycles) {
	goto L1203;
	}
	ntriedb = 0.f;
	if (rmax >= rmax0[1]) {
	goto L1217;
	}
	if (rmax2 >= .15f) {
	goto L1217;
	}
	ipen = cal_1.ndat - llhkl;
	if (ipen > cal_1.nind) {
	goto L1217;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.f;
	if (nr == 1) {
	if (igt > 50.f) {
		if (ntried / igt < 1e3f) {
		if (rmax0[1] > .1f) {
			rmax0[1] -= rmax0[1] * .05f;
			s_wsle(&io___538);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[1], (ftnlen)sizeof(
				real));
			e_wsle();
			s_wsle(&io___539);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[1], (ftnlen)sizeof(
				real));
			e_wsle();
		}
		}
	}
	}

	if (igc > 10000) {
	s_wsle(&io___540);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___541);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	--igc;
	goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = c__;
	cel[igc * 6 - 3] = 90.f;
	cel[igc * 6 - 2] = 90.f;
	cel[igc * 6 - 1] = 120.f;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = c__;
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1240: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	km[igc - 1] = llhkl;
	km2[igc - 1] = cal_1.lhkl;
	ifi[igc - 1] = ifile;
	nsol[igc - 1] = 1;
	vgc[igc - 1] = v1;
	rp[igc - 1] = rmax;
	rp2[igc - 1] = diff;
	if (rp[igc - 1] < rpsmall) {
	rpsmall = rp[igc - 1];
	isee = 1;
	} else {
	isee = 0;
	}
	supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__2);
	brav_(&cal_1.lhkl, ihkl, &ibr);
	ib[igc - 1] = ibr;
	a = cel[igc * 6 - 6];
	cel[igc * 6 - 5] = a;
	c__ = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
	s_wsfe(&io___542);
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___543);
	e_wsle();
	s_wsle(&io___544);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___545);
	e_wsle();
	s_wsle(&io___546);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();
	s_wsle(&io___547);
	e_wsle();
	s_wsle(&io___548);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();

/* ... Refine that cell */

	indic = 2;
	bb[2] = a;
	bb[3] = a;
	bb[4] = c__;
	bb[5] = 90.f;
	bb[6] = 90.f;
	bb[7] = 120.f;
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 0.f;
	afi[6] = 0.f;
	afi[7] = 0.f;
	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = c__;
	for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L1210: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
		cncalc[igc - 1] = (real) ncalc;
		fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
		ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		s_wsfe(&io___549);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___550);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___551);
		e_wsle();
		s_wsfe(&io___552);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___553);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___554);
		e_wsle();
	}
	iref = 1;
	goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
	i__3 = igc - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (ifi[i__ - 1] != ifile) {
		goto L1218;
		}
		vdelt = vgc[igc - 1] / 300.f;
		vp = vgc[igc - 1] + vdelt;
		vm = vgc[igc - 1] - vdelt;
		if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
		goto L1218;
		}
		adelt = cel[igc * 6 - 6] / 500.f;
		ap = cel[igc * 6 - 6] + adelt;
		am = cel[igc * 6 - 6] - adelt;
		if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
		goto L1218;
		}
		++nsol[i__ - 1];
		if (rp[igc - 1] < rp[i__ - 1]) {
		if (isee == 1) {
			s_wsfe(&io___555);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		km[i__ - 1] = km[igc - 1];
		vgc[i__ - 1] = vgc[igc - 1];
		rp[i__ - 1] = rp[igc - 1];
		cel[i__ * 6 - 6] = cel[igc * 6 - 6];
		cel[i__ * 6 - 5] = cel[igc * 6 - 5];
		cel[i__ * 6 - 4] = cel[igc * 6 - 4];
		}
		--igc;
		goto L1219;
L1218:
		;
	}
	if (iverb == 1) {
		s_wsfe(&io___556);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___557);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
L1219:
	;
	} else {
	if (iverb == 1) {
		s_wsfe(&io___558);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___559);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	}


L1217:
	rmax = rmaxref;
	celpre[0] = celold[0];
	celpre[2] = celold[2];

/* ... Stop if max limit of grid Carlo tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (celpre[0] > pma[0]) {
	goto L1216;
	}
	tkill = d_mod(&ntried, &c_b1413);
	if (tkill >= 0.f) {
	killk_(&pressedk);
	if (pressedk) {
		goto L1216;
	}
	}
	goto L1202;
L1216:
	if (rmin == rmax) {
	goto L1298;
	}
	if (iverb == 1) {
	s_wsle(&io___560);
	e_wsle();
	s_wsle(&io___561);
	do_lio(&c__9, &c__1, "Best result : a = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___562);
	do_lio(&c__9, &c__1, "Best result : c = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, "V = ", (ftnlen)4);
	do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___563);
	e_wsle();
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L1298:
	if (pressedk) {
	goto L5000;
	}


	++ihr;
	if (ihr == 2) {
	goto L1290;
	}
L1300:
	if (nsys[2] == 0) {
	goto L1400;
	}

/*    Tetragonal case */


	rpsmall = 1.f;
	s_wsle(&io___564);
	do_lio(&c__9, &c__1, "Tetragonal:   Rp     a       c        V     Nind", (
		ftnlen)48);
	e_wsle();
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 3;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.f;
	ncycles = 500.f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 90.f;
	if (iverb == 1) {
	s_wsle(&io___565);
	e_wsle();
	s_wsle(&io___566);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___567);
	do_lio(&c__9, &c__1, "Grid search :", (ftnlen)13);
	e_wsle();
	s_wsle(&io___568);
	do_lio(&c__9, &c__1, "   Results in tetragonal :", (ftnlen)26);
	e_wsle();
	s_wsle(&io___569);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	}

	if (ngrid == 3) {
	pmin = 2.f;
	pmax = 30.f;
	pmi[0] = pmin;
	pmi[2] = pmin;
	pma[0] = dmax1 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	pma[2] = dmax1 * 4.f;
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	vmin = 8.f;
	vmax = pma[0] * pma[1] * pma[2];
	if (vmax > 4e3f) {
		vmax = 4e3f;
	}
	if (iverb == 1) {
		s_wsle(&io___570);
		do_lio(&c__9, &c__1, " Max(a,c), V ", (ftnlen)13);
		do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
		e_wsle();
	}
	}

	if (iverb == 1) {
	s_wsle(&io___571);
	do_lio(&c__9, &c__1, " Rp  Trial number    a      c        V  Nind I"
		"cod", (ftnlen)49);
	e_wsle();
	s_wsle(&io___572);
	e_wsle();
	}

/*     READ hkl Miller indices in tet.hkl */

	s_copy(tempo, "tet.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___573);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1301: */
	s_rsle(&io___574);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = celpre[0];
	celpre[2] = pmi[2] - spar;
	ifin = 1;

L1302:

/*     Which parameter to vary ? a or c ? */

	ntriedb = 0.f;
	if (ifin == 1) {
	celpre[0] += spar;
	if (celpre[0] > pma[0]) {
		goto L1316;
	}
	ifin = 0;
	ntried += 1.f;
	}
	celpre[2] += spar;
	if (celpre[2] > pma[2]) {
	celpre[2] = pmi[2] - spar;
	ifin = 1;
	goto L1302;
	}
	ntried += 1.f;
	goto L1304;
L1303:
	del = deltab * (1.f - ntriedb / cy);
	i__ = 3;
	if (randi_(&iseed) > .5f) {
	i__ = 1;
	}
	celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - .5f) * 2.f;
	ntriedb += 1.f;
L1304:
	celpre[1] = celpre[0];
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1305: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.f) {
	goto L1316;
	}
	if (ntriedb != 0.f) {
	goto L1306;
	}
	if (v1 > vmax || v1 < vmin) {
	ntried += -1;
	goto L1302;
	}

L1306:
	calcul1_(&diff, &diff2);
	if (cal_1.nmx > cal_1.ndat10) {
	ntried += -1;
	goto L1302;
	}
	if (ntriedb != 0.f) {
	goto L1314;
	}

/* ... Rp value satisfying ??? */

	if (cal_1.lhkl >= nmax) {
	rmax = diff;
	icode = 2;
	if (diff <= rmaxref) {
		icode = 1;
	}
	} else {
	icode = 1;
	}
	celold[0] = celpre[0];
	celold[2] = celpre[2];
	if (diff > rmax) {
	goto L1317;
	}
	if (cal_1.lhkl < nmax) {
	goto L1317;
	}
L1314:
	if (diff <= rmax) {
	llhkl = cal_1.lhkl;
	rmax = diff;
	rmax2 = diff2;
	a = celpre[0];
	c__ = celpre[2];
	v2 = v1;
	if (diff < rmin) {
		rmin = diff;
		bpar[0] = a;
		bpar[2] = c__;
		v3 = v1;
	}

/* ... "Refine" that cell (by Monte Carlo too...) */

	pstartb[0] = celpre[0];
	pstartb[2] = celpre[2];
	}
	if (ntriedb <= ncycles) {
	goto L1303;
	}
	ntriedb = 0.f;
	if (rmax >= rmax0[2]) {
	goto L1317;
	}
	if (rmax2 >= .15f) {
	goto L1317;
	}
	ipen = cal_1.ndat - llhkl;
	if (ipen > cal_1.nind) {
	goto L1317;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.f;
	if (nr == 1) {
	if (igt > 50.f) {
		if (ntried / igt < 1e3f) {
		if (rmax0[2] > .1f) {
			rmax0[2] -= rmax0[2] * .05f;
			s_wsle(&io___575);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[2], (ftnlen)sizeof(
				real));
			e_wsle();
			s_wsle(&io___576);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[2], (ftnlen)sizeof(
				real));
			e_wsle();
		}
		}
	}
	}

	if (igc > 10000) {
	s_wsle(&io___577);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___578);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	--igc;
	goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = c__;
	cel[igc * 6 - 3] = 90.f;
	cel[igc * 6 - 2] = 90.f;
	cel[igc * 6 - 1] = 90.f;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = c__;
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1340: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	km[igc - 1] = llhkl;
	km2[igc - 1] = cal_1.lhkl;
	ifi[igc - 1] = ifile;
	nsol[igc - 1] = 1;
	vgc[igc - 1] = v1;
	rp[igc - 1] = rmax;
	rp2[igc - 1] = diff;
	if (rp[igc - 1] < rpsmall) {
	rpsmall = rp[igc - 1];
	isee = 1;
	} else {
	isee = 0;
	}
	supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__4);
	brav_(&cal_1.lhkl, ihkl, &ibr);
	ib[igc - 1] = ibr;
	a = cel[igc * 6 - 6];
	cel[igc * 6 - 5] = a;
	c__ = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
	s_wsfe(&io___579);
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___580);
	e_wsle();
	s_wsle(&io___581);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___582);
	e_wsle();
	s_wsle(&io___583);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();
	s_wsle(&io___584);
	e_wsle();
	s_wsle(&io___585);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();

/* ... Refine that cell */

	indic = 2;
	bb[2] = a;
	bb[3] = a;
	bb[4] = c__;
	bb[5] = 90.f;
	bb[6] = 90.f;
	bb[7] = 90.f;
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 0.f;
	afi[6] = 0.f;
	afi[7] = 0.f;
	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = c__;
	for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L1310: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
		cncalc[igc - 1] = (real) ncalc;
		fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
		ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		s_wsfe(&io___586);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___587);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___588);
		e_wsle();
		s_wsfe(&io___589);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___590);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___591);
		e_wsle();
	}
	iref = 1;
	goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
	i__3 = igc - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (ifi[i__ - 1] != ifile) {
		goto L1318;
		}
		vdelt = vgc[igc - 1] / 300.f;
		vp = vgc[igc - 1] + vdelt;
		vm = vgc[igc - 1] - vdelt;
		if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
		goto L1318;
		}
		adelt = cel[igc * 6 - 6] / 500.f;
		ap = cel[igc * 6 - 6] + adelt;
		am = cel[igc * 6 - 6] - adelt;
		if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
		goto L1318;
		}
		++nsol[i__ - 1];
		if (rp[igc - 1] < rp[i__ - 1]) {
		if (isee == 1) {
			s_wsfe(&io___592);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		km[i__ - 1] = km[igc - 1];
		vgc[i__ - 1] = vgc[igc - 1];
		rp[i__ - 1] = rp[igc - 1];
		cel[i__ * 6 - 6] = cel[igc * 6 - 6];
		cel[i__ * 6 - 5] = cel[igc * 6 - 5];
		cel[i__ * 6 - 4] = cel[igc * 6 - 4];
		}
		--igc;
		goto L1319;
L1318:
		;
	}
	if (iverb == 1) {
		s_wsfe(&io___593);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___594);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
L1319:
	;
	} else {
	if (iverb == 1) {
		s_wsfe(&io___595);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___596);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	}


L1317:
	rmax = rmaxref;
	celpre[0] = celold[0];
	celpre[2] = celold[2];

/* ... Stop if max limit of grid tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (celpre[0] > pma[0]) {
	goto L1316;
	}
	tkill = d_mod(&ntried, &c_b1413);
	if (tkill >= 0.f) {
	killk_(&pressedk);
	if (pressedk) {
		goto L1316;
	}
	}
	goto L1302;
L1316:
	if (rmin == rmax) {
	goto L1398;
	}
	if (iverb == 1) {
	s_wsle(&io___597);
	e_wsle();
	s_wsle(&io___598);
	do_lio(&c__9, &c__1, "Best result : a = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___599);
	do_lio(&c__9, &c__1, "Best result : c = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, "V = ", (ftnlen)4);
	do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___600);
	e_wsle();
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L1398:
	if (pressedk) {
	goto L5000;
	}


L1400:
	if (nsys[3] == 0) {
	goto L1500;
	}

/*    Orthorhombic case */


	rpsmall = 1.f;
	s_wsle(&io___601);
	do_lio(&c__9, &c__1, "Orthorhombic: Rp     a       b       c        V   "
		"  Nind", (ftnlen)56);
	e_wsle();
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 4;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.f;
	ncycles = 1e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[4] = 90.f;
	celpre[5] = 90.f;
	if (iverb == 1) {
	s_wsle(&io___602);
	e_wsle();
	s_wsle(&io___603);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___604);
	do_lio(&c__9, &c__1, "Grid search :", (ftnlen)13);
	e_wsle();
	s_wsle(&io___605);
	do_lio(&c__9, &c__1, "   Results in orthorhombic :", (ftnlen)28);
	e_wsle();
	s_wsle(&io___606);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	}

	if (ngrid == 3) {
	pmin = 2.f;
	pmax = 20.f;
	pmi[0] = pmin;
	pma[0] = dmax1 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	pmi[1] = pmin;
	pma[1] = dmax2 * 2.1f;
	if (pma[1] > pmax) {
		pma[1] = pmax;
	}
	pmi[2] = pmin;
	pma[2] = dmax3 * 2.1f;
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	vmin = 8.f;
	vmax = pma[0] * pma[1] * pma[2];
	if (vmax > 2e3f) {
		vmax = 2e3f;
	}
	if (iverb == 1) {
		s_wsle(&io___607);
		do_lio(&c__9, &c__1, " Max (a,b,c), V ", (ftnlen)16);
		do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&pma[1], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
		e_wsle();
	}
	}

	if (iverb == 1) {
	s_wsle(&io___608);
	do_lio(&c__9, &c__1, " Rp  Trial number    a   b   c        V  Nind "
		"Icod", (ftnlen)50);
	e_wsle();
	s_wsle(&io___609);
	e_wsle();
	}

/*     READ hkl Miller indices in ort.hkl */

	s_copy(tempo, "ort.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___610);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1401: */
	s_rsle(&io___611);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = pmi[1] - spar;
	celpre[2] = pmi[2] - spar;
	ifin = 1;
	ifin2 = 1;

L1402:

/*     Which parameter to vary ? a or b or c ? */

	ntriedb = 0.f;
	if (ifin == 1) {
	celpre[0] += spar;
	s_wsle(&io___613);
	do_lio(&c__9, &c__1, "  a = ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&celpre[0], (ftnlen)sizeof(real));
	e_wsle();
	if (celpre[0] > pma[0]) {
		goto L1416;
	}
	ifin = 0;
	ntried += 1.f;
	}
	if (ifin2 == 1) {
	celpre[1] += spar;
	ifin2 = 0;
	}
	if (celpre[1] > pma[1]) {
	celpre[1] = pmi[1] - spar;
	ifin = 1;
	goto L1402;
	}
	ntried += 1.f;
	celpre[2] += spar;
	if (celpre[2] > pma[2]) {
	celpre[2] = pmi[2] - spar;
	ifin2 = 1;
	goto L1402;
	}
	ntried += 1.f;
	goto L1404;
L1403:
	del = deltab * (1.f - ntriedb / cy);
	x = randi_(&iseed);
	if (x >= 0.f && x < .33333f) {
	i__ = 1;
	}
	if (x >= .33333f && x < .66666f) {
	i__ = 2;
	}
	if (x >= .66666f && x <= 1.f) {
	i__ = 3;
	}
	celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - .5f) * 2.f;
	ntriedb += 1.f;
L1404:
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1405: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.f) {
	goto L1416;
	}
	if (ntriedb != 0.f) {
	goto L1406;
	}
	if (v1 > vmax || v1 < vmin) {
	ntried += -1.f;
	goto L1402;
	}

L1406:
	calcul1_(&diff, &diff2);
	if (cal_1.nmx > cal_1.ndat10) {
	ntried += -1;
	goto L1402;
	}
	if (ntriedb != 0.f) {
	goto L1414;
	}

/* ... Rp value satisfying ??? */

	if (cal_1.lhkl >= nmax) {
	rmax = diff;
	icode = 2;
	if (diff <= rmaxref) {
		icode = 1;
	}
	} else {
	icode = 1;
	}
	celold[0] = celpre[0];
	celold[1] = celpre[1];
	celold[2] = celpre[2];
	if (diff > rmax) {
	goto L1417;
	}
	if (cal_1.lhkl < nmax) {
	goto L1417;
	}
L1414:
	if (diff <= rmax) {
	llhkl = cal_1.lhkl;
	rmax = diff;
	rmax2 = diff2;
	a = celpre[0];
	b = celpre[1];
	c__ = celpre[2];
	v2 = v1;
	if (diff < rmin) {
		rmin = diff;
		bpar[0] = a;
		bpar[1] = b;
		bpar[2] = c__;
		v3 = v1;
	}

/* ... "Refine" that cell (by Monte Carlo too...) */

	pstartb[0] = celpre[0];
	pstartb[1] = celpre[1];
	pstartb[2] = celpre[2];
	}
	if (ntriedb <= ncycles) {
	goto L1403;
	}
	ntriedb = 0.f;
	if (rmax >= rmax0[3]) {
	goto L1417;
	}
	if (rmax2 >= .15f) {
	goto L1417;
	}
	ipen = cal_1.ndat - llhkl;
	if (ipen > cal_1.nind) {
	goto L1417;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.f;
	if (nr == 1) {
	if (igt > 50.f) {
		if (ntried / igt < 1e4f) {
		if (rmax0[3] > .1f) {
			rmax0[3] -= rmax0[3] * .05f;
			s_wsle(&io___614);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[3], (ftnlen)sizeof(
				real));
			e_wsle();
			s_wsle(&io___615);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[3], (ftnlen)sizeof(
				real));
			e_wsle();
		}
		}
	}
	}

	if (igc > 10000) {
	s_wsle(&io___616);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___617);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	--igc;
	goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = b;
	cel[igc * 6 - 4] = c__;
	cel[igc * 6 - 3] = 90.f;
	cel[igc * 6 - 2] = 90.f;
	cel[igc * 6 - 1] = 90.f;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = b;
	celpre[2] = c__;
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1440: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	km[igc - 1] = llhkl;
	km2[igc - 1] = cal_1.lhkl;
	ifi[igc - 1] = ifile;
	nsol[igc - 1] = 1;
	vgc[igc - 1] = v1;
	rp[igc - 1] = rmax;
	rp2[igc - 1] = diff;
	if (rp[igc - 1] < rpsmall) {
	rpsmall = rp[igc - 1];
	isee = 1;
	} else {
	isee = 0;
	}
	supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__1);
	brav_(&cal_1.lhkl, ihkl, &ibr);
	ib[igc - 1] = ibr;
	a = cel[igc * 6 - 6];
	b = cel[igc * 6 - 5];
	c__ = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
	s_wsfe(&io___618);
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___619);
	e_wsle();
	s_wsle(&io___620);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___621);
	e_wsle();
	s_wsle(&io___622);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();
	s_wsle(&io___623);
	e_wsle();
	s_wsle(&io___624);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();

/* ... Refine that cell */

	indic = 0;
	bb[2] = a;
	bb[3] = b;
	bb[4] = c__;
	bb[5] = 90.f;
	bb[6] = 90.f;
	bb[7] = 90.f;
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 0.f;
	afi[6] = 0.f;
	afi[7] = 0.f;
	celpre[0] = a;
	celpre[1] = b;
	celpre[2] = c__;
	for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L1410: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
		cncalc[igc - 1] = (real) ncalc;
		fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
		ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		s_wsfe(&io___625);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___626);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___627);
		e_wsle();
		s_wsfe(&io___628);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___629);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___630);
		e_wsle();
	}
	iref = 1;
	goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
	i__3 = igc - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (ifi[i__ - 1] != ifile) {
		goto L1418;
		}
		vdelt = vgc[igc - 1] / 300.f;
		vp = vgc[igc - 1] + vdelt;
		vm = vgc[igc - 1] - vdelt;
		if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
		goto L1418;
		}
		adelt = cel[igc * 6 - 6] / 500.f;
		ap = cel[igc * 6 - 6] + adelt;
		am = cel[igc * 6 - 6] - adelt;
		na = 0;
		if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
		na = 1;
		}
		bdelt = cel[igc * 6 - 5] / 500.f;
		bp = cel[igc * 6 - 5] + bdelt;
		bm = cel[igc * 6 - 5] - bdelt;
		nb = 0;
		if (cel[i__ * 6 - 6] > bp || cel[i__ * 6 - 6] < bm) {
		nb = 1;
		}
		cdelt = cel[igc * 6 - 4] / 500.f;
		cp = cel[igc * 6 - 4] + cdelt;
		cm = cel[igc * 6 - 4] - cdelt;
		nc = 0;
		if (cel[i__ * 6 - 6] > cp || cel[i__ * 6 - 6] < cm) {
		nc = 1;
		}
		if (na == 1 && nb == 1 && nc == 1) {
		goto L1418;
		}
		++nsol[i__ - 1];
		if (rp[igc - 1] < rp[i__ - 1]) {
		if (isee == 1) {
			s_wsfe(&io___631);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		km[i__ - 1] = km[igc - 1];
		vgc[i__ - 1] = vgc[igc - 1];
		rp[i__ - 1] = rp[igc - 1];
		cel[i__ * 6 - 6] = cel[igc * 6 - 6];
		cel[i__ * 6 - 5] = cel[igc * 6 - 5];
		cel[i__ * 6 - 4] = cel[igc * 6 - 4];
		}
		--igc;
		goto L1419;
L1418:
		;
	}
	if (iverb == 1) {
		s_wsfe(&io___632);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___633);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
L1419:
	;
	} else {
	if (iverb == 1) {
		s_wsfe(&io___634);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___635);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	}


L1417:
	rmax = rmaxref;
	celpre[0] = celold[0];
	celpre[1] = celold[1];
	celpre[2] = celold[2];

/* ... Stop if max limit of grid tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (celpre[0] > pma[0]) {
	goto L1416;
	}
	tkill = d_mod(&ntried, &c_b1413);
	if (tkill >= 0.f) {
	killk_(&pressedk);
	if (pressedk) {
		goto L1416;
	}
	}
	goto L1402;
L1416:
	if (rmin == rmax) {
	goto L1498;
	}
	if (iverb == 1) {
	s_wsle(&io___636);
	e_wsle();
	s_wsle(&io___637);
	do_lio(&c__9, &c__1, "Best result : a = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, " Rp = ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___638);
	do_lio(&c__9, &c__1, "Best result : b = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[1], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___639);
	do_lio(&c__9, &c__1, "Best result : c = ", (ftnlen)18);
	do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, "V = ", (ftnlen)4);
	do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___640);
	e_wsle();
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L1498:
	if (pressedk) {
	goto L5000;
	}


L1500:
	if (nsys[4] == 0) {
	goto L1600;
	}

/*    Monoclinic case - would be too long in grid search, but... */


	rpsmall = 1.f;
	s_wsle(&io___641);
	do_lio(&c__9, &c__1, "Monoclinic:   Rp     a       b       c       bet  "
		"   V     Nind", (ftnlen)63);
	e_wsle();
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 5;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.f;
	ncycles = 2e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.f;
	celpre[5] = 90.f;
	if (iverb == 1) {
	s_wsle(&io___642);
	e_wsle();
	s_wsle(&io___643);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___644);
	do_lio(&c__9, &c__1, "Grid search :", (ftnlen)13);
	e_wsle();
	s_wsle(&io___645);
	do_lio(&c__9, &c__1, "   Results in monoclinic :", (ftnlen)26);
	e_wsle();
	s_wsle(&io___646);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	}

	if (ngrid == 3) {
	pmin = 2.f;
	pmax = 20.f;
	pmi[0] = pmin;
	pma[0] = dmax1 * 2.1f;
	if (pma[0] > pmax) {
		pma[0] = pmax;
	}
	pmi[1] = pmin;
	pma[1] = dmax1 * 2.1f;
	if (pma[1] > pmax) {
		pma[1] = pmax;
	}
	pmi[2] = pmin;
	pma[2] = dmax2 * 2.1f;
	if (pma[2] > pmax) {
		pma[2] = pmax;
	}
	vmin = 8.f;
	vmax = pma[0] * pma[1] * pma[2];
	if (vmax > 2e3f) {
		vmax = 2e3f;
	}
	if (iverb == 1) {
		s_wsle(&io___647);
		do_lio(&c__9, &c__1, " Max (a,b,c), V ", (ftnlen)16);
		do_lio(&c__4, &c__1, (char *)&pma[0], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&pma[1], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&pma[2], (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
		e_wsle();
	}
	}

	if (iverb == 1) {
	s_wsle(&io___648);
	do_lio(&c__9, &c__1, " Rp  Trial number    a   b   c     bet  V  Nin"
		"d Icod", (ftnlen)52);
	e_wsle();
	s_wsle(&io___649);
	e_wsle();
	}

/*     READ hkl Miller indices in mon.hkl */

	s_copy(tempo, "mon.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___650);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1501: */
	s_rsle(&io___651);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = pmi[1] - spar;
	celpre[2] = pmi[2] - spar;
	celpre[4] = pmi[4] - sang;
	ifin1 = 1;
	ifin2 = 1;
	ifin3 = 1;

L1502:

/*     Which parameter to vary ? a or b or c or bet ? */

	ntriedb = 0.f;
	if (ifin1 == 1) {
	celpre[0] += spar;
	s_wsle(&io___654);
	do_lio(&c__9, &c__1, "  a = ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&celpre[0], (ftnlen)sizeof(real));
	e_wsle();
	if (celpre[0] > pma[0]) {
		goto L1516;
	}
	ifin1 = 0;
	ntried += 1.f;
	}
	if (ifin2 == 1) {
	celpre[1] += spar;
	ifin2 = 0;
	}
	if (celpre[1] > pma[1]) {
	celpre[1] = pmi[1] - spar;
	ifin1 = 1;
	goto L1502;
	}
	ntried += 1.f;
	if (ifin3 == 1) {
	celpre[2] += spar;
	ifin3 = 0;
	}
	if (celpre[2] > pma[2]) {
	celpre[2] = pmi[2] - spar;
	ifin2 = 1;
	goto L1502;
	}
	ntried += 1.f;
	celpre[4] += sang;
	if (celpre[4] > pma[4]) {
	celpre[4] = pmi[4] - sang;
	ifin3 = 1;
	goto L1502;
	}
	ntried += 1.f;
	goto L1504;
L1503:
	del = deltab * (1.f - ntriedb / cy);
	deld = deltad * (1.f - ntriedb / cy);
	x = randi_(&iseed);
	if (x >= 0.f && x < .25f) {
	i__ = 1;
	}
	if (x >= .25f && x < .5f) {
	i__ = 2;
	}
	if (x >= .5f && x < .75f) {
	i__ = 3;
	}
	if (x >= .75f && x <= 1.f) {
	i__ = 5;
	}
	if (i__ != 5) {
	celpre[i__ - 1] = pstartb[i__ - 1] + del * (randi_(&iseed) - .5f) * 
		2.f;
	} else {
	celpre[i__ - 1] = pstartb[i__ - 1] + deld * (randi_(&iseed) - .5f) * 
		2.f;
	}
	ntriedb += 1.f;
L1504:
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1505: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.f) {
	goto L1516;
	}
	if (ntriedb != 0.f) {
	goto L1506;
	}
	if (v1 > vmax || v1 < vmin) {
	ntried += -1.f;
	goto L1502;
	}

L1506:
	calcul1_(&diff, &diff2);
	if (cal_1.nmx > cal_1.ndat10) {
	ntried += -1;
	goto L1502;
	}
	if (ntriedb != 0.f) {
	goto L1514;
	}

/* ... Rp value satisfying ??? */

	if (cal_1.lhkl >= nmax) {
	rmax = diff;
	icode = 2;
	if (diff <= rmaxref) {
		icode = 1;
	}
	} else {
	icode = 1;
	}
	celold[0] = celpre[0];
	celold[1] = celpre[1];
	celold[2] = celpre[2];
	celold[4] = celpre[4];
	if (diff > rmax) {
	goto L1517;
	}
	if (cal_1.lhkl < nmax) {
	goto L1517;
	}
L1514:
	if (diff <= rmax) {
	llhkl = cal_1.lhkl;
	rmax = diff;
	rmax2 = diff2;
	a = celpre[0];
	b = celpre[1];
	c__ = celpre[2];
	bet = celpre[4];
	v2 = v1;
	if (diff < rmin) {
		rmin = diff;
		bpar[0] = a;
		bpar[1] = b;
		bpar[2] = c__;
		bpar[4] = bet;
		v3 = v1;
	}

/* ... "Refine" that cell (by Monte Carlo too...) */

	pstartb[0] = celpre[0];
	pstartb[1] = celpre[1];
	pstartb[2] = celpre[2];
	pstartb[4] = celpre[4];
	}
	if (ntriedb <= ncycles) {
	goto L1503;
	}
	ntriedb = 0.f;
	if (rmax >= rmax0[4]) {
	goto L1517;
	}
	if (rmax2 >= .15f) {
	goto L1517;
	}
	ipen = cal_1.ndat - llhkl;
	if (ipen > cal_1.nind) {
	goto L1517;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.f;
	if (nr == 1) {
	if (igt > 50.f) {
		if (ntried / igt < 1e5f) {
		if (rmax0[4] > .1f) {
			rmax0[4] -= rmax0[4] * .05f;
			s_wsle(&io___655);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[4], (ftnlen)sizeof(
				real));
			e_wsle();
			s_wsle(&io___656);
			do_lio(&c__9, &c__1, "  Rmax reduced by 5%, now Rmax = ", 
				(ftnlen)33);
			do_lio(&c__4, &c__1, (char *)&rmax0[4], (ftnlen)sizeof(
				real));
			e_wsle();
		}
		}
	}
	}

	if (igc > 10000) {
	s_wsle(&io___657);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	s_wsle(&io___658);
	do_lio(&c__9, &c__1, "   More than 10000 good cells = STOP", (ftnlen)
		36);
	e_wsle();
	--igc;
	goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = b;
	cel[igc * 6 - 4] = c__;
	cel[igc * 6 - 3] = 90.f;
	cel[igc * 6 - 2] = bet;
	cel[igc * 6 - 1] = 90.f;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = b;
	celpre[2] = c__;
	celpre[4] = bet;
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L1540: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	km[igc - 1] = llhkl;
	km2[igc - 1] = cal_1.lhkl;
	ifi[igc - 1] = ifile;
	nsol[igc - 1] = 1;
	vgc[igc - 1] = v1;
	rp[igc - 1] = rmax;
	rp2[igc - 1] = diff;
	if (rp[igc - 1] < rpsmall) {
	rpsmall = rp[igc - 1];
	isee = 1;
	} else {
	isee = 0;
	}
	supcel_(&cal_1.lhkl, ihkl, cel, &igc, vgc, &c__1);
	a = cel[igc * 6 - 6];
	b = cel[igc * 6 - 5];
	c__ = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
	s_wsfe(&io___659);
	do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___660);
	e_wsle();
	s_wsle(&io___661);
	do_lio(&c__9, &c__1, "=============================================="
		"======\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___662);
	e_wsle();
	s_wsle(&io___663);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();
	s_wsle(&io___664);
	e_wsle();
	s_wsle(&io___665);
	do_lio(&c__9, &c__1, " YOU HAVE FOUND AN INTERESTING RESULT : Rp < R"
		"min!", (ftnlen)50);
	e_wsle();

/* ... Refine that cell */

	indic = 0;
	bb[2] = a;
	bb[3] = b;
	bb[4] = c__;
	bb[5] = 90.f;
	bb[6] = bet;
	bb[7] = 90.f;
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 0.f;
	afi[6] = 1.f;
	afi[7] = 0.f;
	celpre[0] = a;
	celpre[1] = b;
	celpre[2] = c__;
	celpre[4] = bet;
	for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
/* L1510: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
		}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &igc);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
		cncalc[igc - 1] = (real) ncalc;
		fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
		ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
		s_wsfe(&io___666);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___667);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___668);
		e_wsle();
		s_wsfe(&io___669);
		do_fio(&c__1, (char *)&fm20[igc - 1], (ftnlen)sizeof(real));
		e_wsfe();
		s_wsfe(&io___670);
		do_fio(&c__1, (char *)&ff20[igc - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
		e_wsfe();
		s_wsle(&io___671);
		e_wsle();
	}
	iref = 1;
	goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
	i__3 = igc - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (ifi[i__ - 1] != ifile) {
		goto L1518;
		}
		vdelt = vgc[igc - 1] / 300.f;
		vp = vgc[igc - 1] + vdelt;
		vm = vgc[igc - 1] - vdelt;
		if (vgc[i__ - 1] > vp || vgc[i__ - 1] < vm) {
		goto L1518;
		}
		bdelt = cel[igc * 6 - 5] / 500.f;
		bp = cel[igc * 6 - 5] + bdelt;
		bm = cel[igc * 6 - 5] - bdelt;
		if (cel[i__ * 6 - 5] > bp || cel[i__ * 6 - 5] < bm) {
		goto L1518;
		}
		betdelt = cel[igc * 6 - 2] / 500.f;
		betp = cel[igc * 6 - 2] + betdelt;
		betm = cel[igc * 6 - 2] - betdelt;
		if (cel[i__ * 6 - 2] > betp || cel[i__ * 6 - 2] < betm) {
		goto L1518;
		}
		adelt = cel[igc * 6 - 6] / 500.f;
		ap = cel[igc * 6 - 6] + adelt;
		am = cel[igc * 6 - 6] - adelt;
		na = 0;
		if (cel[i__ * 6 - 6] > ap || cel[i__ * 6 - 6] < am) {
		na = 1;
		}
		cdelt = cel[igc * 6 - 4] / 500.f;
		cp = cel[igc * 6 - 4] + cdelt;
		cm = cel[igc * 6 - 4] - cdelt;
		nc = 0;
		if (cel[i__ * 6 - 6] > cp || cel[i__ * 6 - 6] < cm) {
		nc = 1;
		}
		if (na == 1 && nc == 1) {
		goto L1518;
		}
		++nsol[i__ - 1];
		if (rp[igc - 1] < rp[i__ - 1]) {
		if (isee == 1) {
			s_wsfe(&io___672);
			do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
			do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
			e_wsfe();
		}
		km[i__ - 1] = km[igc - 1];
		vgc[i__ - 1] = vgc[igc - 1];
		rp[i__ - 1] = rp[igc - 1];
		cel[i__ * 6 - 6] = cel[igc * 6 - 6];
		cel[i__ * 6 - 5] = cel[igc * 6 - 5];
		cel[i__ * 6 - 4] = cel[igc * 6 - 4];
		cel[i__ * 6 - 2] = cel[igc * 6 - 2];
		}
		--igc;
		goto L1519;
L1518:
		;
	}
	if (iverb == 1) {
		s_wsfe(&io___673);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___674);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
L1519:
	;
	} else {
	if (iverb == 1) {
		s_wsfe(&io___675);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ntried, (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&icode, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	if (isee == 1) {
		s_wsfe(&io___676);
		do_fio(&c__1, (char *)&rmax, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&a, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&b, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&c__, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&bet, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&v2, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ipen, (ftnlen)sizeof(integer));
		e_wsfe();
	}
	}


L1517:
	rmax = rmaxref;
	celpre[0] = celold[0];
	celpre[1] = celold[1];
	celpre[2] = celold[2];
	celpre[4] = celold[4];

/* ... Stop if max limit of grid tests outpassed */
/*         or if K is pressed (tested every 30000 MC event) */

	if (celpre[0] > pma[0]) {
	goto L1516;
	}
	tkill = d_mod(&ntried, &c_b1413);
	if (tkill >= 0.f) {
	killk_(&pressedk);
	if (pressedk) {
		goto L1516;
	}
	}
	goto L1502;
L1516:
	if (rmin == rmax) {
	goto L1598;
	}
	if (iverb == 1) {
	s_wsle(&io___677);
	e_wsle();
	s_wsle(&io___678);
	do_lio(&c__9, &c__1, "Best result : a =   ", (ftnlen)20);
	do_lio(&c__4, &c__1, (char *)&bpar[0], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, "  Rp = ", (ftnlen)7);
	do_lio(&c__4, &c__1, (char *)&rmin, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___679);
	do_lio(&c__9, &c__1, "Best result : b =   ", (ftnlen)20);
	do_lio(&c__4, &c__1, (char *)&bpar[1], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___680);
	do_lio(&c__9, &c__1, "Best result : c =   ", (ftnlen)20);
	do_lio(&c__4, &c__1, (char *)&bpar[2], (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___681);
	do_lio(&c__9, &c__1, "Best result : bet = ", (ftnlen)20);
	do_lio(&c__4, &c__1, (char *)&bpar[4], (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, "  V = ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&v3, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___682);
	e_wsle();
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	}
	rmin = rmax;
L1598:
	if (pressedk) {
	goto L5000;
	}


L1600:
	if (nsys[5] == 0) {
	goto L5000;
	}

/*    Triclinic case - would be too long in grid search... */


L5000:
	if (igc == 0) {
	goto L6000;
	}
	imem = 0;

/*   Prepare sorted output for CHEKCELL */

	s_wsle(&io___684);
	e_wsle();
	s_wsle(&io___685);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___686);
	do_lio(&c__9, &c__1, "  FINAL LIST OF CELL PROPOSALS, sorted by McM20 :", 
		(ftnlen)49);
	e_wsle();
	s_wsle(&io___687);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___688);
	e_wsle();
	s_wsle(&io___689);
	do_lio(&c__9, &c__1, "      Global list as produced in the .ckm file", (
		ftnlen)46);
	e_wsle();
	s_wsle(&io___690);
	do_lio(&c__9, &c__1, "            (IN=number of indexed lines)", (ftnlen)
		40);
	e_wsle();
	s_wsle(&io___691);
	do_lio(&c__9, &c__1, "  The correct cell has some chances to be just bel"
		"ow", (ftnlen)52);
	e_wsle();
	s_wsle(&io___692);
	e_wsle();
	s_wsfe(&io___693);
	e_wsfe();
	s_wsle(&io___694);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
	goto L1998;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L1998:
	open_write1__(&c__24, tempo, (ftnlen)80);
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".mcm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
	goto L2066;
	}
	filedel_(&c__25, tempo, (ftnlen)80);
L2066:
	open_write1__(&c__25, tempo, (ftnlen)80);

/*  Calculate new F.o.M */

	i__3 = igc;
	for (j = 1; j <= i__3; ++j) {
	if (rp[j - 1] < .001f) {
		rp[j - 1] = .001f;
	}
	x2 = 1.f;
	if (ib[j - 1] == 1) {
		x2 = 4.f;
	}
	if (ib[j - 1] == 2) {
		x2 = 2.f;
	}
	if (ib[j - 1] == 3) {
		x2 = 2.f;
	}
	if (ib[j - 1] == 4) {
		x2 = 2.f;
	}
	if (ib[j - 1] == 5) {
		x2 = 6.f;
	}
	if (ifi[j - 1] == 7) {
		x2 = 6.f;
	}
	x3 = 1.f;
	if (ifi[j - 1] == 1) {
		x3 = 6.f;
	}
	if (ifi[j - 1] == 2) {
		x3 = 4.f;
	}
	if (ifi[j - 1] == 3) {
		x3 = 4.f;
	}
	if (ifi[j - 1] == 4) {
		x3 = 2.f;
	}
	if (ifi[j - 1] == 7) {
		x3 = 6.f;
	}
	xfom[j - 1] = 1.f / rp[j - 1] * 100.f / cncalc[j - 1] * x2 * x3;
/* L18000: */
	}

	sort_(&igc, xfom, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
	s_copy(indxprog, "McMaille4.00", (ftnlen)12, (ftnlen)12);
	ibest = ll[igc - 1];
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 20) {
		goto L20000;
	}
	j = ll[igc + 1 - i__ - 1];
	ipedig = j;
	for (k = 1; k <= 6; ++k) {
/* L2067: */
		celpre[k - 1] = cel[k + j * 6 - 7];
	}
	dcell_(celpre, cal_1.al, &v1);
	ql[0] = cal_1.al[0] * 1e4f;
	ql[1] = cal_1.al[4] * 1e4f;
	ql[2] = cal_1.al[8] * 1e4f;
	ql[3] = cal_1.al[7] * 1e4f;
	ql[4] = cal_1.al[6] * 1e4f;
	ql[5] = cal_1.al[3] * 1e4f;
	if (ib[j - 1] == 1) {
		*(unsigned char *)bl = 'I';
	}
	if (ib[j - 1] == 2) {
		*(unsigned char *)bl = 'A';
	}
	if (ib[j - 1] == 3) {
		*(unsigned char *)bl = 'B';
	}
	if (ib[j - 1] == 4) {
		*(unsigned char *)bl = 'C';
	}
	if (ib[j - 1] == 5) {
		*(unsigned char *)bl = 'F';
	}
	if (ib[j - 1] == 6) {
		*(unsigned char *)bl = 'P';
	}
	if (ifi[j - 1] == 7) {
		*(unsigned char *)bl = 'R';
	}
	if (ifi[j - 1] == 1) {
		s_copy(more, "Cubic *****", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 2) {
		s_copy(more, "Hexag **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 3) {
		s_copy(more, "Tetra **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 4) {
		s_copy(more, "Ortho ***  ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 5) {
		s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 6) {
		s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 7) {
		s_copy(more, "Rhomb **** ", (ftnlen)11, (ftnlen)11);
	}
	vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];
	s_wsfe(&io___708);
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xfom[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___709);
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xfom[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	do_fio(&c__1, bl, (ftnlen)1);
	do_fio(&c__1, more, (ftnlen)11);
	e_wsfe();
	s_wsfe(&io___710);
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xfom[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	do_fio(&c__1, bl, (ftnlen)1);
	do_fio(&c__1, indxprog, (ftnlen)12);
	do_fio(&c__1, datenow, (ftnlen)7);
	do_fio(&c__1, timenow, (ftnlen)8);
	do_fio(&c__1, (char *)&ipedig, (ftnlen)sizeof(integer));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&ql[k - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L2000: */
	}
L20000:
/* L2001: */
	s_wsle(&io___711);
	e_wsle();
	s_wsle(&io___712);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___713);
	do_lio(&c__9, &c__1, "    FINAL LIST OF CELL PROPOSALS, sorted by F(20) :"
		, (ftnlen)51);
	e_wsle();
	s_wsle(&io___714);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___715);
	e_wsle();
	s_wsle(&io___716);
	do_lio(&c__9, &c__1, "      Global list (IN=number of indexed lines)", (
		ftnlen)46);
	e_wsle();
	s_wsle(&io___717);
	do_lio(&c__9, &c__1, "  The correct cell has some chances to be just bel"
		"ow", (ftnlen)52);
	e_wsle();
	s_wsle(&io___718);
	e_wsle();
	s_wsfe(&io___719);
	e_wsfe();
	s_wsle(&io___720);
	e_wsle();
	sort_(&igc, ff20, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 20) {
		goto L20002;
	}
	j = ll[igc + 1 - i__ - 1];
	ipedig = j;
	for (k = 1; k <= 6; ++k) {
/* L20671: */
		celpre[k - 1] = cel[k + j * 6 - 7];
	}
	dcell_(celpre, cal_1.al, &v1);
	if (ib[j - 1] == 1) {
		*(unsigned char *)bl = 'I';
	}
	if (ib[j - 1] == 2) {
		*(unsigned char *)bl = 'A';
	}
	if (ib[j - 1] == 3) {
		*(unsigned char *)bl = 'B';
	}
	if (ib[j - 1] == 4) {
		*(unsigned char *)bl = 'C';
	}
	if (ib[j - 1] == 5) {
		*(unsigned char *)bl = 'F';
	}
	if (ib[j - 1] == 6) {
		*(unsigned char *)bl = 'P';
	}
	if (ifi[j - 1] == 7) {
		*(unsigned char *)bl = 'R';
	}
	if (ifi[j - 1] == 1) {
		s_copy(more, "Cubic *****", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 2) {
		s_copy(more, "Hexag **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 3) {
		s_copy(more, "Tetra **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 4) {
		s_copy(more, "Ortho ***  ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 5) {
		s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 6) {
		s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 7) {
		s_copy(more, "Rhomb **** ", (ftnlen)11, (ftnlen)11);
	}
	vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];
	s_wsfe(&io___721);
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ff20[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	do_fio(&c__1, bl, (ftnlen)1);
	do_fio(&c__1, more, (ftnlen)11);
	e_wsfe();
/* L20001: */
	}
L20002:
	s_wsle(&io___722);
	e_wsle();
	s_wsle(&io___723);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___724);
	do_lio(&c__9, &c__1, "    FINAL LIST OF CELL PROPOSALS, sorted by M(20) :"
		, (ftnlen)51);
	e_wsle();
	s_wsle(&io___725);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___726);
	e_wsle();
	s_wsle(&io___727);
	do_lio(&c__9, &c__1, "      Global list   (IN=number of indexed lines)", (
		ftnlen)48);
	e_wsle();
	s_wsle(&io___728);
	do_lio(&c__9, &c__1, "  The correct cell has some chances to be just bel"
		"ow", (ftnlen)52);
	e_wsle();
	s_wsle(&io___729);
	e_wsle();
	s_wsfe(&io___730);
	e_wsfe();
	s_wsle(&io___731);
	e_wsle();
	sort_(&igc, fm20, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 20) {
		goto L20004;
	}
	j = ll[igc + 1 - i__ - 1];
	ipedig = j;
	for (k = 1; k <= 6; ++k) {
/* L20672: */
		celpre[k - 1] = cel[k + j * 6 - 7];
	}
	dcell_(celpre, cal_1.al, &v1);
	if (ib[j - 1] == 1) {
		*(unsigned char *)bl = 'I';
	}
	if (ib[j - 1] == 2) {
		*(unsigned char *)bl = 'A';
	}
	if (ib[j - 1] == 3) {
		*(unsigned char *)bl = 'B';
	}
	if (ib[j - 1] == 4) {
		*(unsigned char *)bl = 'C';
	}
	if (ib[j - 1] == 5) {
		*(unsigned char *)bl = 'F';
	}
	if (ib[j - 1] == 6) {
		*(unsigned char *)bl = 'P';
	}
	if (ifi[j - 1] == 7) {
		*(unsigned char *)bl = 'R';
	}
	if (ifi[j - 1] == 1) {
		s_copy(more, "Cubic *****", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 2) {
		s_copy(more, "Hexag **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 3) {
		s_copy(more, "Tetra **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 4) {
		s_copy(more, "Ortho ***  ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 5) {
		s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 6) {
		s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 7) {
		s_copy(more, "Rhomb **** ", (ftnlen)11, (ftnlen)11);
	}
	vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];
	s_wsfe(&io___732);
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&fm20[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	do_fio(&c__1, bl, (ftnlen)1);
	do_fio(&c__1, more, (ftnlen)11);
	e_wsfe();
/* L20003: */
	}
L20004:
	s_wsle(&io___733);
	e_wsle();
	s_wsle(&io___734);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___735);
	do_lio(&c__9, &c__1, "  See for the highest F.o.M. above the cell(s) wit"
		"h \r highest symmetry, if any", (ftnlen)78);
	e_wsle();
	s_wsle(&io___736);
	do_lio(&c__9, &c__1, "  (Cubic, hexagonal, etc), they could correspond t"
		"o \r the the right solution", (ftnlen)76);
	e_wsle();
	s_wsle(&io___737);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___738);
	e_wsle();
	cl__1.cerr = 0;
	cl__1.cunit = 24;
	cl__1.csta = 0;
	f_clos(&cl__1);
	cl__1.cerr = 0;
	cl__1.cunit = 25;
	cl__1.csta = 0;
	f_clos(&cl__1);

/*   Make output for cells sorted by symmetry */

	sort_(&igc, rp, ll);
	s_wsle(&io___739);
	e_wsle();
	s_wsle(&io___740);
	do_lio(&c__9, &c__1, "  Cells sorted by symmetry", (ftnlen)26);
	e_wsle();
	s_wsle(&io___741);
	e_wsle();
	s_wsfe(&io___742);
	e_wsfe();
/* L19992: */
	for (i__ = 1; i__ <= 7; ++i__) {
/* L2018: */
	isyst[i__ - 1] = 0;
	}
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (ifi[i__ - 1] == 1) {
		isyst[0] = 1;
	}
	if (ifi[i__ - 1] == 2) {
		isyst[1] = 1;
	}
	if (ifi[i__ - 1] == 3) {
		isyst[2] = 1;
	}
	if (ifi[i__ - 1] == 4) {
		isyst[3] = 1;
	}
	if (ifi[i__ - 1] == 5) {
		isyst[4] = 1;
	}
	if (ifi[i__ - 1] == 6) {
		isyst[5] = 1;
	}
	if (ifi[i__ - 1] == 7) {
		isyst[6] = 1;
	}
/* L2019: */
	}
	for (jifi = 1; jifi <= 7; ++jifi) {
	if (isyst[jifi - 1] == 0) {
		goto L2020;
	}
	s_wsle(&io___745);
	e_wsle();
	switch (jifi) {
		case 1:  goto L2021;
		case 2:  goto L2022;
		case 3:  goto L2023;
		case 4:  goto L2024;
		case 5:  goto L2025;
		case 6:  goto L2026;
		case 7:  goto L2027;
	}
L2021:
	s_wsle(&io___746);
	do_lio(&c__9, &c__1, "   Cubic cells", (ftnlen)14);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_cub.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2521;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2521:
	open_write1__(&c__24, tempo, (ftnlen)80);
	goto L2028;
L2022:
	s_wsle(&io___747);
	do_lio(&c__9, &c__1, "   Hexagonal/trigonal cells", (ftnlen)27);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_hex.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2522;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2522:
	open_write1__(&c__24, tempo, (ftnlen)80);
	goto L2028;
L2023:
	s_wsle(&io___748);
	do_lio(&c__9, &c__1, "   Tetragonal cells", (ftnlen)19);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_tet.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2523;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2523:
	open_write1__(&c__24, tempo, (ftnlen)80);
	goto L2028;
L2024:
	s_wsle(&io___749);
	do_lio(&c__9, &c__1, "   Orthorhombic cells", (ftnlen)21);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_ort.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2524;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2524:
	open_write1__(&c__24, tempo, (ftnlen)80);
	goto L2028;
L2025:
	s_wsle(&io___750);
	do_lio(&c__9, &c__1, "   Monoclinic cells", (ftnlen)19);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_mon.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2525;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2525:
	open_write1__(&c__24, tempo, (ftnlen)80);
	goto L2028;
L2026:
	s_wsle(&io___751);
	do_lio(&c__9, &c__1, "   Triclinic cells", (ftnlen)18);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_tri.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2526;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2526:
	open_write1__(&c__24, tempo, (ftnlen)80);
	goto L2028;
L2027:
	s_wsle(&io___752);
	do_lio(&c__9, &c__1, "   Rhombohedral cells", (ftnlen)21);
	e_wsle();
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_rho.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
		goto L2527;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2527:
	open_write1__(&c__24, tempo, (ftnlen)80);
L2028:
	s_wsle(&io___753);
	e_wsle();
	jjj = 0;
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (i__ > 20) {
		goto L20060;
		}
		j = ll[i__ - 1];
		if (ifi[j - 1] != jifi) {
		goto L2006;
		}
		++jjj;
		lll[jjj - 1] = j;
		if (rp[j - 1] < .001f) {
		rp[j - 1] = .001f;
		}
		x = 1.f / rp[j - 1] * 5.f;
		vr = vgc[j - 1] / vgc[lll[0] - 1];
		s_wsfe(&io___756);
		do_fio(&c__1, (char *)&rp[j - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&nsol[j - 1], (ftnlen)sizeof(integer));
		for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(
			real));
		}
		e_wsfe();
		s_wsfe(&io___757);
		do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
		for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(
			real));
		}
		e_wsfe();
L2006:
		;
	}
L20060:
	nsolmax = nsol[lll[0] - 1];
	if (nsolmax > 5) {
		s_wsle(&io___759);
		e_wsle();
		s_wsle(&io___760);
		do_lio(&c__9, &c__1, "WARNING - WARNING - WARNING :", (ftnlen)29);
		e_wsle();
		s_wsle(&io___761);
		do_lio(&c__9, &c__1, "Same solution found Nsol = ", (ftnlen)27);
		do_lio(&c__3, &c__1, (char *)&nsolmax, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, " times,", (ftnlen)7);
		e_wsle();
		s_wsle(&io___762);
		do_lio(&c__9, &c__1, "you should probably reduce the test number"
			"s...", (ftnlen)46);
		e_wsle();
		s_wsle(&io___763);
		e_wsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 24;
	cl__1.csta = 0;
	f_clos(&cl__1);
L2020:
	;
	}

/* ... Refine the "best" cell if this was not already done */


/*     READ hkl Miller indices in *.hkl */

	ifile = ifi[ibest - 1];
	switch (ifile) {
	case 1:  goto L31;
	case 2:  goto L32;
	case 3:  goto L33;
	case 4:  goto L34;
	case 5:  goto L35;
	case 6:  goto L36;
	case 7:  goto L37;
	}
L31:
	s_copy(tempo, "cub.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___764);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 6;
	if (cal_1.nhkl0 > 400) {
	cal_1.nhkl0 = 400;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L41: */
	s_rsle(&io___765);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L50;
L32:
	s_copy(tempo, "hex.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___766);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L42: */
	s_rsle(&io___767);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L50;
L33:
	s_copy(tempo, "tet.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___768);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L43: */
	s_rsle(&io___769);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L50;
L34:
	s_copy(tempo, "ort.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___770);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L44: */
	s_rsle(&io___771);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L50;
L35:
	s_copy(tempo, "mon.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___772);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L45: */
	s_rsle(&io___773);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L50;
L36:
	s_copy(tempo, "tri.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___774);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L46: */
	s_rsle(&io___775);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L50;
L37:
	s_copy(tempo, "rho.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___776);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 600) {
	cal_1.nhkl0 = 600;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L47: */
	s_rsle(&io___777);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);

L50:

/*      IF(IREF.EQ.1)GO TO 5900 */
	s_wsle(&io___778);
	e_wsle();
	s_wsle(&io___779);
	do_lio(&c__9, &c__1, "    \"Best\" cell with largest McM20 :", (ftnlen)36)
		;
	e_wsle();
	s_wsle(&io___780);
	do_lio(&c__9, &c__1, "    --------------------------------", (ftnlen)36);
	e_wsle();
	s_wsle(&io___781);
	e_wsle();
	j = ibest;
	ifile = ifi[j - 1];
	indic = 0;
	if (ifile == 1) {
	indic = 1;
	}
	if (ifile == 2) {
	indic = 2;
	}
	if (ifile == 3) {
	indic = 2;
	}
	if (ifile == 7) {
	indic = 2;
	}
	bb[2] = cel[j * 6 - 6];
	bb[3] = cel[j * 6 - 5];
	bb[4] = cel[j * 6 - 4];
	bb[5] = cel[j * 6 - 3];
	bb[6] = cel[j * 6 - 2];
	bb[7] = cel[j * 6 - 1];
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 1.f;
	afi[6] = 1.f;
	afi[7] = 1.f;
	if (ifile == 1) {
	afi[5] = 0.f;
	}
	if (ifile == 1) {
	afi[6] = 0.f;
	}
	if (ifile == 1) {
	afi[7] = 0.f;
	}
	if (ifile == 2) {
	afi[5] = 0.f;
	}
	if (ifile == 2) {
	afi[6] = 0.f;
	}
	if (ifile == 2) {
	afi[7] = 0.f;
	}
	if (ifile == 3) {
	afi[5] = 0.f;
	}
	if (ifile == 3) {
	afi[6] = 0.f;
	}
	if (ifile == 3) {
	afi[7] = 0.f;
	}
	if (ifile == 4) {
	afi[5] = 0.f;
	}
	if (ifile == 4) {
	afi[6] = 0.f;
	}
	if (ifile == 4) {
	afi[7] = 0.f;
	}
	if (ifile == 5) {
	afi[5] = 0.f;
	}
	if (ifile == 5) {
	afi[7] = 0.f;
	}
	if (ifile == 7) {
	afi[5] = 0.f;
	}
	if (ifile == 7) {
	afi[6] = 0.f;
	}
	if (ifile == 7) {
	afi[7] = 0.f;
	}
	celpre[0] = bb[2];
	celpre[1] = bb[3];
	celpre[2] = bb[4];
	celpre[3] = bb[5];
	celpre[4] = bb[6];
	celpre[5] = bb[7];
	bb[0] = 0.f;
/* zeropoint after correction... = 0. */
	for (i__ = 1; i__ <= 3; ++i__) {
	for (jj = 1; jj <= 3; ++jj) {
/* L6010: */
		cal_1.al[i__ + jj * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &j);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
	cncalc[igc - 1] = (real) ncalc;
	fm20[j - 1] = qo[19] / (cncalc[j - 1] * 2.f * ddq);
	ff20[j - 1] = 20.f / (cncalc[j - 1] * ddt);
	s_wsfe(&io___783);
	do_fio(&c__1, (char *)&fm20[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___784);
	do_fio(&c__1, (char *)&ff20[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___785);
	e_wsle();
	s_wsfe(&io___786);
	do_fio(&c__1, (char *)&fm20[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___787);
	do_fio(&c__1, (char *)&ff20[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___788);
	e_wsle();
	}
/* L5900: */

/*   Make a .prf with the best solution */


/* .....START TO GENERATE THE "OBSERVED" PROFILE */

/*   Sum on observed peak positions */

	xnhkl = (real) cal_1.ndat;
	hmin = 100.f;
	i__3 = cal_1.ndat;
	for (i__ = 1; i__ <= i__3; ++i__) {
	cal_1.th2[i__ - 1] += bb[0];
/* L1800: */
	cal_1.fobs[i__ - 1] = cal_1.fobs[i__ - 1] / cal_1.sum_f__ * xnhkl * 
		50.f;
	}
	cal_1.sum_f__ = xnhkl * 50.f;
	flog = log(2.f) * 4.f;
	slog = 1.f / sqrt(flog);
	if (w < 0.f) {
	w = 0.f;
	x = 0.f;
	i__3 = nhkl;
	for (i__ = 1; i__ <= i__3; ++i__) {
		w += w1[i__ - 1];
		x += 1.f;
/* L1801: */
	}
	w /= x;
	}

/*  Use FWHM 1/2 the entered value... */

/* Computing 2nd power */
	r__1 = w / 2.f;
	w = r__1 * r__1;


	i__3 = nhkl;
	for (nm = 1; nm <= i__3; ++nm) {
	dd = slabda / (sin(cal_1.th2[nm - 1] / cal_1.pi) * 2.f);
/* Computing 2nd power */
	r__1 = dd;
	sss = 1 / (r__1 * r__1);

/*     Needs to calculate halfwidths and positions */

	d__[nm - 1] = sqrt(1.f / sss);
	sinth = slabda * slabda * sss / 4.f;
	costh = 1.f - sinth;
	tanth = sqrt(sinth / costh);
	theta[nm - 1] = atan(tanth) * cal_1.pi;
	hw[nm - 1] = sqrt(w);
	if (hmin > hw[nm - 1]) {
		hmin = hw[nm - 1];
	}
	hw[nm - 1] *= slog;
	hw4[nm - 1] = hw[nm - 1] / slog * 2.f;
	bbb[nm - 1] = hw[nm - 1] * hw[nm - 1];
/* L1802: */
	}
	dmi = d__[nhkl - 1] - .02f;

/*  Step and positions */
/*  The minimum step is FWHMmin/STEPN */

	step = hmin / 4.f;
	astep = step * 1e3f;
	istep = astep;
	astep = (real) istep;
	step = astep / 1e3f;
	pos[0] = step;
	po = step;
	for (i__ = 2; i__ <= 16000; ++i__) {
	po += step;
	if (po > 160.f) {
		goto L1804;
	}
/* L1803: */
	pos[i__ - 1] = po;
	}
L1804:
	npts = i__ - 1;
/*  Sum at each point */
	i__3 = npts;
	for (k = 1; k <= i__3; ++k) {
	yobs[k - 1] = 0.f;
/*  SUM on all hkl */
	i__1 = nhkl;
	for (j = 1; j <= i__1; ++j) {
		deltap = pos[k - 1] - theta[j - 1];
		if (dabs(deltap) > hw4[j - 1]) {
		goto L1806;
		}
		delt = deltap * deltap;
		omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
		yobs[k - 1] += omeg * cal_1.fobs[j - 1];
L1806:
		;
	}
	if (yobs[k - 1] < .01f) {
		yobs[k - 1] = 0.f;
	}
/* L1805: */
	}
	l = npts + 1;
	i__3 = npts;
	for (k = 1; k <= i__3; ++k) {
	--l;
/* L1807: */
	if (yobs[l - 1] != 0.f) {
		goto L1808;
	}
	}
L1808:
	n2 = l;
	n1 = 1;

	for (i__ = 1; i__ <= 6; ++i__) {
	k = i__ + 2;
/* L1809: */
	celpre[i__ - 1] = bb[k - 1];
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 3; ++j) {
/* L55: */
		cal_1.al[i__ + j * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);

/* ...  Keep only the hkl for d > dmin */

	jh = 0;
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
	for (kk = 1; kk <= 3; ++kk) {
/* L56: */
		hh[kk - 1] = (real) cal_1.ihh[kk + i__ * 3 - 4];
	}
/*     HH IS INDICES OF REFLECTION */
	x = 0.f;
	for (j = 1; j <= 3; ++j) {
		for (k = j; k <= 3; ++k) {
/* L57: */
		x = cal_1.al[j + k * 3 - 4] * hh[j - 1] * hh[k - 1] + x;
		}
	}
	x = 1 / sqrt(x);
	if (x < dmi) {
		goto L59;
	}
	++jh;
	for (kk = 1; kk <= 3; ++kk) {
/* L58: */
		jhh[kk + jh * 3 - 4] = cal_1.ihh[kk + i__ * 3 - 4];
	}
	fcal[jh - 1] = 50.f;
	d__[jh - 1] = x;
/*     X IS D(hkl) FOR REFLECTION HH */
L59:
	;
	}
	nhkl = jh;

/*   Again, calculate 2-theta, etc. */

	i__3 = nhkl;
	for (nm = 1; nm <= i__3; ++nm) {

	dd = d__[nm - 1];
/* Computing 2nd power */
	r__1 = dd;
	sss = 1 / (r__1 * r__1);

/*     Needs to calculate halfwidths and positions again */

	sinth = slabda * slabda * sss / 4.f;
	costh = 1.f - sinth;
	tanth = sqrt(sinth / costh);
	theta[nm - 1] = atan(tanth) * cal_1.pi;
	hw[nm - 1] = sqrt(w);
	if (hmin > hw[nm - 1]) {
		hmin = hw[nm - 1];
	}
	hw[nm - 1] *= slog;
	hw4[nm - 1] = hw[nm - 1] / slog * 2.f;
	bbb[nm - 1] = hw[nm - 1] * hw[nm - 1];
/* L60: */
	}

/* ...  Calculate best Yobs */

/*  Sum at each point */

	i__3 = nhkl;
	for (k = 1; k <= i__3; ++k) {
	somega[k - 1] = 0.f;
/* L1810: */
	cal_1.fobs[k - 1] = 0.f;
	}
	i__3 = n2;
	for (k = 1; k <= i__3; ++k) {
	ycalc[k - 1] = 0.f;
/*  SUM on all hkl */
	kpos = 0;
	omegt = 0.f;
	i__1 = nhkl;
	for (j = 1; j <= i__1; ++j) {
		deltap = pos[k - 1] - theta[j - 1];
		if (dabs(deltap) > hw4[j - 1]) {
		goto L1812;
		}
		if (kpos == 0) {
		nha[k - 1] = j;
		}
		kpos = 1;
		nhb[k - 1] = j;
		delt = deltap * deltap;
		omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
		somega[j - 1] += omeg;
		dump[j - 1] = fcal[j - 1] * omeg;
		ycalc[k - 1] += dump[j - 1];
L1812:
		;
	}
	if (ycalc[k - 1] == 0.f) {
		goto L1813;
	}
	yoy = yobs[k - 1] / ycalc[k - 1];
	i__1 = nhb[k - 1];
	for (j = nha[k - 1]; j <= i__1; ++j) {
		cal_1.fobs[j - 1] += dump[j - 1] * yoy;
/* L1814: */
	}
L1813:
/* L1811: */
	;
	}
	i__3 = nhkl;
	for (k = 1; k <= i__3; ++k) {
	if (somega[k - 1] == 0.f) {
		goto L1816;
	}
	cal_1.fobs[k - 1] /= somega[k - 1];
	goto L1815;
L1816:
	cal_1.fobs[k - 1] = 0.f;
L1815:
	;
	}
	i__3 = nhkl;
	for (j = 1; j <= i__3; ++j) {
	fcal[j - 1] = cal_1.fobs[j - 1];
/* L61: */
	}

/* ...  2 more iterations by Le Bail fit */

	for (kiter = 1; kiter <= 2; ++kiter) {
/*  Sum at each point */
	i__3 = nhkl;
	for (k = 1; k <= i__3; ++k) {
		somega[k - 1] = 0.f;
/* L1817: */
		cal_1.fobs[k - 1] = 0.f;
	}
	i__3 = n2;
	for (k = 1; k <= i__3; ++k) {
		ycalc[k - 1] = 0.f;
/*  SUM on all hkl */
		omegt = 0.f;
		i__1 = nhb[k - 1];
		for (j = nha[k - 1]; j <= i__1; ++j) {
		deltap = pos[k - 1] - theta[j - 1];
		if (dabs(deltap) > hw4[j - 1]) {
			goto L1819;
		}
		delt = deltap * deltap;
		omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
		somega[j - 1] += omeg;
		dump[j - 1] = fcal[j - 1] * omeg;
		ycalc[k - 1] += dump[j - 1];
L1819:
		;
		}
		if (ycalc[k - 1] == 0.f) {
		goto L1820;
		}
		yoy = yobs[k - 1] / ycalc[k - 1];
		i__1 = nhb[k - 1];
		for (j = nha[k - 1]; j <= i__1; ++j) {
		cal_1.fobs[j - 1] += dump[j - 1] * yoy;
/* L1821: */
		}
L1820:
/* L1818: */
		;
	}
	i__3 = nhkl;
	for (k = 1; k <= i__3; ++k) {
		if (somega[k - 1] == 0.f) {
		goto L1823;
		}
		cal_1.fobs[k - 1] /= somega[k - 1];
		goto L1822;
L1823:
		cal_1.fobs[k - 1] = 0.f;
L1822:
		;
	}
	i__3 = nhkl;
	for (j = 1; j <= i__3; ++j) {
		fcal[j - 1] = cal_1.fobs[j - 1];
/* L62: */
	}
/* L1850: */
	}



	diff = 0.f;
/*  Sum at each point */
	i__3 = n2;
	for (k = 1; k <= i__3; ++k) {
	ycalc[k - 1] = 0.f;
	if (yobs[k - 1] == 0.f) {
		goto L1824;
	}
/*  SUM on all hkl */
	omegt = 0.f;
	i__1 = nhb[k - 1];
	for (j = nha[k - 1]; j <= i__1; ++j) {
		deltap = pos[k - 1] - theta[j - 1];
		if (dabs(deltap) > hw4[j - 1]) {
		goto L1825;
		}
		delt = deltap * deltap;
		omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
		ycalc[k - 1] += fcal[j - 1] * omeg;
L1825:
		;
	}
L1824:
	;
	}
	sum_y__ = 0.f;
	i__3 = n2;
	for (k = 1; k <= i__3; ++k) {
	diff += (r__1 = yobs[k - 1] - ycalc[k - 1], dabs(r__1));
	sum_y__ += yobs[k - 1];
/* L1826: */
	}
	diff /= sum_y__;
	s_wsle(&io___828);
	e_wsle();
	s_wsle(&io___829);
	do_lio(&c__9, &c__1, " Final Rp on the .prf = ", (ftnlen)24);
	do_lio(&c__4, &c__1, (char *)&diff, (ftnlen)sizeof(real));
	e_wsle();
	s_wsle(&io___830);
	e_wsle();

/*     Make the .prf */

/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 4, a__1[1] = ".prf";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
	goto L1827;
	}
	filedel_(&c__12, tempo, (ftnlen)80);
L1827:
	open_write1__(&c__12, tempo, (ftnlen)80);
	ivers = 8;
	zero = 0.f;
	s_wsfe(&io___832);
	for (i__ = 1; i__ <= 20; ++i__) {
	do_fio(&c__1, (char *)&text[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___833);
	do_fio(&c__1, "3111   1.0000", (ftnlen)13);
	do_fio(&c__1, (char *)&zero, (ftnlen)sizeof(real));
	e_wsfe();
	thmax = pos[n2 - 1];
	thmin = pos[n1 - 1];
	icn = nhkl;
	i__3 = nhkl;
	for (i__ = 1; i__ <= i__3; ++i__) {
	i1 = jhh[i__ * 3 - 3];
	i2 = jhh[i__ * 3 - 2];
	i3 = jhh[i__ * 3 - 1];
/* L1828: */
	irefs[i__ - 1] = ((i1 + 2432 << 8) + 128 + i2 << 8) + 128 + i3;
	}
	icz = nhkl;
	npts = n2 - n1 + 1;
	s_wsfe(&io___842);
	do_fio(&c__1, (char *)&thmax, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&thmin, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&step, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ivers, (ftnlen)sizeof(integer));
	e_wsfe();
	amda1 = slabda;
	amda2 = slabda;
	npat1 = 1;
	nvk = 0;
	s_wsfe(&io___847);
	do_fio(&c__1, (char *)&npat1, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&npts, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&amda1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&amda2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&zero, (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___848);
	do_fio(&c__1, (char *)&icn, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nvk, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___849);
	i__3 = n2;
	for (j = n1; j <= i__3; ++j) {
	do_fio(&c__1, (char *)&yobs[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___850);
	i__3 = n2;
	for (j = n1; j <= i__3; ++j) {
	do_fio(&c__1, (char *)&ycalc[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___851);
	i__3 = icz;
	for (j = 1; j <= i__3; ++j) {
	do_fio(&c__1, (char *)&irefs[j - 1], (ftnlen)sizeof(integer));
	}
	e_wsfe();
	s_wsfe(&io___852);
	i__3 = icz;
	for (j = 1; j <= i__3; ++j) {
	do_fio(&c__1, (char *)&theta[j - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	nexcrg = 0;
	excrg = 0.f;
	s_wsfe(&io___855);
	do_fio(&c__1, (char *)&nexcrg, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsfe(&io___856);
	do_fio(&c__1, (char *)&excrg, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&excrg, (ftnlen)sizeof(real));
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = 12;
	cl__1.csta = 0;
	f_clos(&cl__1);



	if (igc == 1) {
	goto L6000;
	}

/*   Sort cells by volume */

	s_wsle(&io___857);
	e_wsle();
	s_wsle(&io___858);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___859);
	do_lio(&c__9, &c__1, "       CELL PROPOSALS sorted by increasing volume :"
		, (ftnlen)51);
	e_wsle();
	s_wsle(&io___860);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___861);
	e_wsle();
	s_wsle(&io___862);
	e_wsle();
	sort_(&igc, vgc, ll);
	++imem;
	im[imem - 1] = ll[0];
	s_wsfe(&io___863);
	e_wsfe();
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 20) {
		goto L20030;
	}
	j = ll[i__ - 1];
	if (rp[j - 1] < .001f) {
		rp[j - 1] = .001f;
	}
	vr = vgc[j - 1] / vgc[ll[0] - 1];
	s_wsfe(&io___864);
	do_fio(&c__1, (char *)&rp[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nsol[j - 1], (ftnlen)sizeof(integer));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
/* L2003: */
	}
L20030:
	nsolmax = nsol[ll[0] - 1];
	if (nsolmax > 5) {
	s_wsle(&io___865);
	e_wsle();
	s_wsle(&io___866);
	do_lio(&c__9, &c__1, "WARNING - WARNING - WARNING :", (ftnlen)29);
	e_wsle();
	s_wsle(&io___867);
	do_lio(&c__9, &c__1, "Same solution found Nsol = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&nsolmax, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " times,", (ftnlen)7);
	e_wsle();
	s_wsle(&io___868);
	do_lio(&c__9, &c__1, "you should probably reduce the test numbers...",
		 (ftnlen)46);
	e_wsle();
	s_wsle(&io___869);
	e_wsle();
	}

/* ... Refine the "best" cell if this was not already done */
/*          this time with smallest volume and if Rp < 10% */

	if (iref == 1) {
	goto L5901;
	}
	j = ll[0];
	if (rp[j - 1] > .1f) {
	goto L5901;
	}

/*     READ hkl Miller indices in *.hkl */

	ifile = ifi[ll[0] - 1];
	switch (ifile) {
	case 1:  goto L1731;
	case 2:  goto L1732;
	case 3:  goto L1733;
	case 4:  goto L1734;
	case 5:  goto L1735;
	case 6:  goto L1736;
	case 7:  goto L1737;
	}
L1731:
	s_copy(tempo, "cub.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___870);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 6;
	if (cal_1.nhkl0 > 400) {
	cal_1.nhkl0 = 400;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1741: */
	s_rsle(&io___871);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L1750;
L1732:
	s_copy(tempo, "hex.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___872);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1742: */
	s_rsle(&io___873);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L1750;
L1733:
	s_copy(tempo, "tet.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___874);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 800) {
	cal_1.nhkl0 = 800;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1743: */
	s_rsle(&io___875);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L1750;
L1734:
	s_copy(tempo, "ort.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___876);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1744: */
	s_rsle(&io___877);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L1750;
L1735:
	s_copy(tempo, "mon.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___878);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1745: */
	s_rsle(&io___879);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L1750;
L1736:
	s_copy(tempo, "tri.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___880);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 20;
	if (cal_1.nhkl0 > 1000) {
	cal_1.nhkl0 = 1000;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1746: */
	s_rsle(&io___881);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);
	goto L1750;
L1737:
	s_copy(tempo, "rho.hkl", (ftnlen)80, (ftnlen)7);
	open_read1__(&c__35, tempo, (ftnlen)80);
	s_rsle(&io___882);
	do_lio(&c__3, &c__1, (char *)&cal_1.nhkl0, (ftnlen)sizeof(integer));
	e_rsle();
	cal_1.nhkl0 = cal_1.ndat * 12;
	if (cal_1.nhkl0 > 600) {
	cal_1.nhkl0 = 600;
	}
	i__3 = cal_1.nhkl0;
	for (i__ = 1; i__ <= i__3; ++i__) {
/* L1747: */
	s_rsle(&io___883);
	for (kk = 1; kk <= 3; ++kk) {
		do_lio(&c__3, &c__1, (char *)&cal_1.ihh[kk + i__ * 3 - 4], (
			ftnlen)sizeof(integer));
	}
	e_rsle();
	}
	cl__1.cerr = 0;
	cl__1.cunit = 35;
	cl__1.csta = 0;
	f_clos(&cl__1);


L1750:
	s_wsle(&io___884);
	e_wsle();
	s_wsle(&io___885);
	do_lio(&c__9, &c__1, "    \"Best\" cell with smallest volume :", (ftnlen)
		38);
	e_wsle();
	s_wsle(&io___886);
	do_lio(&c__9, &c__1, "    ---------------------------------", (ftnlen)37);
	e_wsle();
	s_wsle(&io___887);
	e_wsle();
	ifile = ifi[j - 1];
	indic = 0;
	if (ifile == 1) {
	indic = 1;
	}
	if (ifile == 2) {
	indic = 2;
	}
	if (ifile == 3) {
	indic = 2;
	}
	if (ifile == 7) {
	indic = 2;
	}
	bb[2] = cel[j * 6 - 6];
	bb[3] = cel[j * 6 - 5];
	bb[4] = cel[j * 6 - 4];
	bb[5] = cel[j * 6 - 3];
	bb[6] = cel[j * 6 - 2];
	bb[7] = cel[j * 6 - 1];
	afi[2] = 1.f;
	afi[3] = 1.f;
	afi[4] = 1.f;
	afi[5] = 1.f;
	afi[6] = 1.f;
	afi[7] = 1.f;
	if (ifile == 1) {
	afi[5] = 0.f;
	}
	if (ifile == 1) {
	afi[6] = 0.f;
	}
	if (ifile == 1) {
	afi[7] = 0.f;
	}
	if (ifile == 2) {
	afi[5] = 0.f;
	}
	if (ifile == 2) {
	afi[6] = 0.f;
	}
	if (ifile == 2) {
	afi[7] = 0.f;
	}
	if (ifile == 3) {
	afi[5] = 0.f;
	}
	if (ifile == 3) {
	afi[6] = 0.f;
	}
	if (ifile == 3) {
	afi[7] = 0.f;
	}
	if (ifile == 4) {
	afi[5] = 0.f;
	}
	if (ifile == 4) {
	afi[6] = 0.f;
	}
	if (ifile == 4) {
	afi[7] = 0.f;
	}
	if (ifile == 5) {
	afi[5] = 0.f;
	}
	if (ifile == 5) {
	afi[7] = 0.f;
	}
	if (ifile == 7) {
	afi[5] = 0.f;
	}
	if (ifile == 7) {
	afi[6] = 0.f;
	}
	if (ifile == 7) {
	afi[7] = 0.f;
	}
	celpre[0] = bb[2];
	celpre[1] = bb[3];
	celpre[2] = bb[4];
	celpre[3] = bb[5];
	celpre[4] = bb[6];
	celpre[5] = bb[7];
	bb[0] = 0.f;
/* zeropoint after correction... = 0. */
	for (i__ = 1; i__ <= 3; ++i__) {
	for (jj = 1; jj <= 3; ++jj) {
/* L6011: */
		cal_1.al[i__ + jj * 3 - 4] = 0.f;
	}
	}
	dcell_(celpre, cal_1.al, &v1);
	calcul2_(&diff, ihkl, th3, &ncalc, &j);
	celref_(&indic, bb, afi, &cal_1.lhkl, th3, ihkl, &ddt, &ddq);
	if (cal_1.ndat >= 20) {
	cncalc[j - 1] = (real) ncalc;
	fm20[j - 1] = qo[19] / (cncalc[j - 1] * 2.f * ddq);
	ff20[j - 1] = 20.f / (cncalc[j - 1] * ddt);
	s_wsfe(&io___888);
	do_fio(&c__1, (char *)&fm20[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___889);
	do_fio(&c__1, (char *)&ff20[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___890);
	e_wsle();
	s_wsfe(&io___891);
	do_fio(&c__1, (char *)&fm20[j - 1], (ftnlen)sizeof(real));
	e_wsfe();
	s_wsfe(&io___892);
	do_fio(&c__1, (char *)&ff20[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ddt, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&ncalc, (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___893);
	e_wsle();
	}
L5901:



	if (igc == 1) {
	goto L6000;
	}

/*   Sort cells, the most frequently found first */

	s_wsle(&io___894);
	e_wsle();
	s_wsle(&io___895);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___896);
	do_lio(&c__9, &c__1, "       CELL PROPOSALS most frequently found :", (
		ftnlen)45);
	e_wsle();
	s_wsle(&io___897);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___898);
	e_wsle();
	s_wsle(&io___899);
	e_wsle();
	sort2_(&igc, nsol, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
	s_wsfe(&io___900);
	e_wsfe();
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 20) {
		goto L20040;
	}
	j = ll[igc + 1 - i__ - 1];
	if (nsol[j - 1] < 2) {
		goto L2004;
	}
	if (rp[j - 1] < .001f) {
		rp[j - 1] = .001f;
	}
	vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];
	s_wsfe(&io___901);
	do_fio(&c__1, (char *)&rp[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nsol[j - 1], (ftnlen)sizeof(integer));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
L2004:
	;
	}
L20040:


/*   Sort cells with largest number of peak indexed + small Rp2 */

	s_wsle(&io___902);
	e_wsle();
	s_wsle(&io___903);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___904);
	do_lio(&c__9, &c__1, "CELLS with small Rp2 + largest number of peak inde"
		"xed", (ftnlen)53);
	e_wsle();
	s_wsle(&io___905);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r===========================", (ftnlen)80);
	e_wsle();
	s_wsle(&io___906);
	do_lio(&c__9, &c__1, "Rp2 is on peak indexed only, and width divided by "
		"2,", (ftnlen)52);
	e_wsle();
	s_wsle(&io___907);
	do_lio(&c__9, &c__1, "while Rp is on all peaks, and large width.", (
		ftnlen)42);
	e_wsle();
	s_wsle(&io___908);
	e_wsle();
/*      CALL SORT2(IGC,RP2,LL) */

	sort2_(&igc, rp2, ll);

	++imem;
	im[imem - 1] = ll[i__ - 1];
	s_wsfe(&io___909);
	e_wsfe();
	i__3 = igc;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 20) {
		goto L20050;
	}
	j = ll[i__ - 1];
	if (rp2[j - 1] < .001f) {
		rp2[j - 1] = .001f;
	}
	if (rp2[j - 1] > .3f) {
		goto L2005;
	}
	if (km2[j - 1] < nmax) {
		goto L2005;
	}
	vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];
	s_wsfe(&io___910);
	do_fio(&c__1, (char *)&rp2[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&km2[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nsol[j - 1], (ftnlen)sizeof(integer));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
L2005:
	;
	}
L20050:


/*   Sort associations of two cells with largest number of peak indexed */

	if (rmax0[0] < .5f) {
	goto L6000;
	}
	s_wsle(&io___911);
	e_wsle();
	s_wsle(&io___912);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___913);
	do_lio(&c__9, &c__1, "Double cells with largest number of peak indexed", (
		ftnlen)48);
	e_wsle();
	s_wsle(&io___914);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r===========================", (ftnlen)80);
	e_wsle();
	s_wsle(&io___915);
	e_wsle();
	s_wsle(&io___916);
	do_lio(&c__9, &c__1, "WARNING - WARNING - WARNING - WARNING - WARNING", (
		ftnlen)47);
	e_wsle();
	s_wsle(&io___917);
	do_lio(&c__9, &c__1, "           This is the two-phase mode", (ftnlen)37);
	e_wsle();
	s_wsle(&io___918);
	do_lio(&c__9, &c__1, "    It could be better to go back to the lab", (
		ftnlen)44);
	e_wsle();
	s_wsle(&io___919);
	do_lio(&c__9, &c__1, "         and try and make a pure sample", (ftnlen)
		39);
	e_wsle();
	s_wsle(&io___920);
	e_wsle();
	ii = 0;
	i__3 = igc - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (rp2[i__ - 1] > .3f) {
		goto L2009;
	}
	if (km2[i__ - 1] < nmax) {
		goto L2009;
	}
	i__1 = igc;
	for (j = i__ + 1; j <= i__1; ++j) {
		if (rp2[j - 1] > .3f) {
		goto L2008;
		}
		if (rp2[i__ - 1] + rp2[j - 1] > .4f) {
		goto L2008;
		}
		if (km2[j - 1] < nmax) {
		goto L2008;
		}
		++ii;
		if (ii > 99999) {
		goto L2011;
		}
		km3[ii - 1] = 0;
		i__4 = cal_1.ndat;
		for (k = 1; k <= i__4; ++k) {
		km3[ii - 1] = km3[ii - 1] + cal2_1.ind[k + i__ * 100 - 101] + 
			cal2_1.ind[k + j * 100 - 101];
		if (cal2_1.ind[k + i__ * 100 - 101] * cal2_1.ind[k + j * 100 
			- 101] == 1) {
			--km3[ii - 1];
		}
/* L2007: */
		}
		if (km3[ii - 1] < cal_1.ndat - 10) {
		--ii;
		goto L2008;
		}
		id1[ii - 1] = i__;
		id2[ii - 1] = j;
L2008:
		;
	}
L2009:
	;
	}
L2011:
	igc2 = ii;
	sort3_(&igc2, km3, ll2);
	s_wsfe(&io___927);
	e_wsfe();
	vr = 1.f;
/* Writing concatenation */
	i__2[0] = lfile, a__1[0] = file;
	i__2[1] = 8, a__1[1] = "_two.ckm";
	s_cat(tempo, a__1, i__2, &c__2, (ftnlen)80);
	ioin__1.inerr = 0;
	ioin__1.infilen = 80;
	ioin__1.infile = tempo;
	ioin__1.inex = &qex;
	ioin__1.inopen = 0;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (! qex) {
	goto L2530;
	}
	filedel_(&c__24, tempo, (ftnlen)80);
L2530:
	open_write1__(&c__24, tempo, (ftnlen)80);
	i__3 = igc2;
	for (i__ = 1; i__ <= i__3; ++i__) {
	if (i__ > 1000) {
		goto L2012;
	}
	jj = ll2[igc2 + 1 - i__ - 1];
	j = id1[ll2[igc2 + 1 - i__ - 1] - 1];
	if (rp2[j - 1] < .001f) {
		rp2[j - 1] = .001f;
	}
	x = 1.f / rp2[j - 1] * 5.f;
	s_wsfe(&io___928);
	do_fio(&c__1, (char *)&rp2[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&km3[jj - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nsol[j - 1], (ftnlen)sizeof(integer));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___929);
	do_fio(&c__1, (char *)&km3[jj - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
	j = id2[ll2[igc2 + 1 - i__ - 1] - 1];
	if (rp2[j - 1] < .001f) {
		rp2[j - 1] = .001f;
	}
	x = 1.f / rp2[j - 1] * 5.f;
	s_wsfe(&io___930);
	do_fio(&c__1, (char *)&rp2[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&km2[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nsol[j - 1], (ftnlen)sizeof(integer));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___931);
	do_fio(&c__1, (char *)&km2[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vr, (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
		do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsle(&io___932);
	e_wsle();
/* L2010: */
	}
L2012:
	cl__1.cerr = 0;
	cl__1.cunit = 24;
	cl__1.csta = 0;
	f_clos(&cl__1);


L6000:
	s_wsle(&io___933);
	e_wsle();
	s_wsle(&io___934);
	e_wsle();
	s_wsle(&io___935);
	e_wsle();
	s_wsle(&io___936);
	e_wsle();
	s_wsle(&io___937);
	e_wsle();
	s_wsle(&io___938);
	e_wsle();
	s_wsle(&io___939);
	e_wsle();
	s_wsle(&io___940);
	e_wsle();
	s_wsle(&io___941);
	e_wsle();
	s_wsle(&io___942);
	e_wsle();
	s_wsle(&io___943);
	e_wsle();
	s_wsle(&io___944);
	e_wsle();
	s_wsle(&io___945);
	e_wsle();
	s_wsle(&io___946);
	e_wsle();
	s_wsle(&io___947);
	e_wsle();
	s_wsle(&io___948);
	e_wsle();
	s_wsle(&io___949);
	e_wsle();
	s_wsle(&io___950);
	e_wsle();
	s_wsle(&io___951);
	e_wsle();
	s_wsle(&io___952);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___953);
	do_lio(&c__9, &c__1, "         THE SELECTION OF THE \"BEST\" CELL", (
		ftnlen)41);
	e_wsle();
	s_wsle(&io___954);
	do_lio(&c__9, &c__1, "based on McM20, Rp, F(20), M(20), V, high symmetry"
		" ?", (ftnlen)52);
	e_wsle();
	s_wsle(&io___955);
	do_lio(&c__9, &c__1, "           DEPENDS ON YOU, EXCLUSIVELY.", (ftnlen)
		39);
	e_wsle();
	s_wsle(&io___956);
	e_wsle();
	s_wsle(&io___957);
	do_lio(&c__9, &c__1, "                    However...", (ftnlen)30);
	e_wsle();
	s_wsle(&io___958);
	do_lio(&c__9, &c__1, "  It is suggested that the correct cell could be :",
		 (ftnlen)50);
	e_wsle();
	s_wsle(&io___959);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___960);
	e_wsle();
	s_wsfe(&io___961);
	e_wsfe();

/*  Ultimate analysis... */

	i__3 = imem;
	for (i__ = 1; i__ <= i__3; ++i__) {
	imn[i__ - 1] = 1;
/* L20010: */
	}
	if (imem > 1) {
	i__3 = imem - 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
		if (imn[i__ - 1] == 0) {
		goto L20008;
		}
		i__1 = imem;
		for (k = i__ + 1; k <= i__1; ++k) {
		if (imn[k - 1] == 0) {
			goto L30008;
		}
		if (im[k - 1] == im[i__ - 1]) {
			++imn[i__ - 1];
			imn[k - 1] = 0;
		}
L30008:
		;
		}
L20008:
		;
	}
	}
	j = im[0];
	if (ib[j - 1] == 1) {
	*(unsigned char *)bl = 'I';
	}
	if (ib[j - 1] == 2) {
	*(unsigned char *)bl = 'A';
	}
	if (ib[j - 1] == 3) {
	*(unsigned char *)bl = 'B';
	}
	if (ib[j - 1] == 4) {
	*(unsigned char *)bl = 'C';
	}
	if (ib[j - 1] == 5) {
	*(unsigned char *)bl = 'F';
	}
	if (ib[j - 1] == 6) {
	*(unsigned char *)bl = 'P';
	}
	if (ifi[j - 1] == 7) {
	*(unsigned char *)bl = 'R';
	}
	if (ifi[j - 1] == 1) {
	s_copy(more, "Cubic *****", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 2) {
	s_copy(more, "Hexag **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 3) {
	s_copy(more, "Tetra **** ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 4) {
	s_copy(more, "Ortho ***  ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 5) {
	s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 6) {
	s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
	}
	if (ifi[j - 1] == 7) {
	s_copy(more, "Rhomb **** ", (ftnlen)11, (ftnlen)11);
	}
	s_wsfe(&io___963);
	do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&xfom[j - 1], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
	for (k = 1; k <= 6; ++k) {
	do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(real));
	}
	do_fio(&c__1, bl, (ftnlen)1);
	do_fio(&c__1, more, (ftnlen)11);
	e_wsfe();
	s_wsle(&io___964);
	do_lio(&c__9, &c__1, "   Found ", (ftnlen)9);
	do_lio(&c__3, &c__1, (char *)&imn[0], (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " time(s) head of the best lists", (ftnlen)31);
	e_wsle();

	if (imem > 1) {
	xf1 = xfom[im[0] - 1] / 2.f;
	imemt = imem;
	i__3 = imem;
	for (i__ = 2; i__ <= i__3; ++i__) {
		if (imn[i__ - 1] == 0) {
		--imemt;
		goto L20011;
		}
		j = im[i__ - 1];
		if (xfom[j - 1] < xf1) {
		--imemt;
		imn[i__ - 1] = 0;
		}
L20011:
		;
	}
	if (imemt > 1) {
		s_wsle(&io___967);
		e_wsle();
		s_wsle(&io___968);
		e_wsle();
		s_wsle(&io___969);
		do_lio(&c__9, &c__1, "   Other(s) having some chance :", (ftnlen)
			32);
		e_wsle();
		s_wsle(&io___970);
		e_wsle();
		i__3 = imem;
		for (i__ = 2; i__ <= i__3; ++i__) {
		if (imn[i__ - 1] == 0) {
			goto L20009;
		}
		j = im[i__ - 1];
		if (ib[j - 1] == 1) {
			*(unsigned char *)bl = 'I';
		}
		if (ib[j - 1] == 2) {
			*(unsigned char *)bl = 'A';
		}
		if (ib[j - 1] == 3) {
			*(unsigned char *)bl = 'B';
		}
		if (ib[j - 1] == 4) {
			*(unsigned char *)bl = 'C';
		}
		if (ib[j - 1] == 5) {
			*(unsigned char *)bl = 'F';
		}
		if (ib[j - 1] == 6) {
			*(unsigned char *)bl = 'P';
		}
		if (ifi[j - 1] == 7) {
			*(unsigned char *)bl = 'R';
		}
		if (ifi[j - 1] == 1) {
			s_copy(more, "Cubic *****", (ftnlen)11, (ftnlen)11);
		}
		if (ifi[j - 1] == 2) {
			s_copy(more, "Hexag **** ", (ftnlen)11, (ftnlen)11);
		}
		if (ifi[j - 1] == 3) {
			s_copy(more, "Tetra **** ", (ftnlen)11, (ftnlen)11);
		}
		if (ifi[j - 1] == 4) {
			s_copy(more, "Ortho ***  ", (ftnlen)11, (ftnlen)11);
		}
		if (ifi[j - 1] == 5) {
			s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
		}
		if (ifi[j - 1] == 6) {
			s_copy(more, "           ", (ftnlen)11, (ftnlen)11);
		}
		if (ifi[j - 1] == 7) {
			s_copy(more, "Rhomb **** ", (ftnlen)11, (ftnlen)11);
		}
		s_wsfe(&io___971);
		do_fio(&c__1, (char *)&km[j - 1], (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&xfom[j - 1], (ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&vgc[j - 1], (ftnlen)sizeof(real));
		for (k = 1; k <= 6; ++k) {
			do_fio(&c__1, (char *)&cel[k + j * 6 - 7], (ftnlen)sizeof(
				real));
		}
		do_fio(&c__1, bl, (ftnlen)1);
		do_fio(&c__1, more, (ftnlen)11);
		e_wsfe();
		s_wsle(&io___972);
		do_lio(&c__9, &c__1, "   Found ", (ftnlen)9);
		do_lio(&c__3, &c__1, (char *)&imn[i__ - 1], (ftnlen)sizeof(
			integer));
		do_lio(&c__9, &c__1, " time(s) head of the best lists", (
			ftnlen)31);
		e_wsle();
L20009:
		;
		}
	}
	}
	s_wsle(&io___973);
	e_wsle();
	s_wsle(&io___974);
	do_lio(&c__9, &c__1, "=================================================="
		"==\r ===========================", (ftnlen)81);
	e_wsle();
	s_wsle(&io___975);
	e_wsle();
	s_wsle(&io___976);
	e_wsle();
	s_wsle(&io___977);
	e_wsle();
	s_wsle(&io___978);
	e_wsle();
	s_wsle(&io___979);
	e_wsle();
	s_wsle(&io___980);
	e_wsle();



/* L10115: */
/* 1115  FORMAT(I3,F10.0,1X,F5.3,F8.4,F9.1,I3) */
/* 1215  FORMAT(I3,F10.0,1X,F5.3,2F8.4,F9.1,I3) */
/* 1415  FORMAT(I3,F10.0,1X,F5.3,3F8.4,F9.1,I3) */
/* 1515  FORMAT(I3,F10.0,1X,F5.3,3F8.4,F7.2,F9.1,I3) */
/* 1615  FORMAT(I3,F10.0,1X,F5.3,3F8.4,3F7.2,F9.1,I3) */
	goto L3600;
L3500:
	s_wsle(&io___981);
	do_lio(&c__9, &c__1, "  Error reading the first line : TEXT", (ftnlen)37);
	e_wsle();
	s_wsli(&io___982);
	do_lio(&c__9, &c__1, "  Error reading the first line : TEXT", (ftnlen)37);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3501:
	s_wsle(&io___983);
	do_lio(&c__9, &c__1, "  Error reading lambda, etc, line 2", (ftnlen)35);
	e_wsle();
	s_wsli(&io___984);
	do_lio(&c__9, &c__1, "  Error reading lambda, etc, line 2", (ftnlen)35);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3502:
	s_wsle(&io___985);
	do_lio(&c__9, &c__1, "  Error reading symmetry code, line 3", (ftnlen)37);
	e_wsle();
	s_wsli(&io___986);
	do_lio(&c__9, &c__1, "  Error reading symmetry code, line 3", (ftnlen)37);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3503:
	s_wsle(&io___987);
	do_lio(&c__9, &c__1, "  Error reading U, V, W, step", (ftnlen)29);
	e_wsle();
	s_wsli(&io___988);
	do_lio(&c__9, &c__1, "  Error reading U, V, W, Step", (ftnlen)29);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3504:
	s_wsle(&io___989);
	do_lio(&c__9, &c__1, "  Error reading Pmin, Pmax, Vmin", (ftnlen)32);
	e_wsle();
	s_wsli(&io___990);
	do_lio(&c__9, &c__1, "  Error reading Pmin, Pmax, Vmin", (ftnlen)32);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3505:
	s_wsle(&io___991);
	do_lio(&c__9, &c__1, "  Error reading grid steps", (ftnlen)26);
	e_wsle();
	s_wsli(&io___992);
	do_lio(&c__9, &c__1, "  Error reading grid steps", (ftnlen)26);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3506:
	s_wsle(&io___993);
	do_lio(&c__9, &c__1, "  Error reading NTIMELIM", (ftnlen)24);
	e_wsle();
	s_wsli(&io___994);
	do_lio(&c__9, &c__1, "  Error reading NTIMELIM", (ftnlen)24);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
/* L3507: */
	s_wsle(&io___995);
	do_lio(&c__9, &c__1, "  Error reading nstart, rmax, test", (ftnlen)34);
	e_wsle();
	s_wsli(&io___996);
	do_lio(&c__9, &c__1, "  Error reading nstart, rmax, test", (ftnlen)34);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
	goto L3600;
L3508:
	s_wsle(&io___997);
	do_lio(&c__9, &c__1, "  Error reading data angle > 180", (ftnlen)32);
	e_wsle();
	s_wsli(&io___998);
	do_lio(&c__9, &c__1, "  Error reading data angle > 180", (ftnlen)32);
	e_wsli();
	progressview_(buffer, (ftnlen)79);
L3600:
	datn_(datenow, timenow, (ftnlen)7, (ftnlen)8);
	time(&time_end__);
	totaltime = time_end__ - time_begin__;
	s_wsle(&io___1001);
	do_lio(&c__9, &c__1, " Total CPU time elapsed in seconds : ", (ftnlen)37);
	do_lio(&c__4, &c__1, (char *)&totaltime, (ftnlen)sizeof(real));
	e_wsle();
	s_wsfe(&io___1002);
	e_wsfe();
	s_rsle(&io___1003);
	do_lio(&c__9, &c__1, fend, (ftnlen)1);
	e_rsle();
	cl__1.cerr = 0;
	cl__1.cunit = 20;
	cl__1.csta = 0;
	f_clos(&cl__1);
	cl__1.cerr = 0;
	cl__1.cunit = 28;
	cl__1.csta = 0;
	f_clos(&cl__1);
	s_stop("", (ftnlen)0);
	return 0;
} /* MAIN__ */

/* ----------------------------------------------------------------------- */
int progressview_(char *buffer, ftnlen buffer_len)
{
	/* Format strings */
	static char fmt_1[] = "(a79)";

	/* Builtin functions */
	integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

	/* Fortran I/O blocks */
	static cilist io___1005 = { 0, 6, 0, fmt_1, 0 };


	s_wsfe(&io___1005);
	do_fio(&c__1, buffer, (ftnlen)79);
	e_wsfe();
	return 0;
} /* progressview_ */

/* ----------------------------------------------------------------------- */

int datToInt(char *format, struct tm *tm_info)
{
	char temp[10];
	strftime(temp, 10, format, tm_info);    
	return (int) strtol(temp, (char **)NULL, 10);
}

/* ----------------------------------------------------------------------- */

long tz_offset(time_t t) {
	struct tm local = *localtime(&t);
	struct tm utc = *gmtime(&t);
	long diff = ((local.tm_hour - utc.tm_hour) * 60 + (local.tm_min - utc.tm_min));
	int delta_day = local.tm_mday - utc.tm_mday;
	if ((delta_day == 1) || (delta_day < -1)) {
		diff += 24L * 60;
	} else if ((delta_day == -1) || (delta_day > 1)) {
		diff -= 24L * 60;
	}
	return diff;
}

/* ----------------------------------------------------------------------- */

int datn_(char *datenow, char *timenow, ftnlen datenow_len, 
	ftnlen timenow_len)
{
	/* Initialized data */

	static char months[3*12] = "Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" 
		"Aug" "Sep" "Oct" "Nov" "Dec";

	/* Format strings */
	static char fmt_10[] = "(1x,i2,\002-\002,a3,\002-\002,i4,5x,i2,\002 ho"
		"ur \002,i2,\002 min \002,i2,\002 Sec \002)";

	/* System generated locals */
	icilist ici__1;

	/* Builtin functions */
	integer s_wsle(cilist *), e_wsle(void), s_wsfe(cilist *), do_fio(integer *
		, char *, ftnlen), e_wsfe(void);
	int s_copy(char *, char *, ftnlen, ftnlen);
	integer s_wsfi(icilist *), e_wsfi(void);

	/* Local variables */
	static integer dt_values__[8];
	static char year[2];
	// static char da[8], time[10], year[2], zone[5];
	extern int date_and_time__(char *, char *, char *, 
		integer *, ftnlen, ftnlen, ftnlen);

	/* Fortran I/O blocks */
	static cilist io___1011 = { 0, 20, 0, 0, 0 };
	static cilist io___1012 = { 0, 20, 0, fmt_10, 0 };
	static cilist io___1013 = { 0, 20, 0, 0, 0 };
	static time_t timer;


/* Obviously, this subroutine is compiler specific... */

/* In case of problem, remove that subroutine and all */
/* the CALL DATN inside the program. */

/* Date will be printed as follows : */

/*     10-Apr-2000   Hour:  0 Min: 29 Sec: 20 */

/*     dd-mmm-yyyy          hh     mm      ss */

/*     DT_VALUES(3),MONTHS(DT_VALUES(2)),DT_VALUES(1), */
/*     DT_VALUES(5),DT_VALUES(6),DT_VALUES(7) */


/* Return values from DATE_AND_TIME */
/* Month names */

/* Get current date values */

	// date_and_time__(da, time, zone, dt_values__, (ftnlen)8, (ftnlen)10, (	    ftnlen)5);

	struct tm* tm_info;
	char temp[5];

	time(&timer);
	tm_info = localtime(&timer);

	dt_values__[0] = datToInt("%Y", tm_info);
	dt_values__[1] = datToInt("%m", tm_info);
	dt_values__[2] = datToInt("%d", tm_info);
	dt_values__[3] = tz_offset(timer);
	dt_values__[4] = datToInt("%H", tm_info);
	dt_values__[5] = datToInt("%M", tm_info);
	dt_values__[6] = datToInt("%S", tm_info);

/* Format date */

	s_wsle(&io___1011);
	e_wsle();
	s_wsfe(&io___1012);
	do_fio(&c__1, (char *)&dt_values__[2], (ftnlen)sizeof(integer));
	do_fio(&c__1, months + (dt_values__[1] - 1) * 3, (ftnlen)3);
	do_fio(&c__1, (char *)&dt_values__[0], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dt_values__[4], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dt_values__[5], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&dt_values__[6], (ftnlen)sizeof(integer));
	e_wsfe();
	s_wsle(&io___1013);
	e_wsle();
	if (dt_values__[0] == 2006) {
		s_copy(year, "06", (ftnlen)2, (ftnlen)2);
	}
	if (dt_values__[0] == 2007) {
		s_copy(year, "07", (ftnlen)2, (ftnlen)2);
	}
	if (dt_values__[0] == 2008) {
		s_copy(year, "08", (ftnlen)2, (ftnlen)2);
	}
	if (dt_values__[0] == 2009) {
		s_copy(year, "09", (ftnlen)2, (ftnlen)2);
	}
	if (dt_values__[0] == 2010) {
		s_copy(year, "10", (ftnlen)2, (ftnlen)2);
	}
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 7;
	ici__1.iciunit = datenow;
	ici__1.icifmt = "(I2,A3,A2)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&dt_values__[2], (ftnlen)sizeof(integer));
	do_fio(&c__1, months + (dt_values__[1] - 1) * 3, (ftnlen)3);
	do_fio(&c__1, year, (ftnlen)2);
	e_wsfi();
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 8;
	ici__1.iciunit = timenow;
	ici__1.icifmt = "(I2,A,I2,A,I2)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&dt_values__[4], (ftnlen)sizeof(integer));
	do_fio(&c__1, ":", (ftnlen)1);
	do_fio(&c__1, (char *)&dt_values__[5], (ftnlen)sizeof(integer));
	do_fio(&c__1, ":", (ftnlen)1);
	do_fio(&c__1, (char *)&dt_values__[6], (ftnlen)sizeof(integer));
	e_wsfi();

/*  End of compiler specific subroutine */

	return 0;
} /* datn_ */

int open_read1__(integer *unit, char *file, ftnlen file_len)
{
	/* System generated locals */
	olist o__1;

	/* Builtin functions */
	integer f_open(olist *);

	o__1.oerr = 0;
	o__1.ounit = *unit;
	o__1.ofnmlen = file_len;
	o__1.ofnm = file;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	return 0;
} /* open_read1__ */


int open_write1__(integer *unit, char *file, ftnlen file_len)
{
	/* System generated locals */
	olist o__1;

	/* Builtin functions */
	integer f_open(olist *);

	o__1.oerr = 0;
	o__1.ounit = *unit;
	o__1.ofnmlen = file_len;
	o__1.ofnm = file;
	o__1.orl = 0;
	o__1.osta = "NEW";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	return 0;
} /* open_write1__ */


int filedel_(integer *unit, char *file, ftnlen file_len)
{
	/* System generated locals */
	olist o__1;
	cllist cl__1;

	/* Builtin functions */
	integer f_open(olist *), f_clos(cllist *);

	o__1.oerr = 0;
	o__1.ounit = *unit;
	o__1.ofnmlen = file_len;
	o__1.ofnm = file;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	cl__1.cerr = 0;
	cl__1.cunit = *unit;
	cl__1.csta = "delete";
	f_clos(&cl__1);
	return 0;
} /* filedel_ */


/* ----------------------------------------------------------------------- */

/*     SUBROUTINE OPEN_WRITE OPENs a FILE for writing, filling in the */
/*     supplied extension IF none is supplied with the FILE name. */
/*     Allows use of system specIFic facilities of OPEN statement. */

int open_write__(integer *unit, char *file, char *extension, 
	ftnlen file_len, ftnlen extension_len)
{
	/* System generated locals */
	address a__1[2];
	integer i__1[2], i__2;
	olist o__1;
	cllist cl__1;

	/* Builtin functions */
	int s_copy(char *, char *, ftnlen, ftnlen);
	integer i_indx(char *, char *, ftnlen, ftnlen), i_len(char *, ftnlen);
	int s_cat(char *, char **, integer *, integer *, ftnlen);
	integer f_open(olist *), f_clos(cllist *);

	/* Local variables */
	static integer i__, l;
	static char temp[80];

	s_copy(temp, file, (ftnlen)80, file_len);
	i__ = i_indx(file, ".", file_len, (ftnlen)1);
	l = i_len(file, file_len);
	while(*(unsigned char *)&file[l - 1] == ' ') {
	--l;
	}
	if (i__ == 0) {
/* Writing concatenation */
	i__1[0] = l, a__1[0] = file;
	i__1[1] = extension_len, a__1[1] = extension;
	s_cat(temp, a__1, i__1, &c__2, (ftnlen)80);
	}
	o__1.oerr = 1;
	o__1.ounit = *unit;
	o__1.ofnmlen = 80;
	o__1.ofnm = temp;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__2 = f_open(&o__1);
	if (i__2 != 0) {
	goto L10;
	}
	cl__1.cerr = 0;
	cl__1.cunit = *unit;
	cl__1.csta = "delete";
	f_clos(&cl__1);
L10:
	o__1.oerr = 0;
	o__1.ounit = *unit;
	o__1.ofnmlen = 80;
	o__1.ofnm = temp;
	o__1.orl = 0;
	o__1.osta = "NEW";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	return 0;
} /* open_write__ */


/* ----------------------------------------------------------------------- */

/*     SUBROUTINE ESP_INIT must set the seed for the random number */
/*     generator and obtain the current cpu clock reading in seconds */

int esp_init__(integer *iseed)
{
/* PORTLIB/DFPORT is compiler spcecific part, introduced for */
/* using the intrinsic function SECNDS(X) which returns the */
/*  (time in seconds since midnight - X) */

	static cilist io___1004 = { 0, 20, 0, 0, 0 };
	
	// built in functions
	integer s_wsle(cilist *); 
	integer e_wsle(void);
	integer do_lio(integer *, integer *, char *, ftnlen); 
	doublereal randi_(integer *);
	//------------------------------------------

	srand(time(NULL));
	*iseed = (rand() * 100) + 1;
 
/*  run the random number generator N times */
/*  for avoiding effects of starting value */

	int n = *iseed / 2000;
	if (n <= 0) n = 0;
	if (n >= 1000) n = 1000;
	for (int i = 0; i < n; ++i)
	{
		randi_(iseed);
	}
	*iseed = *iseed / 3;
	*iseed = (*iseed * 2) + 1;
	s_wsle(&io___1004);
	do_lio(&c__9, &c__1, " ISEED = ", (ftnlen)9);
	do_lio(&c__3, &c__1, (char *)&iseed, (ftnlen)sizeof(integer));
	e_wsle();
/*      N=ISEED/2000 */
/*      IF(N.LE.0)N=100 */
/*      IF(N.GE.1000)N=1000 */
/*      DO 10 I=1,N */
/*      BIDON=RANDI(ISEED) */
/* 10    ISEED=ISEED/3 */
/*      ISEED=ISEED*2+1 */
/*      WRITE(20,*)' ISEED = ',ISEED */
	return 0;
} /* esp_init__ */


/* ********************************************************************** */

int dcell_(real *celln, real *al, real *v)
{
	/* Builtin functions */
	double cos(doublereal);

	/* Local variables */
	static integer i__, j, k, l;
	static real cell[6];
	extern int perm_(integer *, integer *, integer *), trcl_(
		real *, real *, real *);
	static real degrad, rcelln[6];

	/* Parameter adjustments */
	al -= 4;
	--celln;

	/* Function Body */
	degrad = .017453277777777776f;
	for (i__ = 1; i__ <= 6; ++i__) {
/* L30: */
	cell[i__ - 1] = celln[i__];
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	l = i__ + 3;
	if (cell[l - 1] - 90.f != 0.f) {
		goto L32;
	} else {
		goto L33;
	}
L33:
	cell[l - 1] = 0.f;
	goto L31;
L32:
	cell[l - 1] = cos(degrad * cell[l - 1]);
L31:
	;
	}
	trcl_(cell, rcelln, v);
/*     RCELL IS THE RECIPROCAL CELL CONSTANTS */
	for (i__ = 1; i__ <= 3; ++i__) {
	al[i__ + i__ * 3] = rcelln[i__ - 1] * rcelln[i__ - 1];
	perm_(&i__, &j, &k);
	if (j - k >= 0) {
		goto L36;
	} else {
		goto L35;
	}
L35:
	al[j + k * 3] = rcelln[j - 1] * 2.f * rcelln[k - 1] * rcelln[i__ + 2];
	goto L34;
L36:
	al[k + j * 3] = rcelln[j - 1] * 2.f * rcelln[k - 1] * rcelln[i__ + 2];
L34:
	;
	}
/*      DO 37 I=4,6 */
/* 37    RCELLN(I)=ACOS(RCELLN(I))/DEGRAD */
	return 0;
} /* dcell_ */


/* ----------------------------------------------------------------------- */

int trcl_(real *celln, real *rcelln, real *v)
{
	/* System generated locals */
	real r__1;

	/* Builtin functions */
	double sqrt(doublereal);

	/* Local variables */
	static integer i__, j, k, l;
	static real abc, sina[3];
	extern int perm_(integer *, integer *, integer *);
	static real prod;

/*     TRANSFORMS REAL CELL TO RECIPROCAL OR VICE VERSA */
/*     INPUT CELL IS IN ARRAY CELL AS LENGTHS AND COSINES */
	/* Parameter adjustments */
	--rcelln;
	--celln;

	/* Function Body */
	abc = 1.f;
	prod = 2.f;
	*v = -2.f;
	for (i__ = 1; i__ <= 3; ++i__) {
	l = i__ + 3;
/* Computing 2nd power */
	r__1 = celln[l];
	sina[i__ - 1] = 1.f - r__1 * r__1;
	*v += sina[i__ - 1];
	sina[i__ - 1] = sqrt(sina[i__ - 1]);
	prod *= celln[l];
/* L10: */
	abc *= celln[i__];
	}
	*v = abc * sqrt(*v + prod);
/*      V IS CELL VOLUME */
/*     PUT INVERTED CELL INTO RCELL */
	for (i__ = 1; i__ <= 3; ++i__) {
	perm_(&i__, &j, &k);
	rcelln[i__] = celln[j] * celln[k] * sina[i__ - 1] / *v;
	l = i__ + 3;
/* L20: */
	rcelln[l] = (celln[j + 3] * celln[k + 3] - celln[l]) / (sina[j - 1] * 
		sina[k - 1]);
	}
	return 0;
} /* trcl_ */


/* ----------------------------------------------------------------------- */

int perm_(integer *i__, integer *j, integer *k)
{
	/* System generated locals */
	integer i__1;

/*     PERMS USEFUL COMBINATIONS OF INTEGERS IN THE RANGE 1 TO 3 */
	if ((i__1 = *i__ - 2) < 0) {
	goto L10;
	} else if (i__1 == 0) {
	goto L20;
	} else {
	goto L30;
	}
L10:
	*j = 2;
	*k = 3;
	return 0;
L20:
	*j = 3;
	*k = 1;
	return 0;
L30:
	*j = 1;
	*k = 2;
	return 0;
} /* perm_ */


/* ******************************************************** */

int sort_(integer *n, real *a, integer *l)
{
	/* System generated locals */
	integer i__1;

	/* Local variables */
	static integer i__, j, k, m;
	static real t;
	static integer li, lj, lk, ip, iq, lq, ix, ilt[10], itt, iut[10];

/*     ******************************************************** */
/*     THE SUBROUTINE SORT APPLIES TO THE ARRAY A WITH J ELEMENTS. */
/*     THE INDICES OF THE ORDERED ARRAY ARE PLACED IN ARRAY L */
/*     TO OBTAIN THE ORDERED ARRAY REPLACE THE INDICES I OF A(I) */
/*     WITH I=L(I)---REF.CACM 271 */
/*     ********************************************************* */

	/* Parameter adjustments */
	--l;
	--a;

	/* Function Body */
	j = *n;
	i__ = 1;
	m = 1;
	i__1 = j;
	for (k = 1; k <= i__1; ++k) {
/* L10: */
	l[k] = k;
	}
L20:
	if (j - i__ - 1 <= 0) {
	goto L140;
	} else {
	goto L30;
	}
L30:
	ip = (j + i__) / 2;
	itt = l[ip];
	t = a[itt];
	l[ip] = l[i__];
	iq = j;
	k = i__ + 1;
L40:
	if (k > iq) {
	goto L90;
	}
	lk = l[k];
	if (a[lk] <= t) {
	goto L80;
	}
	iq = iq;
L50:
	if (iq < k) {
	goto L70;
	}
	lq = l[iq];
	if (a[lq] >= t) {
	goto L60;
	}
	ix = l[k];
	l[k] = l[iq];
	l[iq] = ix;
	--iq;
	goto L80;
L60:
	--iq;
	goto L50;
L70:
	iq = k - 1;
	goto L100;
L80:
	++k;
	goto L40;
L90:
L100:
	l[i__] = l[iq];
	l[iq] = itt;
	if ((iq << 1) - i__ - j <= 0) {
	goto L120;
	} else {
	goto L110;
	}
L110:
	ilt[m - 1] = i__;
	iut[m - 1] = iq - 1;
	i__ = iq + 1;
	goto L130;
L120:
	ilt[m - 1] = iq + 1;
	iut[m - 1] = j;
	j = iq - 1;
L130:
	++m;
	goto L20;
L140:
	if (i__ - j >= 0) {
	goto L170;
	} else {
	goto L150;
	}
L150:
	li = l[i__];
	lj = l[j];
	if (a[li] - a[lj] <= 0.f) {
	goto L170;
	} else {
	goto L160;
	}
L160:
	ix = l[i__];
	l[i__] = l[j];
	l[j] = ix;
L170:
	--m;
	if (m <= 0) {
	goto L190;
	} else {
	goto L180;
	}
L180:
	i__ = ilt[m - 1];
	j = iut[m - 1];
	goto L20;
L190:
	return 0;
} /* sort_ */




/* ******************************************************** */

int sort2_(integer *n, integer *na, integer *l)
{
	/* System generated locals */
	integer i__1;

	/* Local variables */
	static integer i__, j, k, m, li, lj, lk, ip, iq, lq, ix, nt, ilt[10], itt,
		 iut[10];

/*     ******************************************************** */
/*     THE SUBROUTINE SORT APPLIES TO THE ARRAY A WITH J ELEMENTS. */
/*     THE INDICES OF THE ORDERED ARRAY ARE PLACED IN ARRAY L */
/*     TO OBTAIN THE ORDERED ARRAY REPLACE THE INDICES I OF A(I) */
/*     WITH I=L(I)---REF.CACM 271 */
/*     ********************************************************* */

	/* Parameter adjustments */
	--l;
	--na;

	/* Function Body */
	j = *n;
	i__ = 1;
	m = 1;
	i__1 = j;
	for (k = 1; k <= i__1; ++k) {
/* L10: */
	l[k] = k;
	}
L20:
	if (j - i__ - 1 <= 0) {
	goto L140;
	} else {
	goto L30;
	}
L30:
	ip = (j + i__) / 2;
	itt = l[ip];
	nt = na[itt];
	l[ip] = l[i__];
	iq = j;
	k = i__ + 1;
L40:
	if (k > iq) {
	goto L90;
	}
	lk = l[k];
	if (na[lk] <= nt) {
	goto L80;
	}
	iq = iq;
L50:
	if (iq < k) {
	goto L70;
	}
	lq = l[iq];
	if (na[lq] >= nt) {
	goto L60;
	}
	ix = l[k];
	l[k] = l[iq];
	l[iq] = ix;
	--iq;
	goto L80;
L60:
	--iq;
	goto L50;
L70:
	iq = k - 1;
	goto L100;
L80:
	++k;
	goto L40;
L90:
L100:
	l[i__] = l[iq];
	l[iq] = itt;
	if ((iq << 1) - i__ - j <= 0) {
	goto L120;
	} else {
	goto L110;
	}
L110:
	ilt[m - 1] = i__;
	iut[m - 1] = iq - 1;
	i__ = iq + 1;
	goto L130;
L120:
	ilt[m - 1] = iq + 1;
	iut[m - 1] = j;
	j = iq - 1;
L130:
	++m;
	goto L20;
L140:
	if (i__ - j >= 0) {
	goto L170;
	} else {
	goto L150;
	}
L150:
	li = l[i__];
	lj = l[j];
	if (na[li] - na[lj] <= 0) {
	goto L170;
	} else {
	goto L160;
	}
L160:
	ix = l[i__];
	l[i__] = l[j];
	l[j] = ix;
L170:
	--m;
	if (m <= 0) {
	goto L190;
	} else {
	goto L180;
	}
L180:
	i__ = ilt[m - 1];
	j = iut[m - 1];
	goto L20;
L190:
	return 0;
} /* sort2_ */




/* ******************************************************** */

int sort3_(integer *n, integer *na, integer *l)
{
	/* System generated locals */
	integer i__1;

	/* Local variables */
	static integer i__, j, k, m, li, lj, lk, ip, iq, lq, ix, nt, ilt[10], itt,
		 iut[10];

/*     ******************************************************** */
/*     THE SUBROUTINE SORT APPLIES TO THE ARRAY A WITH J ELEMENTS. */
/*     THE INDICES OF THE ORDERED ARRAY ARE PLACED IN ARRAY L */
/*     TO OBTAIN THE ORDERED ARRAY REPLACE THE INDICES I OF A(I) */
/*     WITH I=L(I)---REF.CACM 271 */
/*     ********************************************************* */

	/* Parameter adjustments */
	--l;
	--na;

	/* Function Body */
	j = *n;
	i__ = 1;
	m = 1;
	i__1 = j;
	for (k = 1; k <= i__1; ++k) {
/* L10: */
	l[k] = k;
	}
L20:
	if (j - i__ - 1 <= 0) {
	goto L140;
	} else {
	goto L30;
	}
L30:
	ip = (j + i__) / 2;
	itt = l[ip];
	nt = na[itt];
	l[ip] = l[i__];
	iq = j;
	k = i__ + 1;
L40:
	if (k > iq) {
	goto L90;
	}
	lk = l[k];
	if (na[lk] <= nt) {
	goto L80;
	}
	iq = iq;
L50:
	if (iq < k) {
	goto L70;
	}
	lq = l[iq];
	if (na[lq] >= nt) {
	goto L60;
	}
	ix = l[k];
	l[k] = l[iq];
	l[iq] = ix;
	--iq;
	goto L80;
L60:
	--iq;
	goto L50;
L70:
	iq = k - 1;
	goto L100;
L80:
	++k;
	goto L40;
L90:
L100:
	l[i__] = l[iq];
	l[iq] = itt;
	if ((iq << 1) - i__ - j <= 0) {
	goto L120;
	} else {
	goto L110;
	}
L110:
	ilt[m - 1] = i__;
	iut[m - 1] = iq - 1;
	i__ = iq + 1;
	goto L130;
L120:
	ilt[m - 1] = iq + 1;
	iut[m - 1] = j;
	j = iq - 1;
L130:
	++m;
	goto L20;
L140:
	if (i__ - j >= 0) {
	goto L170;
	} else {
	goto L150;
	}
L150:
	li = l[i__];
	lj = l[j];
	if (na[li] - na[lj] <= 0) {
	goto L170;
	} else {
	goto L160;
	}
L160:
	ix = l[i__];
	l[i__] = l[j];
	l[j] = ix;
L170:
	--m;
	if (m <= 0) {
	goto L190;
	} else {
	goto L180;
	}
L180:
	i__ = ilt[m - 1];
	j = iut[m - 1];
	goto L20;
L190:
	return 0;
} /* sort3_ */



int celref_(integer *indi, real *bbb, real *afin, integer *
	nhkl, real *theta, integer *jhkl, real *ddt, real *ddq)
{
	/* Initialized data */

	static real sig[8] = { 0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f };

	/* Format strings */
	static char fmt_10[] = "(\002 PROGRAM *** CELREF ***  (J.LAUGIER & A.FIL"
		"HOL 10/78)\002/)";
	static char fmt_60[] = "(\002 \002,/\002 OBSERVABLE NUMBER    : \002,i5"
		"/\002 ITERATION NUMBER : \002,i5/\002 REFINEMENT CONSTRAINTS :"
		" \002,a5)";
	static char fmt_70[] = "(\002 NUMBER OF INDEPENDENT PARAMETERS : \002,i5"
		"/)";
	static char fmt_80[] = "(\002 INITIAL VALUES :\002)";
	static char fmt_90[] = "(\002 FINAL VALUES   : (STANDARD DEVIATIONS : 2n"
		"d LINE)\002)";
	static char fmt_100[] = "(4x,\002ZERO\002,4x,\002LAMBDA\002,6x,\002A\002"
		",8x,\002B\002,8x,\002C\002,6x,\002ALPHA\002,5x,\002BETA\002,4x"
		",\002GAMMA\002)";
	static char fmt_110[] = "(2x,f6.3,3x,f7.4,3(2x,f7.4),3(2x,f7.3))";
	static char fmt_120[] = "(\002 RECIPROCAL CELL : \002,3(2x,f7.5),3(2x,f7"
		".3)/\002 VOLUME (A**3)  : \002,f12.3/)";
	static char fmt_130[] = "(6x,f2.0,7(7x,f2.0))";
	static char fmt_140[] = "(\002 \002,3x,\002H\002,5x,\002K\002,5x,\002"
		"L\002,2x,\002TH(OBS)\002,4x,\002TH-ZERO\002,4x,\002TH(CALC)\002,"
		"5x,\002DIFF.\002/)";
	static char fmt_150[] = "(2x,3(i3,3x),4(f7.3,4x))";
	static char fmt_170[] = "(\002 ##### ERROR REFLEXION : \002,3i4,f8.3)";
	static char fmt_180[] = "(\002 ##### DATA NUMBER GREATER THAN \002,i4"
		",\002 #####\002)";
	static char fmt_190[] = "(\002 ##### IMPOSSIBLE TO REFINE ALL PARAMETERS"
		" TOGETHER##### THINK, PLEASE ! #####\002)";
	static char fmt_366[] = "(\002 \002,3x,\002H\002,5x,\002K\002,5x,\002"
		"L\002,5x,\002D(OBS)\002,4x,\002D(CALC)\002/)";
	static char fmt_368[] = "(2x,3(i3,3x),2(f7.4,5x))";

	/* System generated locals */
	integer i__1, i__2, i__3, i__4;
	real r__1;

	/* Builtin functions */
	double acos(doublereal);
	integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
	double sqrt(doublereal);
	integer s_wsle(cilist *), e_wsle(void);
	double sin(doublereal);

	/* Local variables */
	static integer i__, j;
	static real r__, y0, y1, y2, y3, y4;
	static integer ik, jj;
	static real rd, qc, qo, rr, dum[3];
	extern int calc_();
	static char icle[3][5] = { "A=B=C", "A=B  ", "NO   "};
	static integer iffi, ifin, ihkl;
	extern int fonc_(real *, real *, real *);
	static integer nrep, ndmax;
	extern int mcrnl_(real *, integer *, real *, real *, 
		integer *, integer *, real *, integer *, U_fp), inver_(real *, 
		real *, real *, integer *);
	static real volum;
	static integer npour;

	/* Fortran I/O blocks */
	static cilist io___1080 = { 0, 0, 0, fmt_10, 0 };
	static cilist io___1086 = { 0, 0, 0, fmt_60, 0 };
	static cilist io___1088 = { 0, 0, 0, fmt_80, 0 };
	static cilist io___1089 = { 0, 0, 0, fmt_100, 0 };
	static cilist io___1090 = { 0, 0, 0, fmt_130, 0 };
	static cilist io___1091 = { 0, 0, 0, fmt_110, 0 };
	static cilist io___1094 = { 0, 0, 0, fmt_120, 0 };
	static cilist io___1096 = { 0, 0, 0, fmt_70, 0 };
	static cilist io___1100 = { 0, 0, 0, fmt_90, 0 };
	static cilist io___1101 = { 0, 0, 0, 0, 0 };
	static cilist io___1102 = { 0, 0, 0, fmt_100, 0 };
	static cilist io___1103 = { 0, 0, 0, fmt_110, 0 };
	static cilist io___1104 = { 0, 0, 0, fmt_110, 0 };
	static cilist io___1105 = { 0, 0, 0, fmt_120, 0 };
	static cilist io___1106 = { 0, 0, 0, fmt_140, 0 };
	static cilist io___1107 = { 0, 6, 0, 0, 0 };
	static cilist io___1108 = { 0, 6, 0, fmt_90, 0 };
	static cilist io___1109 = { 0, 6, 0, fmt_100, 0 };
	static cilist io___1110 = { 0, 6, 0, fmt_110, 0 };
	static cilist io___1111 = { 0, 6, 0, fmt_110, 0 };
	static cilist io___1112 = { 0, 6, 0, fmt_120, 0 };
	static cilist io___1120 = { 0, 0, 0, fmt_150, 0 };
	static cilist io___1121 = { 0, 0, 0, 0, 0 };
	static cilist io___1123 = { 0, 0, 0, fmt_366, 0 };
	static cilist io___1125 = { 0, 0, 0, fmt_368, 0 };
	static cilist io___1126 = { 0, 0, 0, fmt_170, 0 };
	static cilist io___1127 = { 0, 0, 0, fmt_180, 0 };
	static cilist io___1128 = { 0, 0, 0, fmt_190, 0 };



/* .....****************************************************************** */
/* ..... */
/* .....     PROGRAMME *** CELREF ***                         7/10/78 */
/* ..... */
/* .....     AUTEURS : JEAN LAUGIER & ALAIN FILHOL    20/10/78 */
/* ..... */
/* .....     AFFINEMENT LES PARAMETRES DE MAILLE */
/* .....                DU DECALAGE DE ZERO */
/* .....                DE LA LONGUEUR D'ONDE */
/* .....     A PARTIR DES ANGLES THETA DE BRAGG OBSERVES */
/* ..... */
/* .....     METHODE : MOINDRES CARRES NON LINEAIRES */
/* ..... */
/* .....     DONNEES : */
/* .....        1- CARTE COMMENTAIRE          FORMAT(16A5) */
/* .....        2- INDIC,IFIN                 FORMAT(2I) */
/* .....           INDIC  : CONTRAINTE D'AFFINEMENT */
/* .....                    0/1/2 POUR (A,B,C INDEPENDANTS)/(A=B=C)/(A=B) */
/* .....           IFIN   : NOMBRE MAXIMUM DE CYCLES D'AFFINEMENT */
/* ..... */
/* .....        3- B(2),B(1)                       FORMAT(2F) */
/* .....           B(2)   : LONGUEUR D'ONDE */
/* .....           B(1)   : DECALAGE SYSTEMATIQUE DE ZERO */
/* ..... */
/* .....        4- AFI(2),AFI(1)                   FORMAT(2F) */
/* .....           AFI(2) : 0/1 AFFINER (NON)/(OUI) LA LONGUEUR D'ONDE */
/* .....           AFI(1) :  "    "         "       LE DECALAGE DE ZERO */
/* ..... */
/* .....        5- B(3 A 8)                        FORMAT(6F) */
/* .....           B(3)   : PARAMETRE A */
/* .....           B(4)   : PARAMETRE B */
/* .....           B(5)   : PARAMETRE C */
/* .....           B(6)   : ANGLE ALPHA */
/* .....           B(7)   : ANGLE BETA */
/* .....           B(8)   : ANGLE GAMMA */
/* ..... */
/* .....        6- AFI(3 A 8)                      FORMAT(6F) */
/* .....           AFI(3 A 8) : 0/1 POUR AFFINER (NON)/(OUI) */
/* .....                        LES PARAMETRES DE MAILLE */
/* ..... */
/* .....        7 A 7+NR- H,K,L,THETA              FORMAT(3I,F) */
/* .....          ("NR" NOMBRE D'OBSERVATIONS) */
/* .....          (DERNIERE CARTE : H,K,L=0 0 0) */
/* .....           H,K,L  : INDICES DE MILLER */
/* .....           THETA  : ANGLE DE BRAGG OBSERVE */
/* ..... */
/* .....***************************************************************** */


/* $OMP THREADPRIVATE(/TROC/,/TRUC/) */
/*      DATA IWR/20/,ICLE/'A=B=C','A=B  ','NO   '/ */

	troc_1.iwr = 20;

	/* Parameter adjustments */
	jhkl -= 4;
	--theta;
	--afin;
	--bbb;

	/* Function Body */

/* ----- A MODIFIER EN CAS DE CHANGEMENT DES DIMENSIONS */
	ndmax = 200;
/* ----- */
	rd = 180.f / acos(-1.f);
/* ..... */
/* .....ENTREE DES DONNEES */
/* L10: */
/* L20: */
/* L30: */
/* L40: */
/* L50: */
/* L60: */
/* L70: */
/* L80: */
/* L90: */
/* L100: */
/* L110: */
/* L120: */
/* L130: */
/* L140: */
/* L150: */
/* L160: */
/* L170: */
/* L180: */
/* L190: */

/* .....ENTREE DES DONNEES */
/*      call open_read1(15,tmp) */
/* 11    call open_write1(16,tmp) */
/* L11: */
	io___1080.ciunit = troc_1.iwr;
	s_wsfe(&io___1080);
	e_wsfe();
/*      READ(IRID,20)ITITR */
/*      READ(IRID,*)INDIC,IFIN */
	truc_1.indic = *indi;
	*ddt = 0.f;
	*ddq = 0.f;
	ifin = 10;
	if (truc_1.indic == 0 || truc_1.indic > 3) {
	truc_1.indic = 3;
	}
/*      READ(IRID,*)B(2),B(1) */
	truc_1.b[0] = 0.f;
/*      READ(IRID,*)AFI(2),AFI(1) */
	for (i__ = 1; i__ <= 8; ++i__) {
		truc_1.afi[i__ - 1] = afin[i__];
	/* L5500: */
		truc_1.b[i__ - 1] = bbb[i__];
	}
/*      READ(IRID,*)(B(I),I=3,8),(AFI(I),I=3,8) */
	iffi = (integer) (truc_1.afi[2] + truc_1.afi[3] + truc_1.afi[4] + .1f);
	if (iffi == 0 || truc_1.indic == 3) {
	goto L230;
	}
	ik = 3 - truc_1.indic;
	i__1 = ik;
	for (i__ = 1; i__ <= i__1; ++i__) {
	if (truc_1.indic - 2 >= 0) {
		goto L210;
	} else {
		goto L200;
	}
L200:
	truc_1.afi[truc_1.indic + 2 + i__ - 1] = 0.f;
	goto L220;
L210:
	truc_1.afi[truc_1.indic + 1 + i__ - 1] = 0.f;
L220:
	;
	}
	truc_1.afi[2] = 1.f;

L230:
	truc_1.nr = 0;
	i__1 = *nhkl;
	for (truc_1.nr = 1; truc_1.nr <= i__1; ++truc_1.nr) {
		if (truc_1.nr > ndmax) {
			goto L380;
		}
		truc_1.h__[truc_1.nr - 1] = jhkl[truc_1.nr * 3 + 1];
		truc_1.k[truc_1.nr - 1] = jhkl[truc_1.nr * 3 + 2];
		truc_1.l[truc_1.nr - 1] = jhkl[truc_1.nr * 3 + 3];
		ihkl = (i__2 = truc_1.h__[truc_1.nr - 1], abs(i__2)) + (i__3 = 
			truc_1.k[truc_1.nr - 1], abs(i__3)) + (i__4 = truc_1.l[
			truc_1.nr - 1], abs(i__4));
		if (ihkl == 0) {
			goto L260;
		}
		if (theta[truc_1.nr] <= 0.f) {
			goto L370;
		}
	/* L250: */
		truc_1.pds[truc_1.nr - 1] = 1.f;
	/* L240: */
		theta[truc_1.nr] = theta[truc_1.nr] / rd / 2.f;
	}
/* ..... */
/* !!!! 2*THETA EN THETA */
L260:
	truc_1.nr = *nhkl;
	io___1086.ciunit = troc_1.iwr;
	s_wsfe(&io___1086);
	do_fio(&c__1, (char *)&truc_1.nr, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ifin, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&icle[truc_1.indic - 1], (ftnlen)sizeof(icle[truc_1.indic - 1]))
		;
	e_wsfe();
	io___1088.ciunit = troc_1.iwr;
	s_wsfe(&io___1088);
	e_wsfe();
	io___1089.ciunit = troc_1.iwr;
	s_wsfe(&io___1089);
	e_wsfe();
	io___1090.ciunit = troc_1.iwr;
	s_wsfe(&io___1090);
	for (i__ = 1; i__ <= 8; ++i__) {
		do_fio(&c__1, (char *)&truc_1.afi[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	io___1091.ciunit = troc_1.iwr;
	s_wsfe(&io___1091);
	for (i__ = 1; i__ <= 8; ++i__) {
		do_fio(&c__1, (char *)&truc_1.b[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	truc_1.b[0] /= rd;
	truc_1.b[1] *= .5f;
	for (i__ = 6; i__ <= 8; ++i__) {
/* L270: */
	truc_1.b[i__ - 1] /= rd;
	}
/*   ...... CALCUL DES PARAMETRES MAILLE RECIPROQUE */
	inver_(truc_1.b, dum, &volum, &c__0);
/*   ...... */
	for (i__ = 1; i__ <= 3; ++i__) {
/* L280: */
	dum[i__ - 1] = truc_1.b[i__ + 4] * rd;
	}
	io___1094.ciunit = troc_1.iwr;
	s_wsfe(&io___1094);
	for (i__ = 3; i__ <= 5; ++i__) {
	do_fio(&c__1, (char *)&truc_1.b[i__ - 1], (ftnlen)sizeof(real));
	}
	do_fio(&c__3, (char *)&dum[0], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&volum, (ftnlen)sizeof(real));
	e_wsfe();

/* ....."NPAF" : NOMBRE DE PARAMETRES A AFFINER */
/* ....."BB()" : TABLEAU DES PARAMETRES A AFFINER */
	j = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] == 0.f) {
		goto L290;
	}
	++j;
	truc_1.bb[j - 1] = truc_1.b[i__ - 1];
L290:
	;
	}
	truc_1.npaf = j;
	io___1096.ciunit = troc_1.iwr;
	s_wsfe(&io___1096);
	do_fio(&c__1, (char *)&truc_1.npaf, (ftnlen)sizeof(integer));
	e_wsfe();
	if (truc_1.npaf == 8) {
	goto L390;
	}
	truc_1.npaf2 = truc_1.npaf + 2;
/*   ......AFFINEMENT   (PDS() : POIDS (NON-UTILISE POUR CETTE VERSION)) */
	mcrnl_(truc_1.qq, &ndmax, &theta[1], truc_1.bb, &truc_1.npaf, &truc_1.nr, 
		truc_1.pds, &ifin, (U_fp)calc_);
/*   ...... */

/* .....NOUVELLES VALEURS DES PARAMETRES */
	j = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] == 0.f) {
		goto L300;
	}
	++j;
	truc_1.b[i__ - 1] = truc_1.bb[j - 1];
L300:
	;
	}
/*   ......VALEURS DES ANGLES THETA CALCULES */
	fonc_(&theta[1], &r__, &rr);
/*   ...... */

/* .....CALCUL DES ECARTS TYPE */
/* ....."SIG()" LES ECARTS TYPE DES PARAMETRES DE MAILLE QU'IL */
/* .....        CONTIENT SONT CEUX DES PARAMETRES RECIPROQUES. */
	jj = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] != 0.f) {
		goto L320;
	} else {
		goto L310;
	}
L310:
	sig[i__ - 1] = 0.f;
	goto L330;
L320:
	++jj;
	sig[i__ - 1] = sqrt(truc_1.qq[jj + jj * 200 - 201] * r__);
L330:
	;
	}
	if (truc_1.indic == 1 || truc_1.indic == 2) {
	sig[3] = sig[2];
	}
	if (truc_1.indic == 1) {
	sig[4] = sig[2];
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	truc_1.bb[i__ - 1] = truc_1.b[i__ + 1];
/* L340: */
	truc_1.bb[i__ + 2] = truc_1.b[i__ + 4] * rd;
	}

/*   ......RETOUR A LA MAILLE DIRECTE (ET ECARTS TYPE CORRESPONDANTS) */
	inver_(truc_1.b, sig, &volum, &c__1);
/*   ...... */
	sig[0] *= rd;
	sig[1] *= 2.f;

/* .....SORTIE DES RESULTATS */
	volum = 1.f / volum;

/*  Zeropoint in 2-theta to be added (same sense as TREOR, ITO, etc) */

	truc_1.b[0] = -truc_1.b[0] * rd * 2.f;
	sig[0] *= 2.f;

	truc_1.b[1] *= 2.f;
	for (i__ = 6; i__ <= 8; ++i__) {
	sig[i__ - 1] *= rd;
/* L350: */
	truc_1.b[i__ - 1] *= rd;
	}
	io___1100.ciunit = troc_1.iwr;
	s_wsfe(&io___1100);
	e_wsfe();
	io___1101.ciunit = troc_1.iwr;
	s_wsle(&io___1101);
	e_wsle();
	io___1102.ciunit = troc_1.iwr;
	s_wsfe(&io___1102);
	e_wsfe();
	io___1103.ciunit = troc_1.iwr;
	s_wsfe(&io___1103);
	for (i__ = 1; i__ <= 8; ++i__) {
	do_fio(&c__1, (char *)&truc_1.b[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	io___1104.ciunit = troc_1.iwr;
	s_wsfe(&io___1104);
	for (i__ = 1; i__ <= 8; ++i__) {
	do_fio(&c__1, (char *)&sig[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	io___1105.ciunit = troc_1.iwr;
	s_wsfe(&io___1105);
	for (i__ = 1; i__ <= 6; ++i__) {
	do_fio(&c__1, (char *)&truc_1.bb[i__ - 1], (ftnlen)sizeof(real));
	}
	do_fio(&c__1, (char *)&volum, (ftnlen)sizeof(real));
	e_wsfe();
	io___1106.ciunit = troc_1.iwr;
	s_wsfe(&io___1106);
	e_wsfe();
	s_wsle(&io___1107);
	e_wsle();
	s_wsfe(&io___1108);
	e_wsfe();
	s_wsfe(&io___1109);
	e_wsfe();
	s_wsfe(&io___1110);
	for (i__ = 1; i__ <= 8; ++i__) {
	do_fio(&c__1, (char *)&truc_1.b[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___1111);
	for (i__ = 1; i__ <= 8; ++i__) {
	do_fio(&c__1, (char *)&sig[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
	s_wsfe(&io___1112);
	for (i__ = 1; i__ <= 6; ++i__) {
	do_fio(&c__1, (char *)&truc_1.bb[i__ - 1], (ftnlen)sizeof(real));
	}
	do_fio(&c__1, (char *)&volum, (ftnlen)sizeof(real));
	e_wsfe();
	for (i__ = 1; i__ <= 8; ++i__) {
/* L351: */
	bbb[i__] = truc_1.b[i__ - 1];
	}
	npour = 0;
	i__1 = truc_1.nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	y1 = theta[i__] * rd;
	y2 = y1 + truc_1.b[0] / 2.f;
	y3 = truc_1.qq[i__ + truc_1.npaf2 * 200 - 201] * rd + truc_1.b[0] / 
		2.f;
	y4 = y2 - y3;
	y1 *= 2.f;
	y2 *= 2.f;
	y3 *= 2.f;
	y4 *= 2.f;
	if (i__ <= 20) {
		*ddt += dabs(y4);
/* Computing 2nd power */
		r__1 = sin(y2 / 2.f * 3.141593f / 180.f) * 2.f / bbb[2];
		qo = r__1 * r__1;
/* Computing 2nd power */
		r__1 = sin(y3 / 2.f * 3.141593f / 180.f) * 2.f / bbb[2];
		qc = r__1 * r__1;
		*ddq += (r__1 = qo - qc, dabs(r__1));
	}
/* L360: */
	io___1120.ciunit = troc_1.iwr;
	s_wsfe(&io___1120);
	do_fio(&c__1, (char *)&truc_1.h__[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&truc_1.k[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&truc_1.l[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&y1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y2, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y3, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y4, (ftnlen)sizeof(real));
	e_wsfe();
	}
	io___1121.ciunit = troc_1.iwr;
	s_wsle(&io___1121);
	e_wsle();
	*ddt /= 20.f;
	*ddq /= 20.f;
	r__ = sqrt(r__) * 1e3f;
/*      WRITE(IWR,160)R,RR */
	if (npour > 0) {
	goto L400;
	}
/*      print 365 */
/* L365: */
/*      READ(5,*)NREP */
	nrep = 1;
	if (nrep == 1) {
	goto L400;
	}
	io___1123.ciunit = troc_1.iwr;
	s_wsfe(&io___1123);
	e_wsfe();
	i__1 = truc_1.nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	y0 = truc_1.b[1] / 2.f;
	y1 = y0 / sin(theta[i__] - truc_1.b[0] / rd);
	y2 = y0 / sin(truc_1.qq[i__ + truc_1.npaf2 * 200 - 201] - truc_1.b[0] 
		/ rd);
/* L367: */
	io___1125.ciunit = troc_1.iwr;
	s_wsfe(&io___1125);
	do_fio(&c__1, (char *)&truc_1.h__[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&truc_1.k[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&truc_1.l[i__ - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&y1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&y2, (ftnlen)sizeof(real));
	e_wsfe();
	}
/*      print 369 */
/* L369: */
	goto L400;

/* .....MESSAGES D'ERREUR */
L370:
	io___1126.ciunit = troc_1.iwr;
	s_wsfe(&io___1126);
	do_fio(&c__1, (char *)&truc_1.h__[truc_1.nr - 1], (ftnlen)sizeof(integer))
		;
	do_fio(&c__1, (char *)&truc_1.k[truc_1.nr - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&truc_1.l[truc_1.nr - 1], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&theta[truc_1.nr], (ftnlen)sizeof(real));
	e_wsfe();
	--truc_1.nr;
/*      GOTO 240 */
L380:
	io___1127.ciunit = troc_1.iwr;
	s_wsfe(&io___1127);
	do_fio(&c__1, (char *)&ndmax, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L400;
L390:
	io___1128.ciunit = troc_1.iwr;
	s_wsfe(&io___1128);
	e_wsfe();
L400:
	return 0;
} /* celref_ */



int celref2_(integer *indi, real *bbb, real *afin, integer *
	nhkl, real *theta, integer *jhkl, real *ddt, real *ddq)
{
	/* Initialized data */

	static real sig[8] = { 0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f };

	/* System generated locals */
	integer i__1, i__2, i__3, i__4;
	real r__1;

	/* Builtin functions */
	double acos(doublereal), sqrt(doublereal), sin(doublereal);

	/* Local variables */
	static integer i__, j;
	static real r__, y1, y2, y3, y4;
	static integer ik, jj;
	static real rd, qc, qo, rr, dum[3], pip;
	extern int calc_();
	static integer iffi, ifin, ihkl;
	extern int fonc_(real *, real *, real *);
	static integer ndmax;
	extern int mcrnl_(real *, integer *, real *, real *, 
		integer *, integer *, real *, integer *, U_fp), inver_(real *, 
		real *, real *, integer *);
	static real volum;
	static integer npour;


/* $OMP THREADPRIVATE(/TROC/,/TRUC/) */
	/* Parameter adjustments */
	jhkl -= 4;
	--theta;
	--afin;
	--bbb;

	/* Function Body */
	ndmax = 200;
	rd = 180.f / acos(-1.f);
	pip = .0087266472222222221f;
	truc_1.indic = *indi;
	*ddt = 0.f;
	*ddq = 0.f;
	ifin = 10;
	if (truc_1.indic == 0 || truc_1.indic > 3) {
	truc_1.indic = 3;
	}
	truc_1.b[0] = 0.f;
	for (i__ = 1; i__ <= 8; ++i__) {
	truc_1.afi[i__ - 1] = afin[i__];
/* L5500: */
	truc_1.b[i__ - 1] = bbb[i__];
	}
	iffi = (integer) (truc_1.afi[2] + truc_1.afi[3] + truc_1.afi[4] + .1f);
	if (iffi == 0 || truc_1.indic == 3) {
	goto L230;
	}
	ik = 3 - truc_1.indic;
	i__1 = ik;
	for (i__ = 1; i__ <= i__1; ++i__) {
	if (truc_1.indic - 2 >= 0) {
		goto L210;
	} else {
		goto L200;
	}
L200:
	truc_1.afi[truc_1.indic + 2 + i__ - 1] = 0.f;
	goto L220;
L210:
	truc_1.afi[truc_1.indic + 1 + i__ - 1] = 0.f;
L220:
	;
	}
	truc_1.afi[2] = 1.f;

L230:
	truc_1.nr = 0;
	i__1 = *nhkl;
	for (truc_1.nr = 1; truc_1.nr <= i__1; ++truc_1.nr) {
	if (truc_1.nr > ndmax) {
		goto L400;
	}
	truc_1.h__[truc_1.nr - 1] = jhkl[truc_1.nr * 3 + 1];
	truc_1.k[truc_1.nr - 1] = jhkl[truc_1.nr * 3 + 2];
	truc_1.l[truc_1.nr - 1] = jhkl[truc_1.nr * 3 + 3];
	ihkl = (i__2 = truc_1.h__[truc_1.nr - 1], abs(i__2)) + (i__3 = 
		truc_1.k[truc_1.nr - 1], abs(i__3)) + (i__4 = truc_1.l[
		truc_1.nr - 1], abs(i__4));
	if (ihkl == 0) {
		goto L260;
	}
	if (theta[truc_1.nr] <= 0.f) {
		goto L400;
	}
/* L250: */
	truc_1.pds[truc_1.nr - 1] = 1.f;
/* L240: */
	theta[truc_1.nr] = theta[truc_1.nr] / rd / 2.f;
	}
/* ..... */
L260:
	truc_1.nr = *nhkl;
	truc_1.b[0] /= rd;
	truc_1.b[1] *= .5f;
	for (i__ = 6; i__ <= 8; ++i__) {
/* L270: */
	truc_1.b[i__ - 1] /= rd;
	}
	inver_(truc_1.b, dum, &volum, &c__0);
	for (i__ = 1; i__ <= 3; ++i__) {
/* L280: */
	dum[i__ - 1] = truc_1.b[i__ + 4] * rd;
	}
	j = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] == 0.f) {
		goto L290;
	}
	++j;
	truc_1.bb[j - 1] = truc_1.b[i__ - 1];
L290:
	;
	}
	truc_1.npaf = j;
	if (truc_1.npaf == 8) {
	goto L400;
	}
	truc_1.npaf2 = truc_1.npaf + 2;
	mcrnl_(truc_1.qq, &ndmax, &theta[1], truc_1.bb, &truc_1.npaf, &truc_1.nr, 
		truc_1.pds, &ifin, (U_fp)calc_);
	j = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] == 0.f) {
		goto L300;
	}
	++j;
	truc_1.b[i__ - 1] = truc_1.bb[j - 1];
L300:
	;
	}
	fonc_(&theta[1], &r__, &rr);
	jj = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] != 0.f) {
		goto L320;
	} else {
		goto L310;
	}
L310:
	sig[i__ - 1] = 0.f;
	goto L330;
L320:
	++jj;
	sig[i__ - 1] = sqrt(truc_1.qq[jj + jj * 200 - 201] * r__);
L330:
	;
	}
	if (truc_1.indic == 1 || truc_1.indic == 2) {
	sig[3] = sig[2];
	}
	if (truc_1.indic == 1) {
	sig[4] = sig[2];
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	truc_1.bb[i__ - 1] = truc_1.b[i__ + 1];
/* L340: */
	truc_1.bb[i__ + 2] = truc_1.b[i__ + 4] * rd;
	}
	sig[0] *= rd;
	sig[1] *= 2.f;
	truc_1.b[0] = -truc_1.b[0] * rd * 2.f;
	sig[0] *= 2.f;
	truc_1.b[1] *= 2.f;
	for (i__ = 6; i__ <= 8; ++i__) {
	sig[i__ - 1] *= rd;
/* L350: */
	truc_1.b[i__ - 1] *= rd;
	}
	for (i__ = 1; i__ <= 8; ++i__) {
/* L351: */
	bbb[i__] = truc_1.b[i__ - 1];
	}
	npour = 0;
	i__1 = truc_1.nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	y1 = theta[i__] * rd;
	y2 = y1 + truc_1.b[0] / 2.f;
	y3 = truc_1.qq[i__ + truc_1.npaf2 * 200 - 201] * rd + truc_1.b[0] / 
		2.f;
	y4 = y2 - y3;
	y1 *= 2.f;
	y2 *= 2.f;
	y3 *= 2.f;
	y4 *= 2.f;
	if (i__ <= 20) {
		*ddt += dabs(y4);
/* Computing 2nd power */
		r__1 = sin(y2 * pip) * 2.f / bbb[2];
		qo = r__1 * r__1;
/* Computing 2nd power */
		r__1 = sin(y3 * pip) * 2.f / bbb[2];
		qc = r__1 * r__1;
		*ddq += (r__1 = qo - qc, dabs(r__1));
	}
/* L360: */
	}
	*ddt /= 20.f;
	*ddq /= 20.f;
L400:
	return 0;
} /* celref2_ */



int calc_(void)
{
	/* System generated locals */
	integer i__1;
	real r__1, r__2, r__3;

	/* Builtin functions */
	double cos(doublereal), sqrt(doublereal), sin(doublereal), asin(
		doublereal);

	/* Local variables */
	static real d__, f;
	static integer i__, j;
	static real q[2000]	/* was [200][10] */, ae, be, ce, dd;
	static integer ir;
	static real cae, cbe, cce, rad;

/* .....-------------- */
/* $OMP THREADPRIVATE(/TROC/,/TRUC/) */
	j = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
	if (truc_1.afi[i__ - 1] == 0.f) {
		goto L1;
	}
	++j;
	truc_1.b[i__ - 1] = truc_1.bb[j - 1];
L1:
	;
	}
	if (truc_1.indic == 1 || truc_1.indic == 2) {
	truc_1.b[3] = truc_1.b[2];
	}
	if (truc_1.indic == 1) {
	truc_1.b[4] = truc_1.b[2];
	}
	ae = truc_1.b[2];
	be = truc_1.b[3];
	ce = truc_1.b[4];
	cae = cos(truc_1.b[5]);
	cbe = cos(truc_1.b[6]);
	cce = cos(truc_1.b[7]);
	i__1 = truc_1.nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = ae * truc_1.h__[i__ - 1];
/* Computing 2nd power */
	r__2 = be * truc_1.k[i__ - 1];
/* Computing 2nd power */
	r__3 = ce * truc_1.l[i__ - 1];
	dd = r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + (truc_1.h__[i__ - 1] * 
		truc_1.k[i__ - 1] * ae * be * cce + truc_1.k[i__ - 1] * 
		truc_1.l[i__ - 1] * be * ce * cae + truc_1.l[i__ - 1] * 
		truc_1.h__[i__ - 1] * ce * ae * cbe) * 2.f;
	d__ = 1.f / sqrt(dd);
/* Computing 2nd power */
	r__1 = truc_1.b[1];
	rad = sqrt(1.f - r__1 * r__1 * dd);
	f = truc_1.b[1] * d__ / rad;
	q[i__ - 1] = 1.f;
	q[i__ + 199] = 1.f / (d__ * rad);
	q[i__ + 399] = f * truc_1.h__[i__ - 1] * (truc_1.h__[i__ - 1] * ae + 
		truc_1.k[i__ - 1] * be * cce + truc_1.l[i__ - 1] * ce * cbe);
	q[i__ + 599] = f * truc_1.k[i__ - 1] * (truc_1.k[i__ - 1] * be + 
		truc_1.l[i__ - 1] * ce * cae + truc_1.h__[i__ - 1] * ae * cce)
		;
	q[i__ + 799] = f * truc_1.l[i__ - 1] * (truc_1.l[i__ - 1] * ce + 
		truc_1.h__[i__ - 1] * ae * cbe + truc_1.k[i__ - 1] * be * cae)
		;
	q[i__ + 999] = -f * truc_1.k[i__ - 1] * truc_1.l[i__ - 1] * be * ce * 
		sin(truc_1.b[5]);
	q[i__ + 1199] = -f * truc_1.l[i__ - 1] * truc_1.h__[i__ - 1] * ce * 
		ae * sin(truc_1.b[6]);
	q[i__ + 1399] = -f * truc_1.h__[i__ - 1] * truc_1.k[i__ - 1] * ae * 
		be * sin(truc_1.b[7]);
/* L2: */
	truc_1.qq[i__ + truc_1.npaf2 * 200 - 201] = truc_1.b[0] + asin(
		truc_1.b[1] / d__);
	}
	i__1 = truc_1.nr;
	for (ir = 1; ir <= i__1; ++ir) {
	j = 0;
	for (i__ = 1; i__ <= 8; ++i__) {
		if (truc_1.afi[i__ - 1] == 0.f) {
		goto L3;
		}
		++j;
		truc_1.qq[ir + j * 200 - 201] = q[ir + i__ * 200 - 201];
L3:
		;
	}
	}
/*      WRITE(IWR,5)(BB(I),I=1,NPAF) */
/* L5: */
	return 0;
} /* calc_ */



int fonc_(real *theta, real *r__, real *rr)
{
	/* System generated locals */
	integer i__1;
	real r__1, r__2, r__3;

	/* Builtin functions */
	double cos(doublereal), sqrt(doublereal), asin(doublereal);

	/* Local variables */
	static real d__;
	static integer i__;
	static real r1, r2, ae, be, ce, dd, yc, cae, cbe, cce;

/* .....----------------------- */
/* $OMP THREADPRIVATE(/TRUC/) */
	/* Parameter adjustments */
	--theta;

	/* Function Body */
	r1 = 0.f;
	r2 = 0.f;
	ae = truc_1.b[2];
	be = truc_1.b[3];
	ce = truc_1.b[4];
	cae = cos(truc_1.b[5]);
	cbe = cos(truc_1.b[6]);
	cce = cos(truc_1.b[7]);
	i__1 = truc_1.nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	r__1 = ae * truc_1.h__[i__ - 1];
/* Computing 2nd power */
	r__2 = be * truc_1.k[i__ - 1];
/* Computing 2nd power */
	r__3 = ce * truc_1.l[i__ - 1];
	dd = r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + (truc_1.h__[i__ - 1] * 
		truc_1.k[i__ - 1] * ae * be * cce + truc_1.k[i__ - 1] * 
		truc_1.l[i__ - 1] * be * ce * cae + truc_1.l[i__ - 1] * 
		truc_1.h__[i__ - 1] * ce * ae * cbe) * 2.f;
	d__ = 1.f / sqrt(dd);
	yc = truc_1.b[0] + asin(truc_1.b[1] / d__);
/* Computing 2nd power */
	r__1 = yc - theta[i__];
	r1 += r__1 * r__1;
/* Computing 2nd power */
	r__1 = theta[i__];
	r2 += r__1 * r__1;
/* L2: */
	truc_1.qq[i__ + truc_1.npaf2 * 200 - 201] = yc;
	}
	*r__ = r1 / (truc_1.nr - truc_1.npaf);
	*rr = sqrt(r1 / r2);
	return 0;
} /* fonc_ */



int mcrnl_(real *q, integer *id, real *y, real *b, integer *
	m, integer *n, real *p, integer *ifin, S_fp calc)
{
	/* System generated locals */
	integer q_dim1, q_offset, i__1, i__2, i__3;

	/* Local variables */
	static real a[60];
	static integer i__, j, l;
	static real r__;
	static integer m1, m2, il, mm, iq, no, ier, imm;
	extern int matinv_(real *, integer *, integer *);

/* .....------------------------------------------ */

/* .....PROGRAMME DE MOINDRES CARRES NON LINEAIRE */
/* .....UTILISANT L'INVERSION DE MATRICE TABLEAU A UNE DIMENSION */
/* .....LIMITE A 15 PARAMETRES */
/* ..... M = NBRE DE PARAMETRES */
/* ..... N = NBRE DE DONNEES */
/* ..... ID = DIMENSION DU TABLEAU Q(ID,M+1) */
/* ..... A((M*(M-1)/2 + M) = ZONE DE TRAVAIL */

	/* Parameter adjustments */
	q_dim1 = *id;
	q_offset = 1 + q_dim1;
	q -= q_offset;
	--y;
	--b;
	--p;

	/* Function Body */
	mm = *m << 1;
	m1 = *m + 1;
	m2 = *m + 2;
L10:
	--(*ifin);

/* -----CALCUL DES DERIVEES PARTIELLES ET DES Y */
	(*calc)();

/* -----LES Y CALCULES SONT DANS LA COLONNE M+2 */
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	r__ = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
		r__ += (y[i__] - q[i__ + m2 * q_dim1]) * q[i__ + j * q_dim1] * p[
			i__];
	}
/* L30: */
	q[j + m1 * q_dim1] = r__;
	}

/* -----CONSTRUCTION DE LA MATRICE SYMETRIQUE A=Q*QT */
	no = 0;
	i__1 = *m;
	for (iq = 1; iq <= i__1; ++iq) {
	i__2 = *m;
	for (il = iq; il <= i__2; ++il) {
		++no;
		r__ = 0.f;
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L40: */
		r__ += q[i__ + iq * q_dim1] * q[i__ + il * q_dim1] * p[i__];
		}
/* L50: */
		a[no - 1] = r__;
	}
	}

/* -----INVERSION DE LA MATRICE A */
	matinv_(a, m, &ier);
	if (*ifin <= 0) {
	goto L110;
	} else {
	goto L60;
	}

/* -----CALCUL DES NOUVEAUX PARAMETRES */
L60:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	r__ = 0.f;
	imm = (i__ - 1) * (mm - i__) / 2;
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
		if (j - i__ >= 0) {
		goto L70;
		} else {
		goto L80;
		}
L70:
		l = imm + j;
		goto L90;
L80:
		l = (j - 1) * (mm - j) / 2 + i__;
L90:
		r__ += a[l - 1] * q[j + m1 * q_dim1];
	}
/* L100: */
	b[i__] += r__;
	}
	goto L10;

/* -----REMISE DES ELEMENTS DIAGONAUX DANS Q(I,I),APRES LE DERNIER CYCLE */
L110:
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	l = (i__ - 1) * (mm - i__) / 2 + i__;
/* L120: */
	q[i__ + i__ * q_dim1] = a[l - 1];
	}
	return 0;
} /* mcrnl_ */



int matinv_(real *am, integer *n, integer *nfail)
{
	/* System generated locals */
	integer i__1, i__2, i__3;

	/* Builtin functions */
	double sqrt(doublereal);

	/* Local variables */
	static integer i__, j, k, l, m, ii, kdm, kli, kmi, imax;
	static doublereal suma;
	static real term, denom;

/* .....----------------------------- */
/*     ********** SEGMENT 1 OF CHOLESKI INVERSION ********** */
/*     ***** FACTOR MATRIX INTO LOWER TRIANGLE X TRANSPOSE ***** */
	/* Parameter adjustments */
	--am;

	/* Function Body */
	k = 1;
	if ((i__1 = *n - 1) < 0) {
	goto L8;
	} else if (i__1 == 0) {
	goto L10;
	} else {
	goto L20;
	}
L8:
	*nfail = k;
	goto L210;
L10:
	am[1] = 1.f / am[1];
	goto L200;
/*     ***** LOOP M OF A(L,M) ***** */
L20:
	i__1 = *n;
	for (m = 1; m <= i__1; ++m) {
	imax = m - 1;
/*     ***** LOOP L OF A(L,M) ***** */
	i__2 = *n;
	for (l = m; l <= i__2; ++l) {
		suma = 0.f;
		kli = l;
		kmi = m;
		if (imax <= 0) {
		goto L50;
		} else {
		goto L30;
		}
/*     *****SUM OVER I=1,M-1 A(L,I)*A(M,I) ***** */
L30:
		i__3 = imax;
		for (i__ = 1; i__ <= i__3; ++i__) {
		suma += am[kli] * am[kmi];
		j = *n - i__;
		kli += j;
/* L40: */
		kmi += j;
		}
/*     *****TERM=C(L,M)-SUM ***** */
L50:
		term = am[k] - suma;
		if (l - m <= 0) {
		goto L60;
		} else {
		goto L90;
		}
L60:
		if (term <= 0.f) {
		goto L80;
		} else {
		goto L70;
		}
/*     ***** A(M,M)=SQRT(TERM) ***** */
L70:
		denom = sqrt(term);
		am[k] = denom;
		goto L100;
L80:
		*nfail = k;
		goto L210;
/*     ***** A(L,M)=TERM/A(M,M) ***** */
L90:
		am[k] = term / denom;
L100:
		++k;
	}
/* L110: */
	}
/*     ********** SEGMENT 2 OF CHOLESKI INVERSION ********** */
/*     *****INVERSION OF TRIANGULAR MATRIX***** */
/* L120: */
	am[1] = 1.f / am[1];
	kdm = 1;
/*     ***** STEP L OF B(L,M) ***** */
	i__1 = *n;
	for (l = 2; l <= i__1; ++l) {
	kdm = kdm + *n - l + 2;
/*     ***** RECIPROCAL OF DIAGONAL TERM ***** */
	term = 1.f / am[kdm];
	am[kdm] = term;
	kmi = 0;
	kli = l;
	imax = l - 1;
/*     ***** STEP M OF B(L,M) ***** */
	i__2 = imax;
	for (m = 1; m <= i__2; ++m) {
		k = kli;
/*     ***** SUM TERMS ***** */
		suma = 0.f;
		i__3 = imax;
		for (i__ = m; i__ <= i__3; ++i__) {
		ii = kmi + i__;
		suma -= am[kli] * am[ii];
/* L130: */
		kli = kli + *n - i__;
		}
/*     ***** MULT SUM * RECIP OF DIAGONAL ***** */
		am[k] = suma * term;
		j = *n - m;
		kli = k + j;
/* L140: */
		kmi += j;
	}
/* L150: */
	}
/*     ********** SEGMENT 3 OF CHOLESKI INVERSION ********** */
/*     *****PREMULTIPLY LOWER TRIANGLE BY TRANSPOSE***** */
/* L160: */
	k = 1;
	i__1 = *n;
	for (m = 1; m <= i__1; ++m) {
	kli = k;
	i__2 = *n;
	for (l = m; l <= i__2; ++l) {
		kmi = k;
		imax = *n - l + 1;
		suma = 0.f;
		i__3 = imax;
		for (i__ = 1; i__ <= i__3; ++i__) {
		suma += am[kli] * am[kmi];
		++kli;
/* L170: */
		++kmi;
		}
		am[k] = suma;
/* L180: */
		++k;
	}
/* L190: */
	}
L200:
	*nfail = 0;
L210:
	return 0;
} /* matinv_ */



int inver_(real *b, real *db, real *volum, integer *iv)
{
	/* Builtin functions */
	double cos(doublereal), sin(doublereal), sqrt(doublereal), acos(
		doublereal);

	/* Local variables */
	static real d__[6];
	static integer i__, j, k;
	static real q, r__, q2, ad[6], cc[3];
	static integer jj;
	static real sp[3], ss[3], dqd[3], sig[6], cosp[3], sinp[3], cabc2;
	extern int sigma_(real *, real *, real *);

/* ..... */
/* .....  1-CALCUL LES PARAMETRES MAILLE INVERSE */
/* .....  2-CALCULE LES ECARTS TYPE DES PARAMETRES */
/* .....   IV=0 OPTION 1  IV=1 OPTIONS 1 & 2 */
/* ..... */

	/* Parameter adjustments */
	--db;
	--b;

	/* Function Body */
	cabc2 = 0.f;
	for (i__ = 1; i__ <= 3; ++i__) {
	ad[i__ - 1] = b[i__ + 2];
	cosp[i__ - 1] = cos(b[i__ + 5]);
	sinp[i__ - 1] = sin(b[i__ + 5]);
/* L10: */
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	j = i__ % 3 + 1;
	k = (i__ + 1) % 3 + 1;
	dqd[i__ - 1] = cosp[i__ - 1] - cosp[j - 1] * cosp[k - 1];
	ss[i__ - 1] = sinp[j - 1] * sinp[k - 1];
	cc[i__ - 1] = -dqd[i__ - 1] / ss[i__ - 1];
	cabc2 += cosp[i__ - 1] * cosp[i__ - 1];
/* L15: */
	}

	q2 = 1.f - cabc2 + cosp[0] * 2.f * cosp[1] * cosp[2];
	q = sqrt(q2);
	*volum = ad[0] * ad[1] * ad[2] * q;

	for (i__ = 1; i__ <= 3; ++i__) {
	b[i__ + 2] = sinp[i__ - 1] / (ad[i__ - 1] * q);
	b[i__ + 5] = acos(cc[i__ - 1]);
	sp[i__ - 1] = sin(b[i__ + 5]);
/* L20: */
	}

	if (*iv == 0) {
	goto L70;
	}

/* .....DERIVEES DES PARAMETRES A , B , C */
	for (i__ = 1; i__ <= 3; ++i__) {
	j = i__ % 3 + 1;
	k = (i__ + 1) % 3 + 1;
	d__[i__ - 1] = -sinp[i__ - 1] / ad[i__ - 1];
	d__[j - 1] = 0.f;
	d__[k - 1] = 0.f;
	d__[i__ + 2] = cosp[i__ - 1] - sinp[i__ - 1] * sinp[i__ - 1] * dqd[
		i__ - 1] / (q * q);
	d__[j + 2] = -ss[k - 1] * dqd[j - 1] / q2;
	d__[k + 2] = -ss[j - 1] * dqd[k - 1] / q2;
	sigma_(d__, &db[1], &r__);
/* L30: */
	sig[i__ - 1] = sqrt(r__) / (q * ad[i__ - 1]);
	}

/* .....DERIVEES DES ANGLES DE LA MAILLE */
	for (i__ = 1; i__ <= 3; ++i__) {
	for (jj = 1; jj <= 3; ++jj) {
/* L40: */
		d__[jj - 1] = 0.f;
	}
	j = i__ % 3 + 1;
	k = (i__ + 1) % 3 + 1;
	d__[i__ + 2] = sinp[i__ - 1] / ss[i__ - 1];
	d__[j + 2] = cosp[k - 1] / sinp[k - 1] + cosp[j - 1] * cc[i__ - 1] / 
		sinp[j - 1];
	d__[k + 2] = cosp[j - 1] / sinp[j - 1] + cosp[k - 1] * cc[i__ - 1] / 
		sinp[k - 1];
	sigma_(d__, &db[1], &r__);
/* L50: */
	sig[i__ + 2] = sqrt(r__) / sp[i__ - 1];
	}

	for (i__ = 1; i__ <= 6; ++i__) {
/* L60: */
	db[i__ + 2] = sig[i__ - 1];
	}

L70:
	return 0;
} /* inver_ */



int sigma_(real *d__, real *db, real *r__)
{
	/* System generated locals */
	real r__1;

	/* Local variables */
	static integer i__;

/* ..... */
	/* Parameter adjustments */
	--db;
	--d__;

	/* Function Body */
	*r__ = 0.f;
	for (i__ = 1; i__ <= 6; ++i__) {
/* L10: */
/* Computing 2nd power */
	r__1 = d__[i__] * db[i__ + 2];
	*r__ += r__1 * r__1;
	}
	return 0;
} /* sigma_ */


/* ... Here is the program heart... */

int calcul1_(real *diff, real *diff2)
{
	/* System generated locals */
	integer i__1, i__2;
	real r__1;

	/* Builtin functions */
	double sqrt(doublereal), asin(doublereal);

	/* Local variables */
	static integer i__, j, k, l;
	static real x, de[10000];
	static integer jh;
	static real cr[10000], perc[10000], diff1, demax, theta[10000], sinth, 
		sum_f2__;





/* $OMP THREADPRIVATE(/cal/) */

/* ...  Keep only the hkl for d > dmin */

	jh = 0;
	i__1 = cal_2.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	x = 0.f;
	for (j = 1; j <= 3; ++j) {
		for (k = j; k <= 3; ++k) {
/* L107: */
		x = cal_2.al[j + k * 3 - 4] * cal_2.ihh[j + i__ * 3 - 4] * 
			cal_2.ihh[k + i__ * 3 - 4] + x;
		}
	}
	if (x > cal_2.dmin__) {
		goto L109;
	}
	++jh;
/*     X IS 1/D(hkl)**2 FOR REFLECTION IHH */

/*     This should be optimized for speed : */
/*     working only on X, not calculating 2-theta... */
/*     - in fact, tests show that only 10-15% is gained - */

	sinth = cal_2.slabda2 * x;
	sinth = sqrt(sinth);
	theta[jh - 1] = asin(sinth) * cal_2.pi;
	cr[jh - 1] = 0.f;
L109:
	;
	}
	cal_2.nhkl = jh;
	if (cal_2.nhkl > cal_2.ndat10) {
	return 0;
	}

/* ...  Comparison with the data */

	cal_2.lhkl = 0;
	i__1 = cal_2.ndat;
	for (j = 1; j <= i__1; ++j) {
	cal_2.cri[j - 1] = 0.f;
	perc[j - 1] = 0.f;
	demax = cal_2.w2[j - 1];
	i__2 = cal_2.nhkl;
	for (k = 1; k <= i__2; ++k) {
		if (cr[k - 1] == 1.f) {
		goto L111;
		}
		if (theta[k - 1] <= cal_2.difp[j - 1] && theta[k - 1] >= 
			cal_2.difm[j - 1]) {
		de[k - 1] = (r__1 = theta[k - 1] - cal_2.th2[j - 1], dabs(
			r__1));
		if (de[k - 1] <= demax) {
			l = k;
			demax = de[k - 1];
			cal_2.cri[j - 1] = 1.f;
		}
		}
L111:
		;
	}

/*  PERC = percentage of columnar overlap for that peak */

/*  Potential problem here because only one reflection */
/*  overlapping the most closely with the column is */
/*  included (if CRI =1)... */

	if (cal_2.cri[j - 1] == 1.f) {
		perc[j - 1] = 1.f - de[l - 1] / cal_2.w2[j - 1];
		++cal_2.lhkl;
		cr[l - 1] = 1.f;
	}
/* L113: */
	}

/* ...  Calculate "R" */

	diff1 = 0.f;
	sum_f2__ = 0.f;
	i__1 = cal_2.ndat;
	for (k = 1; k <= i__1; ++k) {
	sum_f2__ += cal_2.fobs[k - 1] * cal_2.cri[k - 1];
/* L1122: */
	diff1 += cal_2.cri[k - 1] * cal_2.fobs[k - 1] * perc[k - 1];
	}
	*diff = 1.f - diff1 / cal_2.sum_f__;
	*diff2 = 1.f - diff1 / sum_f2__;
	return 0;
} /* calcul1_ */


int calcul2_(real *diff, integer *ihkl, real *th3, integer *
	ncalc, integer *igc)
{
	/* System generated locals */
	integer i__1, i__2;
	real r__1;

	/* Builtin functions */
	double sqrt(doublereal), asin(doublereal);

	/* Local variables */
	static integer i__, j, k, l;
	static real x;
	static integer l2;
	static real de[10000];
	static integer jh;
	static real qc[10000], cr[10000];
	static integer jj, jjj, jhkl[30000]	/* was [3][10000] */;
	static real perc[10000], demax, theta[10000], sinth, sum_f2__;





/* $OMP THREADPRIVATE(/cal/,/cal2/) */

/* ...  Keep only the hkl for d > dmin */

	/* Parameter adjustments */
	--th3;
	ihkl -= 4;

	/* Function Body */
	jh = 0;
	i__1 = cal_2.nhkl0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	x = 0.f;
	for (j = 1; j <= 3; ++j) {
		for (k = j; k <= 3; ++k) {
/* L107: */
		x = cal_2.al[j + k * 3 - 4] * cal_2.ihh[j + i__ * 3 - 4] * 
			cal_2.ihh[k + i__ * 3 - 4] + x;
		}
	}
	if (x > cal_2.dmin__) {
		goto L109;
	}
	++jh;
	for (j = 1; j <= 3; ++j) {
/* L108: */
		jhkl[j + jh * 3 - 4] = cal_2.ihh[j + i__ * 3 - 4];
	}
/*     X IS 1/D(hkl)**2 FOR REFLECTION IHH */

/*     This should be optimized for speed : */
/*     working only on X, not calculating 2-theta... */

	qc[jh - 1] = x;
	sinth = cal_2.slabda2 * x;
	sinth = sqrt(sinth);
	theta[jh - 1] = asin(sinth) * cal_2.pi;
	cr[jh - 1] = 0.f;
L109:
	;
	}
	cal_2.nhkl = jh;

/* ...  Comparison with the data */

	cal_2.lhkl = 0;
	*ncalc = 0;
	jj = 0;
	i__1 = cal_2.ndat;
	for (j = 1; j <= i__1; ++j) {
	cal_2.cri[j - 1] = 0.f;
	perc[j - 1] = 0.f;
/* CC */
/* CC  Eliminating too spurious peaks here ??? */
/* CC    tolerance on width decreased by a factor 3 */
/* CC */
	demax = cal_2.w2[j - 1] * .33333f;
/* CC */
	i__2 = cal_2.nhkl;
	for (k = 1; k <= i__2; ++k) {
		if (cr[k - 1] == 1.f) {
		goto L111;
		}
		if (theta[k - 1] <= cal_2.difp[j - 1] && theta[k - 1] >= 
			cal_2.difm[j - 1]) {
		de[k - 1] = (r__1 = theta[k - 1] - cal_2.th2[j - 1], dabs(
			r__1));
		if (de[k - 1] <= demax) {
			l = k;
			demax = de[k - 1];
			cal_2.cri[j - 1] = 1.f;
		}
		}
L111:
		;
	}

/*  PERC = percentage of columnar overlap for that peak */

/*  Potential problem here because only one reflection */
/*  overlapping the most closely with the column is */
/*  included (if CRI =1)... */

	if (cal_2.cri[j - 1] == 1.f) {
		perc[j - 1] = 1.f - de[l - 1] / (cal_2.w2[j - 1] * .33333f);
		++cal_2.lhkl;
		cr[l - 1] = 1.f;
		for (l2 = 1; l2 <= 3; ++l2) {
/* L112: */
		ihkl[l2 + cal_2.lhkl * 3] = jhkl[l2 + l * 3 - 4];
		}
		th3[cal_2.lhkl] = cal_2.th2[j - 1];
		jj += cal_2.cri[j - 1];
		if (jj <= 20) {
		jjj = j;
		}
	}
/* L113: */
	}
/*  NCALC is for M(20) FoM */
	i__1 = cal_2.nhkl;
	for (k = 1; k <= i__1; ++k) {
	if (cal_2.th2[jjj - 1] >= theta[k - 1]) {
		++(*ncalc);
	}
/* L114: */
	}

/* ...  Calculate "R" */

	*diff = 0.f;

/*  Change here with SUM_F2 being only on explained reflections... */

	sum_f2__ = 0.f;
	i__1 = cal_2.ndat;
	for (k = 1; k <= i__1; ++k) {
	cal2_1.ind[k + *igc * 100 - 101] = cal_2.cri[k - 1];
	sum_f2__ += cal_2.fobs[k - 1] * cal_2.cri[k - 1];
/* L1122: */
	*diff += cal_2.cri[k - 1] * cal_2.fobs[k - 1] * perc[k - 1];
	}
	*diff = 1.f - *diff / sum_f2__;
	return 0;
} /* calcul2_ */


int peekchar()
{
	int c;

	c = getchar();
	if(c != EOF) ungetc(c, stdin);      /* puts it back */
	
	return c;
}

int killk_(logical *pressedk)
{
	int pressed = peekchar();
	if (pressed)
	{
		char key = getchar();
		if (key == 'K')
		{
			*pressedk = 1;
		}
	}
	return 0;
/*    Checks if the 'K' keystroke has been pressed */

/*      USE DFLIB */
/*      LOGICAL PRESSED */
/*      LOGICAL PRESSEDK */
/*      CHARACTER(1) KEY */
/*      PRESSED = PEEKCHARQQ ( ) */
/* 	PRESSED = /.TRUE./ */
/*     IF(PRESSED)THEN */
/*      KEY = GETCHARQQ() */
/*      IF(KEY.EQ.'K') */
/* 	PRESSEDK=.TRUE. */
/*      ENDIF */
} /* killk_ */

doublereal randi_(integer *ix)
{
	/* System generated locals */
	real ret_val;

	/* Local variables */
	static integer k;
	static doublereal x;
	static integer xhi, xlo, xahi, xalo, leftlo;


/*     A random number generator using the recursion IX=IX*A MOD P */
/*     where A=7**5 and P=2**31-1. The value returned is in the */
/*     range 0.<= ran <1. */

/*     In this form fairly portable as does not require knowledge */
/*     of data storage, but does not adhere to the standard in two respects: */
/*     1) Assumes an integer word length of at least 32 bits */
/*     2) Assumes that a positive integer less than 2**16 may be */
/*        floated without loss of digits. */

/*     This code is based on code published by Linus Schrage in */
/*     T.O.M.S. vol.5 no.2 June 1979 (pp 132-138) */

/*     The method employed is a multiplicative congruential one using a */
/*     multiplier of 7**5 and taking the modulo to 2**31-1, i.e. the */
/*     generator number, x, is updated on each call to the value */
/*     x*7**5  modulo (2**31-1). The result returned is calculated as a */
/*     real number having the value x/(2**31-1) */

/*     a=7**5, b15=2**15, b16=2**16, p= 2**31-1 */

/*     Get the 15 high order bits and 16 low order bits of ix */

	xhi = *ix / 65536;
	xlo = *ix - (xhi << 16);

/*     Multiply low order part */

	xalo = xlo * 16807;

/*     Get high order part of product to carry */

	leftlo = xalo / 65536;

/*     And so obtain high order part of product */

	xahi = xhi * 16807 + leftlo;

/*     Obtain 32nd bit (overflow) of full product */

	k = xahi / 32768;

/*     Put all the bits of the product together and subtract p */
/*     (Must be in this form to prevent overflow. The ()'s are essential) */

	*ix = xalo - (leftlo << 16) - 2147483647 + (xahi - (k << 15) << 16) + k;

/*     If <0 add p back again */

	if (*ix < 0) {
	*ix += 2147483647;
	}

/*     Finally multiply by 1/(2**31-1) to obtain number in range 0-1 */

	xhi = *ix / 65536;
	x = (doublereal) xhi * 65536. + (doublereal) (*ix - (xhi << 16));
	x *= 4.6566128752457969e-10;
	ret_val = (real) x;
	return ret_val;
} /* randi_ */


/* ... Check for supercell, and reduce to minimal cell */

int supcel_(integer *n, integer *ihkl, real *cel, integer *l,
	 real *vgc, integer *js)
{
	/* System generated locals */
	integer i__1, i__2;

	/* Local variables */
	static integer i__, j, k, ii;
	static real xx;
	static integer id2[3], id3[3], id4[3], id5[3], id6[3], ihh, iab2, iab3, 
		iab4, iab5, iab6, ihmax[3];


/*  Sum on the h, k, l */
	/* Parameter adjustments */
	--vgc;
	cel -= 7;
	ihkl -= 4;

	/* Function Body */
	for (j = 1; j <= 3; ++j) {
/*  At the end, if the following codes are still = 1 */
/*     then there could be a common divider */
	id2[j - 1] = 1;
	id3[j - 1] = 1;
	id4[j - 1] = 1;
	id5[j - 1] = 1;
	id6[j - 1] = 1;
	ihmax[j - 1] = 0;
/*  Sum on the hkl data */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*  Test on dividing h, k, l by 2, 3, 4, 5, 6 */
		k = (i__2 = ihkl[j + i__ * 3], abs(i__2));
		if (k > ihmax[j - 1]) {
		ihmax[j - 1] = k;
		}
		if (k == 0) {
		goto L9;
		}
		if (k == 1 || k == 7 || k == 11 || k == 13 || k == 17 || k == 19 
			|| k == 21) {
		id2[j - 1] = 0;
		id3[j - 1] = 0;
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		goto L10;
		}
		if (k == 2) {
		id3[j - 1] = 0;
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 3) {
		id2[j - 1] = 0;
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 4) {
		id3[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 5) {
		id2[j - 1] = 0;
		id3[j - 1] = 0;
		id4[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 6) {
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		}
		if (k == 8) {
		id3[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 9) {
		id2[j - 1] = 0;
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 10) {
		id3[j - 1] = 0;
		id4[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 12) {
		id5[j - 1] = 0;
		}
		if (k == 14) {
		id3[j - 1] = 0;
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 15) {
		id2[j - 1] = 0;
		id4[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 16) {
		id3[j - 1] = 0;
		id5[j - 1] = 0;
		id6[j - 1] = 0;
		}
		if (k == 18) {
		id4[j - 1] = 0;
		id5[j - 1] = 0;
		}
		if (k == 20) {
		id3[j - 1] = 0;
		id6[j - 1] = 0;
		}
L9:
		;
	}
L10:
	;
	}
	if (*js == 3) {
	if (id2[0] != id2[1]) {
		goto L20;
	}
	if (id2[0] != id2[2]) {
		goto L20;
	}
	if (id3[0] != id3[1]) {
		goto L20;
	}
	if (id3[0] != id3[2]) {
		goto L20;
	}
	if (id4[0] != id4[1]) {
		goto L20;
	}
	if (id4[0] != id4[2]) {
		goto L20;
	}
	if (id5[0] != id5[1]) {
		goto L20;
	}
	if (id5[0] != id5[2]) {
		goto L20;
	}
	if (id6[0] != id6[1]) {
		goto L20;
	}
	if (id6[0] != id6[2]) {
		goto L20;
	}
	}
	iab2 = 1;
	iab3 = 1;
	iab4 = 1;
	iab5 = 1;
	iab6 = 1;
	if (*js == 2) {
	if (id2[0] != id2[1]) {
		iab2 = 0;
	}
	if (id3[0] != id3[1]) {
		iab3 = 0;
	}
	if (id4[0] != id4[1]) {
		iab4 = 0;
	}
	if (id5[0] != id5[1]) {
		iab5 = 0;
	}
	if (id6[0] != id6[1]) {
		iab6 = 0;
	}
	}
	if (*js == 4) {
	if (id2[0] != id2[1]) {
		iab2 = 0;
	}
	if (id3[0] != id3[1]) {
		iab3 = 0;
	}
	if (id4[0] != id4[1]) {
		iab4 = 0;
	}
	if (id5[0] != id5[1]) {
		iab5 = 0;
	}
	if (id6[0] != id6[1]) {
		iab6 = 0;
	}
	ihh = 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
		k = ihkl[j * 3 + 1] + ihkl[j * 3 + 2];
		for (i__ = 1; i__ <= 21; ++i__) {
		ii = (i__ << 1) - 1;
		if (k == ii) {
			ihh = 0;
		}
/* L8: */
		}
	}
	if (ihh == 1) {
		xx = 1.41421f;
		cel[*l * 6 + 1] /= xx;
		cel[*l * 6 + 2] /= xx;
		vgc[*l] = vgc[*l] / xx / xx;
	}
	}
	for (j = 1; j <= 2; ++j) {
	if (id6[j - 1] == 1) {
		id3[j - 1] = 0;
		id2[j - 1] = 0;
	}
	if (id4[j - 1] == 1) {
		id2[j - 1] = 0;
	}
	if (id2[j - 1] == 1 && ihmax[j - 1] >= 2 && iab2 == 1) {
		xx = 2.f;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id3[j - 1] == 1 && ihmax[j - 1] >= 3 && iab3 == 1) {
		xx = 3.f;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id4[j - 1] == 1 && ihmax[j - 1] >= 4 && iab4 == 1) {
		xx = 4.f;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id5[j - 1] == 1 && ihmax[j - 1] >= 5 && iab5 == 1) {
		xx = 5.f;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id6[j - 1] == 1 && ihmax[j - 1] >= 6 && iab6 == 1) {
		xx = 6.f;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
/* L11: */
	}
	j = 3;
	if (id6[j - 1] == 1) {
	id3[j - 1] = 0;
	id2[j - 1] = 0;
	}
	if (id4[j - 1] == 1) {
	id2[j - 1] = 0;
	}
	if (id2[j - 1] == 1 && ihmax[j - 1] >= 2) {
	xx = 2.f;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id3[j - 1] == 1 && ihmax[j - 1] >= 3) {
	xx = 3.f;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id4[j - 1] == 1 && ihmax[j - 1] >= 4) {
	xx = 4.f;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id5[j - 1] == 1 && ihmax[j - 1] >= 5) {
	xx = 5.f;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id6[j - 1] == 1 && ihmax[j - 1] >= 6) {
	xx = 6.f;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
L20:
	return 0;
} /* supcel_ */


/*  Test on Bravais lattice */

int brav_(integer *n, integer *ihkl, integer *ibr)
{
	/* System generated locals */
	integer i__1;

	/* Local variables */
	static integer h__, i__, j, k, l, jj, hsum;

	/* Parameter adjustments */
	ihkl -= 4;

	/* Function Body */
	*ibr = 6;

/*     I lattice */

	j = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*      write(*,*)ihkl(1,i),ihkl(2,i),ihkl(3,i) */
	hsum = ihkl[i__ * 3 + 1] + ihkl[i__ * 3 + 2] + ihkl[i__ * 3 + 3];
	k = 1;
	if (hsum / 2 << 1 == hsum) {
		k = 0;
	}
	j += k;
	}
	if (j == 0) {
	*ibr = 1;
	return 0;
	}

/*     F lattice */

	j = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = 1;
	if (ihkl[i__ * 3 + 1] / 2 << 1 == ihkl[i__ * 3 + 1]) {
		h__ = 0;
	}
	k = 1;
	if (ihkl[i__ * 3 + 2] / 2 << 1 == ihkl[i__ * 3 + 2]) {
		k = 0;
	}
	l = 1;
	if (ihkl[i__ * 3 + 3] / 2 << 1 == ihkl[i__ * 3 + 3]) {
		l = 0;
	}
	jj = h__ + k + l;
	if (jj != 0) {
		h__ = 1;
		if (ihkl[i__ * 3 + 1] / 2 << 1 != ihkl[i__ * 3 + 1]) {
		h__ = 0;
		}
		k = 1;
		if (ihkl[i__ * 3 + 2] / 2 << 1 != ihkl[i__ * 3 + 2]) {
		k = 0;
		}
		l = 1;
		if (ihkl[i__ * 3 + 3] / 2 << 1 != ihkl[i__ * 3 + 3]) {
		l = 0;
		}
		jj = h__ + k + l;
	}
	j += jj;
	}
	if (j == 0) {
	*ibr = 5;
	return 0;
	}

/*     A lattice */

	j = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	hsum = ihkl[i__ * 3 + 2] + ihkl[i__ * 3 + 3];
	k = 1;
	if (hsum / 2 << 1 == hsum) {
		k = 0;
	}
	j += k;
	}
	if (j == 0) {
	*ibr = 2;
	return 0;
	}

/*     B lattice */

	j = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	hsum = ihkl[i__ * 3 + 1] + ihkl[i__ * 3 + 3];
	k = 1;
	if (hsum / 2 << 1 == hsum) {
		k = 0;
	}
	j += k;
	}
	if (j == 0) {
	*ibr = 3;
	return 0;
	}

/*     C lattice */

	j = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	hsum = ihkl[i__ * 3 + 1] + ihkl[i__ * 3 + 2];
	k = 1;
	if (hsum / 2 << 1 == hsum) {
		k = 0;
	}
	j += k;
	}
	if (j == 0) {
	*ibr = 4;
	return 0;
	}
	return 0;
} /* brav_ */


/* ------------------------------------------------------------ */

int mcmnam_(integer *ln, char *nam, ftnlen nam_len)
{
	/* System generated locals */
	integer i__1;

	/* Builtin functions */
	int s_copy(char *, char *, ftnlen, ftnlen);

	/* Local variables */
	static integer i__;
	static char kr[80], ks[1];
	extern integer iargc_(void);
	extern int getarg_(integer *, char *, ftnlen);


/* Get generic filename (NAME) and its length (LN) from the command */
/* line. LN set to 0 if the file names are to be defined externally */

	s_copy(kr, " ", (ftnlen)80, (ftnlen)1);
	i__1 = iargc_();
	getarg_(&i__1, kr, (ftnlen)80);
	*ln = 0;
	s_copy(nam, " ", (ftnlen)80, (ftnlen)1);
	for (i__ = 1; i__ <= 80; ++i__) {
	*(unsigned char *)ks = *(unsigned char *)&kr[i__ - 1];
	if (*(unsigned char *)ks == ' ') {
		goto L2;
	}
	++(*ln);
	*(unsigned char *)&nam[*ln - 1] = *(unsigned char *)ks;
L2:
	;
	}
	return 0;
} /* mcmnam_ */

int main () { MAIN__ (); return 0; }
