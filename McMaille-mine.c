/*
*     Version 4.00 parallelized with OpenMP for multi-core processors
*       but this version is slightly modified for monoprocessors
*
*     MAILLE in french = CELL in english
*     Mc for Monte Carlo
*     Pronounce : MacMy
*
************************************************************************
*
*     A Monte Carlo and grid search code for indexing powder patterns
*
*     For more details see the documentation at
*                   http://www.cristal.org/McMaille/
*             or    http://sdpd.univ-lemans.fr/McMaille/
*
*              by A. Le Bail - September 2002 for version 0.9
*                              as well as for versions 1.0, 2.0 and 3.0
*                              October 2006 for version 4.00
*                        alb@cristal.org
*                        http://www.cristal.org/
*
*                        Résidence Cristal - Appt 213
*                        2, rue de Gasperi
*                        72100 Le Mans
*                        FRANCE
*
*   Versions 0.9 : cubic only
*            1.0 : hexagonal/trigonal/rhombohedral, tetragonal,
*                  orthorhombic added, plus .ckm and .prf files
*            2.0 : monoclinic and triclinic added in MC
*                  but not in grid search (too long)
*            3.0 : columnar peak shapes instead of Gaussian
*                  in versions 0.9-2.0
*                  no Le Bail fit contrarily to versions 0.9-2.0
*                  only fit by percentage of inclusion of the
*                  calculated column into the observed one
*            3.02: black box mode
*            3.03: improved Monte Carlo
*            3.04: two-phases mode
*            4.00: automatisation improved : more chances to identify
*                   the correct cell in "black box" mode
*                  Identification of the Bravais lattice
*                  Parallelization by using OpenMP directives
*                  improving the speed with multicore processors
*                  speed x1.7 to 1.8 with "dual core" or "core duo"
*                  speed x3.6 expected with the quad core in 2007
*                  speed x79 expected with the 80-core in 2012...;-)
*
*
*
************************************************************************
*
*    Copyright (C) 2002-2006 Armel Le Bail
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc.,
* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*
*
***********************************************************************
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#define VERSION "4.00"
#define FILENAME_SIZE 20
#define FILE_WITH_EXTENSION_SIZE FILENAME_SIZE + 4
#define N_HKL 10000

int nhkl0, lhkl, ndat;
double dmin, slabda2;
int ihh[30000];
double al[9], cri[N_HKL], difp[N_HKL], difm[N_HKL], th2[N_HKL], fobs[N_HKL], sum_f;
int nind;
double w2[N_HKL];
int nmx, ndat10;
int nhkl;
int ind[1000000];

const double pi = 114.59156;  // 360./3.1415926

//truc
double qq[2000]	/* was [200][10] */, bb[10], b[10];
int h[200], k[200], l[200], npaf;
double afi[10];
int nr, indic;
double pds[200];
int npaf2;

/* ... Here is the program heart... */

int calcul1(double *diff, double *diff2)
{
	/* System generated locals */
	int i1, i2;
	double r1;

	/* Local variables */
	static int i, j, k, l;
	static double x, de[10000];
	static int jh;
	static double cr[10000], perc[10000], diff1, demax, theta[10000], sinth,
		sum_f2;

/* $OMP THREADPRIVATE(/cal/) */

/* ...  Keep only the hkl for d > dmin */

	jh = 0;
	i1 = nhkl0;
	for (i = 1; i <= i1; ++i) {
	x = 0.f;
	for (j = 1; j <= 3; ++j) {
		for (k = j; k <= 3; ++k) {
			x = al[j + k * 3 - 4] * ihh[j + i * 3 - 4] *
			ihh[k + i * 3 - 4] + x;
		}
	}
	if (x > dmin) {
		goto L109;
	}
	++jh;
/*     X IS 1/D(hkl)**2 FOR REFLECTION IHH */

/*     This should be optimized for speed : */
/*     working only on X, not calculating 2-theta... */
/*     - in fact, tests show that only 10-15% is gained - */

	sinth = slabda2 * x;
	sinth = sqrt(sinth);
	theta[jh - 1] = asin(sinth) * pi;
	cr[jh - 1] = 0.f;
L109:
	;
	}
	nhkl = jh;
	if (nhkl > ndat10) {
	return 0;
	}

/* ...  Comparison with the data */

	lhkl = 0;
	i1 = ndat;
	for (j = 1; j <= i1; ++j) {
	cri[j - 1] = 0.f;
	perc[j - 1] = 0.f;
	demax = w2[j - 1];
	i2 = nhkl;
	for (k = 1; k <= i2; ++k) {
		if (cr[k - 1] == 1.f) {
			goto L111;
		}
		if (theta[k - 1] <= difp[j - 1] && theta[k - 1] >= difm[j - 1]) {
			de[k - 1] = (r1 = theta[k - 1] - th2[j - 1], abs(r1));
			if (de[k - 1] <= demax) {
				l = k;
				demax = de[k - 1];
				cri[j - 1] = 1.f;
			}
		}
L111:
		;
	}

/*  PERC = percentage of columnar overlap for that peak */

/*  Potential problem here because only one reflection */
/*  overlapping the most closely with the column is */
/*  included (if CRI =1)... */

	if (cri[j - 1] == 1.f) {
		perc[j - 1] = 1.f - de[l - 1] / w2[j - 1];
		++lhkl;
		cr[l - 1] = 1.f;
	}
/* L113: */
	}

/* ...  Calculate "R" */

	diff1 = 0.f;
	sum_f2 = 0.f;
	i1 = ndat;
	for (k = 1; k <= i1; ++k) {
	sum_f2 += fobs[k - 1] * cri[k - 1];
/* L1122: */
	diff1 += cri[k - 1] * fobs[k - 1] * perc[k - 1];
	}
	*diff = 1.f - diff1 / sum_f;
	*diff2 = 1.f - diff1 / sum_f2;
	return 0;
} /* calcul1_ */

//==============================================================================

int calcul2(double *diff, int *ihkl, double *th3, int *ncalc, int *igc)
{
	/* System generated locals */
	int i1, i2;
	double r1;

	/* Local variables */
	static int i, j, k, l;
	static double x;
	static int l2;
	static double de[10000];
	static int jh;
	static double qc[10000], cr[10000];
	static int jj, jjj, jhkl[30000]	/* was [3][10000] */;
	static double perc[10000], demax, theta[10000], sinth, sum_f2;

/* $OMP THREADPRIVATE(/cal/,/cal2/) */

/* ...  Keep only the hkl for d > dmin */

	/* Parameter adjustments */
	--th3;
	ihkl -= 4;

	/* Function Body */
	jh = 0;
	i1 = nhkl0;
	for (i = 1; i <= i1; ++i) {
	x = 0.f;
	for (j = 1; j <= 3; ++j) {
		for (k = j; k <= 3; ++k) {
/* L107: */
			x = al[j + k * 3 - 4] * ihh[j + i * 3 - 4] *
			ihh[k + i * 3 - 4] + x;
		}
	}
	if (x > dmin) {
		goto L109;
	}
	++jh;
	for (j = 1; j <= 3; ++j) {
/* L108: */
		jhkl[j + jh * 3 - 4] = ihh[j + i * 3 - 4];
	}
/*     X IS 1/D(hkl)**2 FOR REFLECTION IHH */

/*     This should be optimized for speed : */
/*     working only on X, not calculating 2-theta... */

	qc[jh - 1] = x;
	sinth = slabda2 * x;
	sinth = sqrt(sinth);
	theta[jh - 1] = asin(sinth) * pi;
	cr[jh - 1] = 0.f;
L109:
	;
	}
	nhkl = jh;

/* ...  Comparison with the data */

	lhkl = 0;
	*ncalc = 0;
	jj = 0;
	i1 = ndat;
	for (j = 1; j <= i1; ++j) {
	cri[j - 1] = 0.f;
	perc[j - 1] = 0.f;
/* CC */
/* CC  Eliminating too spurious peaks here ??? */
/* CC    tolerance on width decreased by a factor 3 */
/* CC */
	demax = w2[j - 1] * .33333f;
/* CC */
	i2 = nhkl;
	for (k = 1; k <= i2; ++k) {
		if (cr[k - 1] == 1.f) {
		goto L111;
		}
		if (theta[k - 1] <= difp[j - 1] && theta[k - 1] >=
			difm[j - 1]) {
		de[k - 1] = (r1 = theta[k - 1] - th2[j - 1], abs(
			r1));
		if (de[k - 1] <= demax) {
			l = k;
			demax = de[k - 1];
			cri[j - 1] = 1.f;
		}
		}
L111:
		;
	}

/*  PERC = percentage of columnar overlap for that peak */

/*  Potential problem here because only one reflection */
/*  overlapping the most closely with the column is */
/*  included (if CRI =1)... */

	if (cri[j - 1] == 1.f) {
		perc[j - 1] = 1.f - de[l - 1] / (w2[j - 1] * .33333f);
		++lhkl;
		cr[l - 1] = 1.f;
		for (l2 = 1; l2 <= 3; ++l2) {
/* L112: */
			ihkl[l2 + lhkl * 3] = jhkl[l2 + l * 3 - 4];
		}
		th3[lhkl] = th2[j - 1];
		jj += cri[j - 1];
		if (jj <= 20) {
		jjj = j;
		}
	}
/* L113: */
	}
/*  NCALC is for M(20) FoM */
	i1 = nhkl;
	for (k = 1; k <= i1; ++k) {
	if (th2[jjj - 1] >= theta[k - 1]) {
		++(*ncalc);
	}
/* L114: */
	}

/* ...  Calculate "R" */

	*diff = 0.f;

/*  Change here with SUM_F2 being only on explained reflections... */

	sum_f2 = 0.f;
	i1 = ndat;
	for (k = 1; k <= i1; ++k) {
	ind[k + *igc * 100 - 101] = cri[k - 1];
	sum_f2 += fobs[k - 1] * cri[k - 1];
/* L1122: */
	*diff += cri[k - 1] * fobs[k - 1] * perc[k - 1];
	}
	*diff = 1.f - *diff / sum_f2;
	return 0;
} /* calcul2_ */

//==============================================================================

/*  Test on Bravais lattice */

int brav(int *n, int *ihkl, int *ibr)
{
	/* System generated locals */
	int i1;

	/* Local variables */
	static int h, i, j, k, l, jj, hsum;

	/* Parameter adjustments */
	ihkl -= 4;

	/* Function Body */
	*ibr = 6;

/*     I lattice */

	j = 0;
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
/*      write(*,*)ihkl(1,i),ihkl(2,i),ihkl(3,i) */
	hsum = ihkl[i * 3 + 1] + ihkl[i * 3 + 2] + ihkl[i * 3 + 3];
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
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	h = 1;
	if (ihkl[i * 3 + 1] / 2 << 1 == ihkl[i * 3 + 1]) {
		h = 0;
	}
	k = 1;
	if (ihkl[i * 3 + 2] / 2 << 1 == ihkl[i * 3 + 2]) {
		k = 0;
	}
	l = 1;
	if (ihkl[i * 3 + 3] / 2 << 1 == ihkl[i * 3 + 3]) {
		l = 0;
	}
	jj = h + k + l;
	if (jj != 0) {
		h = 1;
		if (ihkl[i * 3 + 1] / 2 << 1 != ihkl[i * 3 + 1]) {
		h = 0;
		}
		k = 1;
		if (ihkl[i * 3 + 2] / 2 << 1 != ihkl[i * 3 + 2]) {
		k = 0;
		}
		l = 1;
		if (ihkl[i * 3 + 3] / 2 << 1 != ihkl[i * 3 + 3]) {
		l = 0;
		}
		jj = h + k + l;
	}
	j += jj;
	}
	if (j == 0) {
	*ibr = 5;
	return 0;
	}

/*     A lattice */

	j = 0;
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	hsum = ihkl[i * 3 + 2] + ihkl[i * 3 + 3];
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
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	hsum = ihkl[i * 3 + 1] + ihkl[i * 3 + 3];
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
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
	hsum = ihkl[i * 3 + 1] + ihkl[i * 3 + 2];
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

//==============================================================================

int sigma(double *d, double *db, double *r)
{
	/* System generated locals */
	double r1;

	/* Local variables */
	static int i;

/* ..... */
	/* Parameter adjustments */
	--db;
	--d;

	/* Function Body */
	*r = 0.f;
	for (i = 1; i <= 6; ++i) {
/* L10: */
/* Computing 2nd power */
	r1 = d[i] * db[i + 2];
	*r += r1 * r1;
	}
	return 0;
} /* sigma_ */

//==============================================================================

int inver(double *b, double *db, double *volum, int *iv)
{

	/* Local variables */
	static double d[6];
	static int i, j, k;
	static double q, r, q2, ad[6], cc[3];
	static int jj;
	static double sp[3], ss[3], dqd[3], sig[6], cosp[3], sinp[3], cabc2;

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
	for (i = 1; i <= 3; ++i) {
	ad[i - 1] = b[i + 2];
	cosp[i - 1] = cos(b[i + 5]);
	sinp[i - 1] = sin(b[i + 5]);
/* L10: */
	}
	for (i = 1; i <= 3; ++i) {
	j = i % 3 + 1;
	k = (i + 1) % 3 + 1;
	dqd[i - 1] = cosp[i - 1] - cosp[j - 1] * cosp[k - 1];
	ss[i - 1] = sinp[j - 1] * sinp[k - 1];
	cc[i - 1] = -dqd[i - 1] / ss[i - 1];
	cabc2 += cosp[i - 1] * cosp[i - 1];
/* L15: */
	}

	q2 = 1.f - cabc2 + cosp[0] * 2.f * cosp[1] * cosp[2];
	q = sqrt(q2);
	*volum = ad[0] * ad[1] * ad[2] * q;

	for (i = 1; i <= 3; ++i) {
	b[i + 2] = sinp[i - 1] / (ad[i - 1] * q);
	b[i + 5] = acos(cc[i - 1]);
	sp[i - 1] = sin(b[i + 5]);
/* L20: */
	}

	if (*iv == 0) {
	goto L70;
	}

/* .....DERIVEES DES PARAMETRES A , B , C */
	for (i = 1; i <= 3; ++i) {
	j = i % 3 + 1;
	k = (i + 1) % 3 + 1;
	d[i - 1] = -sinp[i - 1] / ad[i - 1];
	d[j - 1] = 0.f;
	d[k - 1] = 0.f;
	d[i + 2] = cosp[i - 1] - sinp[i - 1] * sinp[i - 1] * dqd[
		i - 1] / (q * q);
	d[j + 2] = -ss[k - 1] * dqd[j - 1] / q2;
	d[k + 2] = -ss[j - 1] * dqd[k - 1] / q2;
	sigma(d, &db[1], &r);
/* L30: */
	sig[i - 1] = sqrt(r) / (q * ad[i - 1]);
	}

/* .....DERIVEES DES ANGLES DE LA MAILLE */
	for (i = 1; i <= 3; ++i) {
	for (jj = 1; jj <= 3; ++jj) {
/* L40: */
		d[jj - 1] = 0.f;
	}
	j = i % 3 + 1;
	k = (i + 1) % 3 + 1;
	d[i + 2] = sinp[i - 1] / ss[i - 1];
	d[j + 2] = cosp[k - 1] / sinp[k - 1] + cosp[j - 1] * cc[i - 1] /
		sinp[j - 1];
	d[k + 2] = cosp[j - 1] / sinp[j - 1] + cosp[k - 1] * cc[i - 1] /
		sinp[k - 1];
	sigma(d, &db[1], &r);
/* L50: */
	sig[i + 2] = sqrt(r) / sp[i - 1];
	}

	for (i = 1; i <= 6; ++i) {
/* L60: */
	db[i + 2] = sig[i - 1];
	}

L70:
	return 0;
} /* inver */

//==============================================================================

int matinv(double *am, int *n, int *nfail)
{
	/* System generated locals */
	int i1, i2, i3;

	/* Local variables */
	static int i, j, k, l, m, ii, kdm, kli, kmi, imax;
	static double suma;
	static double term, denom;

/* .....----------------------------- */
/*     ********** SEGMENT 1 OF CHOLESKI INVERSION ********** */
/*     ***** FACTOR MATRIX INTO LOWER TRIANGLE X TRANSPOSE ***** */
	/* Parameter adjustments */
	--am;

	/* Function Body */
	k = 1;
	if ((i1 = *n - 1) < 0) {
	goto L8;
	} else if (i1 == 0) {
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
	i1 = *n;
	for (m = 1; m <= i1; ++m) {
	imax = m - 1;
/*     ***** LOOP L OF A(L,M) ***** */
	i2 = *n;
	for (l = m; l <= i2; ++l) {
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
		i3 = imax;
		for (i = 1; i <= i3; ++i) {
		suma += am[kli] * am[kmi];
		j = *n - i;
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
	i1 = *n;
	for (l = 2; l <= i1; ++l) {
	kdm = kdm + *n - l + 2;
/*     ***** RECIPROCAL OF DIAGONAL TERM ***** */
	term = 1.f / am[kdm];
	am[kdm] = term;
	kmi = 0;
	kli = l;
	imax = l - 1;
/*     ***** STEP M OF B(L,M) ***** */
	i2 = imax;
	for (m = 1; m <= i2; ++m) {
		k = kli;
/*     ***** SUM TERMS ***** */
		suma = 0.f;
		i3 = imax;
		for (i = m; i <= i3; ++i) {
		ii = kmi + i;
		suma -= am[kli] * am[ii];
/* L130: */
		kli = kli + *n - i;
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
	i1 = *n;
	for (m = 1; m <= i1; ++m) {
	kli = k;
	i2 = *n;
	for (l = m; l <= i2; ++l) {
		kmi = k;
		imax = *n - l + 1;
		suma = 0.f;
		i3 = imax;
		for (i = 1; i <= i3; ++i) {
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

//==============================================================================

int calc(void)
{
	/* System generated locals */
	int i1;
	double r1, r2, r3;

	/* Local variables */
	static double d, f;
	static int i, j;
	static double q[2000]	/* was [200][10] */, ae, be, ce, dd;
	static int ir;
	static double cae, cbe, cce, rad;

/* .....-------------- */
/* $OMP THREADPRIVATE(/TROC/,/TRUC/) */
	j = 0;
	for (i = 1; i <= 8; ++i) {
	if (afi[i - 1] == 0.f) {
		goto L1;
	}
	++j;
	b[i - 1] = bb[j - 1];
L1:
	;
	}
	if (indic == 1 || indic == 2) {
	b[3] = b[2];
	}
	if (indic == 1) {
	b[4] = b[2];
	}
	ae = b[2];
	be = b[3];
	ce = b[4];
	cae = cos(b[5]);
	cbe = cos(b[6]);
	cce = cos(b[7]);
	i1 = nr;
	for (i = 1; i <= i1; ++i) {
/* Computing 2nd power */
	r1 = ae * h[i - 1];
/* Computing 2nd power */
	r2 = be * k[i - 1];
/* Computing 2nd power */
	r3 = ce * l[i - 1];
	dd = r1 * r1 + r2 * r2 + r3 * r3 + (h[i - 1] *
		k[i - 1] * ae * be * cce + k[i - 1] *
		l[i - 1] * be * ce * cae + l[i - 1] *
		h[i - 1] * ce * ae * cbe) * 2.f;
	d = 1.f / sqrt(dd);
/* Computing 2nd power */
	r1 = b[1];
	rad = sqrt(1.f - r1 * r1 * dd);
	f = b[1] * d / rad;
	q[i - 1] = 1.f;
	q[i + 199] = 1.f / (d * rad);
	q[i + 399] = f * h[i - 1] * (h[i - 1] * ae +
		k[i - 1] * be * cce + l[i - 1] * ce * cbe);
	q[i + 599] = f * k[i - 1] * (k[i - 1] * be +
		l[i - 1] * ce * cae + h[i - 1] * ae * cce)
		;
	q[i + 799] = f * l[i - 1] * (l[i - 1] * ce +
		h[i - 1] * ae * cbe + k[i - 1] * be * cae)
		;
	q[i + 999] = -f * k[i - 1] * l[i - 1] * be * ce *
		sin(b[5]);
	q[i + 1199] = -f * l[i - 1] * h[i - 1] * ce *
		ae * sin(b[6]);
	q[i + 1399] = -f * h[i - 1] * k[i - 1] * ae *
		be * sin(b[7]);
/* L2: */
	qq[i + npaf2 * 200 - 201] = b[0] + asin(
		b[1] / d);
	}
	i1 = nr;
	for (ir = 1; ir <= i1; ++ir) {
	j = 0;
	for (i = 1; i <= 8; ++i) {
		if (afi[i - 1] == 0.f) {
		goto L3;
		}
		++j;
		qq[ir + j * 200 - 201] = q[ir + i * 200 - 201];
L3:
		;
	}
	}
/*      WRITE(IWR,5)(BB(I),I=1,NPAF) */
/* L5: */
	return 0;
} /* calc_ */

//==============================================================================

int mcrnl(double *q, int *id, double *y, double *b, int *
	m, int *n, double *p, int *ifin)
{
	/* System generated locals */
	int q_dim1, q_offset, i1, i2, i3;

	/* Local variables */
	static double a[60];
	static int i, j, l;
	static double r;
	static int m1, m2, il, mm, iq, no, ier, imm;

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
	calc();

/* -----LES Y CALCULES SONT DANS LA COLONNE M+2 */
	i1 = *m;
	for (j = 1; j <= i1; ++j) {
	r = 0.f;
	i2 = *n;
	for (i = 1; i <= i2; ++i) {
/* L20: */
		r += (y[i] - q[i + m2 * q_dim1]) * q[i + j * q_dim1] * p[
			i];
	}
/* L30: */
	q[j + m1 * q_dim1] = r;
	}

/* -----CONSTRUCTION DE LA MATRICE SYMETRIQUE A=Q*QT */
	no = 0;
	i1 = *m;
	for (iq = 1; iq <= i1; ++iq) {
	i2 = *m;
	for (il = iq; il <= i2; ++il) {
		++no;
		r = 0.f;
		i3 = *n;
		for (i = 1; i <= i3; ++i) {
/* L40: */
		r += q[i + iq * q_dim1] * q[i + il * q_dim1] * p[i];
		}
/* L50: */
		a[no - 1] = r;
	}
	}

/* -----INVERSION DE LA MATRICE A */
	matinv(a, m, &ier);
	if (*ifin <= 0) {
	goto L110;
	} else {
	goto L60;
	}

/* -----CALCUL DES NOUVEAUX PARAMETRES */
L60:
	i2 = *m;
	for (i = 1; i <= i2; ++i) {
	r = 0.f;
	imm = (i - 1) * (mm - i) / 2;
	i1 = *m;
	for (j = 1; j <= i1; ++j) {
		if (j - i >= 0) {
		goto L70;
		} else {
		goto L80;
		}
L70:
		l = imm + j;
		goto L90;
L80:
		l = (j - 1) * (mm - j) / 2 + i;
L90:
		r += a[l - 1] * q[j + m1 * q_dim1];
	}
/* L100: */
	b[i] += r;
	}
	goto L10;

/* -----REMISE DES ELEMENTS DIAGONAUX DANS Q(I,I),APRES LE DERNIER CYCLE */
L110:
	i2 = *m;
	for (i = 1; i <= i2; ++i) {
	l = (i - 1) * (mm - i) / 2 + i;
/* L120: */
	q[i + i * q_dim1] = a[l - 1];
	}
	return 0;
}

//==============================================================================

int fonc(double *theta, double *r, double *rr)
{
	/* System generated locals */
	int i1;
	double r1l, r2l, r3l;

	/* Local variables */
	static double d;
	static int i;
	static double r1, r2, ae, be, ce, dd, yc, cae, cbe, cce;

/* .....----------------------- */
/* $OMP THREADPRIVATE(/TRUC/) */
	/* Parameter adjustments */
	--theta;

	/* Function Body */
	r1 = 0.f;
	r2 = 0.f;
	ae = b[2];
	be = b[3];
	ce = b[4];
	cae = cos(b[5]);
	cbe = cos(b[6]);
	cce = cos(b[7]);
	i1 = nr;
	for (i = 1; i <= i1; ++i) {
/* Computing 2nd power */
	r1l = ae * h[i - 1];
/* Computing 2nd power */
	r2l = be * k[i - 1];
/* Computing 2nd power */
	r3l = ce * l[i - 1];
	dd = r1l * r1l + r2l * r2l + r3l * r3l + (h[i - 1] *
		k[i - 1] * ae * be * cce + k[i - 1] *
		l[i - 1] * be * ce * cae + l[i - 1] *
		h[i - 1] * ce * ae * cbe) * 2.f;
	d = 1.f / sqrt(dd);
	yc = b[0] + asin(b[1] / d);
/* Computing 2nd power */
	r1l = yc - theta[i];
	r1 += r1l * r1l;
/* Computing 2nd power */
	r1l = theta[i];
	r2 += r1l * r1l;
/* L2: */
	qq[i + npaf2 * 200 - 201] = yc;
	}
	*r = r1 / (nr - npaf);
	*rr = sqrt(r1 / r2);
	return 0;
} /* fonc */

//==============================================================================

int celref2(int *indi, double *bbb, double *afin, int *
	nhkl, double *theta, int *jhkl, double *ddt, double *ddq)
{
	/* Initialized data */

	static double sig[8] = { 0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f };

	/* System generated locals */
	int i1, i2, i3, i4;
	double r1;

	/* Local variables */
	static int i, j;
	static double r, y1, y2, y3, y4;
	static int ik, jj;
	static double rd, qc, qo, rr, dum[3], pip;
	static int iffi, ifin, ihkl;
	static int ndmax;
	static double volum;
	static int npour;


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
	indic = *indi;
	*ddt = 0.f;
	*ddq = 0.f;
	ifin = 10;
	if (indic == 0 || indic > 3) {
	indic = 3;
	}
	b[0] = 0.f;
	for (i = 1; i <= 8; ++i) {
	afi[i - 1] = afin[i];
/* L5500: */
	b[i - 1] = bbb[i];
	}
	iffi = (int) (afi[2] + afi[3] + afi[4] + .1f);
	if (iffi == 0 || indic == 3) {
	goto L230;
	}
	ik = 3 - indic;
	i1 = ik;
	for (i = 1; i <= i1; ++i) {
	if (indic - 2 >= 0) {
		goto L210;
	} else {
		goto L200;
	}
L200:
	afi[indic + 2 + i - 1] = 0.f;
	goto L220;
L210:
	afi[indic + 1 + i - 1] = 0.f;
L220:
	;
	}
	afi[2] = 1.f;

L230:
	nr = 0;
	i1 = *nhkl;
	for (nr = 1; nr <= i1; ++nr) {
	if (nr > ndmax) {
		goto L400;
	}
	h[nr - 1] = jhkl[nr * 3 + 1];
	k[nr - 1] = jhkl[nr * 3 + 2];
	l[nr - 1] = jhkl[nr * 3 + 3];
	ihkl = (i2 = h[nr - 1], abs(i2)) + (i3 =
		k[nr - 1], abs(i3)) + (i4 = l[
		nr - 1], abs(i4));
	if (ihkl == 0) {
		goto L260;
	}
	if (theta[nr] <= 0.f) {
		goto L400;
	}
/* L250: */
	pds[nr - 1] = 1.f;
/* L240: */
	theta[nr] = theta[nr] / rd / 2.f;
	}
/* ..... */
L260:
	nr = *nhkl;
	b[0] /= rd;
	b[1] *= .5f;
	for (i = 6; i <= 8; ++i) {
/* L270: */
	b[i - 1] /= rd;
	}
	int c0;
	inver(b, dum, &volum, &c0);
	for (i = 1; i <= 3; ++i) {
/* L280: */
	dum[i - 1] = b[i + 4] * rd;
	}
	j = 0;
	for (i = 1; i <= 8; ++i) {
	if (afi[i - 1] == 0.f) {
		goto L290;
	}
	++j;
	bb[j - 1] = b[i - 1];
L290:
	;
	}
	npaf = j;
	if (npaf == 8) {
	goto L400;
	}
	npaf2 = npaf + 2;
	mcrnl(qq, &ndmax, &theta[1], bb, &npaf, &nr,
		pds, &ifin);
	j = 0;
	for (i = 1; i <= 8; ++i) {
	if (afi[i - 1] == 0.f) {
		goto L300;
	}
	++j;
	b[i - 1] = bb[j - 1];
L300:
	;
	}
	fonc(&theta[1], &r, &rr);
	jj = 0;
	for (i = 1; i <= 8; ++i) {
	if (afi[i - 1] != 0.f) {
		goto L320;
	} else {
		goto L310;
	}
L310:
	sig[i - 1] = 0.f;
	goto L330;
L320:
	++jj;
	sig[i - 1] = sqrt(qq[jj + jj * 200 - 201] * r);
L330:
	;
	}
	if (indic == 1 || indic == 2) {
	sig[3] = sig[2];
	}
	if (indic == 1) {
	sig[4] = sig[2];
	}
	for (i = 1; i <= 3; ++i) {
	bb[i - 1] = b[i + 1];
/* L340: */
	bb[i + 2] = b[i + 4] * rd;
	}
	sig[0] *= rd;
	sig[1] *= 2.f;
	b[0] = -b[0] * rd * 2.f;
	sig[0] *= 2.f;
	b[1] *= 2.f;
	for (i = 6; i <= 8; ++i) {
	sig[i - 1] *= rd;
/* L350: */
	b[i - 1] *= rd;
	}
	for (i = 1; i <= 8; ++i) {
/* L351: */
	bbb[i] = b[i - 1];
	}
	npour = 0;
	i1 = nr;
	for (i = 1; i <= i1; ++i) {
	y1 = theta[i] * rd;
	y2 = y1 + b[0] / 2.f;
	y3 = qq[i + npaf2 * 200 - 201] * rd + b[0] /
		2.f;
	y4 = y2 - y3;
	y1 *= 2.f;
	y2 *= 2.f;
	y3 *= 2.f;
	y4 *= 2.f;
	if (i <= 20) {
		*ddt += abs(y4);
/* Computing 2nd power */
		r1 = sin(y2 * pip) * 2.f / bbb[2];
		qo = r1 * r1;
/* Computing 2nd power */
		r1 = sin(y3 * pip) * 2.f / bbb[2];
		qc = r1 * r1;
		*ddq += (r1 = qo - qc, abs(r1));
	}
/* L360: */
	}
	*ddt /= 20.f;
	*ddq /= 20.f;
L400:
	return 0;
} /* celref2_ */

//==============================================================================

int perm(int *i, int *j, int *k)
{
	/* System generated locals */
	int i1;

/*     PERMS USEFUL COMBINATIONS OF INTEGERS IN THE RANGE 1 TO 3 */
	if ((i1 = *i - 2) < 0) {
	goto L10;
	} else if (i1 == 0) {
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

//==============================================================================

int trcl(double *celln, double *rcelln, double *v)
{
	/* System generated locals */
	double r1;

	/* Local variables */
	static int i, j, k, l;
	static double abc, sina[3];
	static double prod;

/*     TRANSFORMS REAL CELL TO RECIPROCAL OR VICE VERSA */
/*     INPUT CELL IS IN ARRAY CELL AS LENGTHS AND COSINES */
	/* Parameter adjustments */
	--rcelln;
	--celln;

	/* Function Body */
	abc = 1.f;
	prod = 2.f;
	*v = -2.f;
	for (i = 1; i <= 3; ++i) {
	l = i + 3;
/* Computing 2nd power */
	r1 = celln[l];
	sina[i - 1] = 1.f - r1 * r1;
	*v += sina[i - 1];
	sina[i - 1] = sqrt(sina[i - 1]);
	prod *= celln[l];
/* L10: */
	abc *= celln[i];
	}
	*v = abc * sqrt(*v + prod);
/*      V IS CELL VOLUME */
/*     PUT INVERTED CELL INTO RCELL */
	for (i = 1; i <= 3; ++i) {
	perm(&i, &j, &k);
	rcelln[i] = celln[j] * celln[k] * sina[i - 1] / *v;
	l = i + 3;
/* L20: */
	rcelln[l] = (celln[j + 3] * celln[k + 3] - celln[l]) / (sina[j - 1] *
		sina[k - 1]);
	}
	return 0;
}

//==============================================================================

/* ... Check for supercell, and reduce to minimal cell */

int supcel(int *n, int *ihkl, double *cel, int *l, double *vgc, int *js)
{
	/* System generated locals */
	int i1, i2;

	/* Local variables */
	static int i, j, k, ii;
	static double xx;
	static int id2[3], id3[3], id4[3], id5[3], id6[3], ihh, iab2, iab3,
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
	i1 = *n;
	for (i = 1; i <= i1; ++i) {
/*  Test on dividing h, k, l by 2, 3, 4, 5, 6 */
		k = (i2 = ihkl[j + i * 3], abs(i2));
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
	i1 = *n;
	for (j = 1; j <= i1; ++j) {
		k = ihkl[j * 3 + 1] + ihkl[j * 3 + 2];
		for (i = 1; i <= 21; ++i) {
		ii = (i << 1) - 1;
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


//==============================================================================

int dcell(double *celln, double *al, double *v)
{
	/* Local variables */
	static int i, j, k, l;
	static double cell[6];
	static double degrad, rcelln[6];

	/* Parameter adjustments */
	al -= 4;
	--celln;

	/* Function Body */
	degrad = .017453277777777776f;
	for (i = 1; i <= 6; ++i) {
/* L30: */
	cell[i - 1] = celln[i];
	}
	for (i = 1; i <= 3; ++i) {
	l = i + 3;
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
	trcl(cell, rcelln, v);
/*     RCELL IS THE RECIPROCAL CELL CONSTANTS */
	for (i = 1; i <= 3; ++i) {
	al[i + i * 3] = rcelln[i - 1] * rcelln[i - 1];
	perm(&i, &j, &k);
	if (j - k >= 0) {
		goto L36;
	} else {
		goto L35;
	}
L35:
	al[j + k * 3] = rcelln[j - 1] * 2.f * rcelln[k - 1] * rcelln[i + 2];
	goto L34;
L36:
	al[k + j * 3] = rcelln[j - 1] * 2.f * rcelln[k - 1] * rcelln[i + 2];
L34:
	;
	}
/*      DO 37 I=4,6 */
/* 37    RCELLN(I)=ACOS(RCELLN(I))/DEGRAD */
	return 0;
}

//==============================================================================

int peekchar()
{
	int c;

	c = getchar();
	if(c != EOF) ungetc(c, stdin);      /* puts it back */

	return c;
}

int killk(int *pressedk)
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
}

//==============================================================================

double randi(int *ix)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	static int k;
	static double x;
	static int xhi, xlo, xahi, xalo, leftlo;


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
/*     double number having the value x/(2**31-1) */

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
	x = (double) xhi * 65536. + (double) (*ix - (xhi << 16));
	x *= 4.6566128752457969e-10;
	ret_val = (double) x;
	return ret_val;
} /* randi_ */

//==============================================================================

int espInit(FILE *file)
{
/* PORTLIB/DFPORT is compiler spcecific part, introduced for */
/* using the intrinsic function SECNDS(X) which returns the */
/*  (time in seconds since midnight - X) */

	srand(time(NULL));
	int iseed = (rand() * 100) + 1;

/*  run the random number generator N times */
/*  for avoiding effects of starting value */

	int n = iseed / 2000;
	if (n <= 0) n = 100;
	if (n >= 1000) n = 1000;
	for (int i = 0; i < n; ++i)
	{
		randi(&iseed);
	}
	iseed = iseed / 3;
	iseed = (iseed * 2) + 1;
	char temp[15];
	snprintf(temp, sizeof(temp), " ISEED = %d", iseed);
	// fwrite(temp, strlen(temp), 1, file);
	return iseed;
}

// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char *trimwhitespace(char *str)
{
	char *end;

	// Trim leading space
	while(isspace((unsigned char)*str)) str++;

	if(*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace((unsigned char)*end)) end--;

	// Write new null terminator
	*(end+1) = 0;

	return str;
}

FILE *openFile(const char *file_name, const char *ext, const char *mode)
{
	//TODO - change [FILENAME_SIZE + EXTENSION_SIZE]
	char result_file_name[FILE_WITH_EXTENSION_SIZE];
	snprintf(result_file_name, FILE_WITH_EXTENSION_SIZE, "%s%s", file_name, ext);
	FILE *result = fopen(result_file_name, mode);
	if (result == NULL)
	{
		printf("Falha ao abrir o arquivo %s\n", result_file_name);
		exit(0);
	}
	return result;
}

char *getFormattedDate()
{
	time_t timer;
	char *result;
	size_t size = 40;

	result = (char *)malloc(size * sizeof(char));

	struct tm* tm_info;

	time(&timer);
	tm_info = localtime(&timer);

	strftime(result, size, "\n\n %d-%b-%Y\t\t%H hour %M min %S Sec\n\n", tm_info);

	return result;
}

void writeFormattedDate(FILE *file)
{
	char *date = getFormattedDate();
	fwrite(date, strlen(date), 1, file);
	free(date);
}

void writeTotalTime(FILE *file, const time_t *begin, time_t *end)
{
	writeFormattedDate(file);
	time(end);
	double total_time = difftime(*end, *begin);
	char buffer[50];
	snprintf(buffer, sizeof(buffer), "\n\n Total CPU time elapsed in seconds : %.f", total_time);
	fwrite(buffer, strlen(buffer), 1, file);
	printf("\n\nType any character to continue : \n");
	getchar();
}

void insertHeader(FILE *file, const char *file_name, const int procs)
{
	char buffer[300];
	snprintf(buffer, sizeof(buffer), "\n =================================================================\n McMaille version %s     by A. Le Bail - 2006 - alb@cristal.org\n =================================================================\n\n Using generic filename : %s\n\n   Number of Processors : \t\t%d\n", VERSION, file_name, procs);
	fwrite(buffer, strlen(buffer), 1, file);
}

//==============================================================================
/*
* Open files using name from command line and standard extensions.
* The OPEN statements may need to be changed for some computers.
* Subroutine SXNM gets the generic filename from the command line.
* If nothing is found, then the user is prompted for the filename.
*/
char *programInit(int argc, char *argv[])
{
	char *file_name;
	file_name = (char *)malloc(FILENAME_SIZE * sizeof(char));

	int write_file_name = (argc != 2);
	do
	{
		if (write_file_name)
		{
			printf("\nEntry file (no extension) ??");
			scanf("%s", file_name);
			trimwhitespace(file_name);
		}
		else
		{
			snprintf(file_name, FILENAME_SIZE, "%s", argv[1]);
			write_file_name = 1;
		}

		char file_name_ext[FILE_WITH_EXTENSION_SIZE];
		snprintf(file_name_ext, sizeof(file_name_ext), "%s.dat", file_name);
		if ( access(file_name_ext, F_OK) != -1)
		{
			return file_name;
		}

		printf("\n\nThat file does not exist, try again...");
	} while (1);

}

FILE *createImpFile(char *file_name, const int procs)
{
	FILE *dat_file = openFile(file_name, ".dat", "r");
	FILE *inp_file = openFile(file_name, ".inp", "w");

	char *buffer = NULL;
	size_t buffsize = 0;
	ssize_t nread;
	while ((nread = getline(&buffer, &buffsize, dat_file)) != -1)
	{
		if (buffer[0] != '!')
		{
			fwrite(buffer, nread, 1, inp_file);
		}
	}

	free(buffer);
	fclose(dat_file);

	printf("\n\nMcMaille version %s\nDataFile : %s\n", VERSION, file_name);

	FILE *imp_file = openFile(file_name, ".imp", "w");
	insertHeader(imp_file, file_name, procs);

	fclose(inp_file);

	return imp_file;
}

//==============================================================================

void err(FILE *file, char *msg)
{
	fwrite(msg, strlen(msg), 1, file);
	exit(0);
}

//==============================================================================

void deleteFile(char *file_name, char *ext)
{
	char temp[FILE_WITH_EXTENSION_SIZE];
	snprintf(temp, sizeof(temp), "%s%s", file_name, ext);
	remove(temp);
}

//==============================================================================

int main(int argc, char *argv[])
{

	int nsys[6], ihkl[30000], nrun;
	char text[20];
	double pstartb[6],delta[3],pstart[3];
	double celpre[6],celold[6],w1[N_HKL],fm20[N_HKL],ff20[N_HKL];
	double bpar[6],cel[60000],rp[N_HKL],vgc[N_HKL],d[N_HKL];
	double ifi[N_HKL],ll[N_HKL],qo[N_HKL],ib[N_HKL];
	double afi[8],bb[8],km[N_HKL],th3[N_HKL];
	double pmi[6],pma[6],nsol[N_HKL],cncalc[N_HKL],xfom[N_HKL];
	double deltct[3],astartt[3],imn[10],im[10];
	double hw[N_HKL],hw4[N_HKL],bbb[N_HKL],fcal[N_HKL],hh[3];
	double pos[16000],yobs[16000],ycalc[16000];
	double dump[N_HKL],somega[N_HKL],theta[N_HKL],rmax0[6];
	double nha[16000],nhb[16000],jhh[3][N_HKL],irefs[N_HKL];
	double isyst[7],lll[N_HKL],ql[6],rp2[N_HKL],km2[N_HKL];
	double km3[100000],ll2[100000],id1[100000],id2[100000];

    const int procs = 1;  // 360./3.1415926

    int pressedk = 0;
    double pndat, ddq, ddt, isee, v3, v2, a, rmax2, llhkl;
    double diff2, diff, v1, del, rmin, interest, tmax, ttmax, ncells, iiseed, ntried, ntriedt;
    double nout, ntriedb;
    int ipen, icode, iseed, ncel, ncalc, c3, ibr;

    double am, ap, adelt, vm, vp, vdelt, nglob, rglob, c;
    int c2, ip;

	//TODO #OMP THREADPRIVATE(/cal/,/cal2/)

	//Search for the number of processors available
	//OMP_GET_NUM_PROCS()
	printf("\nNumber of used processors : \t\t%d\n", procs);

	int n1 = 1;
	int n2 = 2;
	double x = 0.0;

	char *file_name = programInit(argc, argv);
	FILE *imp_file = createImpFile(file_name, procs);

	time_t time_begin, time_end;

	writeFormattedDate(imp_file);
	time(&time_begin);

	FILE *inp_file = openFile(file_name, ".inp", "r");

	char *buffer = NULL;
	size_t buffsize = 0;
	ssize_t nread;
	int sread;

//.....READ PROBLEM IDENTIFICATION


	nread = getline(&buffer, &buffsize, inp_file);
	if (nread == -1) err(imp_file, "  Error reading the first line : TEXT\n");
	strcpy(text, buffer);
	fwrite(text, strlen(text), 1, imp_file);


/*
*.....READ wavelength and type of calculation
*    SLABDA = wavelength
*    ZERO   = zeropoint to be added at thr beginning
*    NGRID  if =. 1  : grid cell generation
*           if = 0  : Monte Carlo cell generation
*           if = 2  : both
*           if = 3  : black box Monte Carlo only...
*           if = -3 : black box Monte Carlo only, without triclinic...
*           if = 4  : black box Monte Carlo + grid search...
*/
	int iverb, ngrid, notric;
	double slabda, zero;

	nread = getline(&buffer, &buffsize, inp_file);
	sread = sscanf(buffer, "%lf %lf %d", &slabda, &zero, &ngrid);
	if (nread == -1 || sread == EOF) err(imp_file, "  Error reading lambda, etc, line 2\n");

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

/*C
C.....READ codes for search in crystalline systems
C
C    NSYS(n)
C     n
C     1   Cubic
C     2   Hexagonal/trigonal/Rhombohedral
C     3   Tetragonal
C     4   Orthorhombic
C     5   Monoclinic
C     6   Triclinic
C
C     if NSYS(n)=0 : no search
C     if NSYS(n)=1 : search
C     if NSYS(2)=2 : search in rhombohedral
C*/
	int nblack = 0;
	if (ngrid == 4)
	{
		nblack = 1;
	}
	if (ngrid == 4)
	{
		ngrid = 3;
	}
	FILE *new_dat_file;
	if (ngrid == 3)
	{
		char temp[100];
		new_dat_file = openFile(file_name, "-new.dat", "w+");
		fwrite(text, strlen(text), 1, new_dat_file);
		snprintf(temp, sizeof(temp), "! Wavelength, zeropoint, Ngrid\n%.6lf %.4lf 0\n! Codes for symmetry\n1 0 0 0 0 0\n", slabda, zero);
		fwrite(temp, strlen(temp), 1, new_dat_file);

		for (int i = 0; i < 6; ++i)
		{
			nsys[i] = 1;
		}
		if (notric == 1)
		{
			nsys[5] = 0;
		}
	}
	else
	{
		nread = getline(&buffer, &buffsize, inp_file);
		sread = sscanf(buffer, "%d %d %d %d %d %d", &nsys[0], &nsys[1], &nsys[2], &nsys[3], &nsys[4], &nsys[5]);
		if (nread == -1 || sread == EOF)
		{
			err(imp_file, "  Error reading symmetry code, line 3\n");
		}
	}


/*C
C.....Read the tolerated error on 2-theta values W
C          which is also the column width of the
C          columnar profile shape
C          and how many non-indexed reflections NIND
C
C     If W is given as negative, then individual W
C     values will be read later (triplets : 2-theta, I, Width)
C          moreover, Width will be multiplied by -W
C          (use W = -1 for no change...)
C*/
	nind;
	double w;
	if (ngrid == 3)
	{
		nind = 3;
		w = slabda * 0.3 / 1.54056;
		char temp[30];
		snprintf(temp, sizeof(temp), "! W, Nind\n%.3lf %d\n", w, nind);
		fwrite(temp, strlen(temp), 1, new_dat_file);
	}
	else
	{
		nread = getline(&buffer, &buffsize, inp_file);
		sread = sscanf(buffer, "%lf %d", &w, &nind);
		if (nread == -1 || sread == EOF)
		{
			err(imp_file, "  Error reading U, V, W, step\n");
		}
	}

/*C
C.....Some printing
C*/

	{
		char temp[500] = "";
		char *cur = temp, *const end = temp + sizeof(temp);
		cur += snprintf(cur, end-cur, "\n Wavelength : %.5lf  Zeropoint : %.4lf\n\n\n ===============================================================================\n", slabda, zero);
		if (ngrid == 0)	cur += snprintf(cur, end-cur, "     Monte Carlo cell generation\n");
		else if (ngrid == 1) cur += snprintf(cur, end-cur, "         Grid cell generation\n");
		else if (ngrid == 2) cur += snprintf(cur, end-cur, " Both searches - Monte Carlo AND grid\n");
		else if (ngrid == 3) cur += snprintf(cur, end-cur, "  -- Black box Monte Carlo search --\n");
		if (nblack == 1) cur += snprintf(cur, end-cur, "  Black box Monte Carlo + grid search \n");

		cur += snprintf(cur, end-cur, " ===============================================================================\n\n Width of the columnar profile shape, W  = %.4lf\n\n Max non-indexed reflections, NIND  = %d\n", w, nind);
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*C
C.....READ Min/Max parameters, volume, Rp
C...  RMI : If Rp < RMI, stop, should be the good cell
C...  RMAX : Keep a refined cell if Rp < Rmax
C...  RMAXREF : if Rp < RMAXREF, refine that cell by Monte Carlo
C
C*/

	double pmin, pmax, vmin, vmax, rmi, rmax, rmaxref;

	if (ngrid == 3) {
		pmin = 2.0;
		pmax = 30.0;
		vmin = 8.0;
		vmax = 27000.0;
		rmi = 0.02f;
		rmax = 0.15;
		rmaxref = 0.40;

		char temp[90];
		snprintf(temp, sizeof(temp), "!Pmin, Pmax, Vmin, Vmax, Rmin, Rmax, Rmaxref\n 2. 50. 8. 125000. 0.05 0.15 0.50\n");
		fwrite(temp, strlen(temp), 1, new_dat_file);
	}
	else
	{
		nread = getline(&buffer, &buffsize, inp_file);
		sread = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf ", &pmin, &pmax, &vmin, &vmax, &rmi, &rmax, &rmaxref);
		if (nread == -1 || sread == EOF)
		{
			err(imp_file, "  Error reading Pmin, Pmax, Vmin\n");
		}
	}
	if (pmin < 0.0)
	{
		nread = getline(&buffer, &buffsize, inp_file);
		sread = sscanf(buffer, "%lf %lf %lf %lf %lf %lf", &pmi[0], &pmi[1], &pmi[2], &pmi[3], &pmi[4], &pmi[5]);
		if (nread == -1 || sread == EOF)
		{
			err(imp_file, "  Error reading Pmin, Pmax, Vmin\n");
		}
	}

	{
		char temp[2] = "\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}
	if (pmin > 0.0)
	{
		for (int i = 0; i < 3; ++i) {
			int j = i + 3;
			pmi[i] = pmin;
			pma[i] = pmax;
			pmi[j] = 60.0;
			pma[j] = 120.0;
		}
	}
	else
	{
		char temp[300];
		snprintf(temp, sizeof(temp), " Min/Max a cell parameter %lf %lf\n Min/Max b cell parameter %lf %lf\n Min/Max c cell parameter %lf %lf\n Min/Max alpha cell parameter %lf %lf\n Min/Max beta cell parameter %lf %lf\n Min/Max gamma cell parameter %lf %lf\n", pmi[0], pma[0], pmi[1], pma[1], pmi[2], pma[2], pmi[3], pma[3], pmi[4], pma[4], pmi[5], pma[5]);
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*C
C  NR is test for automatic Rmax decrease
C*/

	nr = 0;
	if (rmax < 0.0) {
		nr = 1;
		rmax = -rmax;
	}

	{
		char temp[70];
		snprintf(temp, sizeof(temp), "\n  Min/Max Rp, Rmaxref %lf %lf %lf\n", rmi, rmax, rmaxref);
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*C
C.....According to NGRID, read either grid steps
C                              or Monte Carlo parameters
C                              or both
*/

	double spar, sang;
	int timlim, nruns, ntimelim[6];

	if (ngrid == 3) {
		spar = 0.02;
		sang = 0.05;

		char temp[60];
		snprintf(temp, sizeof(temp), "\n Steps on (a,b,c) and angles  %lf %lf\n", spar, sang);
		fwrite(temp, strlen(temp), 1, imp_file);

		snprintf(temp, sizeof(temp), "! Ntests, Nruns\n-100 20\n!  2-theta   Intensity\n");
		fwrite(temp, strlen(temp), 1, new_dat_file);
	}
	else
	{
		if (ngrid == 1 || ngrid == 2) {
/*C
C.....READ grid steps
C     SPAR = step on cell parameters
C     SANG = step on angles
C*/
			nread = getline(&buffer, &buffsize, inp_file);
			sread = sscanf(buffer, "%lf %lf", &spar, &sang);
			if (nread == -1 || sread == EOF)
			{
				err(imp_file, "  Error reading grid steps\n");
			}
			char temp[60];
			snprintf(temp, sizeof(temp), "\n Steps on (a,b,c) and angles  %lf %lf\n", spar, sang);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
		if (ngrid == 0 || ngrid == 2)
		{
/*		C
C  Continue up to NTIMELIM tests
C  Save parameters if Rp < Rmax
C  Make NRUNS times those NTIMELIM tests
C  If NTIMELIM is negative, |NTIMELIM| will apply to cubic
C             and |NTIMELIM|*50. for tetragonal, hexagonal,
C                 etc
C*/

			nread = getline(&buffer, &buffsize, inp_file);
			sread = sscanf(buffer, "%d %d", &timlim, &nruns);
			if (nread == -1 || sread == EOF)
			{
				err(imp_file, "  Error reading NTIMELIM\n");
			}

			if (timlim < 0.0)
			{
				timlim = -timlim;
				ntimelim[0] = timlim;
				ntimelim[1] = ntimelim[0] * 20;
				ntimelim[2] = ntimelim[1];
				ntimelim[3] = ntimelim[2] * 20;
				ntimelim[4] = ntimelim[3] * 20;
				ntimelim[5] = ntimelim[4] * 20;

				char temp[230];
				snprintf(temp, sizeof(temp), "\n N of runs  %d\n N of MC events in cubic         %d\n N of MC events in tetra/hexa    %d\n N of MC events in orthorhombic  %d\n N of MC events in monoclinic    %d\n N of MC events in triclinic     %d\n", nruns, ntimelim[0], ntimelim[1], ntimelim[3], ntimelim[4], ntimelim[5]);
				fwrite(temp, strlen(temp), 1, imp_file);
			}
			else
			{
				for (int i = 0; i < 6; ++i) {
					ntimelim[i] = timlim;
				}
				char temp[45];
				snprintf(temp, sizeof(temp), "\n N of MC events, N of runs  %d %d\n", timlim, nruns);
				fwrite(temp, strlen(temp), 1, imp_file);
			}
		}

		if (ngrid != 0 && ngrid != 1 && ngrid != 2)
		{
			char temp[40];
			snprintf(temp, sizeof(temp), "\n UNKNOWN NGRID PARAMETER : STOP\n");
			fwrite(temp, strlen(temp), 1, imp_file);
			exit(0);
		}
	}

/*C
C... Make a WARNING
C*/
	if (ngrid == 3)
	{
		printf("\n      This is the black box mode, can be long...\n");
		if(notric == 1)
		{
			printf("\n               No triclinic search.\n");
		}
	}
	printf("\n  To cancel and save, type K (capital letter) anytime\n\n");

/*C
C.....And now, read couples of 2-theta and I values
C*/

	sum_f = 0.0;
	ndat = 0;
	while ((nread = getline(&buffer, &buffsize, inp_file)) != -1)
	{
		if (w > 0.0)
		{
			sscanf(buffer, "%lf %lf", &th2[ndat], &fobs[ndat]);
		} else
		{
			sscanf(buffer, "%lf %lf %lf", &th2[ndat], &fobs[ndat], &w1[ndat]);
		}
		if (th2[ndat] >= 180.f) {
			err(imp_file, "  Error reading data angle > 180\n");
		}
		++ndat;
		if (ndat > 100) {
			char *temp = "Max data = 100 !\n";
			fwrite(temp, strlen(temp), 1, imp_file);
			exit(0);
		}
	}
	// --ndat;

/*C
C.... Verify if these are d values or 2-theta
C*/

	if (th2[1] < th2[0]) {
		char *temp = "\n    WARNING : DATA were given as d(A) values\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
		printf("%s\n", temp);
		for (int nda = 1; nda <= ndat; ++nda) {
			th2[nda - 1] = asin(slabda / (th2[nda - 1] * 2.f)) * pi;
		}
	}

	for (int nda = 1; nda <= ndat; ++nda) {
		if (w > 0.f) {
			w1[nda - 1] = w;
			if (ngrid == 3) {
				char temp[20];
				snprintf(temp, sizeof(temp), "%lf %lf", th2[nda - 1], fobs[nda - 1]);
				fwrite(temp, strlen(temp), 1, new_dat_file);
			}
		} else {
			w1[nda - 1] *= -w;
		}
		if (th2[ndat - 1] >= 180.f) {
			err(imp_file, "  Error reading data angle > 180\n");
		}

/*     Addition of the Zeropoint */

		th2[nda - 1] += zero;
		d[nda - 1] = slabda / (sin(th2[nda - 1] / pi) * 2.f);
	/* Computing 2nd power */
		double r1 = d[nda - 1];
		qo[nda - 1] = 1 / (r1 * r1);
	}
	if (ngrid == 3 && ndat >= 20) {
		ndat = 20;
	}
	if (ngrid == 3) {
		fclose(new_dat_file);
	}
	int nhkl = ndat;

/* ... NDAT10 is the max limit for the number of calculated */
/*           peak positions = 10 times the number of */
/*           observed peak positions */

	ndat10 = ndat * 10;
	int ndat2 = ndat << 1;

	for (int i = 1; i <= nhkl; ++i) {
		sum_f += fobs[i - 1];
	}
	int nmax = ndat - nind;

/* .....END OF DATA READING */

	fclose(inp_file);
	deleteFile(file_name, ".inp");

/*C
C  Output of some data
C
*/
	{
		char temp[3000] = "";
		char *cur = temp, *const end = temp + sizeof(temp);
		cur += snprintf(temp, sizeof(temp), "\n\t2-THETA\td(A)\t\tIobs\t\t\tW\n");
		for (int i = 0; i < nhkl; ++i)
		{
			cur += snprintf(cur, end-cur, "\t%.3lf \t%.4lf\t%.3lf \t%.3lf\n", th2[i], d[i], fobs[i], w1[i]);
		}
		cur += snprintf(cur, end-cur, "\n");
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/* ...  Various starting values initialized */

	double deltc, astart, deltab;
	slabda2 = slabda * slabda / 4.f;
	if (pmin > 0.f) {
		for (int i = 1; i <= 3; ++i) {
			deltct[i - 1] = 30.f;
			astartt[i - 1] = 60.f;
			delta[i - 1] = (pmax - pmin) / 2.f;
			pstart[i - 1] = pmin;
		}
		deltc = 15.f;
		astart = 90.f;
	} else {
		for (int i = 1; i <= 3; ++i) {
			int j = i + 3;
			deltct[i - 1] = (pma[j - 1] - pmi[j - 1]) / 2.f;
			astartt[i - 1] = pmi[j - 1];
			delta[i - 1] = (pma[i - 1] - pmi[i - 1]) / 2.f;
			pstart[i - 1] = pmi[i - 1];
		}
		deltc = (pma[4] - pmi[4]) / 2.f;
		astart = pmi[4];
	}
	for (int i = 1; i <= 6; ++i) {
		rmax0[i - 1] = rmax;
	}
	for (int j = 1; j <= ndat; ++j) {
		w2[j - 1] = w1[j - 1] / 2.f;
		difp[j - 1] = th2[j - 1] + w2[j - 1];
		difm[j - 1] = th2[j - 1] - w2[j - 1];
	}
/*   DELTAB = zone explored in angstrom around a good cell */
/*   DMIN acts as a lower d limit for keeping reflections */
	deltab = .02f;
/*      DMIN=SLABDA/(2.*SIN(TH2(NHKL)/PI))-DELTAB */
	dmin = d[nhkl - 1] - deltab;
/*  Dmax values will help to determine max cell parameters */
/*    in black box mode */
	double dmax1 = d[0] + deltab * 2.f;
	double dmax2 = d[1] + deltab * 2.f;
	double dmax3 = d[2] + deltab * 2.f;

/*C
C  Warning on the wavelength...
C*/

	{
		char temp[250] = "";
		char *cur = temp, *const end = temp + sizeof(temp);
		cur += snprintf(cur, end-cur, "\n dmax =  %lf\n Be sure that your choice of max cell parameters\n  ensures the exploration of all possibilities.\n\n", d[0]);
		if (dmax1 > 30.0)
		{
			cur += snprintf(cur, end-cur, "\n WARNING : dmax > 30 A.\n Divide the wavelength by 10 and try again...\n and then, multiply the cell parameters by 10.\n\n");
		}
		cur += snprintf(cur, end-cur, "\n   Dmin = %lf\n", dmin);
		fwrite(temp, strlen(temp), 1, imp_file);
		printf("%s\n", temp);
	}

	dmin = 1 / (dmin * dmin);
/*   DELTAD = zone explored in 1/100 of degrees around a good cell */
	double deltad = 0.2;
/*  IGC = number of retained cells, IGM = max IGC */
/*  IREF = code for refining the best cell if it had */
/*         not Rp < Rmin */
/*  IGT = total number of cells satisfying to Rp < Rmax */
/*        including multiple identical cells */
	int igc = 0;
	double igt = 0.0;
	int iref = 0;
	int igm = 10000;
	int ihr = 0;
	int nruns2 = 1;
	double rpsmall = 1.0;

/*  ESCAPE : a value of 0.15 means that in 15% of the tests, */
/*            a parameter change may be accepted even if that */
/*            change does not lead to any Rp or number of */
/*            indexed reflections improvement */

	double escape = 0.15;

	writeFormattedDate(imp_file);

	{
		char temp[380];
		snprintf(temp, sizeof(temp), "\n===============================================================================\n                RESULTS - RESULTS - RESULTS - RESULTS\n===============================================================================\n                EXPLORED CELL PARAMETERS AND VOLUMES:\n===============================================================================\n\n");
		fwrite(temp, strlen(temp), 1, imp_file);
	}

// *-------------------------------------------------------------------------
// *     Initialisation
// *
	iseed = espInit(imp_file);
// *
// *-------------------------------------------------------------------------
// C
// C
// C.....AND NOW : Generate the cell, either by Monte Carlo
// C               or by grid search
// C
	if (ngrid == 1) goto L700;
// C
// C...  Cell generation by Monte Carlo
// C
	printf("\nMonte Carlo search :\n");
// 	PRINT *
// 	PRINT *,'Monte Carlo search :'
// C
	if (nsys[0] == 0) goto L200;
// C
// C    Cubic case
// C
	printf("Cubic:        Rp     a        V     Nind\n");
// C

	int ifile = 1;
	double ncycles = 200.f;
	double cy = ncycles * 1.1f;
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

	{
		char temp[300] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		cur += snprintf(cur, end-cur, "\nCubic Monte Carlo search :\n Max a, V %lf %lf\n\n", pmax, vmax);
		if (iverb == 1)
		{
			cur += snprintf(cur, end-cur, "  Results in cubic, run, tests : %d %d\n===============================================================================\n Rp  Trial number   a        V    Nind Icod\n\n", nrun, ntimelim[0]);
		}
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*C
C     READ hkl Miller indices in cub.hkl
C*/

	FILE *cub_hkl_file = openFile("cub", ".hkl", "r");
	nread = getline(&buffer, &buffsize, cub_hkl_file);
	sread = sscanf(buffer, "%d", &nhkl0);
	nhkl0 = ndat * 6;
	if (nhkl0 > 400) nhkl0 = 400;
	for (int i = 1; i <= nhkl0; ++i)
	{
		getline(&buffer, &buffsize, cub_hkl_file);
		sscanf(buffer, "%d %d %d", &ihh[1 + i * 3 - 4], &ihh[2 + i * 3 - 4], &ihh[3 + i * 3 - 4]);
	}
	fclose(cub_hkl_file);

	for (nrun = 1; nrun <= nruns; ++nrun) {
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
		ncells = (int) ntimelim[0];
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

		for (ncel = 1; ncel <= ncells; ++ncel) {
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
			celpre[0] = pstart[0] + delta[0] * 2.f * randi(&iseed);
			ntried += 1.f;
			goto L104;
	L103:
			del = deltab * (1.f - ntriedb / cy);
			celpre[0] = pstartb[0] + del * (randi(&iseed) - .5f) * 2.f;
			ntriedb += 1.f;
	L104:
			celpre[1] = celpre[0];
			celpre[2] = celpre[0];
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L105: */
					al[i + j * 3 - 4] = 0.f;
				}
			}
			dcell(celpre, al, &v1);
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
			calcul1(&diff, &diff2);
			if (nmx > ndat10) {
				ntried += -1;
				goto L102;
			}
			if (ntriedb != 0.f) {
				goto L114;
			}

	/* ... Here are the 2 criteria for selecting a cell to be */
	/* ...      "refined" by Monte Carlo (NCYCLES events) */
	/* ... Rp value satisfying (<0.5) ??? or enough hkl explained ??? */

			if (lhkl >= nmax) {
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
			if (lhkl < nmax) {
				goto L117;
			}
	L114:
			if (diff <= rmax) {
				llhkl = lhkl;
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
			ipen = ndat - llhkl;
			if (ipen > nind) {
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
							char temp[50];
							snprintf(temp, sizeof(temp), "  Rmax reduced by 5%%, now Rmax = %lf", rmax0[0]);
							fwrite(temp, strlen(temp), 1, imp_file);
							printf("%s\n", temp);
						}
					}
				}
			}

			if (igc > 10000) {
				char temp[40];
				snprintf(temp, sizeof(temp), "   More than 10000 good cells = STOP");
				fwrite(temp, strlen(temp), 1, imp_file);
				printf("%s\n", temp);
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
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L140: */
					al[i + j * 3 - 4] = 0.f;
				}
			}
			dcell(celpre, al, &v1);

	/* $OMP CRITICAL(STORE2) */

			calcul2(&diff, ihkl, th3, &ncalc, &igc);
			km[igc - 1] = llhkl;
			km2[igc - 1] = lhkl;
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
			supcel(&lhkl, ihkl, cel, &igc, vgc, &c3);
			brav(&lhkl, ihkl, &ibr);
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
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L110: */
					al[i + j * 3 - 4] = 0.f;
				}
			}

	/* $OMP CRITICAL(FOUND) */

			if (rp[igc - 1] < rmi) {
				++interest;
				char temp[90];
				char *temp2 = "\n YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin !\n\n";
				snprintf(temp, sizeof(temp), "\n===============================================================================\n\n");
				fwrite(temp, strlen(temp), 1, imp_file);
				fwrite(temp2, strlen(temp2), 1, imp_file);
				//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
				printf("%.3lf %.4lf %.1lf %d\n%s\n", rmax, a, v2, ipen, temp2);

		/* ... Refine that cell */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
					ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.f * ddq);
					ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
				}
				iref = 1;
				goto L197;
			} else {

	/*  Anyway, calculate the M20 and F20 values */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
					ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.f * ddq);
					ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
				}
			}

	/* Test if cell already found */

			if (igc > 1) {
				int i4 = igc - 1;
				for (int i = 1; i <= i4; ++i) {
					if (ifi[i - 1] != ifile) {
						goto L118;
					}
					double vdelt = vgc[igc - 1] / 300.f;
					double vp = vgc[igc - 1] + vdelt;
					double vm = vgc[igc - 1] - vdelt;
					if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
						goto L118;
					}
					++nsol[i - 1];
					if (rp[igc - 1] < rp[i - 1]) {
						if (isee == 1) {
							//1115
							printf("%.3lf %.4lf %.1lf %d\n", rmax, a, v2, ipen);
						}
						km[i - 1] = km[igc - 1];
						vgc[i - 1] = vgc[igc - 1];
						rp[i - 1] = rp[igc - 1];
						cel[i * 6 - 6] = cel[igc * 6 - 6];
						cel[i * 6 - 5] = cel[igc * 6 - 5];
						cel[i * 6 - 4] = cel[igc * 6 - 4];
					}
					--igc;
					if (nsol[i - 1] > 5) {
						ntried = tmax + 1.f;
						++nout;
					}
					goto L119;
		L118:
					;
				}
				if (iverb == 1) {
					char temp[30];
					snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.1lf %d %d\n", rmax, ntried, a, v2, ipen, icode);
					fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printf("%.3lf %.4lf %.1lf %d\n", rmax, a, v2, ipen);
				}
		L119:
				;
			} else {
				if (iverb == 1) {
					char temp[30];
					snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.1lf %d %d\n", rmax, ntried, a, v2, ipen, icode);
					fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printf("%.3lf %.4lf %.1lf %d\n", rmax, a, v2, ipen);
				}
			}
	L197:
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
		killk(&pressedk);
		if (rmin == rmax) {
			goto L198;
		}
		if (iverb == 1) {
			char temp[80];
			snprintf(temp, sizeof(temp), "\nBest result : a=%lf V=%lf Rp=%lf\n\n", bpar[0], v3, rmin);
			fwrite(temp, strlen(temp), 1, imp_file);
			writeFormattedDate(imp_file);
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
		printf("Hexagonal:    Rp     a       c        V     Nind\n");
		rpsmall = 1.f;
	} else {
		printf("Rhombohedral: Rp     a       c        V     Nind\n");
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
		for (int i = 1; i <= 3; ++i) {
			pmi[i - 1] = pmin;
			delta[i - 1] = (pma[i - 1] - pmi[i - 1]) / 2.f;
	/* L223: */
			pstart[i - 1] = pmi[i - 1];
		}
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 4e3f) {
			vmax = 4e3f;
		}
		ntimelim[1] = vmax * 5.f;
	}

	{
		char temp[100];
		snprintf(temp, sizeof(temp), "\nHexagonal/Trigonal/Rhomboedral Monte Carlo search :\n Max(a,c), V  %lf %lf %lf\n\n", pma[0], pma[2], vmax);
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	if (iverb == 1) {
		char temp[250] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		if (nsys[1] == 1) {
			cur += snprintf(cur, end-cur, "  Results in hexagonal, run, tests : %d %d\n", nrun, ntimelim[1]);
		}
		else if (nsys[1] == 2) {
			cur += snprintf(cur, end-cur, "  Results in rhombohedral, run, tests : %d %d\n", nrun, ntimelim[1]);
		}
		cur += snprintf(cur, end-cur, "===============================================================================\n Rp  Trial number    a      c        V  Nind Icod\n\n");
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	FILE *hex_rho_hkl_file;
	if (nsys[1] == 2)
	{
		hex_rho_hkl_file = openFile("rho", ".hkl", "r");
	}
	else
	{
		hex_rho_hkl_file = openFile("hex", ".hkl", "r");
		ifile = 7;
	}
	nread = getline(&buffer, &buffsize, hex_rho_hkl_file);
	sread = sscanf(buffer, "%d", &nhkl0);
	nhkl0 = ndat * 12;
	if (nsys[1] == 2)
	{
		if (nhkl0 > 600) nhkl0 = 600;
	}
	else
	{
		if (nhkl0 > 800) nhkl0 = 800;
	}

	for (int i = 1; i <= nhkl0; ++i)
	{
		getline(&buffer, &buffsize, hex_rho_hkl_file);
		sscanf(buffer, "%d %d %d", &ihh[1 + i * 3 - 4], &ihh[2 + i * 3 - 4], &ihh[3 + i * 3 - 4]);
	}
	fclose(hex_rho_hkl_file);


	for (nrun = 1; nrun <= nruns; ++nrun) {
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
		ncells = (int) ntimelim[1];
		iiseed = 0;
		ntried = 0.;
		ntriedt = 0.;
		nout = 0;
		celpre[0] = pstart[0] + delta[0] * 2.f * randi(&iseed);
		celpre[1] = celpre[0];
		celpre[2] = pstart[2] + delta[2] * 2.f * randi(&iseed);
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

		for (ncel = 1; ncel <= ncells; ++ncel) {
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
			if (randi(&iseed) > .5f) {
				ip = 1;
			}
			celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.f * randi(&iseed);
			ntried += 1.f;
			goto L204;
	L203:
			del = deltab * (1.f - ntriedb / cy);
			int i = 3;
			if (randi(&iseed) > .5f) {
				i = 1;
			}
			celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.f;
			ntriedb += 1.f;
	L204:
			celpre[1] = celpre[0];
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L205: */
					al[i + j * 3 - 4] = 0.f;
				}
			}
			dcell(celpre, al, &v1);
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
			calcul1(&diff, &diff2);
			if (nmx > ndat10) {
				ntried += -1;
				goto L202;
			}
			if (ntriedb != 0.f) {
				goto L214;
			}

	/* ... Rp value satisfying ??? */

			if (diff < rglob || lhkl > nglob) {
				rglob = diff;
				nglob = lhkl;
				celold[ip - 1] = celpre[ip - 1];
			}
			if (lhkl >= nmax) {
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
			if (lhkl < nmax) {
				goto L217;
			}
	L214:
			if (diff <= rmax) {
				llhkl = lhkl;
				rmax = diff;
				rmax2 = diff2;
				a = celpre[0];
				c = celpre[2];
				v2 = v1;
				if (diff < rmin) {
					rmin = diff;
					bpar[0] = a;
					bpar[2] = c;
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
				celold[ip - 1] = c;
			}
			ntriedb = 0.f;
			if (rmax >= rmax0[1]) {
				goto L217;
			}
			if (rmax2 >= .15f) {
				goto L217;
			}
			ipen = ndat - llhkl;
			if (ipen > nind) {
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
							char temp[50];
							snprintf(temp, sizeof(temp), "  Rmax reduced by 5%%, now Rmax = %lf", rmax0[1]);
							fwrite(temp, strlen(temp), 1, imp_file);
							printf("%s\n", temp);
						}
					}
				}
			}

			if (igc > 10000) {
				char temp[37];
				snprintf(temp, sizeof(temp), "   More than 10000 good cells = STOP");
				fwrite(temp, strlen(temp), 1, imp_file);
				printf("%s\n", temp);
				--igc;
				++interest;
	/*      GO TO 5000 */
			}
			cel[igc * 6 - 6] = a;
			cel[igc * 6 - 5] = a;
			cel[igc * 6 - 4] = c;
			cel[igc * 6 - 3] = 90.f;
			cel[igc * 6 - 2] = 90.f;
			cel[igc * 6 - 1] = 120.f;

	/* $OMP END CRITICAL(STORE1) */

	/* ... Check for supercell */

			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = c;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L240: */
					al[i + j * 3 - 4] = 0.f;
				}
			}
			dcell(celpre, al, &v1);

	/* $OMP CRITICAL(STORE2) */

			calcul2(&diff, ihkl, th3, &ncalc, &igc);
			km[igc - 1] = llhkl;
			km2[igc - 1] = lhkl;
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
			supcel(&lhkl, ihkl, cel, &igc, vgc, &c2);
			brav(&lhkl, ihkl, &ibr);
			ib[igc - 1] = ibr;
			a = cel[igc * 6 - 6];
			cel[igc * 6 - 5] = a;
			c = cel[igc * 6 - 4];
			v2 = vgc[igc - 1];

	/* $OMP END CRITICAL(STORE2) */

	/* ... Check for interesting result */

	/*      IF(INTEREST.GE.1)GO TO 296 */
			indic = 2;
			bb[2] = a;
			bb[3] = a;
			bb[4] = c;
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
			celpre[2] = c;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L210: */
					al[i + j * 3 - 4] = 0.f;
				}
			}

	/* $OMP CRITICAL(FOUND) */

			if (rp[igc - 1] < rmi) {
				++interest;
				char temp[90];
				char *temp2 = "\n YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin !\n\n";
				snprintf(temp, sizeof(temp), "\n===============================================================================\n\n");
				fwrite(temp, strlen(temp), 1, imp_file);
				fwrite(temp2, strlen(temp2), 1, imp_file);
				//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
				printf("%.3lf %.4lf %.4lf %.1lf %d\n%s\n", rmax, a, c, v2, ipen, temp2);

		/* ... Refine that cell */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
					ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
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

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.f * ddq);
					ff20[igc - 1] = 20.f / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.f * ddq);
					ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
				}
			}

	/* Test if cell already found */

			if (igc > 1) {
				int i4 = igc - 1;
				for (i = 1; i <= i4; ++i) {
					if (ifi[i - 1] != ifile) {
						goto L218;
					}
					vdelt = vgc[igc - 1] / 300.f;
					vp = vgc[igc - 1] + vdelt;
					vm = vgc[igc - 1] - vdelt;
					if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
						goto L218;
					}
					adelt = cel[igc * 6 - 6] / 500.f;
					ap = cel[igc * 6 - 6] + adelt;
					am = cel[igc * 6 - 6] - adelt;
					if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
						goto L218;
					}
					++nsol[i - 1];
					if (rp[igc - 1] < rp[i - 1]) {
						if (isee == 1) {
							printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
						}
						km[i - 1] = km[igc - 1];
						vgc[i - 1] = vgc[igc - 1];
						rp[i - 1] = rp[igc - 1];
						cel[i * 6 - 6] = cel[igc * 6 - 6];
						cel[i * 6 - 5] = cel[igc * 6 - 5];
						cel[i * 6 - 4] = cel[igc * 6 - 4];
					}
					--igc;
					if (nsol[i - 1] > 5) {
						ntried = tmax + 1.f;
						++nout;
					}
					goto L219;
		L218:
					;
				}
				if (iverb == 1) {
					char temp[40];
					snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, c, v2, ipen, icode);
					fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
				}
		L219:
				;
			} else {
				if (iverb == 1) {
					char temp[40];
					snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, c, v2, ipen, icode);
					fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
				}
			}
	L297:

	/* $OMP END CRITICAL(FOUND) */

	L217:
			rmax = rmaxref;
			if (randi(&iseed) > escape) {
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
		killk(&pressedk);
		if (rmin == rmax) {
			goto L298;
		}
		if (iverb == 1) {
			char temp[80];
			snprintf(temp, sizeof(temp), "\nBest result : a=%lf Rp=%lf\nBest result : c=%lf V=%lf\n\n", bpar[0], rmin, bpar[2], v3);
			fwrite(temp, strlen(temp), 1, imp_file);
			writeFormattedDate(imp_file);
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
L700:
L5000:

//==============================================================================
	free(file_name);
	fclose(imp_file);

	return 0;
}