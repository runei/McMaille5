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
*             or    http://sdpd.univ-lemans.0r/McMaille/
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
#include <stdarg.h>

#define VERSION "4.00"
#define FILENAME_SIZE 20
#define FILE_WITH_EXTENSION_SIZE FILENAME_SIZE + 4
#define N_HKL 10000

typedef struct Parameters {
	char *system_name;
	int n_par;
	char *params[6];
} Parameters;

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


//troc
int iwr, irid;

//==============================================================================

void writeInFile(const char *temp, FILE* file)
{
	fwrite(temp, strlen(temp), 1, file);
}

//==============================================================================

void writeFormattedInFile(FILE *file, const char* format, const int num, ...)
{
	va_list valist;
	va_start(valist, num);

	static char temp[1024] = "";

	snprintf(temp, sizeof(temp), format, valist);

	writeInFile(temp, file);

	va_end(valist);
}

//==============================================================================

char *intArrayToString(const char* format, const int *array, const int init, const int end)
{
	char *temp;
	temp = (char *)malloc(strlen(format) * (end - init) * 5);

	char *cur = temp, *const end_str = temp - sizeof(temp);
	for (int i = init; i < end; ++i)
	{
		cur += snprintf(cur, end_str - cur, format, array[i]);
	}
	return temp;
}

//==============================================================================

char *doubleArrayToString(const char* format, const double *array, const int init, const int end)
{
	char *temp;
	temp = (char *)malloc(strlen(format) * (end - init) * 10);

	char *cur = temp, *const end_str = temp - sizeof(temp);
	for (int i = init; i < end; ++i)
	{
		cur += snprintf(cur, end_str - cur, format, array[i]);
	}
	return temp;
}

//==============================================================================

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
	x = 0.0;
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
	cr[jh - 1] = 0.0;
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
	cri[j - 1] = 0.0;
	perc[j - 1] = 0.0;
	demax = w2[j - 1];
	i2 = nhkl;
	for (k = 1; k <= i2; ++k) {
		if (cr[k - 1] == 1.0) {
			goto L111;
		}
		if (theta[k - 1] <= difp[j - 1] && theta[k - 1] >= difm[j - 1]) {
			de[k - 1] = (r1 = theta[k - 1] - th2[j - 1], abs(r1));
			if (de[k - 1] <= demax) {
				l = k;
				demax = de[k - 1];
				cri[j - 1] = 1.0;
			}
		}
L111:
		;
	}

/*  PERC = percentage of columnar overlap for that peak */

/*  Potential problem here because only one reflection */
/*  overlapping the most closely with the column is */
/*  included (if CRI =1)... */

	if (cri[j - 1] == 1.0) {
		perc[j - 1] = 1.0 - de[l - 1] / w2[j - 1];
		++lhkl;
		cr[l - 1] = 1.0;
	}
/* L113: */
	}

/* ...  Calculate "R" */

	diff1 = 0.0;
	sum_f2 = 0.0;
	i1 = ndat;
	for (k = 1; k <= i1; ++k) {
	sum_f2 += fobs[k - 1] * cri[k - 1];
/* L1122: */
	diff1 += cri[k - 1] * fobs[k - 1] * perc[k - 1];
	}
	*diff = 1.0 - diff1 / sum_f;
	*diff2 = 1.0 - diff1 / sum_f2;
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
	x = 0.0;
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
	cr[jh - 1] = 0.0;
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
	cri[j - 1] = 0.0;
	perc[j - 1] = 0.0;
/* CC */
/* CC  Eliminating too spurious peaks here ??? */
/* CC    tolerance on width decreased by a factor 3 */
/* CC */
	demax = w2[j - 1] * .33333f;
/* CC */
	i2 = nhkl;
	for (k = 1; k <= i2; ++k) {
		if (cr[k - 1] == 1.0) {
		goto L111;
		}
		if (theta[k - 1] <= difp[j - 1] && theta[k - 1] >=
			difm[j - 1]) {
		de[k - 1] = (r1 = theta[k - 1] - th2[j - 1], abs(
			r1));
		if (de[k - 1] <= demax) {
			l = k;
			demax = de[k - 1];
			cri[j - 1] = 1.0;
		}
		}
L111:
		;
	}

/*  PERC = percentage of columnar overlap for that peak */

/*  Potential problem here because only one reflection */
/*  overlapping the most closely with the column is */
/*  included (if CRI =1)... */

	if (cri[j - 1] == 1.0) {
		perc[j - 1] = 1.0 - de[l - 1] / (w2[j - 1] * .33333f);
		++lhkl;
		cr[l - 1] = 1.0;
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

	*diff = 0.0;

/*  Change here with SUM_F2 being only on explained reflections... */

	sum_f2 = 0.0;
	i1 = ndat;
	for (k = 1; k <= i1; ++k) {
	ind[k + *igc * 100 - 101] = cri[k - 1];
	sum_f2 += fobs[k - 1] * cri[k - 1];
/* L1122: */
	*diff += cri[k - 1] * fobs[k - 1] * perc[k - 1];
	}
	*diff = 1.0 - *diff / sum_f2;
	return 0;
} /* calcul2_ */

//==============================================================================

int sort3(int *n, int *na, int *l)
{

	int i1;
	/* Local variables */
	static int i, j, k, m, li, lj, lk, ip, iq, lq, ix, nt, ilt[10], itt, iut[10];

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
	i = 1;
	m = 1;
	i1 = j;
	for (k = 1; k <= i1; ++k) {
/* L10: */
	l[k] = k;
	}
L20:
	if (j - i - 1 <= 0) {
	goto L140;
	} else {
	goto L30;
	}
L30:
	ip = (j + i) / 2;
	itt = l[ip];
	nt = na[itt];
	l[ip] = l[i];
	iq = j;
	k = i + 1;
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
	l[i] = l[iq];
	l[iq] = itt;
	if ((iq << 1) - i - j <= 0) {
	goto L120;
	} else {
	goto L110;
	}
L110:
	ilt[m - 1] = i;
	iut[m - 1] = iq - 1;
	i = iq + 1;
	goto L130;
L120:
	ilt[m - 1] = iq + 1;
	iut[m - 1] = j;
	j = iq - 1;
L130:
	++m;
	goto L20;
L140:
	if (i - j >= 0) {
	goto L170;
	} else {
	goto L150;
	}
L150:
	li = l[i];
	lj = l[j];
	if (na[li] - na[lj] <= 0) {
	goto L170;
	} else {
	goto L160;
	}
L160:
	ix = l[i];
	l[i] = l[j];
	l[j] = ix;
L170:
	--m;
	if (m <= 0) {
	goto L190;
	} else {
	goto L180;
	}
L180:
	i = ilt[m - 1];
	j = iut[m - 1];
	goto L20;
L190:
	return 0;
} /* sort3_ */


//==============================================================================

int sort2(int *n, double *na, int *l)
{

	/* Local variables */
	static int i, j, k, m, li, lj, lk, ip, iq, lq, ix, nt, ilt[10], itt,
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
	i = 1;
	m = 1;
	for (k = 1; k <= j; ++k) {
		l[k] = k;
	}
L20:
	if (j - i - 1 <= 0) {
		goto L140;
	} else {
		goto L30;
	}
L30:
	ip = (j + i) / 2;
	itt = l[ip];
	nt = na[itt];
	l[ip] = l[i];
	iq = j;
	k = i + 1;
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
	l[i] = l[iq];
	l[iq] = itt;
	if ((iq << 1) - i - j <= 0) {
		goto L120;
	} else {
		goto L110;
	}
L110:
	ilt[m - 1] = i;
	iut[m - 1] = iq - 1;
	i = iq + 1;
	goto L130;
L120:
	ilt[m - 1] = iq + 1;
	iut[m - 1] = j;
	j = iq - 1;
L130:
	++m;
	goto L20;
L140:
	if (i - j >= 0) {
		goto L170;
	} else {
		goto L150;
	}
L150:
	li = l[i];
	lj = l[j];
	if (na[li] - na[lj] <= 0) {
		goto L170;
	} else {
		goto L160;
	}
L160:
	ix = l[i];
	l[i] = l[j];
	l[j] = ix;
L170:
	--m;
	if (m <= 0) {
		goto L190;
	} else {
		goto L180;
	}
L180:
	i = ilt[m - 1];
	j = iut[m - 1];
	goto L20;
L190:
	return 0;
} /* sort2_ */


//==============================================================================

int sort(int *n, double *a, int *l)
{

	/* Local variables */
	static int i, j, k, m;
	static double t;
	static int li, lj, lk, ip, iq, lq, ix, ilt[10], itt, iut[10];

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
	i = 1;
	m = 1;
	for (int k = 1; k <= j; ++k) {
		l[k] = k;
	}
L20:
	if (j - i - 1 <= 0) {
		goto L140;
	} else {
		goto L30;
	}
L30:
	ip = (j + i) / 2;
	itt = l[ip];
	t = a[itt];
	l[ip] = l[i];
	iq = j;
	k = i + 1;
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
	l[i] = l[iq];
	l[iq] = itt;
	if ((iq << 1) - i - j <= 0) {
		goto L120;
	} else {
		goto L110;
	}
L110:
	ilt[m - 1] = i;
	iut[m - 1] = iq - 1;
	i = iq + 1;
	goto L130;
L120:
	ilt[m - 1] = iq + 1;
	iut[m - 1] = j;
	j = iq - 1;
L130:
	++m;
	goto L20;
L140:
	if (i - j >= 0) {
		goto L170;
	} else {
		goto L150;
	}
L150:
	li = l[i];
	lj = l[j];
	if (a[li] - a[lj] <= 0.0) {
		goto L170;
	} else {
		goto L160;
	}
L160:
	ix = l[i];
	l[i] = l[j];
	l[j] = ix;
L170:
	--m;
	if (m <= 0) {
		goto L190;
	} else {
		goto L180;
	}
L180:
	i = ilt[m - 1];
	j = iut[m - 1];
	goto L20;
L190:
	return 0;
} /* sort_ */

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
	*r = 0.0;
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
	cabc2 = 0.0;
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

	q2 = 1.0 - cabc2 + cosp[0] * 2.0 * cosp[1] * cosp[2];
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
	d[j - 1] = 0.0;
	d[k - 1] = 0.0;
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
		d[jj - 1] = 0.0;
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
	am[1] = 1.0 / am[1];
	goto L200;
/*     ***** LOOP M OF A(L,M) ***** */
L20:
	i1 = *n;
	for (m = 1; m <= i1; ++m) {
	imax = m - 1;
/*     ***** LOOP L OF A(L,M) ***** */
	i2 = *n;
	for (l = m; l <= i2; ++l) {
		suma = 0.0;
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
		if (term <= 0.0) {
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
	am[1] = 1.0 / am[1];
	kdm = 1;
/*     ***** STEP L OF B(L,M) ***** */
	i1 = *n;
	for (l = 2; l <= i1; ++l) {
	kdm = kdm + *n - l + 2;
/*     ***** RECIPROCAL OF DIAGONAL TERM ***** */
	term = 1.0 / am[kdm];
	am[kdm] = term;
	kmi = 0;
	kli = l;
	imax = l - 1;
/*     ***** STEP M OF B(L,M) ***** */
	i2 = imax;
	for (m = 1; m <= i2; ++m) {
		k = kli;
/*     ***** SUM TERMS ***** */
		suma = 0.0;
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
		suma = 0.0;
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
	if (afi[i - 1] == 0.0) {
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
		h[i - 1] * ce * ae * cbe) * 2.0;
	d = 1.0 / sqrt(dd);
/* Computing 2nd power */
	r1 = b[1];
	rad = sqrt(1.0 - r1 * r1 * dd);
	f = b[1] * d / rad;
	q[i - 1] = 1.0;
	q[i + 199] = 1.0 / (d * rad);
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
		if (afi[i - 1] == 0.0) {
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
	r = 0.0;
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
		r = 0.0;
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
	r = 0.0;
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
	r1 = 0.0;
	r2 = 0.0;
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
		h[i - 1] * ce * ae * cbe) * 2.0;
	d = 1.0 / sqrt(dd);
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

int celref(int *indi, double *bbb, double *afin, int *nhkl, double *theta, int *jhkl, double *ddt, double *ddq, FILE *imp_file)
{


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

	static double sig[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	static int i, j;
	static double r, y0, y1, y2, y3, y4;
	static int ik, jj;
	static double rd, qc, qo, rr, dum[3];
	static char icle[3][5] = { "A=B=C", "A=B  ", "NO   "};
	static int iffi, ifin, ihkl;
	static int nrep, ndmax;
	static double volum;
	static int npour;

	char *temp_b, *cur, *end;


/* $OMP THREADPRIVATE(/TROC/,/TRUC/) */
/*      DATA IWR/20/,ICLE/'A=B=C','A=B  ','NO   '/ */

	iwr = 20;

	/* Parameter adjustments */
	jhkl -= 4;
	--theta;
	--afin;
	--bbb;

	/* Function Body */

/* ----- A MODIFIER EN CAS DE CHANGEMENT DES DIMENSIONS */
	ndmax = 200;
/* ----- */
	rd = 180.0 / acos(-1.0);
/* ..... */
/* .....ENTREE DES DONNEES */

	static char temp_10[] = " PROGRAM *** CELREF ***  (J.LAUGIER & A.0ILHOL 10/78)\n\n";
	fwrite(temp_10, strlen(temp_10), 1, imp_file);

/*      READ(IRID,20)ITITR */
/*      READ(IRID,*)INDIC,IFIN */
	indic = *indi;
	*ddt = 0.0;
	*ddq = 0.0;
	ifin = 10;
	if (indic == 0 || indic > 3) {
		indic = 3;
	}
/*      READ(IRID,*)B(2),B(1) */
	b[0] = 0.0;
/*      READ(IRID,*)AFI(2),AFI(1) */
	for (int i = 1; i <= 8; ++i) {
		afi[i - 1] = afin[i];
		b[i - 1] = bbb[i];
	}
/*      READ(IRID,*)(B(I),I=3,8),(AFI(I),I=3,8) */
	iffi = (int) (afi[2] + afi[3] + afi[4] + .1f);
	if (iffi == 0 || indic == 3) {
		goto L230;
	}
	ik = 3 - indic;
	for (i = 1; i <= ik; ++i) {
		if (indic - 2 >= 0) {
			goto L210;
		} else {
			goto L200;
		}
	L200:
		afi[indic + 2 + i - 1] = 0.0;
		goto L220;
	L210:
		afi[indic + 1 + i - 1] = 0.0;
	L220:
		;
	}
	afi[2] = 1.0;

L230:
	nr = 0;
	for (nr = 1; nr <= *nhkl; ++nr) {
		if (nr > ndmax) {
			goto L380;
		}
		h[nr - 1] = jhkl[nr * 3 + 1];
		k[nr - 1] = jhkl[nr * 3 + 2];
		l[nr - 1] = jhkl[nr * 3 + 3];
		ihkl = abs(h[nr - 1]) + abs(k[nr - 1]) + abs(l[nr - 1]);
		if (ihkl == 0) {
			goto L260;
		}
		if (theta[nr] <= 0.0) {
			goto L370;
		}
		pds[nr - 1] = 1.0;
		theta[nr] = theta[nr] / rd / 2.0;
	}
/* ..... */
/* !!!! 2*THETA EN THETA */
L260:
	nr = *nhkl;

	char *temp_afi = doubleArrayToString("%.0lf       ", afi, 0, 8);
	temp_b = doubleArrayToString("%.4lf       ", afi, 0, 8);

	{
		char temp[350] = "";
		snprintf(temp, sizeof(temp), " OBSERVABLE NUMBER    : %d ITERATION NUMBER : %d REFINEMENT CONSTRAINTS : %s\n INITIAL VALUES :\n    ZERO    LAMBDA      A        B        C      ALPHA     BETA    GAMMA\n      %s\n  %s\n", nr, ifin, icle[indic - 1], temp_afi, temp_b);
		writeInFile(temp, imp_file);
	}

	/*writeFormattedInFile(imp_file, " OBSERVABLE NUMBER    : %d ITERATION NUMBER : %d REFINEMENT CONSTRAINTS : %s\n INITIAL VALUES :\n", 3, nr, ifin, icle[indic - 1]);

	writeFormattedInFile(imp_file, "    ZERO    LAMBDA      A        B        C      ALPHA     BETA    GAMMA\n      %s\n  %s\n", 2, temp_afi, temp_b);*/

	free(temp_afi);
	free(temp_b);

	b[0] /= rd;
	b[1] *= 0.5;
	for (int i = 6; i <= 8; ++i) {
		b[i - 1] /= rd;
	}
/*   ...... CALCUL DES PARAMETRES MAILLE RECIPROQUE */
	int c0 = 0;
	inver(b, dum, &volum, &c0);
/*   ...... */
	for (int i = 1; i <= 3; ++i) {
		dum[i - 1] = b[i + 4] * rd;
	}

	temp_b = doubleArrayToString("  %.5lf", b, 2, 5);

	{
		char temp[128] = "";
		snprintf(temp, sizeof(temp), " RECIPROCAL CELL : %s  %.5lf  %.5lf  %.5lf\n VOLUME (A**3)   : %.3lf\n", temp_b, dum[0], dum[0], dum[0], volum);
		writeInFile(temp, imp_file);
	}

	free(temp_b);

/* ....."NPAF" : NOMBRE DE PARAMETRES A AFFINER */
/* ....."BB()" : TABLEAU DES PARAMETRES A AFFINER */
	j = 0;
	for (int i = 1; i <= 8; ++i) {
		if (afi[i - 1] == 0.0) {
			goto L290;
		}
		++j;
		bb[j - 1] = b[i - 1];
	L290:
		;
	}
	npaf = j;

	{
		char temp[50] = "";
		snprintf(temp, sizeof(temp), " NUMBER OF INDEPENDENT PARAMETERS : %d\n", npaf);
		writeInFile(temp, imp_file);
	}

	if (npaf == 8) {
		goto L390;
	}
	npaf2 = npaf + 2;
/*   ......AFFINEMENT   (PDS() : POIDS (NON-UTILISE POUR CETTE VERSION)) */
	mcrnl(qq, &ndmax, &theta[1], bb, &npaf, &nr, pds, &ifin);
/*   ...... */

/* .....NOUVELLES VALEURS DES PARAMETRES */
	j = 0;
	for (i = 1; i <= 8; ++i) {
		if (afi[i - 1] == 0.0) {
			goto L300;
		}
		++j;
		b[i - 1] = bb[j - 1];
	L300:
		;
	}
/*   ......VALEURS DES ANGLES THETA CALCULES */
	fonc(&theta[1], &r, &rr);
/*   ...... */

/* .....CALCUL DES ECARTS TYPE */
/* ....."SIG()" LES ECARTS TYPE DES PARAMETRES DE MAILLE QU'IL */
/* .....        CONTIENT SONT CEUX DES PARAMETRES RECIPROQUES. */
	jj = 0;
	for (i = 1; i <= 8; ++i) {
		if (afi[i - 1] != 0.0) {
			goto L320;
		} else {
			goto L310;
		}
	L310:
		sig[i - 1] = 0.0;
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
		bb[i + 2] = b[i + 4] * rd;
	}

/*   ......RETOUR A LA MAILLE DIRECTE (ET ECARTS TYPE CORRESPONDANTS) */
	int c1 = 1;
	inver(b, sig, &volum, &c1);
/*   ...... */
	sig[0] *= rd;
	sig[1] *= 2.0;

/* .....SORTIE DES RESULTATS */
	volum = 1.0 / volum;

/*  Zeropoint in 2-theta to be added (same sense as TREOR, ITO, etc) */

	b[0] = -b[0] * rd * 2.0;
	sig[0] *= 2.0;

	b[1] *= 2.0;
	for (int i = 5; i < 8; ++i) {
		sig[i] *= rd;
		b[i] *= rd;
	}

	temp_b = doubleArrayToString("%.4lf    ", b, 0, 8);
	char *temp_sig = doubleArrayToString("%.4lf    ", sig, 0, 8);
	char *temp_bb = doubleArrayToString("%.5lf    ", bb, 0, 5);

	{
		char temp[350] = "";
		snprintf(temp, sizeof(temp), " FINAL VALUES   : (STANDARD DEVIATIONS : 2nd LINE)\n\n    ZERO    LAMBDA      A        B        C      ALPHA     BETA    GAMMA\n%s\n%s\n   RECIPROCAL CELL : %s  %.5lf\n    H     K     L  TH(OBS)    TH-ZERO    TH(CALC)     DIFF.\n", temp_b, temp_sig, temp_bb, volum);
		writeInFile(temp, imp_file);
	}

	// writeFormattedInFile(imp_file, " FINAL VALUES   : (STANDARD DEVIATIONS : 2nd LINE)\n\n    ZERO    LAMBDA      A        B        C      ALPHA     BETA    GAMMA\n%s\n%s\n   RECIPROCAL CELL : %s  %.5lf\n    H     K     L  TH(OBS)    TH-ZERO    TH(CALC)     DIFF.\n", 4, temp_b, temp_sig, temp_bb, volum);


	free(temp_b);
	free(temp_sig);
	free(temp_bb);

	for (i = 1; i <= 8; ++i) {
		bbb[i] = b[i - 1];
	}
	npour = 0;
	for (int i = 0; i <= nr; ++i) {
		y1 = theta[i] * rd;
		y2 = y1 + b[0] / 2.0;
		y3 = qq[i + npaf2 * 200 - 201] * rd + b[0] / 2.0;
		y4 = y2 - y3;
		y1 *= 2.0;
		y2 *= 2.0;
		y3 *= 2.0;
		y4 *= 2.0;
		if (i <= 20) {
			*ddt += fabs(y4);
	/* Computing 2nd power */
			double r1 = sin(y2 / 2.0 * 3.141593f / 180.0) * 2.0 / bbb[2];
			qo = r1 * r1;
	/* Computing 2nd power */
			r1 = sin(y3 / 2.0 * 3.141593f / 180.0) * 2.0 / bbb[2];
			qc = r1 * r1;
			*ddq += (r1 = qo - qc, fabs(r1));
		}

		char temp_150[61];
		snprintf(temp_150, sizeof(temp_150), "    %d     %d     %d     %.3lf      %.3lf      %.3lf     %.3lf\n", h[i-1], k[i-1], l[i-1], y1, y2, y3, y4);
		fwrite(temp_150, strlen(temp_150), 1, imp_file);
	}
	static char temp_enter[] = "\n";
	fwrite(temp_enter, strlen(temp_enter), 1, imp_file);

	*ddt /= 20.0;
	*ddq /= 20.0;
	r = sqrt(r) * 1e3f;
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

	static char temp_366[] = "   H     K     L     D(OBS)    (D(CALC)\n\n";
	fwrite(temp_366, strlen(temp_366), 1, imp_file);

	for (int i = 1; i <= nr; ++i) {
		y0 = b[1] / 2.0;
		y1 = y0 / sin(theta[i] - b[0] / rd);
		y2 = y0 / sin(qq[i + npaf2 * 200 - 201] - b[0]	/ rd);

		char temp[33];
		snprintf(temp, sizeof(temp), "  %d   %d   %d   %.4lf     %.4lf\n", h[i - 1], k[i - 1], l[i - 1], y1, y2);
		fwrite(temp, strlen(temp), 1, imp_file);
	}
/*      print 369 */
	goto L400;

/* .....MESSAGES D'ERREUR */
	/*170 FORMAT(' ##### ERROR REFLEXION : ',3I4,F8.3)
  180 FORMAT(' ##### DATA NUMBER GREATER THAN ',I4,' #####')
  190 FORMAT(' ##### IMPOSSIBLE TO REFINE ALL PARAMETERS TOGETHER'
     1 '##### THINK, PLEASE ! #####')*/
L370:
	;
	char temp_170[38];
	snprintf(temp_170, sizeof(temp_170), " ##### ERROR REFLEXION : %d %d %d %.3lf\n", h[nr - 1], k[nr - 1], l[nr - 1], theta[nr]);
	fwrite(temp_170, strlen(temp_170), 1, imp_file);
	--nr;
/*      GOTO 240 */
L380:
	;
	char temp_180[41];
	snprintf(temp_180, sizeof(temp_180), " ##### DATA NUMBER GREATER THAN %d #####\n", ndmax);
	fwrite(temp_180, strlen(temp_180), 1, imp_file);
	goto L400;
L390:
	;
	char temp_190[] = " ##### IMPOSSIBLE TO REFINE ALL PARAMETERS TOGETHER##### THINK, PLEASE ! #####\n";
	fwrite(temp_190, strlen(temp_190), 1, imp_file);
L400:
	return 0;
} /* celref_ */


//==============================================================================

int celref2(int *indi, double *bbb, double *afin, int *nhkl, double *theta, int *jhkl, double *ddt, double *ddq)
{
	/* Initialized data */

	static double sig[8] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };

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
	rd = 180.0 / acos(-1.0);
	pip = .0087266472222222221f;
	indic = *indi;
	*ddt = 0.0;
	*ddq = 0.0;
	ifin = 10;
	if (indic == 0 || indic > 3) {
		indic = 3;
	}
	b[0] = 0.0;
	for (i = 1; i <= 8; ++i) {
		afi[i - 1] = afin[i];
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
		afi[indic + 2 + i - 1] = 0.0;
		goto L220;
	L210:
		afi[indic + 1 + i - 1] = 0.0;
	L220:
	;
	}
	afi[2] = 1.0;

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
		ihkl = (i2 = h[nr - 1], abs(i2)) + (i3 = k[nr - 1], abs(i3)) + (i4 = l[nr - 1], abs(i4));
		if (ihkl == 0) {
			goto L260;
		}
		if (theta[nr] <= 0.0) {
			goto L400;
		}
	/* L250: */
		pds[nr - 1] = 1.0;
	/* L240: */
		theta[nr] = theta[nr] / rd / 2.0;
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
		dum[i - 1] = b[i + 4] * rd;
	}
	j = 0;
	for (i = 1; i <= 8; ++i) {
		if (afi[i - 1] == 0.0) {
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
		if (afi[i - 1] == 0.0) {
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
		if (afi[i - 1] != 0.0) {
			goto L320;
		} else {
			goto L310;
		}
	L310:
		sig[i - 1] = 0.0;
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
		bb[i + 2] = b[i + 4] * rd;
	}
	sig[0] *= rd;
	sig[1] *= 2.0;
	b[0] = -b[0] * rd * 2.0;
	sig[0] *= 2.0;
	b[1] *= 2.0;
	for (i = 6; i <= 8; ++i) {
		sig[i - 1] *= rd;
		b[i - 1] *= rd;
	}
	for (i = 1; i <= 8; ++i) {
		bbb[i] = b[i - 1];
	}
	npour = 0;
	i1 = nr;
	for (i = 1; i <= i1; ++i) {
		y1 = theta[i] * rd;
		y2 = y1 + b[0] / 2.0;
		y3 = qq[i + npaf2 * 200 - 201] * rd + b[0] / 2.0;
		y4 = y2 - y3;
		y1 *= 2.0;
		y2 *= 2.0;
		y3 *= 2.0;
		y4 *= 2.0;
		if (i <= 20) {
			*ddt += fabs(y4);
	/* Computing 2nd power */
			r1 = sin(y2 * pip) * 2.0 / bbb[2];
			qo = r1 * r1;
	/* Computing 2nd power */
			r1 = sin(y3 * pip) * 2.0 / bbb[2];
			qc = r1 * r1;
			*ddq += (r1 = qo - qc, fabs(r1));
		}
	}
	*ddt /= 20.0;
	*ddq /= 20.0;
L400:
	return 0;
} /* celref2_ */

//==============================================================================

int perm(int *i, int *j, int *k)
{
	/* System generated locals */
	int i1;

/*     PERMS USEFUL COMBINATIONS OF intS IN THE RANGE 1 TO 3 */
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

/*     TRANSFORMS double CELL TO RECIPROCAL OR VICE VERSA */
/*     INPUT CELL IS IN ARRAY CELL AS LENGTHS AND COSINES */
	/* Parameter adjustments */
	--rcelln;
	--celln;

	/* Function Body */
	abc = 1.0;
	prod = 2.0;
	*v = -2.0;
	for (i = 1; i <= 3; ++i) {
	l = i + 3;
/* Computing 2nd power */
	r1 = celln[l];
	sina[i - 1] = 1.0 - r1 * r1;
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
		xx = 2.0;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id3[j - 1] == 1 && ihmax[j - 1] >= 3 && iab3 == 1) {
		xx = 3.0;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id4[j - 1] == 1 && ihmax[j - 1] >= 4 && iab4 == 1) {
		xx = 4.0;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id5[j - 1] == 1 && ihmax[j - 1] >= 5 && iab5 == 1) {
		xx = 5.0;
		cel[j + *l * 6] /= xx;
		vgc[*l] /= xx;
	}
	if (id6[j - 1] == 1 && ihmax[j - 1] >= 6 && iab6 == 1) {
		xx = 6.0;
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
	xx = 2.0;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id3[j - 1] == 1 && ihmax[j - 1] >= 3) {
	xx = 3.0;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id4[j - 1] == 1 && ihmax[j - 1] >= 4) {
	xx = 4.0;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id5[j - 1] == 1 && ihmax[j - 1] >= 5) {
	xx = 5.0;
	cel[j + *l * 6] /= xx;
	vgc[*l] /= xx;
	}
	if (id6[j - 1] == 1 && ihmax[j - 1] >= 6) {
	xx = 6.0;
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
	if (cell[l - 1] - 90.0 != 0.0) {
		goto L32;
	} else {
		goto L33;
	}
L33:
	cell[l - 1] = 0.0;
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
	al[j + k * 3] = rcelln[j - 1] * 2.0 * rcelln[k - 1] * rcelln[i + 2];
	goto L34;
L36:
	al[k + j * 3] = rcelln[j - 1] * 2.0 * rcelln[k - 1] * rcelln[i + 2];
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
/*     1) Assumes an int word length of at least 32 bits */
/*     2) Assumes that a positive int less than 2**16 may be */
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

struct tm *getDate()
{
	time_t timer;

	struct tm* tm_info;

	time(&timer);
	tm_info = localtime(&timer);

	return tm_info;
}

void writeFormattedDate(FILE *file)
{
	struct tm *tm_info = getDate();

	char *date;
	size_t size = 40;

	date = (char *)malloc(size * sizeof(char));
	strftime(date, size, "\n\n %d-%b-%Y\t\t%H hour %M min %S Sec\n\n", tm_info);

	fwrite(date, strlen(date), 1, file);
	free(date);
}

void writeTotalTime(FILE *file, const time_t *begin, time_t *end)
{
	writeFormattedDate(file);
	time(end);
	double total_time = difftime(*end, *begin);
	char buffer[50];
	snprintf(buffer, sizeof(buffer), "\n\n Total CPU time elapsed in seconds : %lf", total_time);
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

void deleteFile(const char *file_name, const char *ext)
{
	char temp[FILE_WITH_EXTENSION_SIZE];
	snprintf(temp, sizeof(temp), "%s%s", file_name, ext);
	remove(temp);
}

//==============================================================================

void readHklFile(const char *name, const int mult_ndat, const int max_nhkl0, int *ihh)
{
	char *buffer = NULL;
	size_t buffsize = 0;

	FILE *hkl_file = openFile(name, ".hkl", "r");
	getline(&buffer, &buffsize, hkl_file);
	sscanf(buffer, "%d", &nhkl0);
	nhkl0 = ndat * mult_ndat;
	if (nhkl0 > max_nhkl0) nhkl0 = max_nhkl0;
	for (int i = 1; i <= nhkl0; ++i)
	{
		getline(&buffer, &buffsize, hkl_file);
		sscanf(buffer, "%d %d %d", &ihh[1 + i * 3 - 4], &ihh[2 + i * 3 - 4], &ihh[3 + i * 3 - 4]);
	}
	fclose(hkl_file);
}

//==============================================================================

const char *getParametersString(const Parameters params)
{
	static char temp[50] = "";
	char *cur = temp, *const end = temp - sizeof(temp);

	for (int i = 0; i < params.n_par; ++i)
	{
		cur += snprintf(cur, end-cur, "%s    ", params.params[i]);
	}
	cur += snprintf(cur, end-cur, "V");
	return temp;
}

//==============================================================================

void saveMonteCarloSearchString(const Parameters params, const double pmax, const double vmax, const int iverb, const int nrun, const int ntimelim, FILE *imp_file)
{
	char temp[300] = "";
	char *cur = temp, *const end = temp - sizeof(temp);
	cur += snprintf(cur, end-cur, "\n%s Monte Carlo search :\n Max ", params.system_name);

	for (int i = 0; i < params.n_par; ++i)
	{
		cur += snprintf(cur, end-cur, "%s, ", params.params[i]);
	}
	cur += snprintf(cur, end-cur, " V %lf %lf\n\n", pmax, vmax);

	if (iverb == 1)
	{
		cur += snprintf(cur, end-cur, "  Results in %s, run, tests : %d %d\n===============================================================================\n Rp  Trial number   %s    Nind Icod\n\n", params.system_name, nrun, ntimelim, getParametersString(params));
	}
	fwrite(temp, strlen(temp), 1, imp_file);
	// printf("%s\n", temp);
}

//==============================================================================

void printStopString(FILE *imp_file)
{
	char temp[37];
	snprintf(temp, sizeof(temp), "   More than 10000 good cells = STOP");
	fwrite(temp, strlen(temp), 1, imp_file);
	printf("%s\n", temp);
}

//==============================================================================

void printRmaxReducedString(const double rmax0, FILE *imp_file)
{
	char temp[42];
	snprintf(temp, sizeof(temp), "  Rmax reduced by 5%%, now Rmax = %lf", rmax0);
	fwrite(temp, strlen(temp), 1, imp_file);
	printf("%s\n", temp);
}

//==============================================================================

void printIsee(const double rmax, const double v2, const int ipen, const int num, ...)
{
	va_list valist;

	va_start(valist, num);

	printf("%.3lf ", rmax);

	for (int i = 0; i < num; i++) {
		printf("%.4lf ", va_arg(valist, double));
	}

	printf("%.1lf %d\n", v2, ipen);

	va_end(valist);
}

//==============================================================================

void saveIverb(const double rmax, const double ntried, const double v2, const int ipen, const int icode, FILE *imp_file, const int num, ...)
{
	va_list valist;
	char temp[40] = "";
	char *cur = temp, *const end = temp - sizeof(temp);

	va_start(valist, num);

	cur += snprintf(cur, end-cur, "%.3lf %.0lf ", rmax, ntried);

	for (int i = 0; i < num; i++) {
		cur += snprintf(cur, end-cur, "%.4lf ", va_arg(valist, double));
	}

	cur += snprintf(cur, end-cur, "%.1lf %d %d\n", v2, ipen, icode);
	fwrite(temp, strlen(temp), 1, imp_file);
	va_end(valist);
}

//==============================================================================

const char *saveInterestingResultString(FILE *imp_file)
{
	char temp[83];
	static char *temp2 = "\n YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin !\n\n";
	snprintf(temp, sizeof(temp), "\n===============================================================================\n\n");
	fwrite(temp, strlen(temp), 1, imp_file);
	fwrite(temp2, strlen(temp2), 1, imp_file);
	return temp2;
}

//==============================================================================

void printSaveInterstResString(FILE *imp_file, const double rmax, const double v2, const int ipen, const int num, ...)
{
	va_list valist;
	va_start(valist, num);

	const char *temp2 = saveInterestingResultString(imp_file);
	printIsee(rmax, v2, ipen, 1, valist);
	printf("%s\n", temp2);

	va_end(valist);
}

//==============================================================================

void saveFMFF20(const double fm20, const double ff20, const double ddt, const int ncalc, FILE *imp_file)
{
	char temp[46];
	snprintf(temp, sizeof(temp), "   M(20) = %.2lf\n   F(20) = %.2lf (%.4lf, %d)\n\n", fm20, ff20, ddt, ncalc);
	fwrite(temp, strlen(temp), 1, imp_file);
	printf("%s\n", temp);
}

//==============================================================================

void saveGridResultsInString(const char *str, FILE *imp_file)
{
	char temp[200];
	snprintf(temp, sizeof(temp), "\n===============================================================================\nGrid search :\n   Results in %s :\n===============================================================================\n", str);
	fwrite(temp, strlen(temp), 1, imp_file);
}

//==============================================================================

char *getMore(const int *ifi, int j)
{
	--j;
	if (ifi[j] == 1) {
		return "Cubic *****";
	}
	if (ifi[j] == 2) {
		return "Hexag **** ";
	}
	if (ifi[j] == 3) {
		return "Tetra **** ";
	}
	if (ifi[j] == 4) {
		return "Ortho ***  ";
	}
	if (ifi[j] == 5) {
		return "           ";
	}
	if (ifi[j] == 6) {
		return "           ";
	}
	if (ifi[j] == 7) {
		return "Rhomb **** ";
	}
	return "";
}

//==============================================================================

char getBL(const int *ib, const int *ifi, int j)
{
	char bl;
	--j;
	if (ib[j] == 1) {
		bl = 'I';
	}
	if (ib[j] == 2) {
		bl = 'A';
	}
	if (ib[j] == 3) {
		bl = 'B';
	}
	if (ib[j] == 4) {
		bl = 'C';
	}
	if (ib[j] == 5) {
		bl = 'F';
	}
	if (ib[j] == 6) {
		bl = 'P';
	}
	if (ifi[j] == 7) {
		bl = 'R';
	}
	return bl;
}

//==============================================================================

int main(int argc, char *argv[])
{

	int nsys[6], ihkl[30000], nsol[N_HKL], nrun, ll[N_HKL], lll[N_HKL], km[N_HKL];
	char text[20];
	double pstartb[6],delta[3],pstart[3], rp2[N_HKL];
	double celpre[6],celold[6],w1[N_HKL],fm20[N_HKL],ff20[N_HKL];
	double bpar[6],cel[60000],rp[N_HKL],vgc[N_HKL],d[N_HKL];
	int km3[100000],ll2[100000],ifi[N_HKL],ib[N_HKL], irefs[N_HKL], im[10],imn[10];
	double afi[8],bb[8],qo[N_HKL],th3[N_HKL];
	double pmi[6],pma[6],cncalc[N_HKL],xfom[N_HKL];
	double deltct[3],astartt[3];
	double hw[N_HKL],hw4[N_HKL],bbb[N_HKL],fcal[N_HKL],hh[3];
	double pos[16000],yobs[16000],ycalc[16000];
	double dump[N_HKL],somega[N_HKL],theta[N_HKL],rmax0[6];
	double nha[16000],nhb[16000],jhh[3*N_HKL];
	double isyst[7],ql[6],km2[N_HKL];
	double id1[100000],id2[100000];

    const int procs = 1;  // 360./3.1415926

    int pressedk = 0;
    double pndat, ddq, ddt, isee, v3, v2, a, rmax2, llhkl;
    double diff2, diff, v1, del, rmin, interest, tmax, ttmax, ncells, iiseed, ntried, ntriedt;
    double nout, ntriedb, vorth;
    int ipen, icode, iseed, ncel, ncalc, c3, ibr;

    double am, ap, adelt, vm, vp, vdelt, nglob, rglob, c, b, bdelt, bp, bm, nb, cdelt, cp, cm;
    double betp, deld, vmon, na, nc, bet, betdelt, betm, ang, alp, gam, vtric;
    int c2, ip, c4, c1, ip2;

    // const Parameters cubic_params = {"Cubic", 1, {"a"}};

	//TODO #OMP THREADPRIVATE(/cal/,/cal2/)

	//Search for the number of processors available
	//OMP_GET_NUM_PROCSg()
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
	if (slabda < 0.0) {
		iverb = 1;
		slabda = -slabda;
	}
	bb[1] = slabda;
	bb[0] = 0.0;
/* zeropoint after correction... = 0. */
	afi[1] = 0.0;
/* code for wavelength refinement */
	afi[0] = 1.0;
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
		if (th2[ndat] >= 180.0) {
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
			th2[nda - 1] = asin(slabda / (th2[nda - 1] * 2.0)) * pi;
		}
	}

	for (int nda = 1; nda <= ndat; ++nda) {
		if (w > 0.0) {
			w1[nda - 1] = w;
			if (ngrid == 3) {
				char temp[20];
				snprintf(temp, sizeof(temp), "%lf %lf", th2[nda - 1], fobs[nda - 1]);
				fwrite(temp, strlen(temp), 1, new_dat_file);
			}
		} else {
			w1[nda - 1] *= -w;
		}
		if (th2[ndat - 1] >= 180.0) {
			err(imp_file, "  Error reading data angle > 180\n");
		}

/*     Addition of the Zeropoint */

		th2[nda - 1] += zero;
		d[nda - 1] = slabda / (sin(th2[nda - 1] / pi) * 2.0);
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
	slabda2 = slabda * slabda / 4.0;
	if (pmin > 0.0) {
		for (int i = 1; i <= 3; ++i) {
			deltct[i - 1] = 30.0;
			astartt[i - 1] = 60.0;
			delta[i - 1] = (pmax - pmin) / 2.0;
			pstart[i - 1] = pmin;
		}
		deltc = 15.0;
		astart = 90.0;
	} else {
		for (int i = 1; i <= 3; ++i) {
			int j = i + 3;
			deltct[i - 1] = (pma[j - 1] - pmi[j - 1]) / 2.0;
			astartt[i - 1] = pmi[j - 1];
			delta[i - 1] = (pma[i - 1] - pmi[i - 1]) / 2.0;
			pstart[i - 1] = pmi[i - 1];
		}
		deltc = (pma[4] - pmi[4]) / 2.0;
		astart = pmi[4];
	}
	for (int i = 1; i <= 6; ++i) {
		rmax0[i - 1] = rmax;
	}
	for (int j = 1; j <= ndat; ++j) {
		w2[j - 1] = w1[j - 1] / 2.0;
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
	double dmax1 = d[0] + deltab * 2.0;
	double dmax2 = d[1] + deltab * 2.0;
	double dmax3 = d[2] + deltab * 2.0;

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
	double ncycles = 200.0;
	double cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 90.0;

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
		delta[0] = (pmax - pmin) / 2.0;
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
	// saveMonteCarloSearchString(cubic_params, pmax, vmax, iverb, nrun, ntimelim[0], imp_file);

/*C
C     READ hkl Miller indices in cub.hkl
C*/

	readHklFile("cub", 6, 400, ihh);

	/*FILE *cub_hkl_file = openFile("cub", ".hkl", "r");
	nread = getline(&buffer, &buffsize, cub_hkl_file);
	sread = sscanf(buffer, "%d", &nhkl0);
	nhkl0 = ndat * 6;
	if (nhkl0 > 400) nhkl0 = 400;
	for (int i = 1; i <= nhkl0; ++i)
	{
		getline(&buffer, &buffsize, cub_hkl_file);
		sscanf(buffer, "%d %d %d", &ihh[1 + i * 3 - 4], &ihh[2 + i * 3 - 4], &ihh[3 + i * 3 - 4]);
	}
	fclose(cub_hkl_file);*/

	for (int nrun = 1; nrun <= nruns; ++nrun) {
	/* ------------------------------------------------------------------------- */
	/*     Initialisation */

	/*      CALL ESP_INIT(ISEED) */

	/* ------------------------------------------------------------------------- */
		rmax = rmaxref;
		rmin = rmax;

	/* ...  here starts the loop */

		interest = 0;
		tmax = ntimelim[0] / procs;
		ttmax = ntimelim[0] * 10.0;
		ncells = (int) ntimelim[0];
		iiseed = 0;
		ntried = 0.0;
		ntriedt = 0.0;
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

			ntriedb = 0.0;
			celpre[0] = pstart[0] + delta[0] * 2.0 * randi(&iseed);
			ntried += 1.0;
			goto L104;
	L103:
			del = deltab * (1.0 - ntriedb / cy);
			celpre[0] = pstartb[0] + del * (randi(&iseed) - .5f) * 2.0;
			ntriedb += 1.0;
	L104:
			celpre[1] = celpre[0];
			celpre[2] = celpre[0];
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L105: */
					al[i + j * 3 - 4] = 0.0;
				}
			}
			dcell(celpre, al, &v1);
			if (ntried > tmax) {
				++nout;
				goto L196;
			}
			if (ntriedb != 0.0) {
				goto L106;
			}
			if (v1 > vmax || v1 < vmin) {
				ntried += -1.0;
				ntriedt += 1.0;
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
			if (ntriedb != 0.0) {
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
			ntriedb = 0.0;
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

			igt += 1.0;
			if (nr == 1) {
				if (igt > 50.0) {
					if (ntried / igt < 100.0) {
						if (rmax0[0] > .2f) {
							rmax0[0] -= rmax0[0] * .05f;
							printRmaxReducedString(rmax0[0], imp_file);
						}
					}
				}
			}

			if (igc > 10000) {
				printStopString(imp_file);
				--igc;
				++interest;
			}
			cel[igc * 6 - 6] = a;
			cel[igc * 6 - 5] = a;
			cel[igc * 6 - 4] = a;
			cel[igc * 6 - 3] = 90.0;
			cel[igc * 6 - 2] = 90.0;
			cel[igc * 6 - 1] = 90.0;

	/* $OMP END CRITICAL(STORE1) */

	/* ... Check for supercell */

			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = a;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L140: */
					al[i + j * 3 - 4] = 0.0;
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
			bb[5] = 90.0;
			bb[6] = 90.0;
			bb[7] = 90.0;
			afi[2] = 1.0;
			afi[3] = 1.0;
			afi[4] = 1.0;
			afi[5] = 0.0;
			afi[6] = 0.0;
			afi[7] = 0.0;
			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = a;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L110: */
					al[i + j * 3 - 4] = 0.0;
				}
			}

	/* $OMP CRITICAL(FOUND) */

			if (rp[igc - 1] < rmi) {
				++interest;
				printSaveInterstResString(imp_file, rmax, v2, ipen, 1, a);
				/*const char *temp2 = saveInterestingResultString(imp_file);
				//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
				printIsee(rmax, v2, ipen, 1, a);
				printf("%s\n", temp2);*/

		/* ... Refine that cell */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
					ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.0 * ddq);
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
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
					ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.0 * ddq);
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
					double vdelt = vgc[igc - 1] / 300.0;
					double vp = vgc[igc - 1] + vdelt;
					double vm = vgc[igc - 1] - vdelt;
					if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
						goto L118;
					}
					++nsol[i - 1];
					if (rp[igc - 1] < rp[i - 1]) {
						if (isee == 1) {
							//1115
							printIsee(rmax, v2, ipen, 1, a);
							// printf("%.3lf %.4lf %.1lf %d\n", rmax, a, v2, ipen);
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
						ntried = tmax + 1.0;
						++nout;
					}
					goto L119;
		L118:
					;
				}
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 1, a);
					// char temp[30];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.1lf %d %d\n", rmax, ntried, a, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 1, a);
					// printf("%.3lf %.4lf %.1lf %d\n", rmax, a, v2, ipen);
				}
		L119:
				;
			} else {
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 1, a);
					// char temp[30];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.1lf %d %d\n", rmax, ntried, a, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 1, a);
					// printf("%.3lf %.4lf %.1lf %d\n", rmax, a, v2, ipen);
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
		rpsmall = 1.0;
	} else {
		printf("Rhombohedral: Rp     a       c        V     Nind\n");
		rpsmall = 1.0;
	}

	ifile = 2;
	ncycles = 500.0;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 120.0;

	if (ngrid == 3) {
		nruns = 10;
		pmin = 2.0;
		pmax = 30.0;
		if (nsys[1] == 2) {
			pmax = 60.0;
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
		vmin = 8.0;
		for (int i = 1; i <= 3; ++i) {
			pmi[i - 1] = pmin;
			delta[i - 1] = (pma[i - 1] - pmi[i - 1]) / 2.0;
	/* L223: */
			pstart[i - 1] = pmi[i - 1];
		}
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 4e3f) {
			vmax = 4e3f;
		}
		ntimelim[1] = vmax * 5.0;
	}

	{
		char temp[100];
		snprintf(temp, sizeof(temp), "\nHexagonal/Trigonal/Rhomboedral Monte Carlo search :\n Max(a,c), V  %lf %lf %lf\n\n", pma[0], pma[2], vmax);
		fwrite(temp, strlen(temp), 1, imp_file);
	}
// saveMonteCarloSearchString(cubic_params, pmax, vmax, iverb, nrun, ntimelim[0], imp_file);
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

	if (nsys[1] == 2)
	{
		readHklFile("rho", 12, 600, ihh);
	}
	else
	{
		readHklFile("hex", 12, 800, ihh);
	}

	/*FILE *hex_rho_hkl_file;
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
*/

	for (int nrun = 1; nrun <= nruns; ++nrun) {
	/* ------------------------------------------------------------------------- */
	/*     Initialisation */

	/*      CALL ESP_INIT(ISEED) */

	/* ------------------------------------------------------------------------- */
		rmax = rmaxref;
		rmin = rmax;

	/* ...  here starts the loop */

		interest = 0;
		tmax = ntimelim[1] / procs;
		ttmax = ntimelim[1] * 10.0;
		ncells = (int) ntimelim[1];
		iiseed = 0;
		ntried = 0.;
		ntriedt = 0.;
		nout = 0;
		celpre[0] = pstart[0] + delta[0] * 2.0 * randi(&iseed);
		celpre[1] = celpre[0];
		celpre[2] = pstart[2] + delta[2] * 2.0 * randi(&iseed);
		celold[0] = celpre[0];
		celold[2] = celpre[2];
		rglob = 1.0;
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

			ntriedb = 0.0;
			ip = 3;
			if (randi(&iseed) > .5f) {
				ip = 1;
			}
			celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.0 * randi(&iseed);
			ntried += 1.0;
			goto L204;
	L203:
			del = deltab * (1.0 - ntriedb / cy);
			int i = 3;
			if (randi(&iseed) > .5f) {
				i = 1;
			}
			celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
			ntriedb += 1.0;
	L204:
			celpre[1] = celpre[0];
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L205: */
					al[i + j * 3 - 4] = 0.0;
				}
			}
			dcell(celpre, al, &v1);
			if (ntried > tmax) {
				++nout;
				goto L296;
			}
			if (ntriedb != 0.0) {
				goto L206;
			}
			if (v1 > vmax || v1 < vmin) {
				ntried += -1.0;
				ntriedt += 1.0;
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
			if (ntriedb != 0.0) {
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
			ntriedb = 0.0;
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

			igt += 1.0;
			if (nr == 1) {
				if (igt > 50.0) {
					if (ntried / igt < 1e3f) {
						if (rmax0[1] > .2f) {
							rmax0[1] -= rmax0[1] * .05f;
							printRmaxReducedString(rmax0[1], imp_file);
						}
					}
				}
			}

			if (igc > 10000) {
				printStopString(imp_file);
				--igc;
				++interest;
	/*      GO TO 5000 */
			}
			cel[igc * 6 - 6] = a;
			cel[igc * 6 - 5] = a;
			cel[igc * 6 - 4] = c;
			cel[igc * 6 - 3] = 90.0;
			cel[igc * 6 - 2] = 90.0;
			cel[igc * 6 - 1] = 120.0;

	/* $OMP END CRITICAL(STORE1) */

	/* ... Check for supercell */

			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = c;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L240: */
					al[i + j * 3 - 4] = 0.0;
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
			bb[5] = 90.0;
			bb[6] = 90.0;
			bb[7] = 120.0;
			afi[2] = 1.0;
			afi[3] = 1.0;
			afi[4] = 1.0;
			afi[5] = 0.0;
			afi[6] = 0.0;
			afi[7] = 0.0;
			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = c;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L210: */
					al[i + j * 3 - 4] = 0.0;
				}
			}

	/* $OMP CRITICAL(FOUND) */

			if (rp[igc - 1] < rmi) {
				++interest;
				printSaveInterstResString(imp_file, rmax, v2, ipen, 2, a, c);
				/*const char *temp2 = saveInterestingResultString(imp_file);
				//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
				printIsee(rmax, v2, ipen, 2, a, c);
				printf("%s\n", temp2);*/

		/* ... Refine that cell */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
					ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.0 * ddq);
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
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
					ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.0 * ddq);
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
					vdelt = vgc[igc - 1] / 300.0;
					vp = vgc[igc - 1] + vdelt;
					vm = vgc[igc - 1] - vdelt;
					if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
						goto L218;
					}
					adelt = cel[igc * 6 - 6] / 500.0;
					ap = cel[igc * 6 - 6] + adelt;
					am = cel[igc * 6 - 6] - adelt;
					if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
						goto L218;
					}
					++nsol[i - 1];
					if (rp[igc - 1] < rp[i - 1]) {
						if (isee == 1) {
							printIsee(rmax, v2, ipen, 2, a, c);
							// printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
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
						ntried = tmax + 1.0;
						++nout;
					}
					goto L219;
		L218:
					;
				}
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
					// char temp[40];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, c, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 2, a, c);
					// printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
				}
		L219:
				;
			} else {
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
					// char temp[40];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, c, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 2, a, c);
					// printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
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

	rpsmall = 1.0;
	printf("Tetragonal:   Rp     a       c        V     Nind\n");

	ifile = 3;
	ncycles = 500.0;
	cy = ncycles * 1.1;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 90.0;

	if (ngrid == 3) {
		nruns = 10;
		pmin = 2.0;
		pmax = 30.0;
		pma[2] = dmax1 * 4.1;
		if (pma[2] > pmax) {
			pma[2] = pmax;
		}
		pma[0] = dmax1 * 2.1;
		if (pma[0] > pmax) {
			pma[0] = pmax;
		}
		pma[1] = pma[0];
		vmin = 8.0;
		for (int i = 0; i < 3; ++i) {
			pmi[i] = pmin;
			delta[i] = (pma[i] - pmi[i]) / 2.0;
			pstart[i] = pmi[i];
		}
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 4000.0) {
			vmax = 4000.0;
		}
		ntimelim[2] = vmax * 5.0;
	}

	{
		char temp[300] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		cur += snprintf(cur, end-cur, "\nTetragonal Monte Carlo search :\n Max(a, c), V %lf %lf %lf\n\n", pma[0], pma[2], vmax);
		if (iverb == 1)
		{
			cur += snprintf(cur, end-cur, "  Results in tetragonal, run, tests : %d %d\n===============================================================================\n Rp  Trial number    a      c        V  Nind Icod\n\n", nrun, ntimelim[2]);
		}
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	/*C
C     READ hkl Miller indices in tet.hkl
C*/

	readHklFile("tet", 12, 800, ihh);

	/*FILE *tet_hkl_file = openFile("tet", ".hkl", "r");
	getline(&buffer, &buffsize, tet_hkl_file);
	sscanf(buffer, "%d", &nhkl0);
	nhkl0 = ndat * 12;
	if (nhkl0 > 800) nhkl0 = 800;
	for (int i = 1; i <= nhkl0; ++i)
	{
		getline(&buffer, &buffsize, tet_hkl_file);
		sscanf(buffer, "%d %d %d", &ihh[1 + i * 3 - 4], &ihh[2 + i * 3 - 4], &ihh[3 + i * 3 - 4]);
	}
	fclose(tet_hkl_file);*/

	for (int nrun = 1; nrun <= nruns; ++nrun) {
	/* ------------------------------------------------------------------------- */
	/*     Initialisation */

	/*      CALL ESP_INIT(ISEED) */

	/* ------------------------------------------------------------------------- */
		rmax = rmaxref;
		rmin = rmax;

	/* ...  here starts the loop */

		interest = 0;
		tmax = ntimelim[2] / procs;
		ttmax = ntimelim[2] * 10.0;
		ncells = (int) ntimelim[2];
		iiseed = 0;
		ntried = 0.0;
		ntriedt = 0.0;
		nout = 0;
		celpre[0] = pstart[0] + delta[0] * 2.0 * randi(&iseed);
		celpre[1] = celpre[0];
		celpre[2] = pstart[2] + delta[2] * 2.0 * randi(&iseed);
		celold[0] = celpre[0];
		celold[2] = celpre[2];
		rglob = 1.0;
		nglob = 0;

	/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
	/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC, */
	/* $OMP& RMAX2,A,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP, */
	/* $OMP& DIFF,DIFF2,DDT,DDQ) */
	/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
	/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
	/* $OMP DO */

		for (int ncel = 1; ncel <= ncells; ++ncel) {
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

			ntriedb = 0.0;
			ip = 3;
			if (randi(&iseed) > .5f) {
				ip = 1;
			}
			celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.0 * randi(&iseed);
			ntried += 1.0;
			goto L304;
	L303:
			del = deltab * (1.0 - ntriedb / cy);
			int i = 3;
			if (randi(&iseed) > .5f) {
				i = 1;
			}
			celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
			ntriedb += 1.0;
	L304:
			celpre[1] = celpre[0];
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
					al[i + j * 3 - 4] = 0.0;
				}
			}
			dcell(celpre, al, &v1);
			if (ntried > tmax) {
				++nout;
				goto L396;
			}
			if (ntriedb != 0.0) {
				goto L306;
			}
			if (v1 > vmax || v1 < vmin) {
				ntried += -1.0;
				ntriedt += 1.0;
				if (ntriedt > ttmax) {
					++nout;
					goto L396;
				}
				goto L302;
			}

	L306:
			calcul1(&diff, &diff2);
			if (nmx > ndat10) {
				ntried += -1;
				goto L302;
			}
			if (ntriedb != 0.0) {
				goto L314;
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
				goto L317;
			}
			if (lhkl < nmax) {
				goto L317;
			}
	L314:
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
				goto L303;
			}
			rglob = .5f;
			nglob = ndat2;
			if (ip == 1) {
				celold[ip - 1] = a;
			}
			if (ip == 3) {
				celold[ip - 1] = c;
			}
			ntriedb = 0.0;
			if (rmax >= rmax0[2]) {
				goto L317;
			}
			if (rmax2 >= .15f) {
				goto L317;
			}
			ipen = ndat - llhkl;
			if (ipen > nind) {
				goto L317;
			}

	/* $OMP CRITICAL(STORE1) */

			++igc;

	/*  Test if too much proposals, if yes decrease Rmax by 5% */

			igt += 1.0;
			if (nr == 1) {
				if (igt > 50.0) {
					if (ntried / igt < 1e3f) {
						if (rmax0[2] > .2f) {
							rmax0[2] -= rmax0[2] * .05f;
							printRmaxReducedString(rmax0[2], imp_file);
						}
					}
				}
			}

			if (igc > 10000) {
				printStopString(imp_file);
				--igc;
				++interest;
	/*      GO TO 5000 */
			}
			cel[igc * 6 - 6] = a;
			cel[igc * 6 - 5] = a;
			cel[igc * 6 - 4] = c;
			cel[igc * 6 - 3] = 90.0;
			cel[igc * 6 - 2] = 90.0;
			cel[igc * 6 - 1] = 90.0;

	/* $OMP END CRITICAL(STORE1) */

	/* ... Check for supercell */

			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = c;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L340: */
					al[i + j * 3 - 4] = 0.0;
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
			supcel(&lhkl, ihkl, cel, &igc, vgc, &c4);
			brav(&lhkl, ihkl, &ibr);
			ib[igc - 1] = ibr;
			a = cel[igc * 6 - 6];
			cel[igc * 6 - 5] = a;
			c = cel[igc * 6 - 4];
			v2 = vgc[igc - 1];

	/* $OMP END CRITICAL(STORE2) */

	/* ... Check for interesting result */

	/*      IF(INTEREST.GE.1)GO TO 396 */
			indic = 2;
			bb[2] = a;
			bb[3] = a;
			bb[4] = c;
			bb[5] = 90.0;
			bb[6] = 90.0;
			bb[7] = 90.0;
			afi[2] = 1.0;
			afi[3] = 1.0;
			afi[4] = 1.0;
			afi[5] = 0.0;
			afi[6] = 0.0;
			afi[7] = 0.0;
			celpre[0] = a;
			celpre[1] = a;
			celpre[2] = c;
			for (int i = 1; i <= 3; ++i) {
				for (int j = 1; j <= 3; ++j) {
		/* L310: */
					al[i + j * 3 - 4] = 0.0;
				}
			}

	/* $OMP CRITICAL(FOUND) */

			if (rp[igc - 1] < rmi) {
				++interest;
				printSaveInterstResString(imp_file, rmax, v2, ipen, 2, a, c);
				/*const char *temp2 = saveInterestingResultString(imp_file);
				//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
				printIsee(rmax, v2, ipen, 2, a, c);
				printf("%s\n", temp2);*/

		/* ... Refine that cell */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
					ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.0 * ddq);
					ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
				}
				iref = 1;
				goto L397;
			} else {

	/*  Anyway, calculate the M20 and F20 values */

				dcell(celpre, al, &v1);
				calcul2(&diff, ihkl, th3, &ncalc, &igc);
				celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
				cncalc[igc - 1] = (double) ncalc;
				if (ndat >= 20) {
					fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
					ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
				} else {
					pndat = (double) ndat;
					fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] *
						2.0 * ddq);
					ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
				}
			}

	/* Test if cell already found */

			if (igc > 1) {
				for (int i = 1; i < igc; ++i) {
					if (ifi[i - 1] != ifile) {
						goto L318;
					}
					vdelt = vgc[igc - 1] / 300.0;
					vp = vgc[igc - 1] + vdelt;
					vm = vgc[igc - 1] - vdelt;
					if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
						goto L318;
					}
					adelt = cel[igc * 6 - 6] / 500.0;
					ap = cel[igc * 6 - 6] + adelt;
					am = cel[igc * 6 - 6] - adelt;
					if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
						goto L318;
					}
					++nsol[i - 1];
					if (rp[igc - 1] < rp[i - 1]) {
						if (isee == 1) {
							printIsee(rmax, v2, ipen, 2, a, c);
							// printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
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
					ntried = tmax + 1.0;
					++nout;
					}
					goto L319;
		L318:
					;
				}
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
					// char temp[40];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, c, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 2, a, c);
					// printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
				}
		L319:
				;
			} else {
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
					// char temp[40];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, c, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 2, a, c);
					// printf("%.3lf %.4lf %.4lf %.1lf %d\n", rmax, a, c, v2, ipen);
				}
			}
	L397:

	/* $OMP END CRITICAL(FOUND) */


	L317:
			rmax = rmaxref;
			if (randi(&iseed) > escape) {
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
		killk(&pressedk);
		if (rmin == rmax) {
			goto L398;
		}
		if (iverb == 1) {
			char temp[84];
			snprintf(temp, sizeof(temp), "\nBest result : a = %lf Rp = %lf\nBest result : c = %lf V = %lf\n\n", bpar[0], rmin, bpar[2], v3);
			fwrite(temp, strlen(temp), 1, imp_file);
			writeFormattedDate(imp_file);
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

	rpsmall = 1.0;
	printf("Orthorhombic: Rp     a       b       c        V     Nind\n");

	ifile = 4;
	ncycles = 1e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 90.0;

	if (ngrid == 3) {
		nruns2 = 6;
		nruns = 10;
		pmin = 2.0;
		pmax = 20.0;
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
		for (int i = 0; i < 3; ++i) {
			pmi[i] = pmin;
			delta[i] = (pma[i] - pmi[i]) / 2.0;
			pstart[i] = pmi[i];
		}
	}

	{
		char temp[300] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		cur += snprintf(cur, end-cur, "\nOrthorhombic Monte Carlo search :\n Max(a,b,c), V %lf %lf %lf %lf\n\n", pma[0], pma[1], pma[2], vmax);
		if (iverb == 1)
		{
			cur += snprintf(cur, end-cur, "  Results in orthorhombic, run, tests : %d %d\n===============================================================================\n Rp  Trial number    a   b   c       V  Nind icod\n\n", nrun, ntimelim[3]);
		}
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	/*C
C     READ hkl Miller indices in ort.hkl
C*/

	readHklFile("ort", 20, 1000, ihh);

	/*FILE *ort_hkl_file = openFile("ort", ".hkl", "r");
	getline(&buffer, &buffsize, ort_hkl_file);
	sscanf(buffer, "%d", &nhkl0);
	nhkl0 = ndat * 20;
	if (nhkl0 > 1000) nhkl0 = 1000;
	for (int i = 1; i <= nhkl0; ++i)
	{
		getline(&buffer, &buffsize, ort_hkl_file);
		sscanf(buffer, "%d %d %d", &ihh[1 + i * 3 - 4], &ihh[2 + i * 3 - 4], &ihh[3 + i * 3 - 4]);
	}
	fclose(ort_hkl_file);*/


	for (int nrun2 = 1; nrun2 <= nruns2; ++nrun2) {
		if (ngrid == 3) {
			if (nrun2 == 1) {
				vmin = 8.0;
				vmax = 500.0;
				if (vorth < 500.0) {
					vmax = vorth;
				}
				ntimelim[3] = (vmax - vmin) * 20.0;
			}
			if (nrun2 == 2) {
				if (vorth < 500.0) {
					goto L500;
				}
				vmin = 500.0;
				vmax = 1e3f;
				if (vorth < 1e3f) {
					vmax = vorth;
				}
				ntimelim[3] = (vmax - vmin) * 20.0;
			}
			if (nrun2 == 3) {
				if (vorth < 1e3f) {
					goto L500;
				}
				vmin = 1e3f;
				vmax = 1500.0;
				if (vorth < 1500.0) {
					vmax = vorth;
				}
				ntimelim[3] = (vmax - vmin) * 20.0;
			}
			if (nrun2 == 4) {
				if (vorth < 1500.0) {
					goto L500;
				}
				vmin = 1500.0;
				vmax = 2e3f;
				if (vorth < 2e3f) {
					vmax = vorth;
				}
				ntimelim[3] = (vmax - vmin) * 20.0;
			}
			if (nrun2 == 5) {
				if (vorth < 2e3f) {
					goto L500;
				}
				vmin = 2e3f;
				vmax = 2500.0;
				if (vorth < 2500.0) {
					vmax = vorth;
				}
				ntimelim[3] = (vmax - vmin) * 20.0;
			}
			if (nrun2 == 6) {
				if (vorth < 2500.0) {
					goto L500;
				}
				vmin = 2500.0;
				vmax = 3e3f;
				if (vorth < 3e3f) {
					vmax = vorth;
				}
				ntimelim[3] = (vmax - vmin) * 20.0;
			}
		}

		for (int nrun = 1; nrun <= nruns; ++nrun) {
	/* ------------------------------------------------------------------------- */
	/*     Initialisation */

	/*      CALL ESP_INIT(ISEED) */

	/* ------------------------------------------------------------------------- */
			rmax = rmaxref;
			rmin = rmax;

	/* ...  here starts the loop */

			interest = 0;
			tmax = ntimelim[3] / procs;
			ttmax = ntimelim[3] * 10.0;
			ncells = (int) ntimelim[3];
			iiseed = 0;
			ntried = 0.0;
			ntriedt = 0.0;
			nout = 0;
			celpre[0] = pstart[0] + delta[0] * 2.0 * randi(&iseed);
			celpre[1] = pstart[1] + delta[1] * 2.0 * randi(&iseed);
			celpre[2] = pstart[2] + delta[2] * 2.0 * randi(&iseed);
			celold[0] = celpre[0];
			celold[1] = celpre[1];
			celold[2] = celpre[2];
			rglob = 1.0;
			nglob = 0;

	/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
	/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC, */
	/* $OMP& RMAX2,A,B,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X, */
	/* $OMP& DIFF,DIFF2,DDT,DDQ) */
	/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
	/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
	/* $OMP DO */

			for (int ncel = 1; ncel <= ncells; ++ncel) {
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

				ntriedb = 0.0;
				x = randi(&iseed);
				if (x >= 0.0 && x < .33333f) {
					ip = 1;
				}
				if (x >= .33333f && x < .66666f) {
					ip = 2;
				}
				if (x >= .66666f && x <= 1.0) {
					ip = 3;
				}
				celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2 * randi(&iseed);
				ntried += 1.0;
				goto L404;
		L403:
				del = deltab * (1.0 - ntriedb / cy);
				x = randi(&iseed);
				int i = 1;
				if (x >= 0.0 && x < .33333f) {
					i = 1;
				}
				if (x >= .33333f && x < .66666f) {
					i = 2;
				}
				if (x >= .66666f && x < 1.0) {
					i = 3;
				}
				celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
				ntriedb += 1.0;
		L404:
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
					}
				}
				dcell(celpre, al, &v1);
				if (ntried > tmax) {
					++nout;
					goto L496;
				}
				if (ntriedb != 0.0) {
					goto L406;
				}
				if (v1 > vmax || v1 < vmin) {
					ntried += -1.0;
					ntriedt += 1.0;
					if (ntriedt > ttmax) {
						++nout;
						goto L496;
					}
					goto L402;
				}

		L406:
				calcul1(&diff, &diff2);
				if (nmx > ndat10) {
					ntried += -1;
					goto L402;
				}
				if (ntriedb != 0.0) {
					goto L414;
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
					goto L417;
				}
				if (lhkl < nmax) {
					goto L417;
				}
		L414:
				if (diff <= rmax) {
					llhkl = lhkl;
					rmax = diff;
					rmax2 = diff2;
					a = celpre[0];
					b = celpre[1];
					c = celpre[2];
					v2 = v1;
					if (diff < rmin) {
						rmin = diff;
						bpar[0] = a;
						bpar[1] = b;
						bpar[2] = c;
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
					celold[ip - 1] = c;
				}
				ntriedb = 0.0;
				if (rmax >= rmax0[3]) {
					goto L417;
				}
				if (rmax2 >= .15f) {
					goto L417;
				}
				ipen = ndat - llhkl;
				if (ipen > nind) {
					goto L417;
				}

	/* $OMP CRITICAL(STORE1) */

				++igc;

		/*  Test if too much proposals, if yes decrease Rmax by 5% */

				igt += 1.0;
				if (nr == 1) {
					if (igt > 50.0) {
						if (ntried / igt < 1e4f) {
							if (rmax0[3] > .2f) {
								rmax0[3] -= rmax0[3] * .05f;
								printRmaxReducedString(rmax0[3], imp_file);
							}
						}
					}
				}

				if (igc > 10000) {
					printStopString(imp_file);
					--igc;
					++interest;
				}
				cel[igc * 6 - 6] = a;
				cel[igc * 6 - 5] = b;
				cel[igc * 6 - 4] = c;
				cel[igc * 6 - 3] = 90.0;
				cel[igc * 6 - 2] = 90.0;
				cel[igc * 6 - 1] = 90.0;

		/* $OMP END CRITICAL(STORE1) */

		/* ... Check for supercell */

				celpre[0] = a;
				celpre[1] = b;
				celpre[2] = c;
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
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
				supcel(&lhkl, ihkl, cel, &igc, vgc, &c1);
				brav(&lhkl, ihkl, &ibr);
				ib[igc - 1] = ibr;
				a = cel[igc * 6 - 6];
				b = cel[igc * 6 - 5];
				c = cel[igc * 6 - 4];
				v2 = vgc[igc - 1];

		/* $OMP END CRITICAL(STORE2) */

		/* ... Check for interesting result */

		/*      IF(INTEREST.GE.1)GO TO 496 */
				indic = 0;
				bb[2] = a;
				bb[3] = b;
				bb[4] = c;
				bb[5] = 90.0;
				bb[6] = 90.0;
				bb[7] = 90.0;
				afi[2] = 1.0;
				afi[3] = 1.0;
				afi[4] = 1.0;
				afi[5] = 0.0;
				afi[6] = 0.0;
				afi[7] = 0.0;
				celpre[0] = a;
				celpre[1] = b;
				celpre[2] = c;
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
					}
				}

	/* $OMP CRITICAL(FOUND) */

				if (rp[igc - 1] < rmi) {
					++interest;
					printSaveInterstResString(imp_file, rmax, v2, ipen, 3, a, b, c);
					/*const char *temp2 = saveInterestingResultString(imp_file);
					//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
					printIsee(rmax, v2, ipen, 3, a, b, c);
					printf("%s\n", temp2);

		/* ... Refine that cell */

					dcell(celpre, al, &v1);
					calcul2(&diff, ihkl, th3, &ncalc, &igc);
					celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
					cncalc[igc - 1] = (double) ncalc;
					if (ndat >= 20) {
						fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
					} else {
						pndat = (double) ndat;
						fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1]	* 2.0 * ddq);
						ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
					}
					iref = 1;
					goto L497;
				} else {

		/*  Anyway, calculate the M20 and F20 values */

					dcell(celpre, al, &v1);
					calcul2(&diff, ihkl, th3, &ncalc, &igc);
					celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
					cncalc[igc - 1] = (double) ncalc;
					if (ndat >= 20) {
						fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
					} else {
						pndat = (double) ndat;
						fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1]
							* 2.0 * ddq);
						ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
					}
				}

		/* Test if cell already found */

				if (igc > 1) {
					for (int i = 1; i < igc; ++i) {
						if (ifi[i - 1] != ifile) {
							goto L418;
						}
						vdelt = vgc[igc - 1] / 300.0;
						vp = vgc[igc - 1] + vdelt;
						vm = vgc[igc - 1] - vdelt;
						if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
							goto L418;
						}
						adelt = cel[igc * 6 - 6] / 500.0;
						ap = cel[igc * 6 - 6] + adelt;
						am = cel[igc * 6 - 6] - adelt;
						int na = 0;
						if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
							na = 1;
						}
						bdelt = cel[igc * 6 - 5] / 500.0;
						bp = cel[igc * 6 - 5] + bdelt;
						bm = cel[igc * 6 - 5] - bdelt;
						nb = 0;
						if (cel[i * 6 - 6] > bp || cel[i * 6 - 6] < bm) {
							nb = 1;
						}
						cdelt = cel[igc * 6 - 4] / 500.0;
						cp = cel[igc * 6 - 4] + cdelt;
						cm = cel[igc * 6 - 4] - cdelt;
						int nc = 0;
						if (cel[i * 6 - 6] > cp || cel[i * 6 - 6] < cm) {
							nc = 1;
						}
						if (na == 1 && nb == 1 && nc == 1) {
							goto L418;
						}
						++nsol[i - 1];
						if (rp[igc - 1] < rp[i - 1]) {
							if (isee == 1) {
								printIsee(rmax, v2, ipen, 3, a, b, c);
								// printf("%.3lf %.4lf %.4lf %.4lf %.1lf %d\n", rmax, a, b, c, v2, ipen);
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
							ntried = tmax + 1.0;
							++nout;
						}
						goto L419;
			L418:
						;
					}
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 3, a, b, c);
					// char temp[40];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, b, c, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 3, a, b, c);
					// printf("%.3lf %.4lf %.4lf %.4lf %.1lf %d\n", rmax, a, b, c, v2, ipen);
				}
	L419:
				;
			} else {
				if (iverb == 1) {
					saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 3, a, b, c);
					// char temp[40];
					// snprintf(temp, sizeof(temp), "%.3lf %.0lf %.4lf %.4lf %.4lf %.1lf %d %d\n", rmax, ntried, a, b, c, v2, ipen, icode);
					// fwrite(temp, strlen(temp), 1, imp_file);
				}
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 3, a, b, c);
					// printf("%.3lf %.4lf %.4lf %.4lf %.1lf %d\n", rmax, a, b, c, v2, ipen);
				}
			}
	L497:

	/* $OMP END CRITICAL(FOUND) */

	L417:
			rmax = rmaxref;
			if (randi(&iseed) > escape) {
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
			killk(&pressedk);
			if (rmin == rmax) {
				goto L498;
			}
			if (iverb == 1) {
				char temp[111];
				snprintf(temp, sizeof(temp), "\nBest result : a = %lf Rp = %lf\nBest result : b = %lf\nBest result : c = %lf V = %lf\n\n", bpar[0], rmin, bpar[1], bpar[2], v3);
				fwrite(temp, strlen(temp), 1, imp_file);
				writeFormattedDate(imp_file);
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

/* Monoclinic case */

L500:
	if (nsys[4] == 0) {
		goto L600;
	}

	rpsmall = 1.0;
	printf("Monoclinic:   Rp     a       b       c       bet     V     Nind\n");
	ifile = 5;
	ncycles = 2e3f;
	cy = ncycles * 1.1;
	celpre[3] = 90.0;
	celpre[5] = 90.0;

	if (ngrid == 3) {
		nruns = 20;
		nruns2 = 6;
		pmin = 2.0;
		pmax = 20.0;
		pma[0] = dmax1 * 2.1;
		pma[1] = dmax1 * 2.1;
		pma[2] = dmax2 * 2.1;
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
		for (int i = 0; i < 3; ++i) {
			pmi[i] = pmin;
			delta[i] = (pma[i] - pmi[i]) / 2.0;
			pstart[i] = pmi[i];
		}
	}

	{
		char temp[300] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		cur += snprintf(cur, end-cur, "\nMonoclinic Monte Carlo search :\n Max(a,b,c), V %lf %lf %lf %lf\n\n", pma[0], pma[1], pma[2], vmax);
		if (iverb == 1)
		{
			cur += snprintf(cur, end-cur, "  Results in monoclinic, run, tests : %d %d\n===============================================================================\n Rp  Trial number    a   b   c    bet   V  Nind Icod\n\n", nrun, ntimelim[4]);
		}
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	/*C
C     READ hkl Miller indices in mon.hkl
C*/

	readHklFile("mon", 20, 1000, ihh);

	for (int nrun2 = 1; nrun2 <= nruns2; ++nrun2) {

		if (ngrid == 3) {
			if (nrun2 == 1) {
				vmin = 8.0;
				vmax = 500.0;
				if (vmon < 500.0) {
					vmax = vmon;
				}
				ntimelim[4] = (vmax - vmin) * 200.0;
			}
			if (nrun2 == 2) {
				if (vmon < 500.0) {
					goto L600;
				}
				vmin = 500.0;
				vmax = 1e3f;
				if (vmon < 1e3f) {
					vmax = vmon;
				}
				ntimelim[4] = (vmax - vmin) * 200.0;
			}
			if (nrun2 == 3) {
				if (vmon < 1e3f) {
					goto L600;
				}
				vmin = 1e3f;
				vmax = 1500.0;
				if (vmon < 1500.0) {
					vmax = vmon;
				}
				ntimelim[4] = (vmax - vmin) * 200.0;
			}
			if (nrun2 == 4) {
				if (vmon < 1500.0) {
					goto L600;
				}
				vmin = 1500.0;
				vmax = 2e3f;
				if (vmon < 2e3f) {
					vmax = vmon;
				}
				ntimelim[4] = (vmax - vmin) * 200.0;
			}
			if (nrun2 == 5) {
				if (vmon < 2e3f) {
					goto L600;
				}
				vmin = 2e3f;
				vmax = 2500.0;
				if (vmon < 2500.0) {
					vmax = vmon;
				}
				ntimelim[4] = (vmax - vmin) * 200.0;
			}
			if (nrun2 == 6) {
				if (vmon < 2500.0) {
					goto L600;
				}
				vmin = 2500.0;
				vmax = 3e3f;
				if (vmon < 3e3f) {
					vmax = vmon;
				}
				ntimelim[4] = (vmax - vmin) * 200.0;
			}
		}

		for (int nrun = 1; nrun <= nruns; ++nrun) {
	/* ------------------------------------------------------------------------- */
	/*     Initialisation */

	/*      CALL ESP_INIT(ISEED) */

	/* ------------------------------------------------------------------------- */
			rmax = rmaxref;
			rmin = rmax;

	/* ...  here starts the loop */

			interest = 0;
			tmax = ntimelim[4] / procs;
			ttmax = ntimelim[4] * 10.0;
			ncells = (int) ntimelim[4];
			iiseed = 0;
			ntried = 0.0;
			ntriedt = 0.0;
			nout = 0;
			celpre[0] = pstart[0] + delta[0] * 2.0 * randi(&iseed);
			celpre[1] = pstart[1] + delta[1] * 2.0 * randi(&iseed);
			celpre[2] = pstart[2] + delta[2] * 2.0 * randi(&iseed);
			celpre[4] = astart + deltc * 2.0 * randi(&iseed);
			celold[0] = celpre[0];
			celold[1] = celpre[1];
			celold[2] = celpre[2];
			celold[4] = celpre[4];
			rglob = 1.0;
			nglob = 0;
	/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
	/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,DELD,V1,ICODE,LLHKL,IHKL,TH3, */
	/* $OMP& RMAX2,A,B,C,BET,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X, */
	/* $OMP& DIFF,DIFF2,NCALC,DDT,DDQ) */
	/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
	/* $OMP& celpre,celold,rglob,nold,rmin,rmax,bb,afi) */
	/* $OMP DO */

			for (int ncel = 1; ncel <= ncells; ++ncel) {
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

				ntriedb = 0.0;
				x = randi(&iseed);
				if (x >= 0.0 && x < .25f) {
					ip = 1;
				}
				if (x >= .25f && x < .5f) {
					ip = 2;
				}
				if (x >= .5f && x < .75f) {
					ip = 3;
				}
				if (x >= .75f && x <= 1.0) {
					ip = 5;
				}
				if (ip != 5) {
					celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.0 *	randi(&iseed);
				} else {
					celpre[ip - 1] = astart + deltc * 2.0 * randi(&iseed);
				}
				ntried += 1.0;
				goto L504;
		L503:
				del = deltab * (1.0 - ntriedb / cy);
				deld = deltad * (1.0 - ntriedb / cy);
				x = randi(&iseed);
				int i = 1;
				if (x >= 0.0 && x < .25f) {
					i = 1;
				}
				if (x >= .25f && x < .5f) {
					i = 2;
				}
				if (x >= .5f && x < .75f) {
					i = 3;
				}
				if (x >= .75f && x <= 1.0) {
					i = 5;
				}
				if (i != 5) {
					celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
				} else {
					celpre[i - 1] = pstartb[i - 1] + deld * (randi(&iseed) - .5f) * 2.0;
				}
				ntriedb += 1.0;
		L504:
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
					}
				}
				dcell(celpre, al, &v1);
				if (ntried > tmax) {
					++nout;
					goto L596;
				}
				if (ntriedb != 0.0) {
					goto L506;
				}
				if (v1 > vmax || v1 < vmin) {
					ntried += -1.0;
					ntriedt += 1.0;
					if (ntriedt > ttmax) {
						++nout;
						goto L596;
					}
					goto L502;
				}

		L506:
				calcul1(&diff, &diff2);
				if (nmx > ndat10) {
					ntried += -1;
					goto L502;
				}
				if (ntriedb != 0.0) {
					goto L514;
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
					goto L517;
				}
				if (lhkl < nmax) {
					goto L517;
				}
		L514:
				if (diff <= rmax) {
					llhkl = lhkl;
					rmax = diff;
					rmax2 = diff2;
					a = celpre[0];
					b = celpre[1];
					c = celpre[2];
					bet = celpre[4];
					v2 = v1;
					if (diff < rmin) {
						rmin = diff;
						bpar[0] = a;
						bpar[1] = b;
						bpar[2] = c;
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
					celold[ip - 1] = c;
				}
				if (ip == 5) {
					celold[ip - 1] = bet;
				}
				ntriedb = 0.0;
				if (rmax >= rmax0[4]) {
					goto L517;
				}
				if (rmax2 >= .15f) {
					goto L517;
				}
				ipen = ndat - llhkl;
				if (ipen > nind) {
					goto L517;
				}

		/* $OMP CRITICAL(STORE1) */

				++igc;

		/*  Test if too much proposals, if yes decrease Rmax by 5% */

				igt += 1.0;
				if (nr == 1) {
					if (igt > 50.0) {
						if (ntried / igt < 1e5f) {
							if (rmax0[4] > 0.2) {
								rmax0[4] -= rmax0[4] * 0.05;
								printRmaxReducedString(rmax0[4], imp_file);
							}
						}
					}
				}

				if (igc > 10000) {
					printStopString(imp_file);
					--igc;
					++interest;
		/*      GO TO 5000 */
				}
				cel[igc * 6 - 6] = a;
				cel[igc * 6 - 5] = b;
				cel[igc * 6 - 4] = c;
				cel[igc * 6 - 3] = 90.0;
				cel[igc * 6 - 2] = bet;
				cel[igc * 6 - 1] = 90.0;

		/* $OMP END CRITICAL(STORE1) */

		/* ... Check for supercell */

				celpre[0] = a;
				celpre[1] = b;
				celpre[2] = c;
				celpre[4] = bet;
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
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
				supcel(&lhkl, ihkl, cel, &igc, vgc, &c1);
				brav(&lhkl, ihkl, &ibr);
				ib[igc - 1] = ibr;
				a = cel[igc * 6 - 6];
				b = cel[igc * 6 - 5];
				c = cel[igc * 6 - 4];
				v2 = vgc[igc - 1];

		/* $OMP END CRITICAL(STORE2) */

		/* ... Check for interesting result */

		/*      IF(INTEREST.GE.1)GO TO 596 */
				indic = 0;
				bb[2] = a;
				bb[3] = b;
				bb[4] = c;
				bb[5] = 90.0;
				bb[6] = bet;
				bb[7] = 90.0;
				afi[2] = 1.0;
				afi[3] = 1.0;
				afi[4] = 1.0;
				afi[5] = 0.0;
				afi[6] = 1.0;
				afi[7] = 0.0;
				celpre[0] = a;
				celpre[1] = b;
				celpre[2] = c;
				celpre[4] = bet;
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
					}
				}

		/* $OMP CRITICAL(FOUND) */

				if (rp[igc - 1] < rmi) {
					++interest;
					printSaveInterstResString(imp_file, rmax, v2, ipen, 4, a, b, c, bet);
					/*const char *temp2 = saveInterestingResultString(imp_file);
					//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
					printIsee(rmax, v2, ipen, 4, a, b, c, bet);
					printf("%s\n", temp2);

		/* ... Refine that cell */

					dcell(celpre, al, &v1);
					calcul2(&diff, ihkl, th3, &ncalc, &igc);
					celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
					cncalc[igc - 1] = (double) ncalc;
					if (ndat >= 20) {
						fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
					} else {
						pndat = (double) ndat;
						fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] * 2.0 * ddq);
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

					dcell(celpre, al, &v1);
					calcul2(&diff, ihkl, th3, &ncalc, &igc);
					celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
					cncalc[igc - 1] = (double) ncalc;
					if (ndat >= 20) {
						fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
					} else {
						pndat = (double) ndat;
						fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
					}
				}

		/* Test if cell already found */


				if (igc > 1) {
					for (i = 1; i < igc; ++i) {
						if (ifi[i - 1] != ifile) {
							goto L518;
						}
						vdelt = vgc[igc - 1] / 300.0;
						vp = vgc[igc - 1] + vdelt;
						vm = vgc[igc - 1] - vdelt;
						if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
							goto L518;
						}
						bdelt = cel[igc * 6 - 5] / 500.0;
						bp = cel[igc * 6 - 5] + bdelt;
						bm = cel[igc * 6 - 5] - bdelt;
						if (cel[i * 6 - 5] > bp || cel[i * 6 - 5] < bm) {
							goto L518;
						}
						betdelt = cel[igc * 6 - 2] / 500.0;
						betp = cel[igc * 6 - 2] + betdelt;
						betm = cel[igc * 6 - 2] - betdelt;
						if (cel[i * 6 - 2] > betp || cel[i * 6 - 2] < betm) {
							goto L518;
						}
						adelt = cel[igc * 6 - 6] / 500.0;
						ap = cel[igc * 6 - 6] + adelt;
						am = cel[igc * 6 - 6] - adelt;
						na = 0;
						if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
							na = 1;
						}
						cdelt = cel[igc * 6 - 4] / 500.0;
						cp = cel[igc * 6 - 4] + cdelt;
						cm = cel[igc * 6 - 4] - cdelt;
						nc = 0;
						if (cel[i * 6 - 6] > cp || cel[i * 6 - 6] < cm) {
							nc = 1;
						}
						if (na == 1 && nc == 1) {
							goto L518;
						}
						++nsol[i - 1];
						if (rp[igc - 1] < rp[i - 1]) {
							if (isee == 1) {
								printIsee(rmax, v2, ipen, 4, a, b, c, bet);
							}
							km[i - 1] = km[igc - 1];
							vgc[i - 1] = vgc[igc - 1];
							rp[i - 1] = rp[igc - 1];
							cel[i * 6 - 6] = cel[igc * 6 - 6];
							cel[i * 6 - 5] = cel[igc * 6 - 5];
							cel[i * 6 - 4] = cel[igc * 6 - 4];
							cel[i * 6 - 2] = cel[igc * 6 - 2];
						}
						--igc;
						if (nsol[i - 1] > 5) {
							ntried = tmax + 1.0;
							++nout;
						}
						goto L519;
			L518:
						;
					}
					if (iverb == 1) {
						saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 4, a, b, c, bet);
					}
					if (isee == 1) {
						printIsee(rmax, v2, ipen, 4, a, b, c, bet);
					}
		L519:
					;
				} else {
					if (iverb == 1) {
						saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 4, a, b, c, bet);
					}
					if (isee == 1) {
						printIsee(rmax, v2, ipen, 4, a, b, c, bet);
					}
				}
		L597:

		/* $OMP END CRITICAL(FOUND) */

		L517:
				rmax = rmaxref;
				if (randi(&iseed) > escape) {
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
			killk(&pressedk);
			if (rmin == rmax) {
				goto L598;
			}
			if (iverb == 1) {
				char temp[151];
				snprintf(temp, sizeof(temp), "\nBest result : a =    %lf Rp = %lf\nBest result : b =    %lf\nBest result : c =    %lf \nBest result : beta = %lf V = %lf\n\n", bpar[0], rmin, bpar[1], bpar[2], bpar[4], v3);
				fwrite(temp, strlen(temp), 1, imp_file);
				writeFormattedDate(imp_file);
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


	rpsmall = 1.0;
	printf("Triclinic:    Rp     a       b       c       alp    bet    gam     V     Nind\n");
	ifile = 6;
	ncycles = 5e3f;
	cy = ncycles * 1.1f;

	if (ngrid == 3) {
		nruns = 20;
		nruns2 = 8;
		pmin = 2.0;
		pmax = 20.0;
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
		for (int i = 0; i < 3; ++i) {
			pmi[i] = pmin;
			delta[i] = (pma[i] - pmi[i]) / 2.0;
	/* L623: */
			pstart[i] = pmi[i];
		}
	}

	{
		char temp[300] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		cur += snprintf(cur, end-cur, "\nTriclinic Monte Carlo search :\n Max(a,b,c), V %lf %lf %lf %lf\n\n", pma[0], pma[1], pma[2], vmax);
		if (iverb == 1)
		{
			cur += snprintf(cur, end-cur, "  Results in triclinic, run, tests : %d %d\n===============================================================================\n Rp Trial number  a  b  c  alp bet gam  V  Nind Icod\n\n", nrun, ntimelim[5]);
		}
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*C
C     READ hkl Miller indices in tri.hkl
C*/

	readHklFile("tri", 20, 1000, ihh);

	for (int nrun2 = 1; nrun2 <= nruns2; ++nrun2) {


		if (ngrid == 3) {
			if (nrun2 == 1) {
				vmin = 8.0;
				vmax = 250.0;
				if (vtric < 250.0) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
			if (nrun2 == 2) {
				if (vtric < 250.0) {
					goto L700;
				}
				vmin = 250.0;
				vmax = 500.0;
				if (vtric < 500.0) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
			if (nrun2 == 3) {
				if (vtric < 500.0) {
					goto L700;
				}
				vmin = 500.0;
				vmax = 750.0;
				if (vtric < 750.0) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
			if (nrun2 == 4) {
				if (vtric < 750.0) {
					goto L700;
				}
				vmin = 750.0;
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
				vmax = 1250.0;
				if (vtric < 1250.0) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
			if (nrun2 == 6) {
				if (vtric < 1250.0) {
					goto L700;
				}
				vmin = 1250.0;
				vmax = 1500.0;
				if (vtric < 1500.0) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
			if (nrun2 == 7) {
				if (vtric < 1500.0) {
					goto L700;
				}
				vmin = 1500.0;
				vmax = 1750.0;
				if (vtric < 1750.0) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
			if (nrun2 == 8) {
				if (vtric < 1750.0) {
					goto L700;
				}
				vmin = 1750.0;
				vmax = 2e3f;
				if (vtric < 2e3f) {
					vmax = vtric;
				}
				ntimelim[5] = (vmax - vmin) * 4e3f;
			}
		}

		for (int nrun = 1; nrun <= nruns; ++nrun) {
	/* ------------------------------------------------------------------------- */
	/*     Initialisation */

	/*      CALL ESP_INIT(ISEED) */

	/* ------------------------------------------------------------------------- */
			rmax = rmaxref;
			rmin = rmax;

	/* ...  here starts the loop */

			interest = 0;
			tmax = ntimelim[5] / procs;
			ttmax = ntimelim[5] * 10.0;
			ncells = (int) ntimelim[5];
			iiseed = 0;
			ntried = 0.0;
			ntriedt = 0.0;
			nout = 0;

	/* ...  here starts the loop */

			celpre[0] = pstart[0] + delta[0] * 2.0 * randi(&iseed);
			celpre[1] = pstart[1] + delta[1] * 2.0 * randi(&iseed);
			celpre[2] = pstart[2] + delta[2] * 2.0 * randi(&iseed);
			celpre[3] = astartt[0] + deltct[0] * 2.0 * randi(&iseed);
			celpre[4] = astartt[1] + deltct[1] * 2.0 * randi(&iseed);
			celpre[5] = astartt[2] + deltct[2] * 2.0 * randi(&iseed);
			celold[0] = celpre[0];
			celold[1] = celpre[1];
			celold[2] = celpre[2];
			celold[3] = celpre[3];
			celold[4] = celpre[4];
			celold[5] = celpre[5];
			rglob = 1.0;
			nglob = 0;

	/* $OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/) */
	/* $OMP& PRIVATE(NCEL,NTRIEDB,DEL,DELD,V1,ICODE,LLHKL,IHKL,TH3, */
	/* $OMP& RMAX2,A,B,C,ALP,BET,GAM,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X, */
	/* $OMP& DIFF,DIFF2,IP2,ANG,NCALC,DDT,DDQ) */
	/* $OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout, */
	/* $OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi) */
	/* $OMP DO */

			for (int ncel = 1; ncel <= ncells; ++ncel) {
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

				ntriedb = 0.0;
				x = randi(&iseed);
				if (x >= 0.0 && x < .16666f) {
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
				if (x >= .83333f && x <= 1.0) {
					ip = 6;
				}
				if (ip != 4 && ip != 5 && ip != 6) {
					celpre[ip - 1] = pstart[ip - 1] + delta[ip - 1] * 2.0 *	randi(&iseed);
				} else {
					celpre[ip - 1] = astartt[ip - 4] + deltct[ip - 4] * 2.0 * randi(&iseed);
				}
				ang = celpre[3] + celpre[4] + celpre[5];
				if (ang >= 360.0 && ang <= 180.0) {
					goto L696;
				}
				ntried += 1.0;
				goto L604;
		L603:
				del = deltab * (1.0 - ntriedb / cy);
				deld = deltad * (1.0 - ntriedb / cy);
				x = randi(&iseed);
				if (x >= 0.0 && x < .16666f) {
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
				if (x >= .83333f && x <= 1.0) {
					ip2 = 6;
				}
				if (ip2 != 4 && ip2 != 5 && ip2 != 6) {
					celpre[ip2 - 1] = pstartb[ip2 - 1] + del * (randi(&iseed) - .5f) * 2.0;
				} else {
					celpre[ip2 - 1] = pstartb[ip2 - 1] + deld * (randi(&iseed) - .5f) * 2.0;
				}
				ang = celpre[3] + celpre[4] + celpre[5];
				if (ang >= 360.0 && ang <= 180.0) {
					goto L603;
				}
				ntriedb += 1.0;
		L604:
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
					}
				}
				dcell(celpre, al, &v1);
				if (ntried > tmax) {
					++nout;
					goto L696;
				}
				if (ntriedb != 0.0) {
					goto L606;
				}
				if (v1 > vmax || v1 < vmin) {
					ntried += -1.0;
					ntriedt += 1.0;
					if (ntriedt > ttmax) {
						++nout;
						goto L696;
					}
					goto L602;
				}

		L606:
				calcul1(&diff, &diff2);
				if (nmx > ndat10) {
					ntried += -1;
					goto L602;
				}
				if (ntriedb != 0.0) {
					goto L614;
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
					goto L617;
				}
				if (lhkl < nmax) {
					goto L617;
				}
		L614:
				if (diff <= rmax) {
					llhkl = lhkl;
					rmax = diff;
					rmax2 = diff2;
					a = celpre[0];
					b = celpre[1];
					c = celpre[2];
					alp = celpre[3];
					bet = celpre[4];
					gam = celpre[5];
					v2 = v1;
					if (diff < rmin) {
						rmin = diff;
						bpar[0] = a;
						bpar[1] = b;
						bpar[2] = c;
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
					celold[ip - 1] = c;
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
				ntriedb = 0.0;
				if (rmax >= rmax0[5]) {
					goto L617;
				}
				if (rmax2 >= .15f) {
					goto L617;
				}
				ipen = ndat - llhkl;
				if (ipen > nind) {
					goto L617;
				}

		/* $OMP CRITICAL(STORE1) */

				++igc;

		/*  Test if too much proposals, if yes decrease Rmax by 5% */

				igt += 1.0;
				if (nr == 1) {
					if (igt > 50.0) {
						if (ntried / igt < 1e5f) {
							if (rmax0[5] > 0.2) {
								rmax0[5] -= rmax0[5] * 0.05;
								printRmaxReducedString(rmax0[5], imp_file);
							}
						}
					}
				}

				if (igc > 10000) {
					printStopString(imp_file);
					--igc;
					++interest;
		/*      GO TO 5000 */
				}
				cel[igc * 6 - 6] = a;
				cel[igc * 6 - 5] = b;
				cel[igc * 6 - 4] = c;
				cel[igc * 6 - 3] = alp;
				cel[igc * 6 - 2] = bet;
				cel[igc * 6 - 1] = gam;

		/* $OMP END CRITICAL(STORE1) */

		/* ... Check for supercell */

				celpre[0] = a;
				celpre[1] = b;
				celpre[2] = c;
				celpre[3] = alp;
				celpre[4] = bet;
				celpre[5] = gam;
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
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
				supcel(&lhkl, ihkl, cel, &igc, vgc, &c1);
				brav(&lhkl, ihkl, &ibr);
				ib[igc - 1] = ibr;
				a = cel[igc * 6 - 6];
				b = cel[igc * 6 - 5];
				c = cel[igc * 6 - 4];
				v2 = vgc[igc - 1];

		/* $OMP END CRITICAL(STORE2) */

		/* ... Check for interesting result */

		/*      IF(INTEREST.GE.1)GO TO 696 */
				indic = 0;
				bb[2] = a;
				bb[3] = b;
				bb[4] = c;
				bb[5] = alp;
				bb[6] = bet;
				bb[7] = gam;
				afi[2] = 1.0;
				afi[3] = 1.0;
				afi[4] = 1.0;
				afi[5] = 1.0;
				afi[6] = 1.0;
				afi[7] = 1.0;
				celpre[0] = a;
				celpre[1] = b;
				celpre[2] = c;
				celpre[3] = alp;
				celpre[4] = bet;
				celpre[5] = gam;
				for (int i = 1; i <= 3; ++i) {
					for (int j = 1; j <= 3; ++j) {
						al[i + j * 3 - 4] = 0.0;
					}
				}

		/* $OMP CRITICAL(FOUND) */

				if (rp[igc - 1] < rmi) {
					++interest;
					printSaveInterstResString(imp_file, rmax, v2, ipen, 6, a, b, c, alp, bet, gam);
					/*const char *temp2 = saveInterestingResultString(imp_file);
					//1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
					printIsee(rmax, v2, ipen, 6, a, b, c, alp, bet, gam);
					printf("%s\n", temp2);

		/* ... Refine that cell */

					dcell(celpre, al, &v1);
					calcul2(&diff, ihkl, th3, &ncalc, &igc);
					celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
					cncalc[igc - 1] = (double) ncalc;
					if (ndat >= 20) {
						fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
					} else {
						pndat = (double) ndat;
						fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1]	* 2.0 * ddq);
						ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
					}

					iref = 1;
					goto L697;
				} else {

		/*  Anyway, calculate the M20 and F20 values */

					dcell(celpre, al, &v1);
					calcul2(&diff, ihkl, th3, &ncalc, &igc);
					celref2(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq);
					cncalc[igc - 1] = (double) ncalc;
					if (ndat >= 20) {
						fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
						ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
					} else {
						pndat = (double) ndat;
						fm20[igc - 1] = qo[ndat - 1] / (cncalc[igc - 1]	* 2.0 * ddq);
						ff20[igc - 1] = pndat / (cncalc[igc - 1] * ddt);
					}
				}

		/* Test if cell already found */


				if (igc > 1) {
					for (int i = 1; i < igc; ++i) {
						if (ifi[i - 1] != ifile) {
							goto L618;
						}
						vdelt = vgc[igc - 1] / 300.0;
						vp = vgc[igc - 1] + vdelt;
						vm = vgc[igc - 1] - vdelt;
						if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
							goto L618;
						}
						adelt = cel[igc * 6 - 6] / 500.0;
						ap = cel[igc * 6 - 6] + adelt;
						am = cel[igc * 6 - 6] - adelt;
						na = 0;
						if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
							na = 1;
						}
						bdelt = cel[igc * 6 - 5] / 500.0;
						bp = cel[igc * 6 - 5] + bdelt;
						bm = cel[igc * 6 - 5] - bdelt;
						nb = 0;
						if (cel[i * 6 - 6] > bp || cel[i * 6 - 6] < bm) {
							nb = 1;
						}
						cdelt = cel[igc * 6 - 4] / 500.0;
						cp = cel[igc * 6 - 4] + cdelt;
						cm = cel[igc * 6 - 4] - cdelt;
						nc = 0;
						if (cel[i * 6 - 6] > cp || cel[i * 6 - 6] < cm) {
							nc = 1;
						}
						if (na == 1 && nb == 1 && nc == 1) {
							goto L618;
						}
						na = 0;
						if (cel[i * 6 - 5] > ap || cel[i * 6 - 5] < am) {
							na = 1;
						}
						nb = 0;
						if (cel[i * 6 - 5] > bp || cel[i * 6 - 5] < bm) {
							nb = 1;
						}
						nc = 0;
						if (cel[i * 6 - 5] > cp || cel[i * 6 - 5] < cm) {
							nc = 1;
						}
						if (na == 1 && nb == 1 && nc == 1) {
							goto L618;
						}
						++nsol[i - 1];
						if (rp[igc - 1] < rp[i - 1]) {
							if (isee == 1) {
								printIsee(rmax, v2, ipen, 6, a, b, c, alp, bet, gam);
							}
			/*     1WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN */
							km[i - 1] = km[igc - 1];
							vgc[i - 1] = vgc[igc - 1];
							rp[i - 1] = rp[igc - 1];
							cel[i * 6 - 6] = cel[igc * 6 - 6];
							cel[i * 6 - 5] = cel[igc * 6 - 5];
							cel[i * 6 - 4] = cel[igc * 6 - 4];
							cel[i * 6 - 3] = cel[igc * 6 - 3];
							cel[i * 6 - 2] = cel[igc * 6 - 2];
							cel[i * 6 - 1] = cel[igc * 6 - 1];
						}
						--igc;
						if (nsol[i - 1] > 5) {
							ntried = tmax + 1.0;
							++nout;
						}
						goto L619;
			L618:
						;
					}
					if (iverb == 1) {
						saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 6, a, b, c, alp, bet, gam);
					}
					if (isee == 1) {
						printIsee(rmax, v2, ipen, 6, a, b, c, alp, bet, gam);
					}
		L619:
					;
				} else {
					if (iverb == 1) {
						saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 6, a, b, c, alp, bet, gam);
					}
					if (isee == 1) {
						printIsee(rmax, v2, ipen, 6, a, b, c, alp, bet, gam);
					}
				}
		L697:

		/* $OMP END CRITICAL(FOUND) */

		L617:
				rmax = rmaxref;
				if (randi(&iseed) > escape) {
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
			killk(&pressedk);
	/* L616: */
			if (rmin == rmax) {
				goto L698;
			}
			if (iverb == 1) {
				char temp[211];
				snprintf(temp, sizeof(temp), "\nBest result : a =    %lf Rp = %lf\nBest result : b =    %lf\nBest result : c =    %lf \nBest result : alph = %lf\nBest result : beta = %lf\nBest result : gamm = %lf V = %lf\n\n", bpar[0], rmin, bpar[1], bpar[2], bpar[3], bpar[4], bpar[5], v3);
				fwrite(temp, strlen(temp), 1, imp_file);
				writeFormattedDate(imp_file);
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

L700:
	if (ngrid == 0) {
		goto L5000;
	}
	if (nblack == 0 && ngrid == 3) {
		goto L5000;
	}
	printf("\nGrid search :\n");

/* ...  Cell generation by systematic grid */


	if (nsys[0] == 0) {
		goto L1200;
	}

/*    Cubic case */

	rpsmall = 1.0;
	printf("Cubic:        Rp     a       V     Nind\n");

		ifile = 1;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.0;
	ncycles = 200.0;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 90.0;
	if (iverb == 1) {
		saveGridResultsInString("cubic", imp_file);
	}

	if (ngrid == 3) {
		pmin = 2.0;
		pmax = dmax1 * 3.1f;
		pmi[0] = pmin;
		pma[0] = pmax;
		vmin = 8.0;
		vmax = pmax * pmax * pmax;
		if (iverb == 1) {
			char temp[30];
			snprintf(temp, sizeof(temp), " Max a, V %lf %lf\n\n", pmax, vmax);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
	}

	if (iverb == 1) {
		char *temp = " Rp  Trial number    a         V  Nind Icod\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	/*C
C     READ hkl Miller indices in cub.hkl
C*/

	readHklFile("cub", 6, 400, ihh);

	/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
L1102:
	ntriedb = 0.0;
	celpre[0] += spar;
	ntried += 1.0;
	goto L1104;
L1103:
	del = deltab * (1.0 - ntriedb / cy);
	celpre[0] = pstartb[0] + del * (randi(&iseed) - .5f) * 2.0;
	ntriedb += 1.0;
L1104:
	celpre[1] = celpre[0];
	celpre[2] = celpre[0];
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.0) {
		goto L1116;
	}
	if (ntriedb != 0.0) {
		goto L1106;
	}
	if (v1 > vmax || v1 < vmin) {
		ntried += -1.0;
		goto L1102;
	}

L1106:
	calcul1(&diff, &diff2);
	if (nmx > ndat10) {
		ntried += -1;
		goto L1102;
	}
	if (ntriedb != 0.0) {
		goto L1114;
	}

	/* ... Rp value satisfying ??? */

	if (lhkl >= nmax) {
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
	if (lhkl < nmax) {
		goto L1117;
	}
L1114:
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
		goto L1103;
	}
	ntriedb = 0.0;
	if (rmax >= rmax0[0]) {
		goto L1117;
	}
	if (rmax2 >= .15f) {
		goto L1117;
	}
	ipen = ndat - llhkl;
	if (ipen > nind) {
		goto L1117;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.0;
	if (nr == 1) {
		if (igt > 50.0) {
			if (ntried / igt < 100.0) {
				if (rmax0[0] > .1f) {
					rmax0[0] -= rmax0[0] * .05f;
					printRmaxReducedString(rmax0[0], imp_file);
				}
			}
		}
	}

	if (igc > 10000) {
		printStopString(imp_file);
		--igc;
		goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = a;
	cel[igc * 6 - 3] = 90.0;
	cel[igc * 6 - 2] = 90.0;
	cel[igc * 6 - 1] = 90.0;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = a;
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
				al[i + j * 3 - 4] = 0.0;
			}
	}
	dcell(celpre, al, &v1);
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

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
		printSaveInterstResString(imp_file, rmax, v2, ipen, 1, a);
		/*const char *temp2 = saveInterestingResultString(imp_file);
		printIsee(rmax, v2, ipen, 1, a);
		printf("%s\n", temp2);
	/* ... Refine that cell */

		indic = 1;
		bb[2] = a;
		bb[3] = a;
		bb[4] = a;
		bb[5] = 90.0;
		bb[6] = 90.0;
		bb[7] = 90.0;
		afi[2] = 1.0;
		afi[3] = 1.0;
		afi[4] = 1.0;
		afi[5] = 0.0;
		afi[6] = 0.0;
		afi[7] = 0.0;
		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = a;
		for (int i = 1; i <= 3; ++i) {
			for (int j = 1; j <= 3; ++j) {
				al[i + j * 3 - 4] = 0.0;
			}
		}
		dcell(celpre, al, &v1);
		calcul2(&diff, ihkl, th3, &ncalc, &igc);
		celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
		if (ndat >= 20) {
			cncalc[igc - 1] = (double) ncalc;
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
			ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);

			saveFMFF20(fm20[igc - 1], ff20[igc - 1], ddt, ncalc, imp_file);
		}
		iref = 1;
		goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
		for (int i = 1; i < igc; ++i) {
			if (ifi[i - 1] != ifile) {
				goto L1118;
			}
			vdelt = vgc[igc - 1] / 300.0;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
				goto L1118;
			}
			++nsol[i - 1];
			if (rp[igc - 1] < rp[i - 1]) {
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 1, a);
				}
				km[i - 1] = km[igc - 1];
				vgc[i - 1] = vgc[igc - 1];
				rp[i - 1] = rp[igc - 1];
				cel[i * 6 - 6] = cel[igc * 6 - 6];
				cel[i * 6 - 5] = cel[igc * 6 - 5];
				cel[i * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			goto L1119;
	L1118:
			;
		}
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 1, a);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 1, a);
		}
	L1119:
		;
	} else {
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 1, a);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 1, a);
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
	// tkill = d_mod(&ntried, &c_b1413);
	if ((int) ntried % 30000 == 0) {
		killk(&pressedk);
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
		char temp[80];
		snprintf(temp, sizeof(temp), "\nBest result : a=%lf V=%lf Rp=%lf\n\n", bpar[0], v3, rmin);
		fwrite(temp, strlen(temp), 1, imp_file);
		writeFormattedDate(imp_file);
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
			printf("Hexagonal:    Rp     a      c       V     Nind\n");
		}
		if (nsys[1] == 2) {
			printf("Rhombohedral: Rp     a      c       V     Nind\n");
		}
	}
	rpsmall = 1.0;

	/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 2;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.0;
	ncycles = 500.0;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 120.0;
	if (iverb == 1) {
		if (nsys[1] == 1) {
			saveGridResultsInString("hexagonal", imp_file);
		}
		if (nsys[1] == 2) {
			saveGridResultsInString("rhombohedral", imp_file);
		}
	}

	if (ngrid == 3) {
		pmin = 2.0;
		pmax = 30.0;
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
		vmin = 8.0;
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 4e3f) {
			vmax = 4e3f;
		}
		if (iverb == 1) {
			char temp[43];
			snprintf(temp, sizeof(temp), "\n Max(a,c), V  %lf %lf %lf\n", pma[0], pma[2], vmax);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
	}

	if (iverb == 1) {
		char temp[] = " Rp  Trial number    a      c        V  Nind Icod\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*     READ hkl Miller indices in hex.hkl */

	if (nsys[1] == 2)
	{
		readHklFile("rho", 12, 600, ihh);
	}
	else
	{
		readHklFile("hex", 12, 800, ihh);
	}

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = celpre[0];
	celpre[2] = pmi[2] - spar;
	int ifin = 1;

L1202:

/*     Which parameter to vary ? a or c ? */

	ntriedb = 0.0;
	if (ifin == 1) {
		celpre[0] += spar;
		if (celpre[0] > pma[0]) {
			goto L1216;
		}
		ifin = 0;
		ntried += 1.0;
	}
	celpre[2] += spar;
	if (celpre[2] > pma[2]) {
		celpre[2] = pmi[2] - spar;
		ifin = 1;
		goto L1202;
	}
	ntried += 1.0;
	goto L1204;
L1203:
	del = deltab * (1.0 - ntriedb / cy);
	int i = 3;
	if (randi(&iseed) > .5f) {
		i = 1;
	}
	celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
	ntriedb += 1.0;
L1204:
	celpre[1] = celpre[0];
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.0) {
		goto L1216;
	}
	if (ntriedb != 0.0) {
		goto L1206;
	}
	if (v1 > vmax || v1 < vmin) {
		ntried += -1.0;
		goto L1202;
	}

L1206:
	calcul1(&diff, &diff2);
	if (nmx > ndat10) {
		ntried += -1;
		goto L1202;
	}
	if (ntriedb != 0.0) {
		goto L1214;
	}

/* ... Rp value satisfying ??? */

	if (lhkl >= nmax) {
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
	if (lhkl < nmax) {
		goto L1217;
	}
L1214:
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
		goto L1203;
	}
	ntriedb = 0.0;
	if (rmax >= rmax0[1]) {
		goto L1217;
	}
	if (rmax2 >= .15f) {
		goto L1217;
	}
	ipen = ndat - llhkl;
	if (ipen > nind) {
		goto L1217;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.0;
	if (nr == 1) {
		if (igt > 50.0) {
			if (ntried / igt < 1e3f) {
				if (rmax0[1] > .1f) {
					rmax0[1] -= rmax0[1] * .05f;
					printRmaxReducedString(rmax0[1], imp_file);
				}
			}
		}
	}

	if (igc > 10000) {
		printStopString(imp_file);
		--igc;
		goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = c;
	cel[igc * 6 - 3] = 90.0;
	cel[igc * 6 - 2] = 90.0;
	cel[igc * 6 - 1] = 120.0;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = c;
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
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

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
		printSaveInterstResString(imp_file, rmax, v2, ipen, 2, a, c);
	/* ... Refine that cell */

		indic = 2;
		bb[2] = a;
		bb[3] = a;
		bb[4] = c;
		bb[5] = 90.0;
		bb[6] = 90.0;
		bb[7] = 120.0;
		afi[2] = 1.0;
		afi[3] = 1.0;
		afi[4] = 1.0;
		afi[5] = 0.0;
		afi[6] = 0.0;
		afi[7] = 0.0;
		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = c;
		for (int i = 1; i <= 3; ++i) {
			for (int j = 1; j <= 3; ++j) {
				al[i + j * 3 - 4] = 0.0;
			}
		}
		dcell(celpre, al, &v1);
		calcul2(&diff, ihkl, th3, &ncalc, &igc);
		celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
		if (ndat >= 20) {
			cncalc[igc - 1] = (double) ncalc;
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
			ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
			saveFMFF20(fm20[igc - 1], ff20[igc - 1], ddt, ncalc, imp_file);
		}
		iref = 1;
		goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
		for (i = 1; i < igc; ++i) {
			if (ifi[i - 1] != ifile) {
				goto L1218;
			}
			vdelt = vgc[igc - 1] / 300.0;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
				goto L1218;
			}
			adelt = cel[igc * 6 - 6] / 500.0;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
				goto L1218;
			}
			++nsol[i - 1];
			if (rp[igc - 1] < rp[i - 1]) {
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 2, a, c);
				}
				km[i - 1] = km[igc - 1];
				vgc[i - 1] = vgc[igc - 1];
				rp[i - 1] = rp[igc - 1];
				cel[i * 6 - 6] = cel[igc * 6 - 6];
				cel[i * 6 - 5] = cel[igc * 6 - 5];
				cel[i * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			goto L1219;
	L1218:
			;
		}
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 2, a, c);
		}
	L1219:
		;
	} else {
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 2, a, c);
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
	int tkill = (int) ntried % 30000;
	if (tkill >= 0.0) {
		killk(&pressedk);
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
		char temp[84];
		snprintf(temp, sizeof(temp), "\nBest result : a = %lf Rp = %lf\nBest result : c = %lf V = %lf\n\n", bpar[0], rmin, bpar[2], v3);
		fwrite(temp, strlen(temp), 1, imp_file);
		writeFormattedDate(imp_file);
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


	rpsmall = 1.0;
	printf("Tetragonal:   Rp     a       c        V     Nind\n");

	/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 3;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.0;
	ncycles = 500.0;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 90.0;
	if (iverb == 1) {
		saveGridResultsInString("tetragonal", imp_file);
	}

	if (ngrid == 3) {
		pmin = 2.0;
		pmax = 30.0;
		pmi[0] = pmin;
		pmi[2] = pmin;
		pma[0] = dmax1 * 2.1f;
		if (pma[0] > pmax) {
			pma[0] = pmax;
		}
		pma[2] = dmax1 * 4.0;
		if (pma[2] > pmax) {
			pma[2] = pmax;
		}
		vmin = 8.0;
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 4e3f) {
			vmax = 4e3f;
		}
		if (iverb == 1) {
			char temp[43];
			snprintf(temp, sizeof(temp), "\n Max(a,c), V  %lf %lf %lf\n", pma[0], pma[2], vmax);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
	}

	if (iverb == 1) {
		char temp[] = " Rp  Trial number    a      c        V  Nind Icod\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*     READ hkl Miller indices in tet.hkl */

	readHklFile("tet", 12, 800, ihh);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = celpre[0];
	celpre[2] = pmi[2] - spar;
	ifin = 1;

L1302:

/*     Which parameter to vary ? a or c ? */

	ntriedb = 0.0;
	if (ifin == 1) {
		celpre[0] += spar;
		if (celpre[0] > pma[0]) {
			goto L1316;
		}
		ifin = 0;
		ntried += 1.0;
	}
	celpre[2] += spar;
	if (celpre[2] > pma[2]) {
		celpre[2] = pmi[2] - spar;
		ifin = 1;
		goto L1302;
	}
	ntried += 1.0;
	goto L1304;
L1303:
	del = deltab * (1.0 - ntriedb / cy);
	i = 3;
	if (randi(&iseed) > .5f) {
		i = 1;
	}
	celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
	ntriedb += 1.0;
L1304:
	celpre[1] = celpre[0];
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.0) {
		goto L1316;
	}
	if (ntriedb != 0.0) {
		goto L1306;
	}
	if (v1 > vmax || v1 < vmin) {
		ntried += -1;
		goto L1302;
	}

L1306:
	calcul1(&diff, &diff2);
	if (nmx > ndat10) {
		ntried += -1;
		goto L1302;
	}
	if (ntriedb != 0.0) {
		goto L1314;
	}

/* ... Rp value satisfying ??? */

	if (lhkl >= nmax) {
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
	if (lhkl < nmax) {
		goto L1317;
	}
L1314:
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
		goto L1303;
	}
	ntriedb = 0.0;
	if (rmax >= rmax0[2]) {
		goto L1317;
	}
	if (rmax2 >= .15f) {
		goto L1317;
	}
	ipen = ndat - llhkl;
	if (ipen > nind) {
		goto L1317;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.0;
	if (nr == 1) {
		if (igt > 50.0) {
			if (ntried / igt < 1e3f) {
				if (rmax0[2] > .1f) {
					rmax0[2] -= rmax0[2] * .05f;
					printRmaxReducedString(rmax0[2], imp_file);
				}
			}
		}
	}

	if (igc > 10000) {
		printStopString(imp_file);
		--igc;
		goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = a;
	cel[igc * 6 - 4] = c;
	cel[igc * 6 - 3] = 90.0;
	cel[igc * 6 - 2] = 90.0;
	cel[igc * 6 - 1] = 90.0;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = a;
	celpre[2] = c;
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
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
	supcel(&lhkl, ihkl, cel, &igc, vgc, &c4);
	brav(&lhkl, ihkl, &ibr);
	ib[igc - 1] = ibr;
	a = cel[igc * 6 - 6];
	cel[igc * 6 - 5] = a;
	c = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
		printSaveInterstResString(imp_file, rmax, v2, ipen, 2, a, c);

	/* ... Refine that cell */

		indic = 2;
		bb[2] = a;
		bb[3] = a;
		bb[4] = c;
		bb[5] = 90.0;
		bb[6] = 90.0;
		bb[7] = 90.0;
		afi[2] = 1.0;
		afi[3] = 1.0;
		afi[4] = 1.0;
		afi[5] = 0.0;
		afi[6] = 0.0;
		afi[7] = 0.0;
		celpre[0] = a;
		celpre[1] = a;
		celpre[2] = c;
		for (int i = 1; i <= 3; ++i) {
			for (int j = 1; j <= 3; ++j) {
				al[i + j * 3 - 4] = 0.0;
			}
		}
		dcell(celpre, al, &v1);
		calcul2(&diff, ihkl, th3, &ncalc, &igc);
		celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
		if (ndat >= 20) {
			cncalc[igc - 1] = (double) ncalc;
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
			ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
			saveFMFF20(fm20[igc - 1], ff20[igc - 1], ddt, ncalc, imp_file);
		}
		iref = 1;
		goto L5000;
	}

/* Test if cell already found */

	if (igc > 1) {
		for (i = 1; i < igc; ++i) {
			if (ifi[i - 1] != ifile) {
				goto L1318;
			}
			vdelt = vgc[igc - 1] / 300.0;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
				goto L1318;
			}
			adelt = cel[igc * 6 - 6] / 500.0;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
				goto L1318;
			}
			++nsol[i - 1];
			if (rp[igc - 1] < rp[i - 1]) {
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 2, a, c);
				}
				km[i - 1] = km[igc - 1];
				vgc[i - 1] = vgc[igc - 1];
				rp[i - 1] = rp[igc - 1];
				cel[i * 6 - 6] = cel[igc * 6 - 6];
				cel[i * 6 - 5] = cel[igc * 6 - 5];
				cel[i * 6 - 4] = cel[igc * 6 - 4];
				}
			--igc;
			goto L1319;
	L1318:
			;
		}
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 2, a, c);
		}
	L1319:
		;
	} else {
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 2, a, c);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 2, a, c);
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
	tkill = (int) ntried % 30000;
	if (tkill >= 0.0) {
		killk(&pressedk);
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
		char temp[84];
		snprintf(temp, sizeof(temp), "\nBest result : a = %lf Rp = %lf\nBest result : c = %lf V = %lf\n\n", bpar[0], rmin, bpar[2], v3);
		fwrite(temp, strlen(temp), 1, imp_file);
		writeFormattedDate(imp_file);
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


	rpsmall = 1.0;
	printf("Orthorhombic: Rp     a       b       c        V     Nind");
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 4;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.0;
	ncycles = 1e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[4] = 90.0;
	celpre[5] = 90.0;
	if (iverb == 1) {
		saveGridResultsInString("orthorhombic", imp_file);
	}

	if (ngrid == 3) {
		pmin = 2.0;
		pmax = 20.0;
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
		vmin = 8.0;
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 2e3f) {
			vmax = 2e3f;
		}
		if (iverb == 1) {
			char temp[53];
			snprintf(temp, sizeof(temp), " Max (a,b,c) V %lf %lf %lf %lf\n\n", pma[0], pma[1], pma[2], vmax);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
	}

	if (iverb == 1) {
		char temp[] = " Rp  Trial number    a   b   c        V  Nind Icod\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*     READ hkl Miller indices in ort.hkl */

	readHklFile("ort", 20, 1000, ihh);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = pmi[1] - spar;
	celpre[2] = pmi[2] - spar;
	ifin = 1;
	int ifin2 = 1;

L1402:

/*     Which parameter to vary ? a or b or c ? */

	ntriedb = 0.0;
	if (ifin == 1) {
		celpre[0] += spar;
		printf("  a = %lf", celpre[0]);
		if (celpre[0] > pma[0]) {
			goto L1416;
		}
		ifin = 0;
		ntried += 1.0;
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
	ntried += 1.0;
	celpre[2] += spar;
	if (celpre[2] > pma[2]) {
		celpre[2] = pmi[2] - spar;
		ifin2 = 1;
		goto L1402;
	}
	ntried += 1.0;
	goto L1404;
L1403:
	del = deltab * (1.0 - ntriedb / cy);
	x = randi(&iseed);
	if (x >= 0.0 && x < .33333f) {
		i = 1;
	}
	if (x >= .33333f && x < .66666f) {
		i = 2;
	}
	if (x >= .66666f && x <= 1.0) {
		i = 3;
	}
	celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
	ntriedb += 1.0;
L1404:
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.0) {
		goto L1416;
	}
	if (ntriedb != 0.0) {
		goto L1406;
	}
	if (v1 > vmax || v1 < vmin) {
		ntried += -1.0;
		goto L1402;
	}

L1406:
	calcul1(&diff, &diff2);
	if (nmx > ndat10) {
		ntried += -1;
		goto L1402;
	}
	if (ntriedb != 0.0) {
		goto L1414;
	}

/* ... Rp value satisfying ??? */

	if (lhkl >= nmax) {
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
	if (lhkl < nmax) {
		goto L1417;
	}
L1414:
	if (diff <= rmax) {
		llhkl = lhkl;
		rmax = diff;
		rmax2 = diff2;
		a = celpre[0];
		b = celpre[1];
		c = celpre[2];
		v2 = v1;
		if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[1] = b;
			bpar[2] = c;
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
	ntriedb = 0.0;
	if (rmax >= rmax0[3]) {
		goto L1417;
	}
	if (rmax2 >= .15f) {
		goto L1417;
	}
	ipen = ndat - llhkl;
	if (ipen > nind) {
		goto L1417;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.0;
	if (nr == 1) {
		if (igt > 50.0) {
			if (ntried / igt < 1e4f) {
				if (rmax0[3] > .1f) {
					rmax0[3] -= rmax0[3] * .05f;
					printRmaxReducedString(rmax0[3], imp_file);
				}
			}
		}
	}

	if (igc > 10000) {
		printStopString(imp_file);
		--igc;
		goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = b;
	cel[igc * 6 - 4] = c;
	cel[igc * 6 - 3] = 90.0;
	cel[igc * 6 - 2] = 90.0;
	cel[igc * 6 - 1] = 90.0;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = b;
	celpre[2] = c;
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
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
	supcel(&lhkl, ihkl, cel, &igc, vgc, &c1);
	brav(&lhkl, ihkl, &ibr);
	ib[igc - 1] = ibr;
	a = cel[igc * 6 - 6];
	b = cel[igc * 6 - 5];
	c = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
		printSaveInterstResString(imp_file, rmax, v2, ipen, 3, a, b, c);

	/* ... Refine that cell */

		indic = 0;
		bb[2] = a;
		bb[3] = b;
		bb[4] = c;
		bb[5] = 90.0;
		bb[6] = 90.0;
		bb[7] = 90.0;
		afi[2] = 1.0;
		afi[3] = 1.0;
		afi[4] = 1.0;
		afi[5] = 0.0;
		afi[6] = 0.0;
		afi[7] = 0.0;
		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c;
		for (int i = 1; i <= 3; ++i) {
			for (int j = 1; j <= 3; ++j) {
	/* L1410: */
			al[i + j * 3 - 4] = 0.0;
			}
		}
		dcell(celpre, al, &v1);
		calcul2(&diff, ihkl, th3, &ncalc, &igc);
		celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
		if (ndat >= 20) {
			cncalc[igc - 1] = (double) ncalc;
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
			ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
			saveFMFF20(fm20[igc - 1], ff20[igc - 1], ddt, ncalc, imp_file);
		}
		iref = 1;
		goto L5000;
		}

/* Test if cell already found */

	if (igc > 1) {
		for (int i = 1; i < igc; ++i) {
			if (ifi[i - 1] != ifile) {
				goto L1418;
			}
			vdelt = vgc[igc - 1] / 300.0;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
				goto L1418;
			}
			adelt = cel[igc * 6 - 6] / 500.0;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			na = 0;
			if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
				na = 1;
			}
			bdelt = cel[igc * 6 - 5] / 500.0;
			bp = cel[igc * 6 - 5] + bdelt;
			bm = cel[igc * 6 - 5] - bdelt;
			nb = 0;
			if (cel[i * 6 - 6] > bp || cel[i * 6 - 6] < bm) {
				nb = 1;
			}
			cdelt = cel[igc * 6 - 4] / 500.0;
			cp = cel[igc * 6 - 4] + cdelt;
			cm = cel[igc * 6 - 4] - cdelt;
			nc = 0;
			if (cel[i * 6 - 6] > cp || cel[i * 6 - 6] < cm) {
				nc = 1;
			}
			if (na == 1 && nb == 1 && nc == 1) {
				goto L1418;
			}
			++nsol[i - 1];
			if (rp[igc - 1] < rp[i - 1]) {
				if (isee == 1) {
				printIsee(rmax, v2, ipen, 3, a, b, c);
				}
				km[i - 1] = km[igc - 1];
				vgc[i - 1] = vgc[igc - 1];
				rp[i - 1] = rp[igc - 1];
				cel[i * 6 - 6] = cel[igc * 6 - 6];
				cel[i * 6 - 5] = cel[igc * 6 - 5];
				cel[i * 6 - 4] = cel[igc * 6 - 4];
			}
			--igc;
			goto L1419;
	L1418:
			;
		}
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 3, a, b, c);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 3, a, b, c);
		}
	L1419:
		;
	} else {
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 3, a, b, c);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 3, a, b, c);
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
	tkill = (int) ntried % 30000;
	if (tkill >= 0.0) {
		killk(&pressedk);
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
		char temp[111];
		snprintf(temp, sizeof(temp), "\nBest result : a = %lf Rp = %lf\nBest result : b = %lf\nBest result : c = %lf V = %lf\n\n", bpar[0], rmin, bpar[1], bpar[2], v3);
		fwrite(temp, strlen(temp), 1, imp_file);
		writeFormattedDate(imp_file);
	}
	rmin = rmax;
L1498:
	if (pressedk) {
		goto L5000;
	}


L1500:

	/*    Monoclinic case - would be too long in grid search, but... */


	rpsmall = 1.0;
	printf("Monoclinic:   Rp     a       b       c       bet     V     Nind\n");
/* ------------------------------------------------------------------------- */
/*     Initialisation */

/*      CALL ESP_INIT(ISEED) */

/* ------------------------------------------------------------------------- */
	ifile = 5;
	rmax = rmaxref;
	rmin = rmax;
	ntried = 0.0;
	ncycles = 2e3f;
	cy = ncycles * 1.1f;
	celpre[3] = 90.0;
	celpre[5] = 90.0;
	if (iverb == 1) {
		saveGridResultsInString("monoclinic", imp_file);
	}

	if (ngrid == 3) {
		pmin = 2.0;
		pmax = 20.0;
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
		vmin = 8.0;
		vmax = pma[0] * pma[1] * pma[2];
		if (vmax > 2e3f) {
			vmax = 2e3f;
		}
		if (iverb == 1) {
			char temp[53];
			snprintf(temp, sizeof(temp), " Max (a,b,c) V %lf %lf %lf %lf\n\n", pma[0], pma[1], pma[2], vmax);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
	}

	if (iverb == 1) {
		char temp[] = " Rp  Trial number    a   b   c     bet  V  Nind Icod\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}

/*     READ hkl Miller indices in mon.hkl */
	readHklFile("mon", 20, 1000, ihh);

/* ...  here starts the loop */

	celpre[0] = pmi[0] - spar;
	celpre[1] = pmi[1] - spar;
	celpre[2] = pmi[2] - spar;
	celpre[4] = pmi[4] - sang;
	int ifin1 = 1;
	ifin2 = 1;
	int ifin3 = 1;

L1502:

/*     Which parameter to vary ? a or b or c or bet ? */

	ntriedb = 0.0;
	if (ifin1 == 1) {
		celpre[0] += spar;
		printf("  a = %lf", celpre[0]);
		if (celpre[0] > pma[0]) {
			goto L1516;
		}
		ifin1 = 0;
		ntried += 1.0;
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
	ntried += 1.0;
	if (ifin3 == 1) {
		celpre[2] += spar;
		ifin3 = 0;
	}
	if (celpre[2] > pma[2]) {
		celpre[2] = pmi[2] - spar;
		ifin2 = 1;
		goto L1502;
	}
	ntried += 1.0;
	celpre[4] += sang;
	if (celpre[4] > pma[4]) {
		celpre[4] = pmi[4] - sang;
		ifin3 = 1;
		goto L1502;
	}
	ntried += 1.0;
	goto L1504;
L1503:
	del = deltab * (1.0 - ntriedb / cy);
	deld = deltad * (1.0 - ntriedb / cy);
	x = randi(&iseed);
	if (x >= 0.0 && x < .25f) {
		i = 1;
	}
	if (x >= .25f && x < .5f) {
		i = 2;
	}
	if (x >= .5f && x < .75f) {
		i = 3;
	}
	if (x >= .75f && x <= 1.0) {
		i = 5;
	}
	if (i != 5) {
		celpre[i - 1] = pstartb[i - 1] + del * (randi(&iseed) - .5f) * 2.0;
	} else {
		celpre[i - 1] = pstartb[i - 1] + deld * (randi(&iseed) - .5f) *	2.0;
	}
	ntriedb += 1.0;
L1504:
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	if (celpre[0] > pma[0] && ntriedb == 0.0) {
		goto L1516;
	}
	if (ntriedb != 0.0) {
		goto L1506;
	}
	if (v1 > vmax || v1 < vmin) {
		ntried += -1.0;
		goto L1502;
	}

L1506:
	calcul1(&diff, &diff2);
	if (nmx > ndat10) {
		ntried += -1;
		goto L1502;
	}
	if (ntriedb != 0.0) {
		goto L1514;
	}

/* ... Rp value satisfying ??? */

	if (lhkl >= nmax) {
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
	if (lhkl < nmax) {
		goto L1517;
	}
L1514:
	if (diff <= rmax) {
		llhkl = lhkl;
		rmax = diff;
		rmax2 = diff2;
		a = celpre[0];
		b = celpre[1];
		c = celpre[2];
		bet = celpre[4];
		v2 = v1;
		if (diff < rmin) {
			rmin = diff;
			bpar[0] = a;
			bpar[1] = b;
			bpar[2] = c;
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
	ntriedb = 0.0;
	if (rmax >= rmax0[4]) {
		goto L1517;
	}
	if (rmax2 >= .15f) {
		goto L1517;
	}
	ipen = ndat - llhkl;
	if (ipen > nind) {
		goto L1517;
	}
	++igc;

/*  Test if too much proposals, if yes decrease Rmax by 5% */

	igt += 1.0;
	if (nr == 1) {
		if (igt > 50.0) {
			if (ntried / igt < 1e5f) {
				if (rmax0[4] > .1f) {
					rmax0[4] -= rmax0[4] * .05f;
					printRmaxReducedString(rmax0[4], imp_file);
				}
			}
		}
	}

	if (igc > 10000) {
		printStopString(imp_file);
		--igc;
		goto L5000;
	}
	cel[igc * 6 - 6] = a;
	cel[igc * 6 - 5] = b;
	cel[igc * 6 - 4] = c;
	cel[igc * 6 - 3] = 90.0;
	cel[igc * 6 - 2] = bet;
	cel[igc * 6 - 1] = 90.0;

/* ... Check for supercell */

	celpre[0] = a;
	celpre[1] = b;
	celpre[2] = c;
	celpre[4] = bet;
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
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
	supcel(&lhkl, ihkl, cel, &igc, vgc, &c1);
	a = cel[igc * 6 - 6];
	b = cel[igc * 6 - 5];
	c = cel[igc * 6 - 4];
	v2 = vgc[igc - 1];

/* ... Check for interesting result */

	if (rp[igc - 1] < rmi) {
		printSaveInterstResString(imp_file, rmax, v2, ipen, 4, a, b, c, bet);

	/* ... Refine that cell */

		indic = 0;
		bb[2] = a;
		bb[3] = b;
		bb[4] = c;
		bb[5] = 90.0;
		bb[6] = bet;
		bb[7] = 90.0;
		afi[2] = 1.0;
		afi[3] = 1.0;
		afi[4] = 1.0;
		afi[5] = 0.0;
		afi[6] = 1.0;
		afi[7] = 0.0;
		celpre[0] = a;
		celpre[1] = b;
		celpre[2] = c;
		celpre[4] = bet;
		for (int i = 1; i <= 3; ++i) {
			for (int j = 1; j <= 3; ++j) {
				al[i + j * 3 - 4] = 0.0;
			}
		}
		dcell(celpre, al, &v1);
		calcul2(&diff, ihkl, th3, &ncalc, &igc);
		celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
		if (ndat >= 20) {
			cncalc[igc - 1] = (double) ncalc;
			fm20[igc - 1] = qo[19] / (cncalc[igc - 1] * 2.0 * ddq);
			ff20[igc - 1] = 20.0 / (cncalc[igc - 1] * ddt);
			saveFMFF20(fm20[igc - 1], ff20[igc - 1], ddt, ncalc, imp_file);
		}
		iref = 1;
		goto L5000;
	}

/* Test if cell already found */
	if (igc > 1) {
		for (i = 1; i < igc; ++i) {
			if (ifi[i - 1] != ifile) {
				goto L1518;
			}
			vdelt = vgc[igc - 1] / 300.0;
			vp = vgc[igc - 1] + vdelt;
			vm = vgc[igc - 1] - vdelt;
			if (vgc[i - 1] > vp || vgc[i - 1] < vm) {
				goto L1518;
			}
			bdelt = cel[igc * 6 - 5] / 500.0;
			bp = cel[igc * 6 - 5] + bdelt;
			bm = cel[igc * 6 - 5] - bdelt;
			if (cel[i * 6 - 5] > bp || cel[i * 6 - 5] < bm) {
				goto L1518;
			}
			betdelt = cel[igc * 6 - 2] / 500.0;
			betp = cel[igc * 6 - 2] + betdelt;
			betm = cel[igc * 6 - 2] - betdelt;
			if (cel[i * 6 - 2] > betp || cel[i * 6 - 2] < betm) {
				goto L1518;
			}
			adelt = cel[igc * 6 - 6] / 500.0;
			ap = cel[igc * 6 - 6] + adelt;
			am = cel[igc * 6 - 6] - adelt;
			na = 0;
			if (cel[i * 6 - 6] > ap || cel[i * 6 - 6] < am) {
				na = 1;
			}
			cdelt = cel[igc * 6 - 4] / 500.0;
			cp = cel[igc * 6 - 4] + cdelt;
			cm = cel[igc * 6 - 4] - cdelt;
			nc = 0;
			if (cel[i * 6 - 6] > cp || cel[i * 6 - 6] < cm) {
				nc = 1;
			}
			if (na == 1 && nc == 1) {
				goto L1518;
			}
			++nsol[i - 1];
			if (rp[igc - 1] < rp[i - 1]) {
				if (isee == 1) {
					printIsee(rmax, v2, ipen, 4, a, b, c, bet);
				}
				km[i - 1] = km[igc - 1];
				vgc[i - 1] = vgc[igc - 1];
				rp[i - 1] = rp[igc - 1];
				cel[i * 6 - 6] = cel[igc * 6 - 6];
				cel[i * 6 - 5] = cel[igc * 6 - 5];
				cel[i * 6 - 4] = cel[igc * 6 - 4];
				cel[i * 6 - 2] = cel[igc * 6 - 2];
			}
			--igc;
			goto L1519;
	L1518:
			;
		}
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 4, a, b, c, bet);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 4, a, b, c, bet);
		}
	L1519:
		;
	} else {
		if (iverb == 1) {
			saveIverb(rmax, ntried, v2, ipen, icode, imp_file, 4, a, b, c, bet);
		}
		if (isee == 1) {
			printIsee(rmax, v2, ipen, 4, a, b, c, bet);
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
	tkill = (int) ntried % 30000;
	if (tkill >= 0.0) {
		killk(&pressedk);
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
		char temp[151];
		snprintf(temp, sizeof(temp), "\nBest result : a =    %lf Rp = %lf\nBest result : b =    %lf\nBest result : c =    %lf \nBest result : beta = %lf V = %lf\n\n", bpar[0], rmin, bpar[1], bpar[2], bpar[4], v3);
		fwrite(temp, strlen(temp), 1, imp_file);
		writeFormattedDate(imp_file);
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
	int imem = 0;

/*   Prepare sorted output for CHEKCELL */

	{
		char temp[] = "\n===============================================================================\n  FINAL LIST OF CELL PROPOSALS, sorted by McM20 :\n===============================================================================\n\n      Global list as produced in the .ckm file\n            (IN=number of indexed lines)\n  The correct cell has some chances to be just below\n\nIN  F.o.M.    Volume   V/V1      a        b        c       alpha   beta    gamma   Bravais lattice\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}

	FILE *ckm_file = openFile(file_name, ".ckm", "w+");
	FILE *mcm_file = openFile(file_name, ".mcm", "w+");

	/*  Calculate new F.o.M */

	for (int j = 1; j <= igc; ++j) {
		if (rp[j - 1] < .001f) {
			rp[j - 1] = .001f;
		}
		double x2 = 1.0;
		if (ib[j - 1] == 1) {
			x2 = 4.0;
		}
		if (ib[j - 1] == 2) {
			x2 = 2.0;
		}
		if (ib[j - 1] == 3) {
			x2 = 2.0;
		}
		if (ib[j - 1] == 4) {
			x2 = 2.0;
		}
		if (ib[j - 1] == 5) {
			x2 = 6.0;
		}
		if (ifi[j - 1] == 7) {
			x2 = 6.0;
		}
		double x3 = 1.0;
		if (ifi[j - 1] == 1) {
			x3 = 6.0;
		}
		if (ifi[j - 1] == 2) {
			x3 = 4.0;
		}
		if (ifi[j - 1] == 3) {
			x3 = 4.0;
		}
		if (ifi[j - 1] == 4) {
			x3 = 2.0;
		}
		if (ifi[j - 1] == 7) {
			x3 = 6.0;
		}
		xfom[j - 1] = 1.0 / rp[j - 1] * 100.0 / cncalc[j - 1] * x2 * x3;
	}

	sort(&igc, xfom, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
	char indxprog[] = "McMaille4.00";
	int ibest = ll[igc - 1];
	for (int i = 1; i <= igc; ++i) {
		if (i > 20) {
			goto L20000;
		}
		int j = ll[igc + 1 - i - 1];
		int ipedig = j;
		for (int k = 1; k <= 6; ++k) {
			celpre[k - 1] = cel[k + j * 6 - 7];
		}
		dcell(celpre, al, &v1);
		ql[0] = al[0] * 1e4f;
		ql[1] = al[4] * 1e4f;
		ql[2] = al[8] * 1e4f;
		ql[3] = al[7] * 1e4f;
		ql[4] = al[6] * 1e4f;
		ql[5] = al[3] * 1e4f;

		char bl = getBL(ib, ifi, j);
		char more[11] = "";
		snprintf(more, sizeof(more), "%s", getMore(ifi, j));
		double vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];

		char temp[256] = "", temp_cel[100] = "", temp_ql[100] = "";
		char *cur = temp, *const end = temp - sizeof(temp);
		char *cur_cel = temp_cel, *const end_cel = temp_cel - sizeof(temp_cel);
		char *cur_ql = temp_ql, *const end_ql = temp_ql - sizeof(temp_ql);

		for (int k = 1; k <= 6; ++k)
		{
			cur_cel += snprintf(cur_cel, end_cel-cur_cel, " %.4lf", cel[k + j * 6 - 7]);
			cur_ql += snprintf(cur_ql, end_ql-cur_ql, " %.4lf", ql[k - 1]);
		}

		struct tm *tm_info = getDate();
		char *date;
		size_t size = 40;

		date = (char *)malloc(size * sizeof(char));
		strftime(date, size, "%d %b %Y %H %M %S", tm_info);

		snprintf(temp, sizeof(temp), "%d %.2lf %.3lf%.2lf                                       %s\n", km[j - 1], xfom[j - 1], vgc[j - 1], vr, temp_cel);
		fwrite(temp, strlen(temp), 1, ckm_file);

		snprintf(temp, sizeof(temp), "%d %.2lf %.3lf %.2lf   %s    %c  %s\n\n", km[j - 1], xfom[j - 1], vgc[j - 1], vr, temp_cel, bl,  more);
		writeInFile(temp, imp_file);

		// writeFormattedInFile(imp_file, "%d %.2lf %.3lf %.2lf   %s    %c  %s\n\n", 7, km[j - 1], xfom[j - 1], vgc[j - 1], vr, temp_cel, bl,  more);

		snprintf(temp, sizeof(temp), "%d %.2lf %.3lf%.2lf %c %s %s %d %s %s\n\n", km[j - 1], xfom[j - 1], vgc[j - 1], vr, bl, indxprog, date, ipedig, temp_cel, temp_ql);
		fwrite(temp, strlen(temp), 1, mcm_file);

			//%d %.2lf %.3lf%.2lf %c %s %s %d %s %s\n\n", , , km[j - 1], xfom[j - 1], vgc[j - 1], vr, bl, indxprog, date, ipedig, temp_cel, temp_ql);

		free(date);break;

//		fwrite(temp, strlen(temp), 1, imp_file);
	}
L20000:
	;

	writeInFile("\n===============================================================================\n  FINAL LIST OF CELL PROPOSALS, sorted by F(20) :\n===============================================================================\n\n      Global list   (IN=number of indexed lines)\n  The correct cell has some chances to be just below\n\nIN  F(20)    Volume   V/V1      a        b        c       alpha   beta    gamma   Bravais lattice\n\n", imp_file);

	sort(&igc, ff20, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
	for (int i = 1; i <= igc; ++i) {
		if (i > 20) {
			goto L20002;
		}
		int j = ll[igc + 1 - i - 1];
		int ipedig = j;
		for (int k = 1; k <= 6; ++k) {
			celpre[k - 1] = cel[k + j * 6 - 7];
		}
		dcell(celpre, al, &v1);
		char bl = getBL(ib, ifi, j);
		char more[11];
		snprintf(more, sizeof(more), "%s", getMore(ifi, j));
		double vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];

		{
			char temp[100] = "", temp_cel[60] = "";
			char *cur_cel = temp_cel, *const end_cel = temp_cel - sizeof(temp_cel);

			for (int k = 1; k <= 6; ++k)
			{
				cur_cel += snprintf(cur_cel, end_cel-cur_cel, " %.4lf", cel[k + j * 6 - 7]);
			}

			snprintf(temp, sizeof(temp), "%d %.2lf %.3lf%.2lf   %s    %c  %s\n", km[j - 1], ff20[j - 1], vgc[j - 1], vr, temp_cel, bl, more);
			fwrite(temp, strlen(temp), 1, imp_file);
		}

	}
L20002:
	{
		char temp[] = "\n===============================================================================\n  FINAL LIST OF CELL PROPOSALS, sorted by M(20) :\n===============================================================================\n\n      Global list   (IN=number of indexed lines)\n  The correct cell has some chances to be just below\n\nIN  M(20)    Volume   V/V1      a        b        c       alpha   beta    gamma   Bravais lattice\n\n";
		fwrite(temp, strlen(temp), 1, imp_file);
	}
	sort(&igc, fm20, ll);
	++imem;
	im[imem - 1] = ll[igc - 1];
		for (i = 1; i <= igc; ++i) {
		if (i > 20) {
			goto L20004;
		}
		int j = ll[igc + 1 - i - 1];
		int ipedig = j;
		for (int k = 1; k <= 6; ++k) {
			celpre[k - 1] = cel[k + j * 6 - 7];
		}
		dcell(celpre, al, &v1);
		char bl = getBL(ib, ifi, j);
		char more[11] = "";
		snprintf(more, sizeof(more), "%s", getMore(ifi, j));
		double vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];

		{
			char temp[100] = "", temp_cel[60] = "";
			char *cur_cel = temp_cel, *const end_cel = temp_cel - sizeof(temp_cel);

			for (int k = 1; k <= 6; ++k)
			{
				cur_cel += snprintf(cur_cel, end_cel-cur_cel, " %.4lf", cel[k + j * 6 - 7]);
			}

			snprintf(temp, sizeof(temp), "%d %.2lf %.3lf%.2lf   %s    %c  %s\n", km[j - 1], ff20[j - 1], vgc[j - 1], vr, temp_cel, bl, more);
			fwrite(temp, strlen(temp), 1, imp_file);
		}
	}
L20004:

	writeInFile("\n===============================================================================\n  See for the highest F.o.M. above the cell(s) with highest symmetry, if any\n  (Cubic, hexagonal, etc), they could correspond to the the right solution\n===============================================================================\n\n", imp_file);

/*   Make output for cells sorted by symmetry */

	fclose(ckm_file);
	fclose(mcm_file);


	/*   Make output for cells sorted by symmetry */

	writeInFile("\n  Cells sorted by symmetry\n\n    Rp     Vol     Vol/V1 Ind Nsol    a        b         c      alpha  beta  gamma\n", imp_file);

	for (int i = 1; i < 7; ++i) {
		isyst[i] = 0;
	}
	for (int i = 0; i < igc; ++i) {
		if (ifi[i] == 1) {
			isyst[0] = 1;
		}
		if (ifi[i] == 2) {
			isyst[1] = 1;
		}
		if (ifi[i] == 3) {
			isyst[2] = 1;
		}
		if (ifi[i] == 4) {
			isyst[3] = 1;
		}
		if (ifi[i] == 5) {
			isyst[4] = 1;
		}
		if (ifi[i] == 6) {
			isyst[5] = 1;
		}
		if (ifi[i] == 7) {
			isyst[6] = 1;
		}
	}
	for (int jifi = 1; jifi <= 7; ++jifi) {
		if (isyst[jifi - 1] == 0) {
			goto L2020;
		}
		FILE *general_ckm_file;
		switch (jifi) {
			case 1:
				writeInFile("\n   Cubic cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_cub.ckm", "w+");
				break;
			case 2:
				writeInFile("\n   Hexagonal/trigonal cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_hex.ckm", "w+");
				break;
			case 3:
				writeInFile("\n   Tetragonal cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_tet.ckm", "w+");
				break;
			case 4:
				writeInFile("\n   Orthorhombic cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_ort.ckm", "w+");
				break;
			case 5:
				writeInFile("\n   Monoclinic cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_mon.ckm", "w+");
				break;
			case 6:
				writeInFile("\n   Triclinic cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_tri.ckm", "w+");
				break;
			case 7:
				writeInFile("\n   Rhombohedral cells\n\n", imp_file);
				general_ckm_file = openFile(file_name, "_rho.ckm", "w+");
				break;
		}

		int jjj = 0;
		for (int i = 1; i <= igc; ++i) {
			if (i > 20) {
				goto L20060;
			}
			int j = ll[i - 1];
			if (ifi[j - 1] != jifi) {
				goto L2006;
			}
			++jjj;
			lll[jjj - 1] = j;
			if (rp[j - 1] < .001f) {
				rp[j - 1] = .001f;
			}
			x = 1.0 / rp[j - 1] * 5.0;
			double vr = vgc[j - 1] / vgc[lll[0] - 1];
			char temp[80], temp_cel[60];
			char *cur = temp_cel, *const end = temp_cel - sizeof(temp_cel);

			for (int k = 1; k <= 6; ++k) {
				cur += snprintf(cur, end-cur, " %.4lf", cel[k + j * 6 - 7]);
			}

			snprintf(temp, sizeof(temp), "%.3lf %.3lf %.2lf %d %d%s\n", rp[j - 1], vgc[j - 1], vr, km[j - 1], nsol[j - 1], temp_cel);

			writeInFile(temp, imp_file);

			snprintf(temp, sizeof(temp), "%d %.2lf %.3lf %.2lf                                       %s\n", km[j - 1], x, vgc[j - 1], vr, temp_cel);

			writeInFile(temp, general_ckm_file);
	L2006:
			;
		}
	L20060:
		;
		int nsolmax = nsol[lll[0] - 1];
		if (nsolmax > 5) {
			writeFormattedInFile(imp_file, "\nWARNING - WARNING - WARNING :\nSame solution found Nsol = %d times,\nyou should probably reduce the test numbers...\n\n", 1, nsolmax);
		}
		fclose(general_ckm_file);
	L2020:
		;
	}

/* ... Refine the "best" cell if this was not already done */

/*     READ hkl Miller indices in *.hkl */

	ifile = ifi[ibest - 1];
	switch (ifile) {
	case 1:  readHklFile("cub", 6, 400, ihh);
	case 2:  readHklFile("hex", 12, 800, ihh);
	case 3:  readHklFile("tet", 12, 800, ihh);
	case 4:  readHklFile("ort", 20, 1000, ihh);
	case 5:  readHklFile("mon", 20, 1000, ihh);
	case 6:  readHklFile("tri", 20, 1000, ihh);
	case 7:  readHklFile("rho", 12, 600, ihh);
	}


L50:

/*      IF(IREF.EQ.1)GO TO 5900 */
	writeInFile("\n    \"Best\" cell with largest McM20 :\n    --------------------------------\n\n", imp_file);
	int j = ibest;
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
	afi[2] = 1.0;
	afi[3] = 1.0;
	afi[4] = 1.0;
	afi[5] = 1.0;
	afi[6] = 1.0;
	afi[7] = 1.0;
	if (ifile == 1) {
		afi[5] = 0.0;
	}
	if (ifile == 1) {
		afi[6] = 0.0;
	}
	if (ifile == 1) {
		afi[7] = 0.0;
	}
	if (ifile == 2) {
		afi[5] = 0.0;
	}
	if (ifile == 2) {
		afi[6] = 0.0;
	}
	if (ifile == 2) {
		afi[7] = 0.0;
	}
	if (ifile == 3) {
		afi[5] = 0.0;
	}
	if (ifile == 3) {
		afi[6] = 0.0;
	}
	if (ifile == 3) {
		afi[7] = 0.0;
	}
	if (ifile == 4) {
		afi[5] = 0.0;
	}
	if (ifile == 4) {
		afi[6] = 0.0;
	}
	if (ifile == 4) {
		afi[7] = 0.0;
	}
	if (ifile == 5) {
		afi[5] = 0.0;
	}
	if (ifile == 5) {
		afi[7] = 0.0;
	}
	if (ifile == 7) {
		afi[5] = 0.0;
	}
	if (ifile == 7) {
		afi[6] = 0.0;
	}
	if (ifile == 7) {
		afi[7] = 0.0;
	}
	celpre[0] = bb[2];
	celpre[1] = bb[3];
	celpre[2] = bb[4];
	celpre[3] = bb[5];
	celpre[4] = bb[6];
	celpre[5] = bb[7];
	bb[0] = 0.0;
/* zeropoint after correction... = 0. */
	for (int i = 1; i <= 3; ++i) {
		for (int jj = 1; jj <= 3; ++jj) {
			al[i + jj * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	calcul2(&diff, ihkl, th3, &ncalc, &j);
	celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
	if (ndat >= 20) {
		cncalc[igc - 1] = (double) ncalc;
		fm20[j - 1] = qo[19] / (cncalc[j - 1] * 2.0 * ddq);
		ff20[j - 1] = 20.0 / (cncalc[j - 1] * ddt);
		saveFMFF20(fm20[j - 1], ff20[j - 1], ddt, ncalc, imp_file);
	}

/*   Make a .prf with the best solution */


/* .....START TO GENERATE THE "OBSERVED" PROFILE */

/*   Sum on observed peak positions */

	double xnhkl = (double) ndat;
	double hmin = 100.0;
	for (int i = 0; i < ndat; ++i) {
		th2[i] += bb[0];
		fobs[i] = fobs[i] / sum_f * xnhkl *	50.0;
	}
	sum_f = xnhkl * 50.0;
	double flog = log(2.0) * 4.0;
	double slog = 1.0 / sqrt(flog);
	if (w < 0.0) {
		w = 0.0;
		x = 0.0;
		for (int i = 0; i < nhkl; ++i) {
			w += w1[i];
			x += 1.0;
		}
		w /= x;
	}

/*  Use FWHM 1/2 the entered value... */

/* Computing 2nd power */
	double r1 = w / 2.0;
	w = r1 * r1;


	for (int nm = 1; nm <= nhkl; ++nm) {
		double dd = slabda / (sin(th2[nm - 1] / pi) * 2.0);
	/* Computing 2nd power */
		r1 = dd;
		double sss = 1 / (r1 * r1);

	/*     Needs to calculate halfwidths and positions */

		d[nm - 1] = sqrt(1.0 / sss);
		double sinth = slabda * slabda * sss / 4.0;
		double costh = 1.0 - sinth;
		double tanth = sqrt(sinth / costh);
		theta[nm - 1] = atan(tanth) * pi;
		hw[nm - 1] = sqrt(w);
		if (hmin > hw[nm - 1]) {
			hmin = hw[nm - 1];
		}
		hw[nm - 1] *= slog;
		hw4[nm - 1] = hw[nm - 1] / slog * 2.0;
		bbb[nm - 1] = hw[nm - 1] * hw[nm - 1];
	/* L1802: */
	}
	double dmi = d[nhkl - 1] - .02f;

/*  Step and positions */
/*  The minimum step is FWHMmin/STEPN */

	double step = hmin / 4.0;
	double astep = step * 1e3f;
	double istep = astep;
	astep = (double) istep;
	step = astep / 1e3f;
	pos[0] = step;
	double po = step;
	for (int i = 2; i <= 16000; ++i) {
		po += step;
		if (po > 160.0) {
			goto L1804;
		}
	/* L1803: */
		pos[i - 1] = po;
	}
L1804:
	;
	int npts = i - 1;
/*  Sum at each point */
	for (int k = 1; k <= npts; ++k) {
		yobs[k - 1] = 0.0;
	/*  SUM on all hkl */
		for (int j = 1; j <= nhkl; ++j) {
			double deltap = pos[k - 1] - theta[j - 1];
			if (fabs(deltap) > hw4[j - 1]) {
				goto L1806;
			}
			double delt = deltap * deltap;
			double omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
			yobs[k - 1] += omeg * fobs[j - 1];
	L1806:
			;
		}
		if (yobs[k - 1] < .01f) {
			yobs[k - 1] = 0.0;
		}
	/* L1805: */
	}
	int l = npts + 1;
	for (int k = 1; k <= npts; ++k) {
		--l;
		if (yobs[l - 1] != 0.0) {
			goto L1808;
		}
	}
L1808:
	n2 = l;
	n1 = 1;

	for (int i = 1; i <= 6; ++i) {
		int k = i + 2;
		celpre[i - 1] = bb[k - 1];
	}
	for (int i = 1; i <= 3; ++i) {
		for (int j = 1; j <= 3; ++j) {
			al[i + j * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);

/* ...  Keep only the hkl for d > dmin */

	int jh = 0;
	for (int i = 1; i <= nhkl0; ++i) {
		for (int kk = 1; kk <= 3; ++kk) {
			hh[kk - 1] = (double) ihh[kk + i * 3 - 4];
		}
	/*     HH IS INDICES OF REFLECTION */
		x = 0.0;
		for (int j = 1; j <= 3; ++j) {
			for (int k = j; k <= 3; ++k) {
				x = al[j + k * 3 - 4] * hh[j - 1] * hh[k - 1] + x;
			}
		}
		x = 1 / sqrt(x);
		if (x < dmi) {
			goto L59;
		}
		++jh;
		for (int kk = 1; kk <= 3; ++kk) {
			jhh[kk + jh * 3 - 4] = ihh[kk + i * 3 - 4];
		}
		fcal[jh - 1] = 50.0;
		d[jh - 1] = x;
	/*     X IS D(hkl) FOR REFLECTION HH */
	L59:
		;
	}
	nhkl = jh;

/*   Again, calculate 2-theta, etc. */

	for (int nm = 1; nm <= nhkl; ++nm) {

		double dd = d[nm - 1];
	/* Computing 2nd power */
		r1 = dd;
		double sss = 1 / (r1 * r1);

	/*     Needs to calculate halfwidths and positions again */

		double sinth = slabda * slabda * sss / 4.0;
		double costh = 1.0 - sinth;
		double tanth = sqrt(sinth / costh);
		theta[nm - 1] = atan(tanth) * pi;
		hw[nm - 1] = sqrt(w);
		if (hmin > hw[nm - 1]) {
			hmin = hw[nm - 1];
		}
		hw[nm - 1] *= slog;
		hw4[nm - 1] = hw[nm - 1] / slog * 2.0;
		bbb[nm - 1] = hw[nm - 1] * hw[nm - 1];
	}

/* ...  Calculate best Yobs */

/*  Sum at each point */

	for (int k = 1; k <= nhkl; ++k) {
		somega[k - 1] = 0.0;
		fobs[k - 1] = 0.0;
	}
	for (int k = 1; k <= n2; ++k) {
		ycalc[k - 1] = 0.0;
	/*  SUM on all hkl */
		int kpos = 0;
		double omegt = 0.0;
		for (int j = 1; j <= nhkl; ++j) {
			double deltap = pos[k - 1] - theta[j - 1];
			if (fabs(deltap) > hw4[j - 1]) {
				goto L1812;
			}
			if (kpos == 0) {
				nha[k - 1] = j;
			}
			kpos = 1;
			nhb[k - 1] = j;
			double delt = deltap * deltap;
			double omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
			somega[j - 1] += omeg;
			dump[j - 1] = fcal[j - 1] * omeg;
			ycalc[k - 1] += dump[j - 1];
	L1812:
			;
		}
		if (ycalc[k - 1] == 0.0) {
			goto L1813;
		}
		double yoy = yobs[k - 1] / ycalc[k - 1];
		for (j = nha[k - 1]; j <= nhb[k - 1]; ++j) {
			fobs[j - 1] += dump[j - 1] * yoy;
		}
	L1813:
		;
	}
	for (int k = 1; k <= nhkl; ++k) {
		if (somega[k - 1] == 0.0) {
			goto L1816;
		}
		fobs[k - 1] /= somega[k - 1];
		goto L1815;
	L1816:
		fobs[k - 1] = 0.0;
	L1815:
		;
	}
	for (j = 1; j <= nhkl; ++j) {
		fcal[j - 1] = fobs[j - 1];
	}

/* ...  2 more iterations by Le Bail fit */

	for (int kiter = 1; kiter <= 2; ++kiter) {
	/*  Sum at each point */
		for (int k = 1; k <= nhkl; ++k) {
			somega[k - 1] = 0.0;
	/* L1817: */
			fobs[k - 1] = 0.0;
		}
		for (int k = 1; k <= n2; ++k) {
			ycalc[k - 1] = 0.0;
	/*  SUM on all hkl */
			double omegt = 0.0;
			for (int j = nha[k - 1]; j <= nhb[k - 1]; ++j) {
				double deltap = pos[k - 1] - theta[j - 1];
				if (fabs(deltap) > hw4[j - 1]) {
					goto L1819;
				}
				double delt = deltap * deltap;
				double omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
				somega[j - 1] += omeg;
				dump[j - 1] = fcal[j - 1] * omeg;
				ycalc[k - 1] += dump[j - 1];
		L1819:
				;
			}
			if (ycalc[k - 1] == 0.0) {
				goto L1820;
			}
			double yoy = yobs[k - 1] / ycalc[k - 1];
			for (int j = nha[k - 1]; j <= nhb[k - 1]; ++j) {
				fobs[j - 1] += dump[j - 1] * yoy;
			}
	L1820:
			;
		}
		for (int k = 1; k <= nhkl; ++k) {
			if (somega[k - 1] == 0.0) {
				goto L1823;
			}
			fobs[k - 1] /= somega[k - 1];
			goto L1822;
	L1823:
			fobs[k - 1] = 0.0;
	L1822:
			;
		}
		for (int j = 1; j <= nhkl; ++j) {
			fcal[j - 1] = fobs[j - 1];
		}
	}


	diff = 0.0;
/*  Sum at each point */
	for (int k = 1; k <= n2; ++k) {
		ycalc[k - 1] = 0.0;
		if (yobs[k - 1] == 0.0) {
			goto L1824;
		}
	/*  SUM on all hkl */
		double omegt = 0.0;
		for (int j = nha[k - 1]; j <= nhb[k - 1]; ++j) {
			double deltap = pos[k - 1] - theta[j - 1];
			if (fabs(deltap) > hw4[j - 1]) {
				goto L1825;
			}
			double delt = deltap * deltap;
			double omeg = exp(-delt / bbb[j - 1]) / hw[j - 1];
			ycalc[k - 1] += fcal[j - 1] * omeg;
	L1825:
			;
		}
	L1824:
		;
	}
	double sum_y = 0.0;
	for (int k = 1; k <= n2; ++k) {
		diff += (r1 = yobs[k - 1] - ycalc[k - 1], fabs(r1));
		sum_y += yobs[k - 1];
	}
	diff /= sum_y;

	writeFormattedInFile(imp_file, "\n Final Rp on the .prf = %d\n\n", 1, diff);

/*     Make the .prf */

/* Writing concatenation */
	FILE *prf_file = openFile(file_name, ".prf", "w+");

	int ivers = 8;
	zero = 0.0;

	writeFormattedInFile(prf_file, "%s\n3111   1.0000    %.4lf", 2, text, zero);

	double thmax = pos[n2 - 1];
	double thmin = pos[n1 - 1];
	int icn = nhkl;
	for (i = 1; i <= nhkl; ++i) {
		int i1 = jhh[i * 3 - 3];
		int i2 = jhh[i * 3 - 2];
		int i3 = jhh[i * 3 - 1];

		irefs[i - 1] = ((((i1 + 2432) << 8) + 128 + i2) << 8) + 128 + i3;
	}
	int icz = nhkl;
	npts = n2 - n1 + 1;

	writeFormattedInFile(prf_file, "%.3lf %.3lf %.5lf %d\n", 4, thmax, thmin, step, ivers);

	double amda1 = slabda;
	double amda2 = slabda;
	int npat1 = 1;
	int nvk = 0;

	int nexcrg = 0;
	double excrg = 0.0;

	{
		char *temp_yobs = doubleArrayToString(" %.0lf", yobs, n1-1, n2);
		char *temp_ycalc = doubleArrayToString(" %.0lf", ycalc, n1-1, n2);
		char *temp_irefs = intArrayToString(" %d", irefs, 0, icz);
		char *temp_theta = doubleArrayToString(" %.3lf", theta, 0, icz);


		char temp[100];
		snprintf(temp, sizeof(temp), "%d %d %.5lf %.5lf %.5lf\n%d %d\n%s\n%s\n%s\n%s\n%d\n%.2lf %.2lf\n", npat1, npts, amda1, amda2, zero, icn, nvk, temp_yobs, temp_ycalc, temp_irefs, temp_theta, nexcrg, excrg, excrg);

		writeInFile(temp, prf_file);

		// writeFormattedInFile(prf_file, "%d %d %.5lf %.5lf %.5lf\n%d %d\n%s\n%s\n%s\n%s\n%d\n%.2lf %.2lf\n", 14, npat1, npts, amda1, amda2, zero, icn, nvk, temp_yobs, temp_ycalc, temp_irefs, temp_theta, nexcrg, excrg, excrg);
		// writeFormattedInFile(prf_file, "%d %d %.5lf %.5lf %.5lf\n%d %d\n%s", 7, npat1, npts, amda1, amda2, zero, icn, nvk, temp_yobs);

		free(temp_yobs);
		free(temp_ycalc);
		free(temp_irefs);
		free(temp_theta);

	}

	fclose(prf_file);



	if (igc == 1) {
		goto L6000;
	}

	/*   Sort cells by volume */

	writeInFile("\n===============================================================================\n       CELL PROPOSALS sorted by increasing volume :\n===============================================================================\n\n\n", imp_file);

	sort(&igc, vgc, ll);
	++imem;
	im[imem - 1] = ll[0];

	writeInFile("    Rp     Vol     Vol/V1 Ind Nsol    a        b         c      alpha  beta  gamma\n", imp_file);

	for (int i = 1; i <= igc; ++i) {
		if (i > 20) {
			goto L20030;
		}
		j = ll[i - 1];
		if (rp[j - 1] < .001f) {
			rp[j - 1] = .001f;
		}
		double vr = vgc[j - 1] / vgc[ll[0] - 1];
		char *temp = doubleArrayToString(" %.4lf", cel, 0, 6);

		writeFormattedInFile(imp_file, "%.3lf %.3lf %.2lf %d %d %s\n", 6, rp[j - 1], vgc[j - 1], vr, km[j - 1], nsol[j - 1], temp);

		free(temp);

	}
L20030:
	;
	int nsolmax = nsol[ll[0] - 1];
	if (nsolmax > 5) {
		writeFormattedInFile(imp_file, "\nWARNING - WARNING - WARNING :\nSame solution found Nsol = %d times,\nyou should probably reduce the test numbers...\n\n", 1, nsolmax);
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
	case 1:  readHklFile("cub", 6, 400, ihh);
	case 2:  readHklFile("hex", 12, 800, ihh);
	case 3:  readHklFile("tet", 12, 800, ihh);
	case 4:  readHklFile("ort", 20, 1000, ihh);
	case 5:  readHklFile("mon", 20, 1000, ihh);
	case 6:  readHklFile("tri", 20, 1000, ihh);
	case 7:  readHklFile("rho", 12, 600, ihh);
	}

L1750:

	writeInFile("\n    \"Best\" cell with smallest volume :\n    ---------------------------------\n\n", imp_file);

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
	afi[2] = 1.0;
	afi[3] = 1.0;
	afi[4] = 1.0;
	afi[5] = 1.0;
	afi[6] = 1.0;
	afi[7] = 1.0;
	if (ifile == 1) {
		afi[5] = 0.0;
	}
	if (ifile == 1) {
		afi[6] = 0.0;
	}
	if (ifile == 1) {
		afi[7] = 0.0;
	}
	if (ifile == 2) {
		afi[5] = 0.0;
	}
	if (ifile == 2) {
		afi[6] = 0.0;
	}
	if (ifile == 2) {
		afi[7] = 0.0;
	}
	if (ifile == 3) {
		afi[5] = 0.0;
	}
	if (ifile == 3) {
		afi[6] = 0.0;
	}
	if (ifile == 3) {
		afi[7] = 0.0;
	}
	if (ifile == 4) {
		afi[5] = 0.0;
	}
	if (ifile == 4) {
		afi[6] = 0.0;
	}
	if (ifile == 4) {
		afi[7] = 0.0;
	}
	if (ifile == 5) {
		afi[5] = 0.0;
	}
	if (ifile == 5) {
		afi[7] = 0.0;
	}
	if (ifile == 7) {
		afi[5] = 0.0;
	}
	if (ifile == 7) {
		afi[6] = 0.0;
	}
	if (ifile == 7) {
		afi[7] = 0.0;
	}
	celpre[0] = bb[2];
	celpre[1] = bb[3];
	celpre[2] = bb[4];
	celpre[3] = bb[5];
	celpre[4] = bb[6];
	celpre[5] = bb[7];
	bb[0] = 0.0;
/* zeropoint after correction... = 0. */
	for (int i = 1; i <= 3; ++i) {
		for (int jj = 1; jj <= 3; ++jj) {
			al[i + jj * 3 - 4] = 0.0;
		}
	}
	dcell(celpre, al, &v1);
	calcul2(&diff, ihkl, th3, &ncalc, &j);
	celref(&indic, bb, afi, &lhkl, th3, ihkl, &ddt, &ddq, imp_file);
	if (ndat >= 20) {
		cncalc[j - 1] = (double) ncalc;
		fm20[j - 1] = qo[19] / (cncalc[j - 1] * 2.0 * ddq);
		ff20[j - 1] = 20.0 / (cncalc[j - 1] * ddt);
		saveFMFF20(fm20[j-1], ff20[j-1], ddt, ncalc, imp_file);
	}
L5901:

	if (igc == 1) {
		goto L6000;
	}

/*   Sort cells, the most frequently found first */

	writeInFile("\n===============================================================================\n       CELL PROPOSALS most frequently found :\n===============================================================================\n\n", imp_file);

	++imem;
	im[imem - 1] = ll[igc - 1];

	writeInFile("    Rp     Vol     Vol/V1 Ind Nsol    a        b         c      alpha  beta  gamma\n", imp_file);

	for (int i = 1; i <= igc; ++i) {
		if (i > 20) {
			goto L20040;
		}
		j = ll[igc + 1 - i - 1];
		if (nsol[j - 1] < 2) {
			goto L2004;
		}
		if (rp[j - 1] < .001f) {
			rp[j - 1] = .001f;
		}
		double vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];

		char *temp = doubleArrayToString(" %.4lf", cel, 0, 6);

		writeFormattedInFile(imp_file, "%.3lf %.3lf %.2lf %d %d %s\n", 6, rp[j - 1], vgc[j - 1], vr, km[j - 1], nsol[j - 1], temp);

		free(temp);
	L2004:
		;
	}
L20040:


/*   Sort cells with largest number of peak indexed + small Rp2 */

	writeInFile("\n===============================================================================\nCELLS with small Rp2 + largest number of peak indexed\n===============================================================================\n\nRp2 is on peak indexed only, and width divided by 2,\nwhile Rp is on all peaks, and large width.\n\n", imp_file);

/*      CALL SORT2(IGC,RP2,LL) */

	sort2(&igc, rp2, ll);

	++imem;
	im[imem - 1] = ll[i - 1];
	writeInFile("    Rp     Vol     Vol/V1 Ind Nsol    a        b         c      alpha  beta  gamma\n", imp_file);

	for (int i = 1; i <= igc; ++i) {
		if (i > 20) {
			goto L20050;
		}
		j = ll[i - 1];
		if (rp2[j - 1] < .001f) {
			rp2[j - 1] = .001f;
		}
		if (rp2[j - 1] > .3f) {
			goto L2005;
		}
		if (km2[j - 1] < nmax) {
			goto L2005;
		}
		double vr = vgc[j - 1] / vgc[ll[igc - 1] - 1];

		char *temp = doubleArrayToString(" %.4lf", cel, 0, 6);

		writeFormattedInFile(imp_file, "%.3lf %.3lf %.2lf %d %d %s\n", 6, rp[j - 1], vgc[j - 1], vr, km[j - 1], nsol[j - 1], temp);

		free(temp);
	L2005:
		;
		}
L20050:

/*   Sort associations of two cells with largest number of peak indexed */

	if (rmax0[0] < .5f) {
		goto L6000;
	}

	writeInFile("\n===============================================================================\nDouble cells with largest number of peak indexed\n===============================================================================\n\nWARNING - WARNING - WARNING - WARNING - WARNING\n           This is the two-phase mode\n    It could be better to go back to the lab\n         and try and make a pure sample\n\n", imp_file);

	int ii = 0;
	for (int i = 1; i < igc; ++i) {
		if (rp2[i - 1] > .3f) {
			goto L2009;
		}
		if (km2[i - 1] < nmax) {
			goto L2009;
		}
		for (j = i + 1; j <= igc; ++j) {
			if (rp2[j - 1] > .3f) {
			goto L2008;
			}
			if (rp2[i - 1] + rp2[j - 1] > .4f) {
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
			for (int k = 1; k <= ndat; ++k) {
					km3[ii - 1] = km3[ii - 1] + ind[k + i * 100 - 101] + ind[k + j * 100 - 101];
				if (ind[k + i * 100 - 101] * ind[k + j * 100 - 101] == 1) {
					--km3[ii - 1];
				}
			}
			if (km3[ii - 1] < ndat - 10) {
				--ii;
				goto L2008;
			}
			id1[ii - 1] = i;
			id2[ii - 1] = j;
	L2008:
			;
		}
	L2009:
		;
		}
L2011:
	;
	int igc2 = ii;
	sort3(&igc2, km3, ll2);

	writeInFile("    Rp     Vol     Vol/V1 Ind Nsol    a        b         c      alpha  beta  gamma\n", imp_file);

	double vr = 1.0;

	FILE *two_ckm_file = openFile("two", ".ckm", "w+");

	for (int i = 1; i <= igc2; ++i) {
		if (i > 1000) {
			goto L2012;
		}
		double jj = ll2[igc2 + 1 - i - 1];
		j = id1[ll2[igc2 + 1 - i - 1] - 1];
		if (rp2[j - 1] < .001f) {
			rp2[j - 1] = .001f;
		}
		x = 1.0 / rp2[j - 1] * 5.0;

		char *temp_cel = doubleArrayToString(" %.4lf", cel, 0, 6);

		writeFormattedInFile(imp_file, "%.3lf %.3lf %.2lf %d %d %s\n", 6, rp2[j - 1], vgc[j - 1], vr, km3[j - 1], nsol[j - 1], temp_cel);

		writeFormattedInFile(two_ckm_file, "%d %.2lf %.3lf %.2lf                                       %s\n", 5, km3[j - 1], x, vgc[j - 1], vr, temp_cel);


		j = id2[ll2[igc2 + 1 - i - 1] - 1];
		if (rp2[j - 1] < .001f) {
			rp2[j - 1] = .001f;
		}
		x = 1.0 / rp2[j - 1] * 5.0;

		writeFormattedInFile(imp_file, "%.3lf %.3lf %.2lf %d %d %s\n", 6, rp2[j - 1], vgc[j - 1], vr, km2[j - 1], nsol[j - 1], temp_cel);

		writeFormattedInFile(two_ckm_file, "%d %.2lf %.3lf %.2lf                                       %s\n", 5, km2[j - 1], x, vgc[j - 1], vr, temp_cel);

		free(temp_cel);

	/* L2010: */
	}
L2012:
	fclose(two_ckm_file);


L6000:
	writeInFile("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n===============================================================================         THE SELECTION OF THE \"BEST\" CELL\nbased on McM20, Rp, F(20), M(20), V, high symmetry ?\n           DEPENDS ON YOU, EXCLUSIVELY.\n\n                    However...\n  It is suggested that the correct cell could be :\n===============================================================================\n\nIN  F.o.M.    Volume             a        b        c       alpha   beta    gamma   Bravais lattice\n", imp_file);

/*  Ultimate analysis... */

	for (i = 1; i <= imem; ++i) {
		imn[i - 1] = 1;
	}
	if (imem > 1) {
		int i3 = imem - 1;
		for (i = 1; i <= i3; ++i) {
			if (imn[i - 1] == 0) {
			goto L20008;
			}
			for (int k = i + 1; k <= imem; ++k) {
			if (imn[k - 1] == 0) {
				goto L30008;
			}
			if (im[k - 1] == im[i - 1]) {
				++imn[i - 1];
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
	char bl = getBL(ib, ifi, j);
	{
		char *temp_cel = doubleArrayToString(" %.4lf", cel, 0, 6);

		char temp[128] = "";
		snprintf(temp, sizeof(temp), "%d %.2lf %.3lf   %s    %c  %s\nFound %d time(s) head of the best lists\n", km[j - 1], xfom[j - 1], vgc[j - 1], temp_cel, bl, getMore(ifi, j), imn[0]);
		writeInFile(temp, imp_file);

		free(temp_cel);
	}

	if (imem > 1) {
		double xf1 = xfom[im[0] - 1] / 2.0;
		double imemt = imem;
		for (i = 2; i <= imem; ++i) {
			if (imn[i - 1] == 0) {
			--imemt;
			goto L20011;
			}
			j = im[i - 1];
			if (xfom[j - 1] < xf1) {
			--imemt;
			imn[i - 1] = 0;
			}
	L20011:
			;
		}
		if (imemt > 1) {
			writeInFile("\n\n   Other(s) having some chance :\n", imp_file);

				for (i = 2; i <= imem; ++i) {
				if (imn[i - 1] == 0) {
					goto L20009;
				}
				j = im[i - 1];

				char bl = getBL(ib, ifi, j);
				{
					char *temp_cel = doubleArrayToString(" %.4lf", cel, 0, 6);

					writeFormattedInFile(imp_file, "%d %.2lf %.3lf   %s    %c  %s\nFound %d time(s) head of the best lists\n", 7, km[j - 1], xfom[j - 1], vgc[j - 1], temp_cel, bl, getMore(ifi, j), imn[0]);

					free(temp_cel);
				}
		L20009:
				;
			}
		}
	}
	writeInFile("\n===============================================================================\n\n\n\n\n\n", imp_file);



//==============================================================================
	free(file_name);
	fclose(imp_file);
	// fclose(new_dat_file);

	return 0;
}