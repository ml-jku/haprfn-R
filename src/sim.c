#include <Rinternals.h>

int lengthIntersect(SEXP xS, SEXP yS) {
  const int *x = INTEGER(xS);
  const int *y = INTEGER(yS);
	const int nx = LENGTH(xS);
	const int ny = LENGTH(yS);
	int len = 0;
	int ix = 0;
	int iy = 0;

	while (ix < nx && iy < ny) {
		if (x[ix] == y[iy]) {
			len = len + 1;
			ix = ix + 1;
			iy = iy + 1;
			while (ix < nx && x[ix - 1] == x[ix]) {
				ix = ix + 1;
			}
			while (iy < ny && y[iy - 1] == y[iy]) {
				iy = iy + 1;
			}
		} else if (x[ix] > y[iy]) {
			iy = iy + 1;
		} else {
			ix = ix + 1;
		}
	}

	return len;
}

int lengthUnion(SEXP xS, SEXP yS) {
  const int nx = LENGTH(xS);
	const int ny = LENGTH(yS);
	const int *x = INTEGER(xS);
  const int *y = INTEGER(yS);	
  int len = 0;
	int ix = 0;
	int iy = 0;

	while (ix < nx && iy < ny) {
		if (x[ix] == y[iy]) {
			len = len + 1;
			ix = ix + 1;
			iy = iy + 1;
			while (ix < nx && x[ix - 1] == x[ix]) {
				ix = ix + 1;
			}
			while (iy < ny && y[iy - 1] == y[iy]) {
				iy = iy + 1;
			}
		} else {
			len = len + 1;
		}
	}

	while (ix < nx) {
		len = len + 1;
		ix = ix + 1;
		while (ix < nx && x[ix - 1] == x[ix]) {
				ix = ix + 1;
			}
	}
	
	while (iy < ny) {
		len = len + 1;
		iy = iy + 1;
		while (iy < ny && y[iy - 1] == y[iy]) {
				iy = iy + 1;
			}
	}

	return len;
}

SEXP similarityMeasure(SEXP xS, SEXP yS, SEXP simvS, SEXP minInterS) {
  const int nx = LENGTH(xS);
  const int ny = LENGTH(yS);
	const int len = lengthIntersect(xS, yS);
  const int minInter = INTEGER(minInterS)[0];
  const int simv = INTEGER(simvS)[0];

  double ret = 0.0;
	if (len >= minInter) {
		switch(simv) {
			case 0: ret = len / lengthUnion(xS, yS); break;
			case 1: ret = 2 * len / (nx + ny); break;
			case 2: ret = len / (nx > ny ? ny : nx); break;
			case 3: ret = len / (nx > ny ? nx : ny); break;
		}
	}
  
  SEXP result = PROTECT(allocVector(REALSXP, 1));
  REAL(result)[0] = ret;
  UNPROTECT(1);

	return result;
}
