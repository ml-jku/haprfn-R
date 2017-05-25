#include <Rinternals.h>
#include <stdio.h>
#include "htslib/vcf.h"

int hello() {
return 1;
}

/*
- samplesPerFeature from the fabia package:
 basically builds up the sparse matrix the other way round.
 The output is a list with elements: sL (List with one element per feature: each element is a vector of samples where the feature is not zero.)
 nsL Vector of feature length containing number of samples having a non-zero feature value.
*/
/*SEXP samplesPerFeature(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS) {
  // open the sparse matrix file
  const char *file_name = CHAR(STRING_ELT(file_nameS, 0));


  FILE* file = fopen(file_name, "r");
  if (file == NULL) {
    Rprintf("File %s not found! Stop.\n", file_name);
    return R_NilValue;
  }

  unsigned int nnz, m;

  int ret = fscanf(file, "%u\n", &nnz);
  if (ret < 1) {
    Rprintf("Wrong file format. Was expecting an unsigned integer (NNZ).\n");
    return R_NilValue;
  }

  ret = fscanf(file, "%u\n", &m);
  if (ret < 1) {
    Rprintf("Wrong file format. Was expecting an unsigned integer (M).\n");
    return R_NilValue;
  }

  SEXP entries = PROTECT(allocVector(REALSXP, nnz));
  SEXP columns = PROTECT(allocVector(REALSXP, nnz));
  SEXP row_indices = PROTECT(allocVector(REALSXP, m + 1));
  SEXP out = PROTECT(allocVector(VECSXP, 3));

  doubl
  while ((ret = f) {

}

SEXP names = PROTECT(allocVector(STRSXP, 3));
             SET_STRING_ELT(names, 0, mkChar("entries"));
             SET_STRING_ELT(names, 1, mkChar("columns"));
             SET_STRING_ELT(names, 2, mkChar("row_indices"));

             SET_VECTOR_ELT(out, 0, entries);
             SET_VECTOR_ELT(out, 1, columns);
             SET_VECTOR_ELT(out, 2, row_indices);
             setAttrib(out, R_NamesSymbol, names);

             UNPROTECT(5);
             return out;
}

/* - readSamplesSpfabia from the fabia package: extract the information for specific samples from the sparse matrix. 
So in our case a simple subsetting. The output originally is a dense matrix. 
But it is up to you whether you prefer a sparse matrix, 
but then you probably would have to re-write a lot of functions that use the output of this function. */

/*SEXP readSamplesSpfabic(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS) {


  FILE *pFile;


  char sst[200];


  int hpp = 0, ig = 0, jg = 0, samp = 0, ret = 0;

  int  i = 0, j = 0, i1 = 0, i2 = 0, n = 0, nn = 0;

  double fs = 0.0;





  const char *file_name = CHAR(STRING_ELT(file_nameS, 0));


  int *xa;
  int **xind;
  double **xval;

  double lowerB = (double)(REAL(lowerBS)[0]);
  double upperB = (double)(REAL(upperBS)[0]);

  int nsamp = length(samplesS);

  int *samples =  INTEGER(samplesS);
  if (samples[0] > 0) {
    samp = 0;

  } else {
    samp = -1;
  }

  sst[0] = 0;
  strcat(sst, file_name);
  strcat(sst, ".txt");

  pFile = fopen (sst, "r");

  if (pFile == NULL) {
    Rprintf("File >%s< not found! Stop.\n", sst);
    return R_NilValue;
  }

  ret = fscanf(pFile, "%d\n", &nn);
  if (ret < 1) {
    Rprintf("Wrong file format.\n");
    return R_NilValue;
  }

  if (!(nn > 0)) {
    fclose (pFile);
    Rprintf("Wrong file format (sparse file format required)! Stop.\n");
    return R_NilValue;
  }

  ret = fscanf(pFile, "%d\n", &n);
  if (ret < 1) {
    Rprintf("Wrong file format.\n");
    return R_NilValue;
  }

  if (!(n > 0)) {
    fclose (pFile);
    Rprintf("Wrong file format (sparse file format required)! Stop.\n");
    return R_NilValue;
  }

  if (samp < 0) {

    xa = (int *) R_Calloc(nn, int);
    xind = (int **) R_Calloc(nn, int *);
    xind[0] = R_Calloc((long) nn * n, int);
    for (i = 0; i < nn ; i++)
    {
      xind[i] =  xind[0] + i * n;
    }
    xval = (double **) R_Calloc(nn, double *);
    xval[0] = R_Calloc((long) nn * n, double);
    for (i = 0; i < nn; i++)
    {
      xval[i] = xval[0] + i * n;
    }


    for (i = 0; i < nn; i ++)
    {
      ret = fscanf(pFile, "%d\n", &ig);
      xa[i] = ig;
      for (j = 0; j <  ig; j ++) {
        ret = fscanf(pFile, "%d", &hpp);
        xind[i][j] = hpp;
      }
      ret = fscanf(pFile, "\n");
      for (j = 0; j < ig; j ++) {
        ret = fscanf(pFile, "%lf", &fs);
        xval[i][j] = fs;
      }
      ret = fscanf(pFile, "\n");
    }

  } else {
    xa = (int *) R_Calloc(nsamp, int);
    xind = (int **) R_Calloc(nsamp, int *);
    xind[0] = R_Calloc((long) nsamp * n, int);
    for (i = 0; i < nsamp ; i++)
    {
      xind[i] =  xind[0] + i * n;
    }
    xval = (double **) R_Calloc(nsamp, double *);
    xval[0] = R_Calloc((long) nsamp * n, double);
    for (i = 0; i < nsamp; i++)
    {
      xval[i] = xval[0] + i * n;
    }

    for (i = 0; i < nn; i ++)
    {
      if ((samples[samp] - 1) > i)
      {
        ret = fscanf(pFile, "%d\n", &ig);
        for (j = 0; j <  ig; j ++) {
          ret = fscanf(pFile, "%d", &hpp);
        }
        ret = fscanf(pFile, "\n");
        for (j = 0; j < ig; j ++) {
          ret = fscanf(pFile, "%lf", &fs);
        }
        ret = fscanf(pFile, "\n");

      } else {
        ret = fscanf(pFile, "%d\n", &ig);
        xa[samp] = ig;
        for (j = 0; j <  ig; j ++) {
          ret = fscanf(pFile, "%d", &hpp);
          xind[samp][j] = hpp;
        }
        ret = fscanf(pFile, "\n");
        for (j = 0; j < ig; j ++) {
          ret = fscanf(pFile, "%lf", &fs);
          xval[samp][j] = fs;
        }
        ret = fscanf(pFile, "\n");
        samp++;
        if (samp == nsamp) break;
      }

    }
    if (samp != nsamp)
    {
      Rprintf("Only %d of %d samples found! Some sample numbers are too large. Continue.\n", samp, nsamp);
    }
    nn = samp;
    Rprintf("Using %d samples!\n", samp);
  }
  fclose (pFile);

  int *Psi = R_Calloc(n, int);
  double *XX = R_Calloc(n, double);

  for (i1 = 0; i1 < n; i1++)
  {
    XX[i1] = 0.0;
  }

  for (i2 = 0; i2 < nn; i2++) {
    for (ig = 0; ig < xa[i2]; ig++) {
      XX[xind[i2][ig]] += xval[i2][ig];
    }
  }


  for (i2 = 0; i2 < n; i2++) {
    if ((XX[i2] > lowerB) && (XX[i2] < upperB)) {
      Psi[i2] = 0;
    } else {
      Psi[i2] = 1;
    }
  }


  for (i2 = 0; i2 < nn; i2++) {
    for (ig = 0, jg = 0; (ig + jg) < xa[i2];) {
      if ( Psi[xind[i2][ig + jg]] == 0 )
      {
        if (jg > 0) {
          xval[i2][ig] = xval[i2][ig + jg];
          xind[i2][ig] = xind[i2][ig + jg];
        }
        ig++;
      } else
      {
        jg++;
      }
    }
    xa[i2] -= jg;
  }






  SEXP X_n;
  PROTECT(X_n = allocMatrix(REALSXP, n, nn));

  for (i1 = 0; i1 < nn; i1++)
  {
    for (i2 = 0; i2 < n; i2++)
      REAL(X_n)[i2 + n * i1] = 0.0;
  }
  for (i2 = 0; i2 < nn; i2++) {
    for (ig = 0; ig < xa[i2]; ig++) {
      REAL(X_n)[xind[i2][ig] + n * i2] = xval[i2][ig];
    }
  }

  R_Free (xa );
  R_Free (xind[0]);
  R_Free (xind );
  R_Free (xval[0]);
  R_Free (xval );

  R_Free (XX);
  R_Free (Psi);


  SEXP namesRET;
  PROTECT(namesRET = allocVector(STRSXP, 1));
  SET_STRING_ELT(namesRET, 0, mkChar("X"));

  SEXP RET;
  PROTECT(RET = allocVector(VECSXP, 1));
  SET_VECTOR_ELT(RET, 0, X_n);
  setAttrib(RET, R_NamesSymbol, namesRET);
  UNPROTECT(3);
  return (RET);

}

/*




  FILE *pFile;


  char sst[200];


  int hpp = 0, ig = 0, samp = 0;



  int  i = 0, j = 0, i1 = 0, i2 = 0, n = 0, nn = 0, ret = 0;

  double fs = 0.0;

  const char *file_name = CHAR(STRING_ELT(file_nameS, 0));


  int *xa;
  int **xind;
  double **xval;

  double lowerB = (double)(REAL(lowerBS)[0]);
  double upperB = (double)(REAL(upperBS)[0]);

  int nsamp = length(samplesS);

  int *samples =  INTEGER(samplesS);
  if (samples[0] > 0) {
    samp = 0;

  } else {
    samp = -1;
  }

  sst[0] = 0;
  strcat(sst, file_name);
  strcat(sst, ".txt");

  pFile = fopen (sst, "r");

  if (pFile == NULL) {
    Rprintf("File >%s< not found! Stop.\n", sst);
    return R_NilValue;
  }

  ret = fscanf(pFile, "%d\n", &nn);
  if (ret < 1) {
    Rprintf("Wrong file format.\n");
    return R_NilValue;
  }

  if (!(nn > 0)) {
    fclose (pFile);
    Rprintf("Wrong file format (sparse file format required)! Stop.\n");
    return R_NilValue;
  }

  ret = fscanf(pFile, "%d\n", &n);
  if (ret < 1) {
    Rprintf("Wrong file format.\n");
    return R_NilValue;
  }

  if (!(n > 0)) {
    fclose (pFile);
    Rprintf("Wrong file format (sparse file format required)! Stop.\n");
    return R_NilValue;
  }

  SEXP xLL, PhiA;

  PROTECT(xLL = allocVector(VECSXP, n));
  PROTECT(PhiA = allocVector(INTSXP, n));
  int *Phi = INTEGER(PhiA);
  int *Psi = R_Calloc(n, int);
  double *XX = R_Calloc(n, double);


  if (samp < 0) {

    xa = (int *) R_Calloc(nn, int);
    xind = (int **) R_Calloc(nn, int *);
    xind[0] = R_Calloc((long) nn * n, int);
    for (i = 0; i < nn ; i++)
    {
      xind[i] =  xind[0] + i * n;
    }
    xval = (double **) R_Calloc(nn, double *);
    xval[0] = R_Calloc((long) nn * n, double);
    for (i = 0; i < nn; i++)
    {
      xval[i] = xval[0] + i * n;
    }


    for (i = 0; i < nn; i ++)
    {
      ret = fscanf(pFile, "%d\n", &ig);
      xa[i] = ig;
      for (j = 0; j <  ig; j ++) {
        ret = fscanf(pFile, "%d", &hpp);
        xind[i][j] = hpp;
      }
      ret = fscanf(pFile, "\n");
      for (j = 0; j < ig; j ++) {
        ret = fscanf(pFile, "%lf", &fs);
        xval[i][j] = fs;
      }
      ret = fscanf(pFile, "\n");
    }

  } else {
    xa = (int *) R_Calloc(nsamp, int);
    xind = (int **) R_Calloc(nsamp, int *);
    xind[0] = R_Calloc((long) nsamp * n, int);
    for (i = 0; i < nsamp ; i++)
    {
      xind[i] =  xind[0] + i * n;
    }
    xval = (double **) R_Calloc(nsamp, double *);
    xval[0] = R_Calloc((long) nsamp * n, double);
    for (i = 0; i < nsamp; i++)
    {
      xval[i] = xval[0] + i * n;
    }

    for (i = 0; i < nn; i ++)
    {
      if ((samples[samp] - 1) > i)
      {
        ret = fscanf(pFile, "%d\n", &ig);
        for (j = 0; j <  ig; j ++) {
          ret = fscanf(pFile, "%d", &hpp);
        }
        ret = fscanf(pFile, "\n");
        for (j = 0; j < ig; j ++) {
          ret = fscanf(pFile, "%lf", &fs);
        }
        ret = fscanf(pFile, "\n");

      } else {
        ret = fscanf(pFile, "%d\n", &ig);
        xa[samp] = ig;
        for (j = 0; j <  ig; j ++) {
          ret = fscanf(pFile, "%d", &hpp);
          xind[samp][j] = hpp;
        }
        ret = fscanf(pFile, "\n");
        for (j = 0; j < ig; j ++) {
          ret = fscanf(pFile, "%lf", &fs);
          xval[samp][j] = fs;
        }
        ret = fscanf(pFile, "\n");
        samp++;
        if (samp == nsamp) break;
      }

    }
    if (samp != nsamp)
    {
      Rprintf("Only %d of %d samples found! Some sample numbers are too large. Continue.\n", samp, nsamp);
    }
    nn = samp;
    Rprintf("Using %d samples!\n", samp);
  }
  fclose (pFile);



  for (i1 = 0; i1 < n; i1++)
  {
    XX[i1] = 0.0;
    Psi[i1] = 0;
    Phi[i1] = 0;
  }


  for (i2 = 0; i2 < nn; i2++) {
    for (ig = 0; ig < xa[i2]; ig++) {
      XX[xind[i2][ig]] += xval[i2][ig];
      Psi[xind[i2][ig]]++;
    }
  }



  for (i2 = 0; i2 < n; i2++) {
    if ((XX[i2] > lowerB) && (XX[i2] < upperB) && (Psi[i2] > 0)) {
      SET_VECTOR_ELT(xLL, i2, allocVector(INTSXP, Psi[i2]));
    } else {
      Psi[i2] = 0;
      SET_VECTOR_ELT(xLL, i2, allocVector(INTSXP, 1));
      INTEGER(VECTOR_ELT(xLL, i2))[0] = 0;

    }
  }



  for (i2 = 0; i2 < nn; i2++) {
    for (ig = 0; ig < xa[i2]; ig++) {
      if (Psi[xind[i2][ig]] > 0) {
        INTEGER(VECTOR_ELT(xLL, xind[i2][ig]))[Phi[xind[i2][ig]]] = i2 + 1;
        Phi[xind[i2][ig]]++;
      }
    }
  }


  R_Free (xa );
  R_Free (xind[0]);
  R_Free (xind );
  R_Free (xval[0]);
  R_Free (xval );

  R_Free (XX);
  R_Free (Psi);


  SEXP namesRET;
  PROTECT(namesRET = allocVector(STRSXP, 2));
  SET_STRING_ELT(namesRET, 0, mkChar("sL"));
  SET_STRING_ELT(namesRET, 1, mkChar("nsL"));

  SEXP RET;
  PROTECT(RET = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(RET, 0, xLL);
  SET_VECTOR_ELT(RET, 1, PhiA);
  setAttrib(RET, R_NamesSymbol, namesRET);
  UNPROTECT(4);
  return (RET);

}*/
