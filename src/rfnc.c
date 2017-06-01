#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#define GET_MACRO(_1,_2,_3,_4,NAME,...) NAME

#define _ERROR(file, str) fclose(file); Rprintf(str); return -1;

#define READ(...) GET_MACRO(__VA_ARGS__, READ4, READ3)(__VA_ARGS__)
#define READ4(from, template, variable, str) \
    if (fscanf(from, template, variable) < 1) { _ERROR(from, str) }
#define READ3(from, template, variable) \
    READ4(from, template, variable, "Wrong file format. \n")

#define CHECK(expression, file) if (!(expression)) { _ERROR(file, "Error: " #expression " is not true!") }


/********************************************
 * COMMON FUNCTIONS                         *
 ********************************************/

void _setMemoryDouble(double* array, const double value, const int n) {
  for (int i = 0; i < n; i++) {
    array[i] = value;
  }
}

void _setMemoryInt(int* array, const int value, const int n) {
  for (int i = 0; i < n; i++) {
    array[i] = value;
  }
}

int _readAllSamplesAndCalculateColSums(FILE* file, int* row_ptr, int* col_ind, double* val, 
    const int nrow, const int nnz, const double lowerB, const double upperB, 
    double** col_sums, int* ncol) {
  // read row pointers
  for (int i = 0; i < nrow + 1; i++) {
    READ(file, "%d", (row_ptr + i));
  }

  // read column indices
  int max_col = -1;
  for (int i = 0; i < nnz; i++) {
    READ(file, "%d", (col_ind + i));
    if (col_ind[i] >= max_col) {
      max_col = col_ind[i];
    }
  }
  *ncol = max_col + 1;

  *col_sums = (double*) R_alloc(*ncol, sizeof(double));
  _setMemoryDouble(*col_sums, 0.0, *ncol);

  // read values and calculate column sums
  for (int i = 0; i < nnz; i++) {
    READ(file, "%lf", (val + i));
    (*col_sums)[col_ind[i]] += val[i];
  }

  fclose(file);

  return 0;
}

int _readFilteredSamplesAndCalculateColSums(FILE* file, int** row_ptr, int* col_ind, double* val, 
    int* nrow, int* nnz, const double lowerB, const double upperB, const int* samples, 
    const int nsamp, double** col_sums, int* ncol) {
  int cur_samp = 0;
  
  int* new_row_ptr = (int*) R_alloc(nsamp + 1, sizeof(int));
  new_row_ptr[0] = 0;
  
  // read row pointers and build the filtered row pointers
  for (int i = 0; i < *nrow + 1; i++) {
    READ(file, "%d", (*row_ptr + i));
    if (cur_samp < nsamp && samples[cur_samp] == i) {
      new_row_ptr[cur_samp + 1] = new_row_ptr[cur_samp] + (*row_ptr)[i] - (*row_ptr)[i - 1];
      cur_samp++;
    }
  }

  int cur_row = 0;
  int max_col = -1;
  cur_samp = 0;
  
  int int_trash;
  double double_trash;

  int it = 0;
  int cur_val = 0;

  // read column indices
  while (cur_samp < nsamp) {
    
    // skip empty rows
    while (it == (*row_ptr)[cur_row + 1]) {
      cur_row++;
    }

    // skip samples with empty rows
    while ((samples[cur_samp] - 1) < cur_row) {
      cur_samp++;
    }

    if (cur_samp >= nsamp || samples[cur_samp] > *nrow) {
      break;
    }

    // here: samples[cur_samp] - 1 is at least cur_row
    // while cur_row is lesser than
    while ((samples[cur_samp] - 1) > cur_row) {
      while (it < (*row_ptr)[cur_row + 1]) {
        READ(file, "%d", &int_trash);
        it++;
      }
      cur_row++;
    }

    // here samples[cur_samp] - 1 == cur_row
    while (it < (*row_ptr)[cur_row + 1]) {
      READ(file, "%d", (col_ind + cur_val));
      if (max_col < col_ind[cur_val]) {
        max_col = col_ind[cur_val];
      }
      it++;
      cur_val++;
    }

    // go onto the next sample
    cur_samp++;
  }
  // read rest
  for (; it < *nnz; it++) {
    READ(file, "%d", &int_trash);
    if (max_col < int_trash) {
        max_col = int_trash;
    }
  }

  *ncol = max_col + 1;

  *col_sums = (double*) R_alloc(*ncol, sizeof(double));
  _setMemoryDouble(*col_sums, 0.0, *ncol);

  cur_row = 0;
  cur_samp = 0;
  it = 0;
  cur_val = 0;
  // read values and calculate column sums
  while (cur_samp < nsamp) {
    // skip empty rows
    while (it == (*row_ptr)[cur_row + 1]) {
      cur_row++;
    }

    // skip samples with empty rows
    while ((samples[cur_samp] - 1) < cur_row) {
      cur_samp++;
    }

    if (cur_samp >= nsamp || samples[cur_samp] > *nrow) {
      break;
    }

    // here: samples[cur_samp] - 1 is at least cur_row
    // while cur_row is lesser than
    while ((samples[cur_samp] - 1) > cur_row) {
      while (it < (*row_ptr)[cur_row + 1]) {
        READ(file, "%lf", &double_trash);
        it++;
      }
      cur_row++;
    }

    // here samples[cur_samp] - 1 == cur_row
    while (it < (*row_ptr)[cur_row + 1]) {
      READ(file, "%lf", (val + cur_val));

      (*col_sums)[col_ind[cur_val]] += val[cur_val];
      it++;
      cur_val++;
    }

    // go onto the next sample
    cur_samp++;
  }

  fclose(file);

  *row_ptr = new_row_ptr;
  *nnz = cur_val;
  *nrow = cur_samp;

  return 0;
}

int _readSamplesAndCalculateColSums(const char* file_name, int** row_ptr, int** col_ind, double** val, 
    int* nrow, int* nnz, const double lowerB, const double upperB, 
    const int* samples, const int nsamp, double** col_sums, int* ncol) {
  FILE* file = fopen(file_name, "r");

  if (file == NULL) {
    Rprintf("File %s not found!\n", file_name);
    return -1;
  }

  READ(file, "%d\n", nrow);
  CHECK(nrow > 0, file);

  READ(file, "%d\n", nnz);
  CHECK(nnz > 0, file);

  *row_ptr = (int*) R_alloc(*nrow + 1, sizeof(int));
  *col_ind = (int*) R_alloc(*nnz, sizeof(int));
  *val = (double*) R_alloc(*nnz, sizeof(double));

  if (samples[0] <= 0) {
    return _readAllSamplesAndCalculateColSums(file, *row_ptr, *col_ind, *val, *nrow, *nnz, lowerB, upperB, 
      col_sums, ncol);
  } else {
    return _readFilteredSamplesAndCalculateColSums(file, row_ptr, *col_ind, *val, nrow, nnz, lowerB,
      upperB, samples, nsamp, col_sums, ncol);
  }
}


/*********************************************
 * READ SAMPLES SPARSE RFN                   *
 ********************************************/

void _sparseToDenseMatrix(double* dense, const int* row_ptr, const int* col_ind, const double* val, 
    const int ncol, const int nnz, const double* col_sums, const double lowerB, 
    const double upperB) {
  int cur_row = 0;
  for (int i = 0; i < nnz; i++) {
    while (i >= row_ptr[cur_row + 1]) {
      cur_row++;
    }
    int col = col_ind[i];
    if (col_sums[col] > lowerB && col_sums[col] < upperB) {
      dense[cur_row * ncol + col] = val[i];
    }
   }
}

SEXP _filterColumnsAndCreateMatrix(const int* row_ptr, const int* col_ind, const double* val, 
    const int nnz, const int nrow, const int ncol, const double* col_sums, double lowerB, 
    double upperB) {
  SEXP XS = PROTECT(allocMatrix(REALSXP, ncol, nrow));
  double* X = REAL(XS);
  _setMemoryDouble(X, 0.0, nrow * ncol);

  _sparseToDenseMatrix(X, row_ptr, col_ind, val, ncol, nnz, col_sums, lowerB, upperB);

  UNPROTECT(1);
  return XS;
}

// sampleS are sorted integers. LENGTH(sampleS) <= nrow
SEXP readSamplesSpRfn(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS) {
  const char* file_name = CHAR(STRING_ELT(file_nameS, 0));
  
  const double lowerB = (double) (REAL(lowerBS)[0]);
  const double upperB = (double) (REAL(upperBS)[0]);
  const int* samples = INTEGER(samplesS);
  const int nsamp = LENGTH(samplesS);

  int* row_ptr;
  int* col_ind;
  double* val;
  int nrow;
  int nnz;

  double* col_sums;
  int ncol;

  if (_readSamplesAndCalculateColSums(file_name, &row_ptr, &col_ind, &val, &nrow, &nnz, lowerB, 
    upperB, samples, nsamp, &col_sums, &ncol) < 0) {
    return R_NilValue;
  }

  return _filterColumnsAndCreateMatrix(row_ptr, col_ind, val, nnz, nrow, ncol, col_sums, lowerB, 
      upperB);
}


/*********************************************
 * SAMPLES PER FEATURE                       *
 *********************************************/

SEXP samplesPerFeature(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS) {
  const char* file_name = CHAR(STRING_ELT(file_nameS, 0));

  const double lowerB = (double) (REAL(lowerBS)[0]);
  const double upperB = (double) (REAL(upperBS)[0]);
  const int* samples = INTEGER(samplesS);
  const int nsamp = LENGTH(samplesS);

  int* row_ptr;
  int* col_ind;
  double* val;
  int nrow;
  int nnz;

  double* col_sums;
  int ncol;

  if (_readSamplesAndCalculateColSums(file_name, &row_ptr, &col_ind, &val, &nrow, &nnz, lowerB, 
    upperB, samples, nsamp, &col_sums, &ncol) < 0) {
    return R_NilValue;
  }

  SEXP nsLS = PROTECT(allocVector(INTSXP, ncol));
  SEXP sLS = PROTECT(allocVector(VECSXP, ncol));

  int* nsL = INTEGER(nsLS);
  _setMemoryInt(nsL, 0, ncol);

  for (int i = 0; i < nnz; i++) {
    int col = col_ind[i];
    if (col_sums[col] > lowerB && col_sums[col] < upperB) {
      nsL[col]++;
    }
  }

  for (int i = 0; i < ncol; i++) {
    if (nsL[i] > 0) {
      SET_VECTOR_ELT(sLS, i, allocVector(INTSXP, nsL[i]));
    } else {
      SET_VECTOR_ELT(sLS, i, allocVector(INTSXP, 1));
      INTEGER(VECTOR_ELT(sLS, i))[0] = 0;
    }
  }

  int* indices = (int*) R_alloc(ncol, sizeof(int));
  _setMemoryInt(indices, 0, ncol);

  int cur_row = 0;
  for (int i = 0; i < nnz; i++) {
    while (i >= row_ptr[cur_row + 1]) {
      cur_row++;
    }
    int col = col_ind[i];
    if (nsL[col] > 0) {
      INTEGER(VECTOR_ELT(sLS, col))[indices[col]] = cur_row + 1;
      
      indices[col]++;
    }
  }

  SEXP namesS = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(namesS, 0, mkChar("sL"));
  SET_STRING_ELT(namesS, 1, mkChar("nsL"));

  SEXP outS = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(outS, 0, sLS);
  SET_VECTOR_ELT(outS, 1, nsLS);

  setAttrib(outS, R_NamesSymbol, namesS);
  UNPROTECT(4);

  return outS;
}


R_CallMethodDef callMethods[] = {
  {"readSamplesSpRfn", (DL_FUNC) &readSamplesSpRfn, 4},
  {"samplesPerFeature", (DL_FUNC) &readSamplesSpRfn, 4},
  {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

int main() {
  return 1;
}
