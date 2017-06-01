#include <R_ext/Rdynload.h>
#include <Rinternals.h>

#define GET_MACRO(_1,_2,_3,_4,NAME,...) NAME

#define _ERROR(file, str) fclose(file); Rprintf(str); return R_NilValue;

#define READ(...) GET_MACRO(__VA_ARGS__, READ4, READ3)(__VA_ARGS__)
#define READ4(from, template, variable, str) \
    if (fscanf(from, template, variable) < 1) { _ERROR(from, str) }
#define READ3(from, template, variable) \
    READ4(from, template, variable, "Wrong file format. \n")

#define CHECK(expression, file) if (!(expression)) { _ERROR(file, "Error: " #expression " is not true!") }


void _calculateColumnDelta(const double* col_sums, const int ncol, const double lowerB, 
    const double upperB, int* col_delta, int* delta) {
  *delta = 0;
  for (int i = 0; i < ncol; i++) {
    if (col_sums[i] > lowerB && col_sums[i] < upperB) {
      col_delta[i] = *delta;
    } else {
      (*delta)++;
      col_delta[i] = -1;
    }
  }
}

void _setMemory(double* array, const double value, const int n) {
  for (int i = 0; i < n; i++) {
    array[i] = 0.0;
  }
}

void _sparseToDenseMatrix(double* dense, const int* row_ptr, const int* col_ind,
    const int* col_delta, const double* val, const int nnz, const int nrow) {
  int cur_row = 0;
  for (int i = 0; i < nnz; i++) {
    while (i >= row_ptr[cur_row + 1]) {
      cur_row++;
    }
    if (col_delta[col_ind[i]] > -1) {
      dense[(col_ind[i] - col_delta[col_ind[i]]) * nrow + cur_row] = val[i];  
    }
  }
}

SEXP _filterColumnsAndCreateMatrix(const int* row_ptr, const int* col_ind, const double* val, 
    const int nnz, const int nrow, const int ncol,
    const double* col_sums, double lowerB, double upperB) {
  // filter by bounds
  int* col_delta = (int*) R_alloc(ncol, sizeof(int));
  int delta;
  _calculateColumnDelta(col_sums, ncol, lowerB, upperB, col_delta, &delta);
  int new_ncol = ncol - delta;

  SEXP XS = PROTECT(allocMatrix(REALSXP, nrow, new_ncol));
  double* X = REAL(XS);
  _setMemory(X, 0.0, nrow * new_ncol);

  _sparseToDenseMatrix(X, row_ptr, col_ind, col_delta, val, nnz, nrow);

  UNPROTECT(1);
  return XS;
}

/**
 * 
 * readSamplesSpRfnAll
 *
 *
 */
SEXP _readSamplesSpRfnAll(FILE* file, int* row_ptr, int* col_ind, double* val, 
    int nrow, int nnz, double lowerB, double upperB) {
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
  int ncol = max_col + 1;

  double* col_sums = (double*) R_alloc(ncol, sizeof(double));
  _setMemory(col_sums, 0.0, ncol);

  // read values and calculate column sums
  for (int i = 0; i < nnz; i++) {
    READ(file, "%lf", (val + i));
    col_sums[col_ind[i]] += val[i];
  }
  fclose(file);

  return _filterColumnsAndCreateMatrix(row_ptr, col_ind, val, nnz, nrow, ncol, col_sums, lowerB, 
      upperB);
}

/**
 *
 *
 *
 */
SEXP _readSamplesSpRfnFilter(FILE* file, const int* samples, int* row_ptr, int* col_ind, 
    double* val, const int nrow, const int nnz, const int nsamp, const double lowerB, 
    const double upperB) {
  int cur_samp = 0;
  
  int* new_row_ptr = (int*) R_alloc(nsamp + 1, sizeof(int));
  new_row_ptr[0] = 0;
  
  // read row pointers and build the filtered row pointers
  for (int i = 0; i < nrow + 1; i++) {
    READ(file, "%d", (row_ptr + i));
    if (cur_samp < nsamp && samples[cur_samp] == i) {
      new_row_ptr[cur_samp + 1] = new_row_ptr[cur_samp] + row_ptr[i] - row_ptr[i - 1];
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
    while (it == row_ptr[cur_row + 1]) {
      cur_row++;
    }

    // skip samples with empty rows
    while ((samples[cur_samp] - 1) < cur_row) {
      cur_samp++;
    }
    if (cur_samp >= nsamp || samples[cur_samp] > nrow) {
      break;
    }

    // here: samples[cur_samp] - 1 is at least cur_row
    // while cur_row is lesser than
    while ((samples[cur_samp] - 1) > cur_row) {
      while (it < row_ptr[cur_row + 1]) {
        READ(file, "%d", &int_trash);
        it++;
      }
      cur_row++;
    }

    // here samples[cur_samp] - 1 == cur_row
    while (it < row_ptr[cur_row + 1]) {
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
  int new_nrow = cur_samp;

  // read rest
  for (; it < nnz; it++) {
    READ(file, "%d", &int_trash);
  }

  int new_nnz = cur_val;
  int ncol = max_col + 1;

  double* col_sums = (double*) R_alloc(ncol, sizeof(double));
  _setMemory(col_sums, 0.0, ncol);

  cur_row = 0;
  cur_samp = 0;
  it = 0;
  cur_val = 0;
  // read values and calculate column sums
  while (cur_samp < nsamp) {
    // skip empty rows
    while (it == row_ptr[cur_row + 1]) {
      cur_row++;
    }

    // skip samples with empty rows
    while ((samples[cur_samp] - 1) < cur_row) {
      cur_samp++;
    }

    if (cur_samp >= nsamp || samples[cur_samp] > nrow) {
      break;
    }


    // here: samples[cur_samp] - 1 is at least cur_row
    // while cur_row is lesser than
    while ((samples[cur_samp] - 1) > cur_row) {
      while (it < row_ptr[cur_row + 1]) {
        READ(file, "%lf", &double_trash);
        it++;
      }
      cur_row++;
    }

    // here samples[cur_samp] - 1 == cur_row
    while (it < row_ptr[cur_row + 1]) {
      READ(file, "%lf", (val + cur_val));

      col_sums[col_ind[cur_val]] += val[cur_val];
      it++;
      cur_val++;
    }

    // go onto the next sample
    cur_samp++;
  }

  fclose(file);

  return _filterColumnsAndCreateMatrix(new_row_ptr, col_ind, val, new_nnz, new_nrow, ncol, col_sums, lowerB, 
      upperB);
}

// sampleS are sorted integers. LENGTH(sampleS) <= nrow
SEXP readSamplesSpRfn(SEXP file_nameS, SEXP samplesS, SEXP lowerBS, SEXP upperBS) {
  const char* file_name = CHAR(STRING_ELT(file_nameS, 0));
  
  const double lowerB = (double) (REAL(lowerBS)[0]);
  const double upperB = (double) (REAL(upperBS)[0]);
  const int* samples = INTEGER(samplesS);
  const int nsamp = LENGTH(samplesS);

  FILE* file = fopen(file_name, "r");

  if (file == NULL) {
    Rprintf("File %s not found!\n", file_name);
    return R_NilValue;
  }

  int nrow, nnz;

  READ(file, "%d\n", &nrow);
  CHECK(nrow > 0, file);

  READ(file, "%d\n", &nnz);
  CHECK(nnz > 0, file);

  int* row_ptr = (int*) R_alloc(nrow + 1, sizeof(int));
  int* col_ind = (int*) R_alloc(nnz, sizeof(int));
  double* val = (double*) R_alloc(nnz, sizeof(double));

  if (samples[0] <= 0) {
    return _readSamplesSpRfnAll(file, row_ptr, col_ind, val, nrow, nnz, lowerB, upperB);
  } else {
    return _readSamplesSpRfnFilter(file, samples, row_ptr, col_ind, val, nrow, nnz, nsamp, 
      lowerB, upperB);
  }
}



R_CallMethodDef callMethods[] = {
  {"readSamplesSpRfn", (DL_FUNC) &readSamplesSpRfn, 4},
  {NULL, NULL, 0}
};

void R_init_myLib(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

int main() {
  return 1;
}
