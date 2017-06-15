#include <R_ext/Boolean.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "interface.h"
#include "vcf.h"
#include "kstring.h"

#define MAX_PLOIDY 2
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef enum {
  kDefault = 0,
  kMajor = 1,
  kMinor = 2,
  kCalculate = 3,
  kInfo = 4,
  kAbort = 5
} missing_t;
/****
 * Sparse matrix structure
 ****/

typedef struct {
  size_t nrow;
  size_t nnz;
  unsigned short *value;
  size_t *column;
  size_t *rows_pointer;
} sparse_matrix_t;

sparse_matrix_t* matrix_init(const size_t nrow, const size_t nnz) {
  sparse_matrix_t *matrix = (sparse_matrix_t*) calloc(1, sizeof(sparse_matrix_t));
  matrix->nrow = nrow;
  matrix->nnz = nnz;
  matrix->value = (unsigned short*) calloc(nnz, sizeof(unsigned short));
  matrix->column = (size_t*) calloc(nnz, sizeof(size_t));
  matrix->rows_pointer = (size_t*) calloc(nrow + 1, sizeof(size_t));

  return matrix;
}

void matrix_destroy(sparse_matrix_t* matrix) {
  free(matrix->value);
  free(matrix->column);
  free(matrix->rows_pointer);
  free(matrix);
}

/****
 * File operations
 ****/


static const char AnnotationPostfix[] = "_annot.txt";
static const char IndividualsPostfix[] = "_individuals.txt";
static const char InfoPostfix[] = "_info.txt";
static const char MatrixPostfix[] = "_mat.txt";
static const char VcfGzPostfix[] = ".vcf.gz";
static const char VcfPostfix[] = ".vcf";
static const size_t IgnoreInterval = -1;


// file_name must be released
char* create_file_name(const char *file_name, const char *prefix, const char *postfix, const size_t lower_interval, const size_t upper_interval) {
  static char buffer[100];
  
  if (lower_interval == IgnoreInterval || upper_interval == IgnoreInterval) {
    buffer[0] = 0;
  } else {
    sprintf(buffer, "_%zd_%zd", lower_interval, upper_interval);  
  }

  char *fn = (char*) calloc(strlen(prefix) + strlen(file_name) + strlen(buffer) + strlen(postfix) + 1, sizeof(char));
  strcat(fn, prefix);
  strcat(fn, file_name);
  strcat(fn, buffer);
  strcat(fn, postfix);

  return fn;
}

FILE* open_file(const char *file_name, const char *prefix, const char *postfix, const size_t lower_interval, const size_t upper_interval) {
  char *fn = create_file_name(file_name, prefix, postfix, lower_interval, upper_interval);
  FILE *file = fopen(fn, "w");
  if (!file) {
    REprintf("Cannot open file %s\n", fn);
  }
  free(fn);
  return file;
}

void write_to_annotation_file(const char *file_name, const char *prefix, const char *string, const size_t lower_interval, const size_t upper_interval) {
  FILE *file = open_file(file_name, prefix, AnnotationPostfix, lower_interval, upper_interval);
  if (!file) {
    REprintf("Cannot write annotation file for interval %zd-%zd\n", lower_interval, upper_interval);
    return;
  }
  fputs(string, file);
  fclose(file);
}

unsigned short** create_matrix(const size_t nrow, const size_t ncol) {
  unsigned short** matrix = (unsigned short**) R_alloc(nrow, sizeof(unsigned short*));
  for (size_t i = 0; i < nrow; i++) {
    matrix[i] = (unsigned short*) R_alloc(ncol, sizeof(unsigned short));
  }
  return matrix;
}

void set_zerom(unsigned short **matrix, const size_t nrow, const size_t ncol) {
  for (size_t i = 0; i < nrow; i++) {
    memset(matrix[i], 0, ncol * sizeof(unsigned short));
  }
}

void set_zerov(unsigned int *vector, const size_t nrow) {
  memset(vector, 0, nrow * sizeof(unsigned int));
}

sparse_matrix_t* dense_to_sparse(unsigned short **de_matrix, unsigned int *nnz, const size_t nrow, const size_t ncol) {
  unsigned int nonzero = 0;
  for (size_t i = 0; i < nrow; i++) {
    nonzero += nnz[i];
  }
  sparse_matrix_t* sp_matrix = matrix_init(nrow, nonzero);
  size_t cur_val = 0;
  for (size_t i = 0; i < nrow; i++) {
    sp_matrix->rows_pointer[i + 1] = sp_matrix->rows_pointer[i];
    for (size_t j = 0; j < ncol; j++) {
      if (de_matrix[i][j] != 0) {
        sp_matrix->value[cur_val] = de_matrix[i][j];
        sp_matrix->column[cur_val] = j;
        sp_matrix->rows_pointer[i + 1]++;
        cur_val++;
      }
    }
  }
  return sp_matrix;
}

void write_sparse_matrix(sparse_matrix_t *matrix, FILE *file) {
  fprintf(file, "%zd\n", matrix->nrow);
  fprintf(file, "%zd\n", matrix->nnz);
  for (size_t i = 0; i < matrix->nrow + 1; i++) {
    fprintf(file, "%zd ", matrix->rows_pointer[i]);
  }
  fprintf(file, "\n");
  for (size_t i = 0; i < matrix->nnz; i++) {
    fprintf(file, "%zd ", matrix->column[i]);
  }
  fprintf(file, "\n");
  for (size_t i = 0; i < matrix->nnz; i++) {
    fprintf(file, "%u ", matrix->value[i]);
  }
}

void write_dense_matrices_as_sparse(unsigned short **matrix, unsigned int *nnz, const size_t nrow, const size_t ncol, 
    const char *file_name, const char *prefix, const size_t lower_interval, const size_t upper_interval) {
  FILE *file = open_file(file_name, prefix, MatrixPostfix, lower_interval, upper_interval);
  if (!file) {
    REprintf("Cannot write sparse matrix to file for interval %zd-%zd\n", lower_interval, upper_interval);
    return;
  }
  
  sparse_matrix_t *sp_matrix = dense_to_sparse(matrix, nnz, nrow, ncol);
  write_sparse_matrix(sp_matrix, file);

  matrix_destroy(sp_matrix);
  fclose(file);
}

void flip_matrix(unsigned short **matrix, unsigned int *nnz, const size_t nrow, const size_t ncol) {
  for (size_t i = 0; i < nrow; i++) {
    if (ncol / 2 < nnz[i]) {
      nnz[i] = ncol - nnz[i];
      for (size_t j = 0; j < ncol; j++) {
        matrix[i][j] = 1 - matrix[i][j];
      }
    }
  }
}

char* bcf_snp_id(bcf1_t *bcf) {
  bcf_unpack(bcf, BCF_UN_STR);
  return bcf->d.id;
}

float bcf_allele_frequency(bcf_hdr_t *hdr, bcf1_t *bcf) {
  float *value;
  int n;
  bcf_unpack(bcf, BCF_UN_INFO);
  int retval = bcf_get_info_float(hdr, bcf, "AF", &value, &n);
  
  return retval >= 0 && n > 0 ? value[0] : -1.0;
}

void print_bcf_error(bcf1_t *bcf, const char *error_message, const char *solution) {
  REprintf("Error (SNP ID: %s)!\n", bcf_snp_id(bcf));
  REprintf("\t%s\n", error_message);
  if (solution != NULL) {
    REprintf("\tSolution: %s\n", solution);
  }
  REprintf("\tAborting.\n");
}

void vcf2sparse(SEXP file_nameS, SEXP prefix_pathS, SEXP interval_sizeS, SEXP shift_sizeS, SEXP annotateS, SEXP genotypesS, SEXP haplotypesS, SEXP missing_valuesS, SEXP output_fileS) {
  const char *file_name = CHAR(STRING_ELT(file_nameS, 0));
  const char *prefix_path = isNull(prefix_pathS) ? "" : CHAR(STRING_ELT(prefix_pathS, 0));
  const char *output_file = isNull(output_fileS) ? file_name : CHAR(STRING_ELT(output_fileS, 0));
  const Rboolean annotate = LOGICAL(annotateS)[0];
  const Rboolean genotypes = LOGICAL(genotypesS)[0];
  const Rboolean haplotypes = LOGICAL(haplotypesS)[0];
  const size_t interval_size = INTEGER(interval_sizeS)[0];
  const size_t shift_size = INTEGER(shift_sizeS)[0];
  const missing_t missing_values = (missing_t) INTEGER(missing_valuesS)[0];

  htsFile *file = NULL;
  bcf_hdr_t *hdr = NULL;
  bcf1_t *bcf = NULL;
  FILE *individuals_file = NULL;
  char *vcfgz_file_name = NULL;
  char *vcf_file_name = NULL;


  vcfgz_file_name = create_file_name(file_name, prefix_path, VcfGzPostfix, IgnoreInterval, IgnoreInterval);
  if (!(file = bcf_open(vcfgz_file_name, "r"))) {

    vcf_file_name = create_file_name(file_name, prefix_path, VcfPostfix, IgnoreInterval, IgnoreInterval);
    if (!(file = bcf_open(vcf_file_name, "r"))) {

      REprintf("Cannot open any: \n\t%s\n\t%s\n", vcfgz_file_name, vcf_file_name);
      goto cleanup;
    }
  }

  // Read header file
  if (!(hdr = bcf_hdr_read(file))) {
    REprintf("Could not read VCF header!\n");
    goto cleanup;
  }

  // Print _individuals.txt
  individuals_file = open_file(output_file, prefix_path, IndividualsPostfix, IgnoreInterval, IgnoreInterval);
  if (!individuals_file) {
    REprintf("Could not open individuals file for writing!\n");
    goto cleanup;
  }

  size_t nsamp = bcf_hdr_nsamples(hdr);
  for (int i = 0; i < nsamp; i++) {
    fprintf(individuals_file, "%d %s\n", i + 1, hdr->samples[i]);
  }
  fclose(individuals_file);
	
  Rprintf("Individuals: %zu\n", nsamp);

  unsigned short **current_matrix = create_matrix(interval_size, nsamp * MAX_PLOIDY);
  unsigned short **next_matrix = create_matrix(interval_size, nsamp * MAX_PLOIDY);
  unsigned short **tmpm;
  set_zerom(current_matrix, interval_size, nsamp * MAX_PLOIDY);
  set_zerom(next_matrix, interval_size, nsamp * MAX_PLOIDY);

  unsigned int *current_nnz = (unsigned int*) R_alloc(interval_size, sizeof(unsigned int));
  unsigned int *next_nnz = (unsigned int*) R_alloc(interval_size, sizeof(unsigned int));
  unsigned int *tmpv;
  set_zerov(current_nnz, interval_size);
  set_zerov(next_nnz, interval_size);

  size_t *missing = (size_t*) R_alloc(nsamp * MAX_PLOIDY, sizeof(size_t));
  memset(missing, 0, nsamp * MAX_PLOIDY * sizeof(size_t));
  size_t missing_ind = 0;

  bcf = bcf_init();

  kstring_t* buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
  kstring_t* current_buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
  kstring_t* next_buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
  kstring_t* tmps;

  size_t current_interval = 0;
  size_t next_interval = 0;
  size_t n_interval = 0;

  size_t nsnp = 0;
  while (bcf_read(file, hdr, bcf) >= 0) {
    nsnp++;

    int32_t *gt_arr = NULL, ngt_arr = 0;
    int ngt = bcf_get_genotypes(hdr, bcf, &gt_arr, &ngt_arr);

    if (ngt > 0) {
      int max_ploidy = ngt / nsamp;

      if (max_ploidy > MAX_PLOIDY) {
        print_bcf_error(bcf, "Cannot handle poliploidy.", NULL);
        goto cleanup;
      }

      for (size_t i = 0; i < nsamp; i++) {
        int32_t *ptr = gt_arr + i * max_ploidy;

        for (size_t j = 0; j < MAX_PLOIDY; j++) {
          int allele_index;
          float frequency;

          if (bcf_gt_is_missing(ptr[j])) {
            switch (missing_values) {
              case kAbort:
                print_bcf_error(bcf, "Missing genotype found.",
                  "Try using the default value for the 'missingValues' parameter.");
                goto cleanup;
              case kMajor:
                allele_index = 0;
                break;
              case kMinor:
                allele_index = 1;
                break;
              case kDefault:
              case kInfo:
                frequency = bcf_allele_frequency(hdr, bcf);
                if (missing_values == kInfo) {
                  if (frequency < 0) {
                    print_bcf_error(bcf, "Missing AF info field.",
                      "Try using the default value for the 'missingValues' parameter.");
                    goto cleanup;
                  } else if (frequency > 1) {
                    print_bcf_error(bcf, "AF info field contains frequency greater than one.", NULL);
                    goto cleanup;
                  }
                } 
                if (frequency >= 0 && frequency <= 1) {
                  allele_index = frequency > 0.5 ? 1 : 0;
                  break;
                }
                // fallthrough for kDefault
              case kCalculate:
                missing[missing_ind++] = i * MAX_PLOIDY + j;
                // calculate
                break;
            }
          } else if (ptr[j] == bcf_int32_vector_end) {
            if (j > 0 && missing_ind > 0 && missing[missing_ind - 1] == i * MAX_PLOIDY + (j - 1)) {
              missing[missing_ind++] = i * MAX_PLOIDY + j;
            } else if (j > 0) {
              allele_index = current_matrix[current_interval][i * MAX_PLOIDY + (j - 1)];
            } else {
              print_bcf_error(bcf, "No genotype found.\n", NULL);
              goto cleanup;
            }
          } else {
            allele_index = bcf_gt_allele(ptr[j]);
          }

          if (missing_ind > 0 && missing[missing_ind - 1] == i * MAX_PLOIDY + j) {
            continue;
          }

          if (allele_index > 1) {
            print_bcf_error(bcf, "Multiallelic SNP found.", "Please split these sites into mutliple rows.");
            goto cleanup;
          }

          if (allele_index < 0) {
            print_bcf_error(bcf, "Negative allele index found.", NULL);
            goto cleanup;
          }

          if (allele_index == 1) {
            size_t sample_haplo_index = i * MAX_PLOIDY + j;
            current_matrix[current_interval][sample_haplo_index] = allele_index;
            current_nnz[current_interval]++;

            if (current_interval >= shift_size) {
              next_matrix[next_interval][sample_haplo_index] = allele_index;
              next_nnz[next_interval]++;
            }
          }
        } // for all alleles
      } // for all samples

      if (missing_ind > 0) {
        // estimate major allele
        // current_nnz is # of minor alleles
        float maf = current_nnz[current_interval] / (nsamp * MAX_PLOIDY - missing_ind);
        if (maf > 0.5) {
          for (size_t i = 0; i < missing_ind; i++) {
            current_matrix[current_interval][missing[i]] = 1;
            if (current_interval >= shift_size) {
              next_matrix[next_interval][missing[i]] = 1;
            }
          }
          current_nnz[current_interval] += missing_ind;
          if (current_interval >= shift_size) {
            next_nnz[next_interval] += missing_ind;
          }
        }
      }

      missing_ind = 0;
    }

    if (annotate) {
      // annotate add columns
      vcf_format(hdr, bcf, buffer, 1);
      kputsn(buffer->s, buffer->l, current_buffer);

      if (current_interval >= shift_size) {
        kputsn(buffer->s, buffer->l, next_buffer);
      }

      buffer->s[0] = 0;
      buffer->l = 0;
    } // annotate

    if (current_interval >= shift_size) {
      next_interval++;
    }

    current_interval++;

    // start a new interval
    if (current_interval % interval_size == 0) {
      flip_matrix(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY);
      size_t lower_interval = n_interval * shift_size;
      write_dense_matrices_as_sparse(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY, output_file, prefix_path, lower_interval, lower_interval + interval_size);

      tmpm = current_matrix;
      current_matrix = next_matrix;
      next_matrix = tmpm;
      set_zerom(next_matrix, interval_size, nsamp * MAX_PLOIDY);

      tmpv = current_nnz;
      current_nnz = next_nnz;
      next_nnz = tmpv;
      set_zerov(next_nnz, interval_size);

      if (annotate) {
        write_to_annotation_file(output_file, prefix_path, current_buffer->s, lower_interval, lower_interval + interval_size);

        tmps = current_buffer;
        current_buffer = next_buffer;
        next_buffer = tmps;

        next_buffer->s[0] = 0;
        next_buffer->l = 0;
      }

      current_interval = next_interval;
      next_interval = 0;

      n_interval++;
    }

    if (nsnp % 1000 == 0) {
      Rprintf("Read SNV: %zu\r", nsnp);
    }
  }

  if (current_interval > 0) {
      flip_matrix(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY);

      size_t lower_interval = n_interval * shift_size;
      write_dense_matrices_as_sparse(current_matrix, current_nnz, current_interval, nsamp * MAX_PLOIDY, output_file, prefix_path, lower_interval, lower_interval + current_interval);

      if (annotate) {
        write_to_annotation_file(output_file, prefix_path, current_buffer->s, lower_interval, lower_interval + current_interval);  
      }
  }

  FILE *info = open_file(output_file, prefix_path, InfoPostfix, IgnoreInterval, IgnoreInterval);
  if (!info) {
    REprintf("Cannot write to info file!\n");
  } else {
    fprintf(info, "%zu\n%zu", nsamp, nsnp);
    fclose(info);
  }


cleanup:
  if (vcfgz_file_name) free(vcfgz_file_name);
  if (vcf_file_name) free(vcf_file_name);
  if (bcf) bcf_destroy(bcf);
  if (hdr) bcf_hdr_destroy(hdr);
}
