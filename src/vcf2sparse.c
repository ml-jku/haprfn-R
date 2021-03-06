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

#define kstring_free(str) { free(str->s); free(str); }

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

static const char VcfGzPostfix[] = ".vcf.gz";
static const char VcfPostfix[] = ".vcf";
static const size_t IgnoreInterval = -1;

#define PATH_SEPARATOR '/'

// file_name must be released
char* create_file_name(const char *file_name, const char *prefix, const char *postfix, const size_t lower_interval, const size_t upper_interval) {
  static char buffer[100];
  static char separator[2];  
  
  if (lower_interval == IgnoreInterval || upper_interval == IgnoreInterval) {
    buffer[0] = 0;
  } else {
    sprintf(buffer, "_%zd_%zd", lower_interval, upper_interval);  
  }

  size_t prefix_len = strlen(prefix);
  if (prefix_len > 0 && prefix[prefix_len - 1] != PATH_SEPARATOR) {
    separator[0] = PATH_SEPARATOR;
    separator[1] = 0;
  } else {
    separator[0] = 0;
  }

  char *fn = (char*) calloc(prefix_len + strlen(separator) + strlen(file_name) + strlen(buffer) + strlen(postfix) + 1, sizeof(char));
  strcat(fn, prefix);
  strcat(fn, separator);
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

void write_to_annotation_file(const char *file_name, const char *prefix, const char *postix, const char *string, const size_t lower_interval, const size_t upper_interval) {
  FILE *file = open_file(file_name, prefix, postix, lower_interval, upper_interval);
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

sparse_matrix_t* dense_to_sparse(unsigned short **de_matrix, const unsigned int *nnz, const size_t nrow, const size_t ncol) {
  unsigned int nonzero = 0;
  for (size_t i = 0; i < nrow; i++) {
    nonzero += nnz[i];
  }
  sparse_matrix_t* sp_matrix = matrix_init(ncol, nonzero);
  size_t cur_val = 0;

  // i and j switched to transpose the matrix
  for (size_t j = 0; j < ncol; j++) {
    sp_matrix->rows_pointer[j + 1] = sp_matrix->rows_pointer[j];
    for (size_t i = 0; i < nrow; i++) {
      if (de_matrix[i][j] != 0) {
        sp_matrix->value[cur_val] = de_matrix[i][j];
        sp_matrix->column[cur_val] = i;
        sp_matrix->rows_pointer[j + 1]++;
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

// This method changes the variable matrix and nnz
void haplotypes_to_genotypes(unsigned short **matrix, unsigned int *nnz, const size_t nrow, const size_t ncol) {
  unsigned int *h = (unsigned int*) R_alloc(MAX_PLOIDY, sizeof(unsigned int));

  for (size_t i = 0; i < nrow; i++) {
    for (size_t j = 0; j < ncol / MAX_PLOIDY; j++) {
      for (size_t p = 0; p < MAX_PLOIDY; p++) {
        h[p] = matrix[i][j * MAX_PLOIDY + p];
      }
      // general solution makes it complicated
      // only check the first two

      if (h[0] == 0 && h[1] == 0) {
        matrix[i][j] = 0;
      } else if (h[0] == 1 && h[1] == 1) {
        matrix[i][j] = 2;
        nnz[i] -= 1;
      } else {
        matrix[i][j] = 1;
      }
    }
  }
}

void write_dense_matrices_as_sparse(unsigned short **matrix, const unsigned int *nnz, const size_t nrow, const size_t ncol, 
    const char *file_name, const char *prefix, const char *postfix, const size_t lower_interval, const size_t upper_interval) {
  FILE *file = open_file(file_name, prefix, postfix, lower_interval, upper_interval);
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

float bcf_allele_frequency(bcf_hdr_t *hdr, bcf1_t *bcf, float **buffer, int *buffer_size) {
  bcf_unpack(bcf, BCF_UN_INFO);
  int retval = bcf_get_info_float(hdr, bcf, "AF", buffer, buffer_size);
  
  return retval >= 0 && *buffer_size > 0 ? *buffer[0] : -1.0;
}

void print_bcf_error(bcf1_t *bcf, const char *error_message, const char *solution) {
  REprintf("Error (SNP ID: %s)!\n", bcf_snp_id(bcf));
  REprintf("\t%s\n", error_message);
  if (solution != NULL) {
    REprintf("\tSolution: %s\n", solution);
  }
  REprintf("\tAborting.\n");
}

void vcf2sparse(SEXP file_nameS, SEXP prefix_pathS, SEXP interval_sizeS, SEXP shift_sizeS,
    SEXP annotateS, SEXP genotypesS, SEXP haplotypesS, SEXP missing_valuesS,
    SEXP annotation_postfixS, SEXP genotypes_postfixS, SEXP haplotypes_postfixS, SEXP info_postfixS,
    SEXP individuals_postfixS, SEXP output_fileS, SEXP output_prefixS) {
  const char *file_name = CHAR(STRING_ELT(file_nameS, 0));
  const char *prefix_path = isNull(prefix_pathS) ? "" : CHAR(STRING_ELT(prefix_pathS, 0));
  const char *output_file = isNull(output_fileS) ? file_name : CHAR(STRING_ELT(output_fileS, 0));
  const char *output_prefix = isNull(output_prefixS) ? prefix_path : CHAR(STRING_ELT(output_prefixS, 0));
  const Rboolean annotate = LOGICAL(annotateS)[0];
  const Rboolean genotypes = LOGICAL(genotypesS)[0];
  const Rboolean haplotypes = LOGICAL(haplotypesS)[0];
  const size_t interval_size = INTEGER(interval_sizeS)[0];
  const size_t shift_size = INTEGER(shift_sizeS)[0];
  const missing_t missing_values = (missing_t) INTEGER(missing_valuesS)[0];

  const char *annotation_postfix = isNull(annotation_postfixS) ? "" : CHAR(STRING_ELT(annotation_postfixS, 0));
  const char *genotypes_postfix = isNull(genotypes_postfixS) ? "" : CHAR(STRING_ELT(genotypes_postfixS, 0));
  const char *haplotypes_postfix = isNull(haplotypes_postfixS) ? "" : CHAR(STRING_ELT(haplotypes_postfixS, 0));
  const char *info_postfix = isNull(info_postfixS) ? "" : CHAR(STRING_ELT(info_postfixS, 0));
  const char *individuals_postfix = isNull(individuals_postfixS) ? "" : CHAR(STRING_ELT(individuals_postfixS, 0));

  htsFile *file = NULL;
  bcf_hdr_t *hdr = NULL;
  bcf1_t *bcf = NULL;
  FILE *individuals_file = NULL;
  char *vcfgz_file_name = NULL;
  char *vcf_file_name = NULL;

  kstring_t* buffer = NULL;
  kstring_t* freq_flip_col = NULL;
  kstring_t* current_buffer = NULL;
  kstring_t* next_buffer = NULL;
  kstring_t* tmps;

  int32_t *gt_arr = NULL, ngt_arr = 0;
  int float_buffer_size = 1;
  float *float_buffer = NULL;

  vcfgz_file_name = create_file_name(file_name, prefix_path, VcfGzPostfix, IgnoreInterval, IgnoreInterval);
  if (!(file = bcf_open(vcfgz_file_name, "r"))) {

    vcf_file_name = create_file_name(file_name, prefix_path, VcfPostfix, IgnoreInterval, IgnoreInterval);
    if (!(file = bcf_open(vcf_file_name, "r"))) {

      REprintf("Cannot open any of: \n\t%s\n\t%s\n", vcfgz_file_name, vcf_file_name);
      goto cleanup;
    }
  }

  // Read header file
  if (!(hdr = bcf_hdr_read(file))) {
    REprintf("Could not read VCF header!\n");
    goto cleanup;
  }

  // Print _individuals.txt
  individuals_file = open_file(output_file, output_prefix, individuals_postfix, IgnoreInterval, IgnoreInterval);
  if (!individuals_file) {
    REprintf("Could not open individuals file for writing!\n");
    goto cleanup;
  }

  size_t nsamp = bcf_hdr_nsamples(hdr);
  for (size_t i = 0; i < nsamp; i++) {
    fprintf(individuals_file, "%zu %s\n", i + 1, hdr->samples[i]);
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

  buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
  freq_flip_col = (kstring_t*) calloc(1, sizeof(kstring_t));
  current_buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
  next_buffer = (kstring_t*) calloc(1, sizeof(kstring_t));

  float_buffer = calloc(float_buffer_size, sizeof(float));

  size_t current_interval = 0;
  size_t next_interval = 0;
  size_t n_interval = 0;

  size_t nsnp = 0;
  while (bcf_read(file, hdr, bcf) >= 0) {
    nsnp++;

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
                frequency = bcf_allele_frequency(hdr, bcf, &float_buffer, &float_buffer_size);
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
          } else if (ptr[j] == bcf_int32_vector_end) { // unexpected end of genotype field
            if (j > 0 && missing_ind > 0 && missing[missing_ind - 1] == i * MAX_PLOIDY + (j - 1)) {
              // if the previous is also missing then set this to missing
              missing[missing_ind++] = i * MAX_PLOIDY + j;
            } else if (j > 0) {
              // if the previous is not missing, use that value for this
              allele_index = current_matrix[current_interval][i * MAX_PLOIDY + (j - 1)];
            } else {
              print_bcf_error(bcf, "No genotype found.\n", NULL);
              goto cleanup;
            }
          } else {
            allele_index = bcf_gt_allele(ptr[j]);
          }

          if (missing_ind > 0 && missing[missing_ind - 1] == i * MAX_PLOIDY + j) {
            // if the current one is missing
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

        float maf = ((float) current_nnz[current_interval]) / (nsamp * MAX_PLOIDY - missing_ind);
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
      // delete newline character
      buffer->l = buffer->l - 1;
      float frequency = ((float) current_nnz[current_interval]) / (nsamp * MAX_PLOIDY);
      ksprintf(freq_flip_col, "\t%.8f\t%d\n", frequency, frequency > 0.5);
      kputsn(freq_flip_col->s, freq_flip_col->l, buffer);

      kputsn(buffer->s, buffer->l, current_buffer);

      if (current_interval >= shift_size) {
        kputsn(buffer->s, buffer->l, next_buffer);
      }

      buffer->s[0] = 0;
      buffer->l = 0;

      freq_flip_col->s[0] = 0;
      freq_flip_col->l = 0;
    } // annotate

    if (current_interval >= shift_size) {
      next_interval++;
    }

    current_interval++;

    // start a new interval
    if (current_interval % interval_size == 0) {
      flip_matrix(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY);
      size_t lower_interval = n_interval * shift_size;
      if (haplotypes) {
        write_dense_matrices_as_sparse(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY, output_file, output_prefix, haplotypes_postfix, lower_interval, lower_interval + interval_size);
      }
      if (genotypes) {
        haplotypes_to_genotypes(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY);
        write_dense_matrices_as_sparse(current_matrix, current_nnz, interval_size, nsamp, output_file, output_prefix, genotypes_postfix, lower_interval, lower_interval + interval_size);
      }

      tmpm = current_matrix;
      current_matrix = next_matrix;
      next_matrix = tmpm;
      set_zerom(next_matrix, interval_size, nsamp * MAX_PLOIDY);

      tmpv = current_nnz;
      current_nnz = next_nnz;
      next_nnz = tmpv;
      set_zerov(next_nnz, interval_size);

      if (annotate) {
        write_to_annotation_file(output_file, output_prefix, annotation_postfix, current_buffer->s, lower_interval, lower_interval + interval_size);

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

    // log message
    if (nsnp % 1000 == 0) {
      Rprintf("Read SNV: %zu\r", nsnp);
    }
  }

  if (current_interval > 0) {
      flip_matrix(current_matrix, current_nnz, interval_size, nsamp * MAX_PLOIDY);

      size_t lower_interval = n_interval * shift_size;
      if (haplotypes) {
        write_dense_matrices_as_sparse(current_matrix, current_nnz, current_interval, nsamp * MAX_PLOIDY, output_file, output_prefix, haplotypes_postfix, lower_interval, lower_interval + current_interval);
      }
      if (genotypes) {
        haplotypes_to_genotypes(current_matrix, current_nnz, current_interval, nsamp * MAX_PLOIDY);
        write_dense_matrices_as_sparse(current_matrix, current_nnz, current_interval, nsamp, output_file, output_prefix, genotypes_postfix, lower_interval, lower_interval + current_interval);
      }

      if (annotate) {
        write_to_annotation_file(output_file, output_prefix, annotation_postfix, current_buffer->s, lower_interval, lower_interval + current_interval);  
      }
  }

  FILE *info = open_file(output_file, output_prefix, info_postfix, IgnoreInterval, IgnoreInterval);
  if (!info) {
    REprintf("Cannot write to info file!\n");
  } else {
    fprintf(info, "%zu\n%zu", nsamp, nsnp);
    fclose(info);
  }


cleanup:
  if (file) hts_close(file);
  if (buffer) kstring_free(buffer);
  if (freq_flip_col) kstring_free(freq_flip_col);
  if (current_buffer) kstring_free(current_buffer);
  if (next_buffer) kstring_free(next_buffer);
  if (float_buffer) free(float_buffer);
  if (vcfgz_file_name) free(vcfgz_file_name);
  if (vcf_file_name) free(vcf_file_name);
  if (gt_arr) free(gt_arr);
  if (bcf) bcf_destroy(bcf);
  if (hdr) bcf_hdr_destroy(hdr);
}
