/* sim_threaded_rng.c
   Reworked:
   - thread-local RNG (splitmix64 + xorshift128+)
   - thread-local normal generator (no static shared state)
   - preallocate arrays outside inner generation loop (no malloc/free inside loop)
   - keep existing functionality otherwise
   Compile:
     gcc -O2 -fopenmp -o sim_threaded_rng sim_threaded_rng.c -lm
   Run:
     OMP_NUM_THREADS=8 ./sim_threaded_rng | tee sim.log
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <stdint.h>
#include <omp.h>

/* Parameters (unchanged) */
const int num_generations = 50000;
const int num_genes = 100;
const int pop_size = 1000;
const int chr_pairs = 4;
const int chr_pair_counts[4] = {30,30,20,20}; /* sum=100 */
const double chr_length = 100000.0;
const double recomb_rate = 1e-8;
const double recomb_chr = recomb_rate * chr_length; /* 0.001 */

const double min_expr = 100000.0;
const double powerlaw_param = 2.0;

const int gene_rows = num_genes * 2;
const int total_rows = pop_size * gene_rows;    /* 200000 */
const int ng = pop_size * num_genes;            /* 100000 */

const double Lnt = 1000.0;
const double Lp = 300.0;
const double e_epistasis = 100.0;
const double min_correct_frac = 0.5;
const double shape_est_TR = 5.86, scale_est_TR = 0.173;
const double shape_est_TL = 5.44, scale_est_TL = 0.197;

/* Parameter combinations */
const double mut_rates[] = {5e-4, 6e-4, 7e-4, 8e-4, 9e-4};
const double param_cs[] = {1e-5, 1e-6, 1e-7, 1e-8};
const int num_mut_rates = 5;
const int num_param_cs = 4;
const int total_combinations = num_mut_rates * num_param_cs;

/* === Thread-local RNG: splitmix64 seed -> xorshift128+ === */

typedef struct {
    uint64_t s0;
    uint64_t s1;
    /* for normal generator spare */
    int has_spare;
    double spare;
} RNG;

/* splitmix64 for seeding */
static inline uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* Initialize RNG with a 64-bit seed; fill s0,s1 */
static void rng_init(RNG *rng, uint64_t seed) {
    uint64_t x = seed;
    rng->s0 = splitmix64(&x);
    rng->s1 = splitmix64(&x);
    /* ensure not both zero */
    if (rng->s0 == 0 && rng->s1 == 0) {
        rng->s0 = 0x9e3779b97f4a7c15ULL;
        rng->s1 = 0x6a09e667f3bcc909ULL;
    }
    rng->has_spare = 0;
    rng->spare = 0.0;
}

/* xorshift128+ next uint64 */
static inline uint64_t rng_next_u64(RNG *rng) {
    uint64_t s1 = rng->s0;
    const uint64_t s0 = rng->s1;
    rng->s0 = s0;
    s1 ^= s1 << 23; // a
    rng->s1 = (s1 ^ s0 ^ (s1 >> 17) ^ (s0 >> 26));
    return rng->s1 + s0;
}

/* uniform [0,1) double */
static inline double rng_uniform(RNG *rng) {
    /* take top 53 bits to produce double in [0,1) */
    uint64_t x = rng_next_u64(rng);
    /* shift to 53 bits */
    const uint64_t two53 = (1ULL << 53);
    double d = (double)(x >> 11) / (double)two53;
    return d;
}

/* normal N(0,1) via Box-Muller, using RNG->has_spare */
static double rng_normal(RNG *rng) {
    if (rng->has_spare) {
        rng->has_spare = 0;
        return rng->spare;
    } else {
        double u, v, s;
        do {
            u = 2.0 * rng_uniform(rng) - 1.0;
            v = 2.0 * rng_uniform(rng) - 1.0;
            s = u*u + v*v;
        } while (s >= 1.0 || s == 0.0);
        double mul = sqrt(-2.0 * log(s) / s);
        rng->spare = v * mul;
        rng->has_spare = 1;
        return u * mul;
    }
}

/* power-law rpar_one(xm,a) = xm / v^(1/a) */
static double rpar_one(RNG *rng, double xm, double a) {
    double v = rng_uniform(rng);
    if (v < 1e-16) v = 1e-16;
    return xm / pow(v, 1.0 / a);
}

/* Poisson sample (Knuth), using uniform from rng */
static int rpois_one(RNG *rng, double lambda) {
    if (lambda <= 0.0) return 0;
    if (lambda < 30.0) { /* Knuth for small lambda */
        double L = exp(-lambda);
        int k = 0;
        double p = 1.0;
        do {
            k++;
            p *= rng_uniform(rng);
        } while (p > L);
        return k - 1;
    } else { /* for larger lambda, use normal approx (fast) */
        double k = floor(lambda + sqrt(lambda) * rng_normal(rng));
        if (k < 0) return 0;
        return (int)k;
    }
}

/* runif_ab */
static inline double runif_ab(RNG *rng, double a, double b) {
    return a + (b - a) * rng_uniform(rng);
}

/* === Gamma sampling: Marsaglia and Tsang, using rng_normal and rng_uniform === */
static double rgamma_one(RNG *rng, double shape, double scale) {
    if (shape <= 0.0) return 0.0;
    if (shape < 1.0) {
        double u = rng_uniform(rng);
        return rgamma_one(rng, shape + 1.0, scale) * pow(u, 1.0 / shape);
    }
    double d = shape - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    while (1) {
        double x = rng_normal(rng);
        double v = 1.0 + c * x;
        if (v <= 0) continue;
        v = v * v * v;
        double u = rng_uniform(rng);
        if (u < 1.0 - 0.0331 * (x * x) * (x * x)) return scale * d * v;
        if (log(u) < 0.5 * x * x + d * (1 - v + log(v))) return scale * d * v;
    }
}

/* === Pearson correlation unchanged (no RNG) === */
static double pearson_corr(const double *x, const double *y, int n) {
    if (n <= 1) return NAN;
    double sx = 0.0, sy = 0.0;
    for (int i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; }
    double mx = sx / n, my = sy / n;
    double num = 0.0, vx = 0.0, vy = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = x[i] - mx, dy = y[i] - my;
        num += dx * dy;
        vx += dx * dx;
        vy += dy * dy;
    }
    double denom = sqrt(vx * vy);
    if (denom == 0.0) return NAN;
    return num / denom;
}

/* CSV saving helpers (unchanged) */
static void save_pop_csv(const char *basename, int gen,
                         const double *beta1, const double *beta2, const double *N_arr,
                         const int *chromosome_of_gene, const char *chr_id_sequence) {
    char fname[300];
    snprintf(fname, sizeof(fname), "%s.pop.csv", basename);
    FILE *f = fopen(fname, "w");
    if (!f) { perror("fopen pop"); return; }
    fprintf(f, "row,ind,gene_id,chr_pair,chr_id,beta1,beta2,N\n");
    for (int ind = 1; ind <= pop_size; ++ind) {
        for (int g = 1; g <= num_genes; ++g) {
            int baseRow = (ind - 1) * gene_rows + (g - 1) * 2;
            for (int a = 0; a < 2; ++a) {
                int idx = baseRow + a;
                fprintf(f, "%d,%d,%d,%d,%c,%.12g,%.12g,%.12g\n",
                        idx + 1, ind, g, chromosome_of_gene[g-1],
                        chr_id_sequence[a], beta1[idx], beta2[idx],
                        N_arr[(ind-1)*num_genes + (g-1)]);
            }
        }
    }
    fclose(f);
}

static void save_track_csv(const char *basename, int gen, const double *trackFitness, const double *trackCor, int upto) {
    char fname[300];
    snprintf(fname, sizeof(fname), "%s.track.csv", basename);
    FILE *f = fopen(fname, "w");
    if (!f) { perror("fopen track"); return; }
    fprintf(f, "generation,fitness,cor\n");
    for (int i = 0; i < upto; ++i) {
        fprintf(f, "%d,%.12g,%.12g\n", i+1, trackFitness[i], trackCor[i]);
    }
    fclose(f);
}

/* compute_mut_effect: modified to accept RNG pointer and use thread-local RNG */
static double compute_mut_effect(RNG *rng, double beta_val, int n_site, int which_beta, double mut_rate) {
    if (beta_val <= 0.0) return 1.0;

    int k = rpois_one(rng, n_site * mut_rate);
    if (k == 0) {
        double mErrGene = 1.0 - pow(1.0 - beta_val, n_site);
        double mAvgErrSite = 1.0 - pow(1.0 - mErrGene, 1.0 / n_site);
        return mAvgErrSite / beta_val;
    } else {
        double mutated_prob = 1.0;
        for (int j = 0; j < k; ++j) {
            double eff;
            if (which_beta == 1) {
                eff = rgamma_one(rng, shape_est_TR, scale_est_TR);
            } else {
                eff = rgamma_one(rng, shape_est_TL, scale_est_TL);
            }

            double mErrSite = beta_val * eff;
            if (mErrSite < 0.0) mErrSite = 0.0;
            if (mErrSite > (1.0 - 1e-6)) mErrSite = 1.0 - 1e-6;
            mutated_prob *= (1.0 - mErrSite);
        }
        double unchanged_prob = pow(1.0 - beta_val, n_site - k);
        double mErrGene = 1.0 - unchanged_prob * mutated_prob;
        double mAvgErrSite = 1.0 - pow(1.0 - mErrGene, 1.0 / n_site);
        return mAvgErrSite / beta_val;
    }
}

/* Single simulation run for a mut_rate & param_c */
static void run_simulation(double mut_rate, double param_c) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int thread_id = omp_get_thread_num();
    uint64_t seed_base = ((uint64_t)tv.tv_sec << 32) ^ (uint64_t)tv.tv_usec ^ (uint64_t)(thread_id * 0x9e3779b1);

    RNG rng;
    rng_init(&rng, seed_base);

    char out_basename[128];
    snprintf(out_basename, sizeof(out_basename), "sim_%.3g_%.1g", mut_rate, param_c);

    printf("Thread %d: Simulation start. mut_rate=%.3g, param_c=%.1g, pop_size=%d, num_genes=%d, generations=%d\n",
           thread_id, mut_rate, param_c, pop_size, num_genes, num_generations);

    /* chromosome mapping */
    int chromosome_of_gene[num_genes];
    {
        int p = 0;
        for (int cp = 0; cp < chr_pairs; ++cp) {
            int cnt = chr_pair_counts[cp];
            for (int j = 0; j < cnt && p < num_genes; ++j) {
                chromosome_of_gene[p++] = cp + 1;
            }
        }
        if (p != num_genes) {
            fprintf(stderr, "Error: chromosome gene count mismatch (expected %d got %d)\n", num_genes, p);
            return;
        }
    }
    char chr_id_sequence[3] = {'F','M','\0'};

    /* allocate main arrays */
    double *beta1 = (double*)malloc(sizeof(double) * total_rows);
    double *beta2 = (double*)malloc(sizeof(double) * total_rows);
    double *N_arr = (double*)malloc(sizeof(double) * ng);
    if (!beta1 || !beta2 || !N_arr) { fprintf(stderr, "alloc fail main\n"); return; }

    /* initialize */
    for (int i = 0; i < total_rows; ++i) {
        beta1[i] = 1e-5;
        beta2[i] = 1e-4;
    }
    for (int i = 0; i < ng; ++i) {
        N_arr[i] = rpar_one(&rng, min_expr, powerlaw_param);
    }

    /* allocate other required arrays ONCE per simulation (reuse inside loop) */
    double *beta1_mean = (double*)malloc(sizeof(double) * ng);
    double *beta2_mean = (double*)malloc(sizeof(double) * ng);
    double *err1 = (double*)malloc(sizeof(double) * ng);
    double *err2 = (double*)malloc(sizeof(double) * ng);
    double *totalCorrect = (double*)malloc(sizeof(double) * ng);
    double *bothErr = (double*)malloc(sizeof(double) * ng);
    double *trackFitness = (double*)malloc(sizeof(double) * num_generations);
    double *trackCor = (double*)malloc(sizeof(double) * num_generations);
    double *fitness = (double*)malloc(sizeof(double) * pop_size);
    double *cdf = (double*)malloc(sizeof(double) * pop_size);
    int *prolif_M = (int*)malloc(sizeof(int) * pop_size);
    int *prolif_F = (int*)malloc(sizeof(int) * pop_size);

    int rows_per_selected = pop_size * num_genes;
    int *selected_indices_F = (int*)malloc(sizeof(int) * rows_per_selected);
    int *selected_indices_M = (int*)malloc(sizeof(int) * rows_per_selected);

    double *popF_beta1 = (double*)malloc(sizeof(double) * rows_per_selected);
    double *popF_beta2 = (double*)malloc(sizeof(double) * rows_per_selected);
    double *popM_beta1 = (double*)malloc(sizeof(double) * rows_per_selected);
    double *popM_beta2 = (double*)malloc(sizeof(double) * rows_per_selected);

    double *new_beta1 = (double*)malloc(sizeof(double) * total_rows);
    double *new_beta2 = (double*)malloc(sizeof(double) * total_rows);

    /* arrays used for per-generation statistics - allocate once */
    double *cellCorr = (double*)malloc(sizeof(double) * pop_size);
    double *x = (double*)malloc(sizeof(double) * num_genes);
    double *y = (double*)malloc(sizeof(double) * num_genes);
    double *tmp = (double*)malloc(sizeof(double) * pop_size);

    if (!beta1_mean || !beta2_mean || !err1 || !err2 || !totalCorrect || !bothErr
        || !trackFitness || !trackCor || !fitness || !cdf || !prolif_M || !prolif_F
        || !selected_indices_F || !selected_indices_M || !popF_beta1 || !popF_beta2 || !popM_beta1 || !popM_beta2
        || !new_beta1 || !new_beta2 || !cellCorr || !x || !y || !tmp) {
        fprintf(stderr, "alloc fail tmp\n");
        /* free what we allocated before returning */
        free(beta1); free(beta2); free(N_arr);
        free(beta1_mean); free(beta2_mean); free(err1); free(err2); free(totalCorrect); free(bothErr);
        free(trackFitness); free(trackCor); free(fitness); free(cdf); free(prolif_M); free(prolif_F);
        free(selected_indices_F); free(selected_indices_M);
        free(popF_beta1); free(popF_beta2); free(popM_beta1); free(popM_beta2);
        free(new_beta1); free(new_beta2);
        free(cellCorr); free(x); free(y); free(tmp);
        return;
    }

    /* main loop */
    for (int thisGen = 1; thisGen <= num_generations; ++thisGen) {
        /* 1 & 2) compute per-allele mut_effect (site-wise independent) and apply */
        for (int r = 0; r < total_rows; ++r) {
            beta1[r] *= compute_mut_effect(&rng, beta1[r], (int)Lnt, 1, mut_rate);
            beta2[r] *= compute_mut_effect(&rng, beta2[r], (int)Lp, 2, mut_rate);
        }

        /* 3) per-(ind,gene) summaries */
        for (int ind = 1; ind <= pop_size; ++ind) {
            for (int g = 1; g <= num_genes; ++g) {
                int alleleA = (ind-1)*gene_rows + (g-1)*2;
                int alleleB = alleleA + 1;
                int flat = (ind-1)*num_genes + (g-1);
                double b1m = 0.5 * (beta1[alleleA] + beta1[alleleB]);
                double b2m = 0.5 * (beta2[alleleA] + beta2[alleleB]);
                beta1_mean[flat] = b1m;
                beta2_mean[flat] = b2m;
                double er1v = 1.0 - pow(1.0 - b1m, Lnt);
                double er2v = 1.0 - pow(1.0 - b2m, Lp);
                err1[flat] = er1v;
                err2[flat] = er2v;
                totalCorrect[flat] = (1.0 - er1v) * (1.0 - er2v);
                bothErr[flat] = er1v * er2v;
            }
        }

        /* 4) fitness per individual */
        for (int ind = 1; ind <= pop_size; ++ind) {
            double totalErrCnt = 0.0;
            int all_enough_correct = 1;
            for (int g = 1; g <= num_genes; ++g) {
                int flat = (ind-1)*num_genes + (g-1);
                double Nval = N_arr[flat];
                double er1v = err1[flat], er2v = err2[flat];
                totalErrCnt += Nval * (er1v*(1.0 - er2v) + er2v*(1.0 - er1v) + e_epistasis * er1v * er2v);
                if (!(totalCorrect[flat] > min_correct_frac)) all_enough_correct = 0;
            }
            double enoughCorrect = all_enough_correct ? 1.0 : 0.0;
            fitness[ind-1] = enoughCorrect * exp(-param_c * totalErrCnt);
        }

        /* 5) sample parents using cdf */
        double total_fitness = 0.0;
        for (int i = 0; i < pop_size; ++i) total_fitness += fitness[i];
        if (total_fitness <= 0.0) {
            double step = 1.0 / pop_size;
            cdf[0] = step;
            for (int i = 1; i < pop_size; ++i) cdf[i] = cdf[i-1] + step;
        } else {
            double cum = 0.0;
            for (int i = 0; i < pop_size; ++i) {
                cum += fitness[i] / total_fitness;
                cdf[i] = cum;
            }
            cdf[pop_size-1] = 1.0;
        }

        for (int i = 0; i < pop_size; ++i) {
            double v = rng_uniform(&rng);
            int lo = 0, hi = pop_size - 1;
            while (lo < hi) {
                int mid = (lo + hi) >> 1;
                if (v <= cdf[mid]) hi = mid; else lo = mid + 1;
            }
            prolif_M[i] = lo + 1;
            v = rng_uniform(&rng); lo = 0; hi = pop_size - 1;
            while (lo < hi) {
                int mid = (lo + hi) >> 1;
                if (v <= cdf[mid]) hi = mid; else lo = mid + 1;
            }
            prolif_F[i] = lo + 1;
        }

        /* 6) select indices for chosen parent alleles depending on p_F / p_M */
        char p_F = (rng_uniform(&rng) < 0.5) ? 'M' : 'F';
        char p_M = (rng_uniform(&rng) < 0.5) ? 'M' : 'F';
        int allele_offset_F = (p_F == 'F') ? 0 : 1;
        int allele_offset_M = (p_M == 'F') ? 0 : 1;

        int selFcount = 0, selMcount = 0;
        for (int i = 0; i < pop_size; ++i) {
            int pid = prolif_F[i];
            int start_row = (pid - 1) * gene_rows;
            for (int g = 0; g < num_genes; ++g) {
                int idx = start_row + g*2 + allele_offset_F;
                selected_indices_F[selFcount++] = idx;
            }
        }
        for (int i = 0; i < pop_size; ++i) {
            int pid = prolif_M[i];
            int start_row = (pid - 1) * gene_rows;
            for (int g = 0; g < num_genes; ++g) {
                int idx = start_row + g*2 + allele_offset_M;
                selected_indices_M[selMcount++] = idx;
            }
        }

        /* 7) build popF/popM arrays (copy values) */
        for (int i = 0; i < rows_per_selected; ++i) {
            int idxF = selected_indices_F[i];
            int idxM = selected_indices_M[i];
            popF_beta1[i] = beta1[idxF];
            popF_beta2[i] = beta2[idxF];
            popM_beta1[i] = beta1[idxM];
            popM_beta2[i] = beta2[idxM];
        }

        /* 8) recombination events if any (k = sum_k_all) */
        int k0 = rpois_one(&rng, recomb_chr);
        int sum_k_all = k0 * pop_size * 4;
        if (sum_k_all > 0) {
            for (int ev = 0; ev < sum_k_all; ++ev) {
                int x = (int)(floor(rng_uniform(&rng) * pop_size)) + 1;
                int y = (int)(floor(rng_uniform(&rng) * 4)) + 1;
                /* find genes range for chromosome y */
                int acc = 0, start_gene = 0, end_gene = -1;
                for (int cp = 0; cp < chr_pairs; ++cp) {
                    int cnt = chr_pair_counts[cp];
                    if (y == cp+1) { start_gene = acc; end_gene = acc + cnt - 1; break; }
                    acc += cnt;
                }
                if (end_gene < start_gene) continue;
                int block_start = (x - 1) * num_genes;
                int F_start = block_start + start_gene;
                int F_end = block_start + end_gene;
                int len = F_end - F_start + 1;
                if (len <= 1) continue;
                int recomb_point = 1 + (int)(floor(rng_uniform(&rng) * (len - 1))); /* 1..len-1 */
                for (int t = recomb_point; t < len; ++t) {
                    int iF = F_start + t;
                    int iM = (x - 1) * num_genes + start_gene + t;
                    double tb1 = popF_beta1[iF];
                    double tb2 = popF_beta2[iF];
                    popF_beta1[iF] = popM_beta1[iM];
                    popF_beta2[iF] = popM_beta2[iM];
                    popM_beta1[iM] = tb1;
                    popM_beta2[iM] = tb2;
                }
            }
        }

        /* 9) interleave popF and popM into new_beta arrays */
        for (int i = 0; i < rows_per_selected; ++i) {
            int g0 = i * 2;
            new_beta1[g0] = popF_beta1[i];
            new_beta2[g0] = popF_beta2[i];
            new_beta1[g0 + 1] = popM_beta1[i];
            new_beta2[g0 + 1] = popM_beta2[i];
        }

        /* 10) replace main beta arrays with new ones */
        memcpy(beta1, new_beta1, sizeof(double) * total_rows);
        memcpy(beta2, new_beta2, sizeof(double) * total_rows);

        /* 11) tracking: per-individual correlation median & mean fitness */
        for (int ind = 0; ind < pop_size; ++ind) {
            int start_row = ind * gene_rows;
            for (int g = 0; g < num_genes; ++g) {
                int alleleA = start_row + g*2;
                int alleleB = alleleA + 1;
                x[g] = 0.5 * (beta1[alleleA] + beta1[alleleB]);
                y[g] = 0.5 * (beta2[alleleA] + beta2[alleleB]);
            }
            cellCorr[ind] = pearson_corr(x, y, num_genes);
        }

        /* median of correlations (simple selection sort to find median) */
        memcpy(tmp, cellCorr, sizeof(double) * pop_size);
        for (int ia = 0; ia < pop_size-1; ++ia) {
            for (int jb = ia+1; jb < pop_size; ++jb) {
                if (tmp[ia] > tmp[jb]) {
                    double tt = tmp[ia]; tmp[ia] = tmp[jb]; tmp[jb] = tt;
                }
            }
        }
        if (pop_size % 2 == 1)
            trackCor[thisGen-1] = tmp[pop_size/2];
        else
            trackCor[thisGen-1] = 0.5 * (tmp[pop_size/2 - 1] + tmp[pop_size/2]);

        /* mean fitness */
        double mean_fitness = 0.0;
        for (int i = 0; i < pop_size; ++i) mean_fitness += fitness[i];
        trackFitness[thisGen-1] = mean_fitness / pop_size;

        if (thisGen % 10000 == 0) {
            printf("%s Generation %d . Average fitness %.12g correlation at %.12g\n", out_basename, thisGen, mean_fitness / pop_size, trackCor[thisGen-1]);
            fflush(stdout);
            save_pop_csv(out_basename, thisGen, beta1, beta2, N_arr, chromosome_of_gene, chr_id_sequence);
            save_track_csv(out_basename, thisGen, trackFitness, trackCor, thisGen);
        }
    } /* end generation loop */

    /* final save */
    save_pop_csv(out_basename, num_generations, beta1, beta2, N_arr, chromosome_of_gene, chr_id_sequence);
    save_track_csv(out_basename, num_generations, trackFitness, trackCor, num_generations);

    /* free memory */
    free(beta1); free(beta2); free(N_arr);
    free(beta1_mean); free(beta2_mean); free(err1); free(err2); free(totalCorrect); free(bothErr);
    free(trackFitness); free(trackCor); free(fitness); free(cdf); free(prolif_M); free(prolif_F);
    free(selected_indices_F); free(selected_indices_M);
    free(popF_beta1); free(popF_beta2); free(popM_beta1); free(popM_beta2);
    free(new_beta1); free(new_beta2);
    free(cellCorr); free(x); free(y); free(tmp);

    printf("Thread %d: Simulation finished. Outputs written with basename %s\n", thread_id, out_basename);
}

/* main: iterate parameter combinations in parallel */
int main(int argc, char **argv) {
    printf("Starting multi-threaded simulation with %d parameter combinations\n", total_combinations);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < total_combinations; ++i) {
        int mut_idx = i / num_param_cs;
        int param_idx = i % num_param_cs;
        double mut_rate = mut_rates[mut_idx];
        double param_c = param_cs[param_idx];
        run_simulation(mut_rate, param_c);
    }

    printf("All simulations finished.\n");
    return 0;
}
