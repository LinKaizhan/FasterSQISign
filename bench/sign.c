#define _XOPEN_SOURCE

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <time.h>
#include <pari/pari.h>

#include "precomputed.h"
#include "sqisign.h"
#include "constants.h"

static __inline__ uint64_t rdtsc(void)
{
    uint32_t hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return lo | (uint64_t) hi << 32;
}

int main(int argc, char **argv) {
  int keys = 10, samples = 100, seed = 1;
  uint64_t c_total = 0;
  float t_total = 0;

  int opt;
  while ((opt = getopt(argc, argv, "k:s:r:h")) != -1) {
    switch (opt) {
    case 'k':
      keys = atoi(optarg);
      break;
    case 's':
      samples = atoi(optarg);
      break;
    case 'r':
      seed = atoi(optarg);
      break;
    default:
      fprintf(stderr,
	      "Usage: %s [-k keys] [-s samples] [-r seed]\n",
	      argv[0]);
      exit(-1);
    }
  }

  pari_init(800000000, 1<<18);
  init_precomputations();

  setrand(mkintn(1, seed));
  srand48(seed);

  printf("### Sign\n");
  printf("# key\tcycles\t\tms\t\tlength\n");
  for (int k = 0; k < keys; k++) {
    public_key pk;
    secret_key sk;
    keygen(&pk, &sk);
    #ifdef PRECOMP
    eichler_package eichler;
    eichler.beta=NULL; eichler.J = NULL; eichler.delta=NULL; eichler.I = NULL; eichler.L=NULL;
    eichler.gamma = NULL;
    GEN x1_gen = gen_0, x2_gen = gen_0;
    proj L_K, KD, P3, P, Q, PQ;
    proj phi_I_basis_ini[7];
    precompute_for_first_sign(&eichler, &x1_gen, &x2_gen, &L_K, &KD, &P3, &P, &Q, &PQ, phi_I_basis_ini, &sk);
    #endif
    for (int i = 0; i < samples; i++) {
      // signature Sigma;
      compressed_signature comp_sigma;
      init_compressed_sig(&comp_sigma);
      uintbig m;
      randombytes(m.c, 32);
      #ifdef PRECOMP
      GEN X1 = x1_gen, X2 = x2_gen;
      proj PP = P, QQ = Q, PPQQ = PQ;
      proj phi_I_basis[7];
      phi_I_basis[0] = phi_I_basis_ini[0];
      phi_I_basis[1] = phi_I_basis_ini[1];
      phi_I_basis[2] = phi_I_basis_ini[2];
      phi_I_basis[3] = phi_I_basis_ini[3];
      phi_I_basis[4] = phi_I_basis_ini[4];
      phi_I_basis[5] = phi_I_basis_ini[5];
      phi_I_basis[6] = phi_I_basis_ini[6];
      #endif
      clock_t t = -clock();
      uint64_t c = -rdtsc();
      #ifdef PRECOMP
      sign_new2(&comp_sigma, &sk, &pk, &m, &X1, &X2, &L_K, &KD, &P3, &PP, &QQ, &PPQQ, phi_I_basis, &eichler);
      #else
      sign_new(&comp_sigma ,&sk,&pk, &m);
      #endif
      // sign(&comp_sigma ,&sk, &pk, &m);
      c += rdtsc();
      t += clock();
      free_compressed_sig(&comp_sigma);
  //     int len = 0;
  //     for (int j = 0; j < Sigma.sigma.len; j++)
	// len += Sigma.sigma.phi[j].len;

      printf("%d\t%" PRIu64 "\t%.3lf\t%ld\n", k, c, 1000. * t / CLOCKS_PER_SEC, signing_length);
      c_total += c;
      t_total += (float)(1000. * t / CLOCKS_PER_SEC);
    }
  }
  printf("\nAvg\t\t%" PRIu64 "\t%.3lf\n", c_total/(keys*samples), t_total/(keys*samples));

  return 0;
}
