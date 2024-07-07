/* Copyright 2022 UCLouvain, Belgium and PQM4 contributors
 *
 * This file is part of pqm4_masked.
 *
 * pqm4_masked is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, version 3.
 *
 * pqm4_masked is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * pqm4_masked. If not, see <https://www.gnu.org/licenses/>.
 */
#include "poly.h"
#include "cbd.h"
#include "ntt.h"
#include "params.h"
#include "symmetric.h"

#include <stdint.h>


/*************************************************
 * Name:        poly_decompress
 *
 * Description: De-serialization and subsequent decompression of a polynomial;
 *              approximate inverse of poly_compress
 *
 * Arguments:   - poly *r:                pointer to output polynomial
 *              - const unsigned char *a: pointer to input byte array (of length
 *KYBER_POLYCOMPRESSEDBYTES bytes)
 **************************************************/
void poly_decompress(poly *r, const unsigned char *a) {
  int i;
#if (KYBER_POLYCOMPRESSEDBYTES == 96)
  for (i = 0; i < KYBER_N; i += 8) {
    r->coeffs[i + 0] = (((a[0] & 7) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 1] = ((((a[0] >> 3) & 7) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 2] = ((((a[0] >> 6) | ((a[1] << 2) & 4)) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 3] = ((((a[1] >> 1) & 7) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 4] = ((((a[1] >> 4) & 7) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 5] = ((((a[1] >> 7) | ((a[2] << 1) & 6)) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 6] = ((((a[2] >> 2) & 7) * KYBER_Q) + 4) >> 3;
    r->coeffs[i + 7] = ((((a[2] >> 5)) * KYBER_Q) + 4) >> 3;
    a += 3;
  }
#elif (KYBER_POLYCOMPRESSEDBYTES == 128)
  for (i = 0; i < KYBER_N; i += 8) {
    r->coeffs[i + 0] = (((a[0] & 15) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 1] = (((a[0] >> 4) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 2] = (((a[1] & 15) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 3] = (((a[1] >> 4) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 4] = (((a[2] & 15) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 5] = (((a[2] >> 4) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 6] = (((a[3] & 15) * KYBER_Q) + 8) >> 4;
    r->coeffs[i + 7] = (((a[3] >> 4) * KYBER_Q) + 8) >> 4;
    a += 4;
  }
#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
  for (i = 0; i < KYBER_N; i += 8) {
    r->coeffs[i + 0] = (((a[0] & 31) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 1] =
        ((((a[0] >> 5) | ((a[1] & 3) << 3)) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 2] = ((((a[1] >> 2) & 31) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 3] =
        ((((a[1] >> 7) | ((a[2] & 15) << 1)) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 4] =
        ((((a[2] >> 4) | ((a[3] & 1) << 4)) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 5] = ((((a[3] >> 1) & 31) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 6] =
        ((((a[3] >> 6) | ((a[4] & 7) << 2)) * KYBER_Q) + 16) >> 5;
    r->coeffs[i + 7] = (((a[4] >> 3) * KYBER_Q) + 16) >> 5;
    a += 5;
  }
#else
#error "KYBER_POLYCOMPRESSEDBYTES needs to be in {96, 128, 160}"
#endif
}


/*************************************************
 * Name:        poly_unpackdecompress
 *
 * Description: Deserialization and subsequent compression of a polynomial of a
 *polyvec, Used to uncompress a polyvec one poly at a time in a loop.
 *
 * Arguments:   - const poly *r:     pointer to output polynomial
 *              - unsigned char *a:  pointer to input byte string representation
 *of a polyvec (of length KYBER_POLYVECCOMPRESSEDBYTES)
 *              - int i:             index of poly in polyvec to decompress
 **************************************************/
void poly_unpackdecompress(poly *r, const unsigned char *a, int i) {
  int j;
#if (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 352))
  for (j = 0; j < KYBER_N / 8; j++) {
    r->coeffs[8 * j + 0] =
        (((a[352 * i + 11 * j + 0] |
           (((uint32_t)a[352 * i + 11 * j + 1] & 0x07) << 8)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 1] =
        ((((a[352 * i + 11 * j + 1] >> 3) |
           (((uint32_t)a[352 * i + 11 * j + 2] & 0x3f) << 5)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 2] =
        ((((a[352 * i + 11 * j + 2] >> 6) |
           (((uint32_t)a[352 * i + 11 * j + 3] & 0xff) << 2) |
           (((uint32_t)a[352 * i + 11 * j + 4] & 0x01) << 10)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 3] =
        ((((a[352 * i + 11 * j + 4] >> 1) |
           (((uint32_t)a[352 * i + 11 * j + 5] & 0x0f) << 7)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 4] =
        ((((a[352 * i + 11 * j + 5] >> 4) |
           (((uint32_t)a[352 * i + 11 * j + 6] & 0x7f) << 4)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 5] =
        ((((a[352 * i + 11 * j + 6] >> 7) |
           (((uint32_t)a[352 * i + 11 * j + 7] & 0xff) << 1) |
           (((uint32_t)a[352 * i + 11 * j + 8] & 0x03) << 9)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 6] =
        ((((a[352 * i + 11 * j + 8] >> 2) |
           (((uint32_t)a[352 * i + 11 * j + 9] & 0x1f) << 6)) *
          KYBER_Q) +
         1024) >>
        11;
    r->coeffs[8 * j + 7] =
        ((((a[352 * i + 11 * j + 9] >> 5) |
           (((uint32_t)a[352 * i + 11 * j + 10] & 0xff) << 3)) *
          KYBER_Q) +
         1024) >>
        11;
  }
#elif (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 320))
  for (j = 0; j < KYBER_N / 4; j++) {
    r->coeffs[4 * j + 0] =
        (((a[320 * i + 5 * j + 0] |
           (((uint32_t)a[320 * i + 5 * j + 1] & 0x03) << 8)) *
          KYBER_Q) +
         512) >>
        10;
    r->coeffs[4 * j + 1] =
        ((((a[320 * i + 5 * j + 1] >> 2) |
           (((uint32_t)a[320 * i + 5 * j + 2] & 0x0f) << 6)) *
          KYBER_Q) +
         512) >>
        10;
    r->coeffs[4 * j + 2] =
        ((((a[320 * i + 5 * j + 2] >> 4) |
           (((uint32_t)a[320 * i + 5 * j + 3] & 0x3f) << 4)) *
          KYBER_Q) +
         512) >>
        10;
    r->coeffs[4 * j + 3] =
        ((((a[320 * i + 5 * j + 3] >> 6) |
           (((uint32_t)a[320 * i + 5 * j + 4] & 0xff) << 2)) *
          KYBER_Q) +
         512) >>
        10;
  }
#else
#error "KYBER_POLYVECCOMPRESSEDBYTES needs to be in {320*KYBER_K, 352*KYBER_K}"
#endif
}


/*************************************************
 * Name:        poly_tobytes
 *
 * Description: Serialization of a polynomial
 *
 * Arguments:   - unsigned char *r: pointer to output byte array (needs space
 *for KYBER_POLYBYTES bytes)
 *              - const poly *a:    pointer to input polynomial
 **************************************************/
void poly_tobytes(unsigned char *r, poly *a) {
  int i;
  uint16_t t0, t1;

  poly_reduce(a);

  for (i = 0; i < KYBER_N / 2; i++) {
    t0 = a->coeffs[2 * i];
    t1 = a->coeffs[2 * i + 1];
    r[3 * i] = t0 & 0xff;
    r[3 * i + 1] = (t0 >> 8) | ((t1 & 0xf) << 4);
    r[3 * i + 2] = (t1 >> 4) & 0xff;
  }
}

/*************************************************
 * Name:        poly_frombytes
 *
 * Description: De-serialization of a polynomial;
 *              inverse of poly_tobytes
 *
 * Arguments:   - poly *r:                pointer to output polynomial
 *              - const unsigned char *a: pointer to input byte array (of
 *KYBER_POLYBYTES bytes)
 **************************************************/
void poly_frombytes(poly *r, const unsigned char *a) {
  int i;

  for (i = 0; i < KYBER_N / 2; i++) {
    r->coeffs[2 * i] = a[3 * i] | ((uint16_t)a[3 * i + 1] & 0x0f) << 8;
    r->coeffs[2 * i + 1] = a[3 * i + 1] >> 4 | ((uint16_t)a[3 * i + 2] & 0xff)
                                                   << 4;
  }
}


/*************************************************
 * Name:        poly_getnoise
 *
 * Description: Sample a polynomial deterministically from a seed and a nonce,
 *              with output polynomial close to centered binomial distribution
 *              with parameter KYBER_ETA
 *
 * Arguments:   - poly *r:                   pointer to output polynomial
 *              - const unsigned char *seed: pointer to input seed (pointing to
 *array of length KYBER_SYMBYTES bytes)
 *              - unsigned char nonce:       one-byte input nonce
 *              - int add:                   boolean to indicate to accumulate
 *into r
 **************************************************/
void poly_noise(poly *r, const unsigned char *seed, unsigned char nonce,
                int add) {
  unsigned char buf[KYBER_ETA * KYBER_N / 4];

  prf(buf, KYBER_ETA * KYBER_N / 4, seed, nonce);
  cbd(r, buf, add);
}

/*************************************************
 * Name:        poly_ntt
 *
 * Description: Computes negacyclic number-theoretic transform (NTT) of
 *              a polynomial in place;
 *              inputs assumed to be in normal order, output in bitreversed
 *order
 *
 * Arguments:   - uint16_t *r: pointer to in/output polynomial
 **************************************************/
void poly_ntt(poly *r) {
  ntt(r->coeffs);
}

/*************************************************
 * Name:        poly_invntt
 *
 * Description: Computes inverse of negacyclic number-theoretic transform (NTT)
 *of a polynomial in place; inputs assumed to be in bitreversed order, output in
 *normal order
 *
 * Arguments:   - uint16_t *a: pointer to in/output polynomial
 **************************************************/
void poly_invntt(poly *r) {
  invntt(r->coeffs);
}


void poly_basemul_i16(int16_t *r, const int16_t *a, const int16_t *b) {
  basemul_asm(r, a, b, zetas);
}


extern void asm_barrett_reduce(int16_t *r);
/*************************************************
 * Name:        poly_reduce
 *
 * Description: Applies Barrett reduction to all coefficients of a polynomial
 *              for details of the Barrett reduction see comments in reduce.c
 *
 * Arguments:   - poly *r:       pointer to input/output polynomial
 **************************************************/
void poly_reduce(poly *r) { asm_barrett_reduce(r->coeffs); }

