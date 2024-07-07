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
#ifndef POLY_H
#define POLY_H

#include "params.h"

#include <stdint.h>

#define poly_getnoise(p, seed, nonce) poly_noise(p, seed, nonce, 0)

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct {
  int16_t coeffs[KYBER_N];
} poly;

void poly_decompress(poly *r, const unsigned char *a);

void poly_unpackdecompress(poly *r, const unsigned char *a, int i);

void poly_tobytes(unsigned char *r, poly *a);
void poly_frombytes(poly *r, const unsigned char *a);

void poly_noise(poly *r, const unsigned char *seed, unsigned char nonce,
                int add);

void poly_ntt(poly *r);
void poly_invntt(poly *r);
void poly_basemul_i16(int16_t *r, const int16_t *a, const int16_t *b);

void poly_reduce(poly *r);

#endif
