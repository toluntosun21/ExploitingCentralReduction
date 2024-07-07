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
#ifndef MASKED_UTILS_H
#define MASKED_UTILS_H
#include "masked.h"
#include "poly.h"
#include "randombytes.h"
#include <stdint.h>


extern int prng_off;


inline uint32_t rand32() {
  uint32_t rand;
  if (prng_off){
    return 0;
  }
  else
  {
    randombytes((uint8_t*) &rand, sizeof(rand));
    return rand;
  }
}

inline void rand_q(int16_t v[2]) {
  uint32_t r;
  do {
    r = rand32();
  } while (r > (387U * KYBER_Q * KYBER_Q));
  r = r % (KYBER_Q * KYBER_Q);
  v[0] = r % (KYBER_Q);
  v[1] = r / KYBER_Q;
}

void masked_poly(StrAPoly mp, const poly *p);
void unmasked_poly(poly *p, const StrAPoly mp);
#endif
