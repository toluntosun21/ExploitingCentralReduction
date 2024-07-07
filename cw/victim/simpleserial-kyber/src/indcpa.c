#include "indcpa.h"
#include "ntt.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"

#include <stdint.h>
#include <string.h>


/*************************************************
* Name:        indcpa_keypair_derand
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
*                             (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
*              - const uint8_t *coins: pointer to input randomness
*                             (of length KYBER_SYMBYTES bytes)
**************************************************/
void indcpa_keypair_derand(unsigned char *pk,
                           unsigned char *sk, 
                           const unsigned char *coins){
  poly skp;
  unsigned char buf[2 * KYBER_SYMBYTES];
  unsigned char *noiseseed = buf + KYBER_SYMBYTES;
  unsigned char nonce = 0;

  hash_g(buf, coins, KYBER_SYMBYTES);

  poly_getnoise(&skp, noiseseed, nonce++);

  poly_ntt(&skp);
  
  poly_tobytes(sk, &skp);
}

