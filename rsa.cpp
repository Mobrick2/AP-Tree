//
//  rsa.cpp
//  IVORL
//
//  Created by Boxiang Dong on 12/1/17.
//  Copyright © 2017 Boxiang Dong. All rights reserved.
//

#include "rsa.hpp"
#include <random>
#include <functional>

#include "gmp.h"
#include "gmpxx.h"

using namespace std;


mpz_class RSA::p;
mpz_class RSA::q;
mpz_class RSA::N;
mpz_class RSA::e;
mpz_class RSA::d;
mpz_class RSA::fin;
std::hash<std::string>  RSA::str_hash;

void
RSA::initialize () {
    gmp_randclass rand (gmp_randinit_mt); //the random generator
    mp_bitcnt_t lambda ;    //number of bits
    default_random_engine generator;
    uniform_int_distribution<int> distribution(1, 10000);
    unsigned long int seed;
    mpz_class one;
    
    //initialize the random number generator
    seed = distribution(generator);
    rand.seed(seed);
    
    //generate p and q
    lambda = 100;
    do {
        mpz_class tmp = rand.get_z_bits(lambda);
        mpz_nextprime(RSA::p.get_mpz_t(), tmp.get_mpz_t());
    }
    while(mpz_probab_prime_p(RSA::p.get_mpz_t(), 30) == 0);
    
    do {
        mpz_class tmp = rand.get_z_bits(lambda);
        mpz_nextprime(RSA::q.get_mpz_t(), tmp.get_mpz_t());
    }
    while(mpz_probab_prime_p(RSA::q.get_mpz_t(), 30) == 0);

    //calculate N
    RSA::N = RSA::p * RSA::q;
    
    //calculate fin
    one = 1;
    RSA::fin = lcm(RSA::p-one, RSA::q-one);

    //generate e
    do {
        do {
            RSA::e = rand.get_z_range(RSA::fin);
        }
        while ((cmp(gcd(RSA::e, RSA::fin), 1) != 0) || (cmp(RSA::e, 1) == 0) || (cmp(RSA::e, 1) == 0));
    }
    while (mpz_invert(RSA::d.get_mpz_t(), RSA::e.get_mpz_t(), RSA::fin.get_mpz_t()) == 0);
    
    //initialize hash
    RSA::str_hash = std::hash<std::string> ();
}


mpz_class
RSA::sign (const std::string &s) {
    size_t h;
    mpz_class sigma, ht;
    
    //generate a hash value h from s
    h = RSA::str_hash(s);
    ht = h;
    
    //generate h^e (mod N)
    mpz_powm(sigma.get_mpz_t(), ht.get_mpz_t(), RSA::e.get_mpz_t(), RSA::N.get_mpz_t());
    
    return sigma;
}


bool
RSA::check (const std::string &s, const mpz_class &sigma) {
    size_t h;
    mpz_class ht;
    
    //generate h from s
    h = RSA::str_hash(s);
    
    //get ht=sigma^d (mod N)
    mpz_powm(ht.get_mpz_t(), sigma.get_mpz_t(), RSA::d.get_mpz_t(), RSA::N.get_mpz_t());
    
    //compare h with ht
    if (cmp(ht, h) == 0)
        return true;
    else
        return false;
}

mpz_class
RSA::condense (const std::vector<mpz_class> &v) {
    mpz_class sigma;
    size_t size;
    
    size = v.size();
    
    sigma = 1;
    
    for (int i = 0; i < v.size(); i++)
        sigma = (sigma * v[i]) % RSA::N;
    
    return sigma;
}

mpz_class
RSA::condense (const std::vector<mpz_class> &v, int l, int r) {
    mpz_class sigma;
    size_t size;
    
    size = v.size();
    
    sigma = 1;
    
    for (int i = l; i <= r ; i++)
        sigma = (sigma * v[i]) % RSA::N;
    
    return sigma;
}

mpz_class
RSA::condenseTwo (mpz_class &v,mpz_class &p ) {
    
    mpz_class sigma;
    sigma = 1;
    sigma = (sigma * v) % RSA::N;
    sigma = (sigma * p) % RSA::N;
    
    return sigma;
}

bool
RSA::check(const std::vector<std::string> &s, const mpz_class &sigma) {
    mpz_class hs;   //ht1*ht2...
    mpz_class ht;   //ht=hash(s[i])
    mpz_class sigma_h;  //sigma^d (mod N)
    
    //calculate sigma^d (mod N)
    mpz_powm(sigma_h.get_mpz_t(), sigma.get_mpz_t(), RSA::d.get_mpz_t(), RSA::N.get_mpz_t());
    
    //calculate h[s[i]] (mod N)
    hs = 1;
    for (int i = 0; i < s.size(); i++) {
        ht = RSA::str_hash(s[i]);
        hs = (hs * ht) % RSA::N;
    }
    
    //compare hs with sigma_h
    if (cmp(sigma_h, hs) == 0)
        return true;
    else
        return false;
}
