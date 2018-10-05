//
//  rsa.hpp
//  IVORL
//
//  Created by Boxiang Dong on 12/1/17.
//  Copyright © 2017 Boxiang Dong. All rights reserved.
//

#ifndef rsa_hpp
#define rsa_hpp

#include <stdio.h>
#include <functional>
#include <string>
#include <vector>
#include "gmpxx.h"

class RSA {
public:
    RSA();
    static void initialize ();  //assign p, q, N and hash
    static mpz_class sign (const std::string &s);  //sign s, include the result in sigma, h^e (mod N)
    static bool check (const std::string &s, const mpz_class &sigma);  //check if sigma is the correct signature of s, sigma^d = h (mod N)
    static mpz_class condense (const std::vector<mpz_class> &v);  //generate a condensed signature from a set of signatures
    static mpz_class condense (const std::vector<mpz_class> &v, int l, int r); // enerate a condensed signature from a range
    static mpz_class condenseTwo (mpz_class &v,mpz_class &p );
    static bool check (const std::vector<std::string> &s, const mpz_class &sigma); //check if the condensed sigma is from the set of strings s
    
    
    static mpz_class p;
    static mpz_class q;
    static mpz_class N; //N=pq
    static mpz_class fin;   //fin=lcm(p-1, q-1)
    static mpz_class d; //private key
    static mpz_class e; //public key, 1<e<fin, gcd(e, fin)=1
    static std::hash<std::string> str_hash;
};

#endif /* rsa_hpp */
