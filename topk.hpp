//
//  topk.hpp
//  pdbnew
//
//  Created by BoZhang on 2/26/18.
//  Copyright © 2018 BoZhang. All rights reserved.
//

#ifndef topk_hpp
#define topk_hpp

#include <stdio.h>
#include "groupfilter.h"
#include <vector>
#include <math.h>
#include <float.h>
#include <set>
using namespace std;

struct Answer{
    
    vector<bool> key;
    double value;
};

struct greater1{
    bool operator()(const Answer& a,const Answer& b) const{
        return a.value > b.value ;
    }
};



double getLowerbound(vector<double> &probServer, int topk, int kind);
void getTopk(map< vector<bool>, double > &result, vector<vector<bool>> &id,
             vector<double> &prob, int topk);


void deterministic(vector<double> &probServer, VPResult & test, int topk);

void brute_force(vector<double> &probServer, VPResult & test, int topk);
void mc_simulationvector(vector<double> &probServer, VPResult & test, double delta, int topk,
                         vector<double> &ran);
int mc_simulationvector2(vector<double> &probServer, VPResult & test,
                          double delta, int topk, vector<double> &ran);
void brute_force2(vector<double> &probServer, vector<double> & allanswer);

unsigned long long choose(unsigned long long n, unsigned long long k) ;
int estimateK (int topk, int m);
int estimateK2 (int k, int m);

vector< Answer > Merge_TopK( vector<double> & hittuples, int left, int right, int topk );
vector< Answer > CombineResult( vector< Answer > &RL, vector< Answer > &RR, int topk);

#endif /* topk_hpp */
