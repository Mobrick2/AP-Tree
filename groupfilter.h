//
//  groupfilter.h
//  pdb
//
//  Created by BoZhang on 2/21/18.
//  Copyright © 2018 BoZhang. All rights reserved.
//

#ifndef groupfilter_h
#define groupfilter_h

#include <map>
#include <bitset>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;


struct AnswerNode{
    vector<bool> key;
    double value;
    vector<AnswerNode *> child;
};

class VPResult{
public:
    VPResult(int n){
        hitnum = n;
    }
    
    VPResult();
    
    
    map < vector<bool>, double> result;
    map < vector<bool>, double> represent;
    //the number of hit tuple
    int hitnum;
    
    void insert_rp( vector<bool> r , double p){
        result.insert(make_pair(r, p));
    }
    
    //generate all possible result
    void create_allrp( vector<double> &tp);
    
    //put 2 (r,p) into one group; we have hitnum of groups
    bool checkallgroupresult(vector<double> &tp);
    bool checkallgroupresult2(vector<double> &tp);
    double checkgroupresult(vector<double> &tp);
    bool checkGroupRatio(vector<double> &tp);
    
    void findRepresentative(vector<double> &tp, map<vector<bool>, double> &vonodes);
    void findRepreForRepre(vector<double> &tp, map<vector<bool>, double> &vonodes);
    
    AnswerNode * createMWT(vector<double> &tp);
    
    int gethitnum(){
        return this->hitnum;
    }
    
    
    void printresult(){
        for ( auto itm : result){
            cout << "(" << itm.second << ")" << endl;
        }
    }
    
    void printTrueProb(){
        for ( auto itm : result){
            cout << "(" <<  pow(10,itm.second) << ")" << endl;
        }
    }
    
    void printRepreProb(){
        for ( auto itm : represent){
            cout << "(" << pow(10,itm.second) << ")" << endl;
        }
    }
    
};


#endif /* groupfilter_h */
