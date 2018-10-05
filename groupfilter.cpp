//
//  groupfilter.cpp
//  pdb
//
//  Created by BoZhang on 2/21/18.
//  Copyright © 2018 BoZhang. All rights reserved.
//

#include <stdio.h>
#include "groupfilter.h"
#include <cmath>        // std::pow


#define HITNUM 10
int ghitnum = 10 ;

struct less_than_key
{
    inline bool operator() (const pair<int, double>& p1, const pair<int, double>& p2){
        
        int c1 = 0, c2 = 0 ;
        for (int i = 0 ; i < HITNUM; ++ i){
            if ( p1.first & (1 << i) ) c1 ++;
            if ( p2.first & (1 << i) ) c2 ++;
        }
        
        return (c1 < c2);
    }
};

void VPResult::create_allrp( vector<double> &tp){
    
    int n = pow(2,tp.size())-1;
    for ( int i = 0 ; i <= n; ++ i ){
        
        double p = 0 ;
        vector<bool> r(hitnum,false);
        
        for ( int j = 0 ; j < hitnum; ++ j ){
            if ( i & (1 << j) ) {
                p+= log10(tp[j]);
                r[j] = true;
            }
            else p+= log10(1-tp[j]);
        }
        result.insert(make_pair(r,p));
        
        //in order to save memory, add the following 
        result.clear();
    }
}

void VPResult:: findRepreForRepre(vector<double> &tp, map<vector<bool>, double> &vonodes){
    
    map < vector<bool>, double> tmp;
    map <vector<bool>, double> nextlevel;
    
    for ( auto itm : vonodes){
        int flag = 0 ;
        int count = 0 ;
        pair<vector<bool>, double> jtm;
        for ( int i = 0 ; i < hitnum; ++ i){
            
            vector<bool> prop = itm.first;
            
            if (  prop[i] == true ){
                count++;
                continue;
            }
            
            else{
                prop[i] = true;
                auto it = result.find(prop);
                
                if ( it != result.end() ){
                    pair< vector<bool>, double> jtm = make_pair(it->first, it->second);
                    nextlevel.insert(jtm);
                    flag = 1;
                    break;
                }
                else{
                    jtm = make_pair(prop, itm.second + log10(tp[i]/(1-tp[i])));
                }
            }
        }
        
        if ( flag == 0 && count != hitnum ) {
            represent.insert(jtm);
            nextlevel.insert(jtm);
        }
    }
    
    vonodes.clear();
    for ( auto itm : tmp){
        represent.insert(itm);
        result.insert(itm);
    }
    vonodes = nextlevel;
    
    //for ( auto itm : vonodes) cout << "(" << itm.first << "," << itm.second << ")";
    //cout << "..." << endl;
    
}

void VPResult::findRepresentative(vector<double> &tp, map<vector<bool>, double> &vonodes){
    
    for ( auto itm : result){
        
        int flag = 0 ;
        int count = 0 ;
        pair<vector<bool>, double> jtm;
        for ( int i = 0 ; i < hitnum; ++ i){
            
            vector<bool> prop = itm.first;
            
            //find the next level node, cause next level the number of item increse 1
            if (  prop[i] == true ){
                count ++;
                continue;
            }
            
            else{
                prop[i] = true;
                auto it = result.find(prop);
                if ( it != result.end() ){
                    flag = 1 ;
                    pair<vector<bool>, double> jtm = make_pair(it->first, it->second);
                    vonodes.insert(jtm);
                    break;
                }
                else{
                    jtm = make_pair(prop, itm.second + log10(tp[i]/(1-tp[i])));
                }
            }
        }
        
        if ( flag == 0 && count != hitnum ){
            represent.insert(jtm);
            vonodes.insert(jtm);
        }
    }
    
    for( auto jtm : represent){
        result.insert(jtm);
    }
    
    
    while( vonodes.size() > 0){
        findRepreForRepre(tp,vonodes);
    }
    
}

bool VPResult::checkallgroupresult2(vector<double> &tp){
    
    int n = pow(2,tp.size())-1;
    double p = 0, p1 = 0  ;
    vector<bool> r(hitnum,false);
    vector<bool> r1(hitnum,false);
    
    int i =  0 ;
    for ( int j = 0 ; j < hitnum; ++ j ){
        if ( i & (1 << j) ) {
            p+= log10(tp[j]);
            p1 += log10(tp[j]);
            r[j] = true;
            r1[j] = true;
        }
        else {
            p+= log10(1-tp[j]);
            p1+= log10(1-tp[j]);
        }
    }
    
    int count = 0 ;
    for ( i = 0 ; i <= n ; ++ i){
        double diff = p - p1;
        double diffreal = log10(tp[0]/(1-tp[0]));
        if ( fabs(diff - diffreal) > 0.00000001  ) count++;
        else count --;
    }
    
    return true;
}
bool VPResult::checkallgroupresult(vector<double> &tp){
    
    for ( auto itm : result){
        
        for ( int i = 0; i < hitnum; ++ i){
            
            if (itm.first[i] == false){
                
                vector<bool> nextkey = itm.first;
                nextkey[i] = true;
                auto jtm = result.find(nextkey);
                
                double diff = jtm->second - itm.second;
                double diffreal = log10(tp[i]/(1-tp[i]));
                if ( fabs(diff - diffreal) > 0.00000001  ) return false;
                break;
            }
        }
    }
    
    return true;
}

AnswerNode * VPResult::createMWT(vector<double> &tp){
    AnswerNode * root = new AnswerNode;
    return root;
}

double VPResult::checkgroupresult(vector<double> &tp){
    
    //cout << "result size =" << result.size() << endl;
    
    double timepassed = 0 ;
    
    for ( auto itm : result){
        
        double diffreal ;
        if (itm.first[0] == true ){
            diffreal = log10((1-tp[0])/tp[0]);
        }
        else{
            diffreal = log10(tp[0]/(1-tp[0]));
        }
        
        double diff = itm.second;
        if ( fabs(diff - diffreal) > 0.00000001  ){
            timepassed += 1;
        }
        
        /*
        for ( int i = 0 ; i < hitnum; ++ i){
            
            vector<bool> prop = itm.first;
            
            if (  prop[i] == true ){
                
                prop[i] = false;
                auto jtm = result.find(prop);
                if ( result.find(prop) != result.end() ){
                    
                    clock_t begin = clock();
                    double diff = jtm->second - itm.second;
                    double diffreal = log10((1-tp[i])/tp[i]);
                    
                    if ( fabs(diff - diffreal) > 0.00000001  ) return false;
                    clock_t end = clock();
                    timepassed += (double(end-begin) / CLOCKS_PER_SEC);
                    break;
                }
            }
            
            else{
                prop[i] = true;
                auto jtm = result.find(prop);
                if ( result.find(prop) != result.end() ){
                    clock_t begin = clock();
                    double diff = jtm->second - itm.second;
                    double diffreal = log10(tp[i]/(1-tp[i]));
                    
                    if ( fabs(diff - diffreal) > 0.00000001  ) return false;
                    clock_t end = clock();
                    timepassed += (double(end-begin) / CLOCKS_PER_SEC);
                    break;
                }
            }
        }
        */
        
        
    }
    
    return timepassed;
}

bool VPResult::checkGroupRatio(vector<double> &tp){
    
    for ( auto itm : result){
        
        for ( int i = 0 ; i < hitnum; ++ i){
            
            
            vector<bool> prop = itm.first;
            
            if ( prop[i] == false ){
                
                prop[i] = true;
                if (result.find(prop) == result.end() ){
                    continue;
                }
                
                else {
                    auto jtm = result.find(prop);
                    double diff = itm.second - jtm->second;
                    double diffreal = log10((1-tp[i])/tp[i]);
                    
                    if ( fabs(diff - diffreal) > 0.00000001  ) return false;;
                    break;
                }
            }
        }
    }
    return true;
}

