//
//  topk.cpp
//  pdbnew
//
//  Created by BoZhang on 2/26/18.
//  Copyright © 2018 BoZhang. All rights reserved.
//

#include "topk.hpp"

void printTopk(vector<int> &id, vector<double> &prob){
    
    int n = int (id.size());
    for ( int i = 0 ; i < n ; ++ i){
        cout << "(" << id[i] << "," << prob[i] << ")\n";
    }
}

void brute_force2(vector<double> &probServer, vector<double> & allanswer){
    
    int m = int (probServer.size());
    int N = int(pow(2,m));
    
    for(int i = 0 ; i < N; ++ i){
        
        double pr = 0 ;
        vector<bool> key(m,false);
        
        for ( int j = 0 ; j < m; ++ j){
            if ( i & (1 << j ) ){
                pr += log10(probServer[j]);
                key[j] = true;
            }
            else{
                pr += log10(1-probServer[j]);
            }
        }
        allanswer.push_back(pr);
    }
    
    
}


void deterministic(vector<double> &probServer, VPResult & test, int topk){
    
    
    vector< Answer > results;
    
    int m = int (probServer.size());
    //find the maximum value
    vector<bool> obj(m, true);
    double prob = 0 ;
    
    for ( int i = 0 ; i < m; ++ i){
        if (probServer[i] > (1-probServer[i]) ){
            prob += log10(probServer[i]);
        }
        else{
            prob += log10(1-probServer[i]);
            obj[i] = false;
        }
    }
    
    Answer ans; ans.key = obj; ans.value = prob;
    results.push_back(ans);
    make_heap(results.begin(), results.end(),greater1());
    map < double, int> countProb;
    countProb.insert(make_pair(prob, 1));
    
    for ( int i = 0 ; i < m ; ++ i){
        
        vector< Answer > results_copy;
        for ( auto itm : results){
            results_copy.push_back(itm);
        }
        
        int k = int (results_copy.size()) ;
        
        for ( int j = 0 ; j < k ; ++ j){
            
            double tmp_value = results_copy[j].value;
            vector<bool> tmp_key = results_copy[j].key;
            
            if ( tmp_key[i] == false){
                tmp_value += (log10( probServer[i]/(1-probServer[i])) );
                tmp_key[i] = true;
            }
            else {
                tmp_key[i] = false;
                tmp_value += (log10((1-probServer[i])/ probServer[i]) );
            }
            
            if ( results.size() < topk ){
                Answer new_element; new_element.value = tmp_value;
                new_element.key = tmp_key;
                results.push_back(new_element);
                make_heap(results.begin(), results.end(),greater1());
            }
            else{
                if (  fabs(results[0].value - tmp_value) < 0.0000001 ){
                    Answer new_element; new_element.value = tmp_value;
                    new_element.key = tmp_key;
                    results.push_back(new_element);
                    make_heap(results.begin(), results.end(),greater1());
                }
                else if ( results[0].value < tmp_value){
                    
                    int mincount = countProb.find(results[0].value)->second;
                    
                    if ( (results.size()+1 - mincount) == topk){
                        for ( int l = 0 ; l < mincount; ++ l){
                            results.erase (results.begin());
                            make_heap(results.begin(), results.end(),greater1());
                        }
                    }
                    
                    Answer new_element; new_element.value = tmp_value;
                    new_element.key = tmp_key;
                    
                    results.push_back(new_element);
                    make_heap(results.begin(), results.end(),greater1());
                }
            }
            
            
            //count the number of each possible answers appeared
            if ( countProb.find(tmp_value) != countProb.end() ){
                auto jtm = countProb.find(tmp_value);
                jtm->second = jtm->second+1;
            }
            else{
                countProb.insert(make_pair(tmp_value, 1));
            }
        }
    }
    
    for ( int i = 0 ; i < results.size(); ++ i){
        test.insert_rp(results[i].key, results[i].value);
    }
}


void brute_force(vector<double> &probServer, VPResult & test, int topk){
    
    
    int m = int (probServer.size());
    int N = int(pow(2,m));
    
    vector< vector<bool> > id(topk,vector<bool>(m,false));
    vector<double> prob(topk,0);
    
    for ( int i = 0 ; i < topk; ++i ){
        
        double pr = 0 ;
        vector<bool> key(m,false);
        
        for ( int j = 0 ; j < m; ++ j){
            if ( i & (1 << j ) ){
                pr += log10(probServer[j]);
                key[j] = true;
            }
            else{
                pr += log10(1-probServer[j]);
            }
        }
        prob[i] = pr;
        id[i] = key;
    }
    
    for ( int i = 0 ; i < topk; ++ i){
        for ( int j = 0 ; j < topk-i-1; ++ j){
            if ( prob[j] > prob[j+1] ){
                swap(prob[j],prob[j+1]);
                swap(id[j],id[j+1]);
            }
        }
    }
            
            
    for(int i = topk ; i < N; ++ i){
        
        double pr = 0 ;
        vector<bool> key(m,false);
        
        for ( int j = 0 ; j < m; ++ j){
            if ( i & (1 << j ) ){
                pr += log10(probServer[j]);
                key[j] = true;
            }
            else{
                pr += log10(1-probServer[j]);
            }
        }
        
        if (pr <= prob[0]) continue;
        
        if ( pr >= prob[topk-1] ){
            for (int k = 0 ; k < topk-1; ++ k){
                prob[k] = prob[k+1];
                id[k] = id[k+1];
            }
            prob[topk-1] = pr;
            id[topk-1] = key;
            continue;
        }
        
        int index = 0 ;
        for (index = 0 ; index < topk ; ++ index){
            if (prob[index] < pr) continue;
            else break;
        }
        
        for (int k = 0 ; k < index-1; ++ k ){
            prob[k] = prob[k+1];
            id[k] = id[k+1];
        }
        
        if ( index == 0 ){
            continue;
        }
        
        prob[index-1] = pr;
        id[index-1] = key;
    }
    
    for ( int i = 0 ; i < topk; ++ i){
        test.insert_rp(id[i], prob[i]);
    }
}

void getTopk(map< vector<bool>, double > &result, vector<vector<bool>> &id,
             vector<double> &prob, int topk){
    
    int i = 0 ;
    for (auto itm : result){
        
        if ( i < topk){
            prob[i] = itm.second;
            id[i] = itm.first;
        }
        else{
            break;
        }
        i++;
    }
    
    for ( i = 0 ; i < topk; ++ i){
        for ( int j = 0 ; j < topk-i-1; ++ j){
            if ( prob[j] > prob[j+1] ){
                swap(prob[j],prob[j+1]);
                swap(id[j],id[j+1]);
            }
        }
    }
    
    i = 0 ;
    for (auto itm : result){
        
        if ( i < topk ) {
            i++;
            continue;
        }
        
        double pr = itm.second;
        if ( pr <= prob[0] ) {
            i++;
            continue;
        }
        
        if ( pr >= prob[topk-1] ){
            for (int k = 0 ; k < topk-1 ; ++ k){
                prob[k] = prob[k+1];
                id[k] = id[k+1];
            }
            
            prob[topk-1] = pr;
            id[topk-1] = itm.first;
            i++;
            continue;
        }
        
        //1 2 3 4 5 (4)
        int index = 0 ;
        for ( index = 0 ; index < topk; ++ index){
            if (prob[index] < pr ) continue;
            else break;
        }
        
        for (int k = 0 ; k < index-1 ; ++ k){
            prob[k] = prob[k+1];
            id[k] = id[k+1];
        }
        
        if ( index == 0 ){
            i++; continue;
        }
        
        prob[index-1] = pr;
        id[index-1] = itm.first;
        i++;
    }
    
}

double getLowerbound( vector<double> &probServer, int topk, int kind){
    
    int m = int (probServer.size());
    int k = ceil(log2(topk));
    double plb = 0;
    
    vector<double> tmp = probServer;
    for ( int i = 0 ; i < m; ++ i ){
        if ( tmp[i] < (1-tmp[i]) ){
            tmp[i] = 1-tmp[i];
        }
    }
    sort(tmp.begin(),tmp.end());
    
    if ( kind == 0 ){
        for ( int i = 0 ; i < m; ++ i){
            if ( i < k ) plb += log10(1-tmp[i]);
            else plb += log10(tmp[i]);
        }
    }
    
    else{
        
        
        map< vector<bool>, double > result;
        
        int count = 0 ;
        
        for ( int j = 0 ; j <= 25 ; ++ j){
            
            
            vector<bool> v(m);
            fill(v.end() - j, v.end(), true);
            
            do {
                double tmpProb = 0 ;
                for (int i = 0; i < m; ++i) {
                    if ( !v[i]) {
                        tmpProb += log10(tmp[i]);
                    }
                    else{
                        tmpProb += log10(1-tmp[i]);
                    }
                }
                
                result.insert(make_pair(v, tmpProb));
                count ++;
            } while (next_permutation(v.begin(), v.end()));
            
            if ( count > topk ) break;
        }
        
        //get the topk answer from map and use the k-th answer as bound.
        vector<vector<bool>> id(topk,vector<bool>(m,false));
        vector<double> prob(topk, -DBL_MAX);
        getTopk(result, id, prob, topk);
        
        plb = prob[topk-1];
        
    }
    
    
    return plb;
    
    
}

unsigned long long
choose(unsigned long long n, unsigned long long k) {
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
    
}

int permulateValues( int n, int m){
    
    int result = 1 ;
    for ( int i = n ; i >= n-m+1 ; -- i){
        result *= i;
    }
    return result;
}


int mc_simulationvector2(vector<double> &probServer, VPResult & test,
                         double delta, int topk, vector<double> &ran){
    
    int k = ceil(log2(topk));
    int m = int (probServer.size());
    double plb = getLowerbound(probServer, topk, 0);
    double plb2 = getLowerbound(probServer, topk, 1);
    if ( plb < plb2 ) plb = plb2;
    
    int x = 0;
    for ( int i = 0 ; i <= k ; ++ i){
        x += choose(m,i);
    }
    
    return x;
}

int estimateK (int topk, int m){
    
    int x = 0 ;
    int k = ceil(log2(topk));
    for ( int i = 0 ; i < k ; ++ i){
        x += choose(m,i);
        
        if ( x >= topk ){
            break;
        }
    }
    
    return x+1;
}

int estimateK2 (int k, int m){
    
    int x = 0 ;
    for ( int i = 0 ; i < k ; ++ i){
        x += choose(m,i);
    }
    
    return x+1;
}

void mc_simulationvector(vector<double> &probServer, VPResult & test,
                         double delta, int topk, vector<double> &ran){
    
    
    int m = int (probServer.size());
    double plb = getLowerbound(probServer, topk, 0);
    //cout << plb << endl;
    double plb2 = getLowerbound(probServer, topk, 1);
    
    //cout  << plb << " " << plb2 << endl;
    //clock_t end = clock();
    //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //cout << elapsed_secs << "\t";
    
    if ( plb < plb2) {
        plb = plb2;
        //cout << 1 ;
    }
    
    //clock_t begin = clock();
    double realplb = pow(10,plb);
    //cout << "**" << realplb << endl;
    //double nthRoot= pow (delta,1.0/topk);
    //int N = ceil(log2(1- nthRoot)/log2(1-realplb)) ;
    int N = ceil(log2((1-delta))/log2(1-realplb));
    
    //double plbreal = 0.02;
    //double delta2 = 1 - pow(1-realplb, m*topk) ;
    cout << ceil(log2((1-delta))/log2(1-realplb)) << " " << m*topk <<  endl;
    //int N = ceil(log2((1-delta))/log2(1-realplb));
    //cout << N << endl;
    
    //int index = 0 ;
 
    plb = 0;
    vector<bool> obj(test.hitnum,true);
    for ( int i = 0 ; i < m; ++ i ){
        
        if ( probServer[i] < (1-probServer[i]) ){
            plb += log10(1-probServer[i]);
            obj[i] = false;
        }
        else{
            plb += log10(probServer[i]);
        }
    }
    
    
    /*
    srand((unsigned)time(NULL));
    vector<bool> obj(test.hitnum,false);
    for (int k = 0; k < m; ++ k){
        
        double num =((double) rand() / (RAND_MAX) );
        //double num = ran[index++] / RAND_MAX;
        if ( probServer[k] >= num ){
            obj[k] = true;
            plb += log10(probServer[k]);
        }
        else{
            plb += log10(1-probServer[k]);
        }
    }
    */
    
    
    
    vector<Answer> results;
    Answer ans; ans.value = plb; ans.key = obj;
    results.push_back(ans);
    set < vector<bool> > keys;
    keys.insert(obj);
    

    double ratio = 1;
    double Pt = 1;
    int bit; double num;
    int k = ceil(log2(topk));
    double tmpplb = plb;
    vector<bool> tmpobj = obj;

    for (int j = 0 ; j < N ; ++ j){
        
        plb = tmpplb;
        obj = tmpobj;
        
        for ( int l = 1 ; l <= k ; ++ l){
        
        bit = rand() % m;
        //bit = int (ran[index++]) % m;
        if ( obj[bit] == true ){
            ratio = (1-probServer[bit])/probServer[bit];
            Pt = fmin(1,ratio);
        }
        else{
            ratio = probServer[bit]/(1-probServer[bit]);
            Pt = fmin(1,ratio);
        }
        
        num = ((double) rand() / (RAND_MAX) );
        //num = ran[index++] / RAND_MAX;
        if (Pt >= num){
            if ( obj[bit] == true ){
                obj[bit] =  false;
            }
            else{
                obj[bit] = true;
            }
            
            plb += log10(ratio);

        }
        }
        
        if ( keys.find(obj) == keys.end() && results.size() < topk ){
            Answer newans; newans.key = obj;
            newans.value = plb;
            results.push_back(newans);
            make_heap(results.begin(),  results.end(), greater1());
            keys.insert(obj);
        }
        else{
            if ( keys.find(obj) == keys.end() && plb > results[0].value ){
                results[0].value = plb;
                results[0].key = obj;
                make_heap(results.begin(),  results.end(), greater1());
                keys.insert(obj);
            }
        }
        
    }
    
    //clock_t end = clock();
    //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //cout << elapsed_secs << "\t----\n";
    
    
    //begin = clock();
    
    for ( int i = 0 ; i < results.size(); ++ i){
        test.insert_rp(results[i].key, results[i].value);
    }
    
    //end = clock();
    //elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //cout << elapsed_secs << "\t";
    //cout << pow(2, m) << "\t" << N << endl;
}


vector< Answer > CombineResult( vector<Answer> & RL, vector< Answer > & RR, int topk){
    
    vector< Answer > R;
    multimap< double, vector<bool>, greater<double> > comR;
    
    int k1 = int (RL.size()), k2 = int (RR.size());
    k1 = min(k1, topk);
    
    for ( int i = 0 ; i < k1 ; ++ i){
        
        int m = (topk % (i+1) == 0) ? topk/(i+1) : topk/(i+1)+1 ;
        k2 = min(k2, m);

        for( int j = 0 ; j < k2 ; ++ j ){
            vector<bool> comkey = RL[i].key;
            comkey.insert(comkey.end(), RR[j].key.begin(), RR[j].key.end());
            double comval = RL[i].value + RR[j].value;
            
            double fakeval = 0 ;
            for ( int j = 0 ; j < k1+k2-1; ++ j){
                fakeval += log10(RL[i].value);
            }
            
            comR.insert(make_pair(comval, comkey));
        }
    }
    
    int k = 0 ;
    double kthvalue = 0;
    for ( auto itm : comR){
        
        if ( k < topk ){
            Answer ans; ans.key = itm.second; ans.value = itm.first;
            R.push_back(ans);
            
            if ( k == topk - 1){
                kthvalue = itm.first;
            }
        }
        else{
            
            if ( fabs(kthvalue - itm.first) < 0.00000001 ){
                Answer ans; ans.key = itm.second; ans.value = itm.first;
                R.push_back(ans);
            }
            else{
                break;
            }
        }
        k++;
    }

    return R;
}

vector< Answer > Merge_TopK( vector<double> & hittuples, int left, int right, int topk ){
    
    vector< Answer > R;
    if ( left < right ){
        int mid = (right-left)/2 + left;
        vector< Answer > RL = Merge_TopK(hittuples, left, mid, topk);
        vector< Answer > RR = Merge_TopK(hittuples, mid+1, right, topk);
        R = CombineResult(RL,RR, topk);
    }
    else{
        vector<bool> key1(1,true);
        double value1 = log10(hittuples[left]);
        vector<bool> key2(1,false);
        double value2 = log10(1-hittuples[left]);
        
        if (value1 < value2 ){
            double tmp = value1; value1 = value2; value2 = tmp;
            key1[0] = false; key2[0] = true;
        }
        
        Answer ans1, ans2;
        ans1.key = key1; ans1.value = value1;
        ans2.key = key2; ans2.value = value2;
        
        
        if ( R.size() < topk ){
            R.push_back(ans1);
        }
        
        if ( R.size() < topk  ||  (fabs(value2-value1) < 0.00000001) ){
            R.push_back(ans2);
        }
        
    }
    
    return R;
}
