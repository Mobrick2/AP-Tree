  //
//  main.cpp
//  pdb
//
//  Created by BoZhang on 2/1/18.
//  Copyright © 2018 BoZhang. All rights reserved.
//

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>     /* atof */
#include <queue>
#include "sha.hpp"
#include <math.h>       /* log10 */
#include <stack>
#include "rsa.hpp"
#include "groupfilter.h"
#include "topk.hpp"
#include <numeric>
#include <random>
#include <float.h>
#include <stdio.h>

#include <queue>


//#include <NTL/ZZ.h>
//#include <NTL/vector.h>
//#include <NTL/ZZ_pXFactoring.h>

using namespace std;
//using namespace NTL;

double seedproServer = 0 ;
double seedproLocal = 0 ;
vector<double> probServer;
double seedtime = 0;
string mutplesize = "" ;

template <typename Map>
bool map_compare (Map const &lhs, Map const &rhs) {
    // No predicate needed because there is operator== for pairs already.
    int count = 0 ;
    //cout << "lhs size" << lhs.size() << " " << "rhs size:" << rhs.size() << endl;
    
    vector< vector<bool> > a;
    vector< vector<bool> > b;
    vector< vector<bool> > c(lhs.size(), vector<bool>(lhs.begin()->first.size()));
    auto jtm = rhs.begin();
    for ( auto itm = lhs.begin() ; itm != lhs.end(); ++ itm){
        
        if ( jtm == rhs.end() ) break;
        a.push_back(itm->first);
        b.push_back(jtm->first);
        jtm++;
    }
    
    sort(a.begin(),a.end());
    sort(b.begin(),b.end());
    auto it = set_intersection (a.begin(), a.end() , b.begin(), b.end(), c.begin());
    c.resize(it-c.begin());
    count = int ( c.size());
    
    if ( count != lhs.size() ) {
        cout << double(count) / double(rhs.size()) << endl;
        return false;
    }
    
    cout << 1 << endl;
    return true;
   // return lhs.size() == rhs.size()
   // && std::equal(lhs.begin(), lhs.end(),
     //             rhs.begin());
    
    
}

void readData(string pathtofile, vector<double> &tuples){
    
    ifstream infile (pathtofile);
    
    while( infile ){
        string s;
        if (!getline( infile, s )) break;
        
        istringstream ss( s );
        vector <string> record;
        
        while (ss)
        {
            string s;
            if (!getline( ss, s, ',' )) break;
            record.push_back( s );
        }
        
        tuples.push_back( stod(record[2]) );
    }
}

void readData3(string pathtofile, vector<double> &tuples){
    
    ifstream infile (pathtofile);
    
    while( infile ){
        string s;
        if (!getline( infile, s )) break;
        
        istringstream ss( s );
        vector <string> record;
        
        while (ss)
        {
            string s;
            if (!getline( ss, s, ',' )) break;
            record.push_back( s );
        }
        
        if ( record[3] != "")
            tuples.push_back( stod(record[3]) );
    }
}

void readData2(string pathtofile, vector<double> &tuples, vector<double> & prob){
    
    ifstream infile (pathtofile);
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,0.1);
    multimap<double , double> tmp;
    int i = 0 ;
    while( infile ){
        string s;
        if (!getline( infile, s )) break;
        
        istringstream ss( s );
        vector <string> record;
        
        while (ss)
        {
            string s;
            if (!getline( ss, s, ',' )) break;
            record.push_back( s );
        }
        
        if( record[3] == "" ) continue;
        
        double r= stod(record[3]);
        double p;
        
        
        
        if( record[6] == "R/V" ){
            p = distribution(generator)+0.8;
        }
        else if( record[6] == "VIS" ){
            p = distribution(generator)+0.7;
        }
        else if( record[6] == "RAD" ){
            p = distribution(generator)+0.6;
        }
        else if( record[6] == "SAT-LOW" ){
            p = distribution(generator)+0.5;
        }
        else if( record[6] == "SAT-MED" ){
            p = distribution(generator)+0.4;
        }
        else{
            p = distribution(generator)+0.3;
        }
        
        tmp.insert(make_pair(r,p));
        i++;
        
    }
    
    for (auto itm : tmp){
        tuples.push_back(itm.first);
        prob.push_back(itm.second);
    }
}

void readData3(string pathtofile, vector<double> &tuples, vector<double> & prob){
    
    ifstream infile (pathtofile);
    default_random_engine generator;
    normal_distribution<double> distribution(0.1,0.1);
    multimap<double , double> tmp;
    
    while( infile ){
        string s;
        if (!getline( infile, s )) break;
        
        istringstream ss( s );
        vector <string> record;
        
        while (ss)
        {
            string s;
            if (!getline( ss, s, ',' )) break;
            record.push_back( s );
        }
        
        if ( record[3] == "" ) continue;
        double r= stod(record[3]);
        double p;
        
        if( record[7] == "R/V" ){
            p = distribution(generator)+0.8;
        }
        else if( record[7] == "VIS" ){
            p = distribution(generator)+0.7;
        }
        else if( record[7] == "RAD" ){
            p = distribution(generator)+0.6;
        }
        else if( record[7] == "SAT-LOW" ){
            p = distribution(generator)+0.5;
        }
        else if( record[7] == "SAT-MED" ){
            p = distribution(generator)+0.4;
        }
        else{
            p = distribution(generator)+0.3;
        }
        
        tmp.insert(make_pair(r,p));
    }
    
    for (auto itm : tmp){
        tuples.push_back(itm.first);
        prob.push_back(itm.second);
    }
}

void assignUniformProb(vector<double> &prob, int n){
    
    double AOD;
    srand((unsigned)time(NULL));
    for (int i = 0; i < n ; ++i){
        AOD=((double) rand() / (RAND_MAX)) ;
        prob.push_back(AOD);
    }
}

void assignNormalProb(vector<double> &prob, int n, double m, double stdev ){
    
    unsigned seed = (unsigned)time(NULL);
    //std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    normal_distribution<double> distribution(m,stdev);
    for ( int i = 0 ; i < n ; ++ i){
        prob.push_back(distribution(generator));
    }
    
}


struct APBTree {
    
    //double p;
    int lindex;
    int rindex;
    
    APBTree * left;
    APBTree * right;
    string h;
    string childh;
    double accV;
    mpz_class tuplesig;
};


APBTree * createRoot(vector <double> &dataset, vector<double> &prob, vector<mpz_class> &sigs, int left, int right){
    
    APBTree * parent = new APBTree();
    parent->lindex = left;
    parent->rindex = right;
    stack<APBTree * > st;
    st.push(parent);
    
    while( !st.empty() ){
        
        APBTree * root = st.top();
        
        //it's a leaf node
        if ( root->lindex == root->rindex){

            string digest = "";
            double accValue = 0 ;
            mpz_class tuplestr;
            string cdigest = "";
            
            //min max values
            digest += to_string(dataset[root->lindex]);
            digest += to_string(dataset[root->rindex]);
            
            //acc value
            accValue += log10(prob[root->lindex]);
            root->accV = accValue;
            digest += to_string(accValue);
            
            //tuple hash
            tuplestr = sigs[root->lindex];
            root->tuplesig = tuplestr;
            
            stringstream buffer;
            buffer << tuplestr;
            digest += buffer.str();
            
            
            clock_t begin = clock();
            root->h = sha224(digest);
            clock_t end = clock();
            double ch = double(end - begin) / CLOCKS_PER_SEC;
            cout << ch << endl;
            
            root->h = sha224(digest);
            root->childh = "";
            st.pop();
        }
        
        else{
            
            string digest = "";
            digest += to_string(dataset[root->lindex]);
            digest += to_string(dataset[root->rindex]);
            
            if (root->left != NULL && root->right == NULL){
                
                root->accV = root->left->accV;
                digest += to_string(root->accV);
                
                root->tuplesig = root->left->tuplesig;
                stringstream buffer;
                buffer << root->tuplesig;
                digest += buffer.str();
                
                root->h = sha224(digest);
                
                string cdigest = root->left->h;
                cdigest = sha224(cdigest);
                root->childh = cdigest;
                
                root->h = sha224(root->h+root->childh);
                st.pop();
                
            }
            else if (root->left == NULL && root->right != NULL){
                
                root->accV = root->right->accV;
                digest += to_string(root->accV);
                
                root->tuplesig = root->right->tuplesig;
                stringstream buffer;
                buffer << root->tuplesig;
                digest += buffer.str();
                
                root->h = sha224(digest);
                
                string cdigest = root->right->h;
                cdigest = sha224(cdigest);
                root->childh = cdigest;
                
                root->h = sha224(root->h+root->childh);
                st.pop();
                
            }
            else if (root->left != NULL && root->right != NULL){
                
                root->accV += (root->left->accV+root->right->accV);
                digest += to_string(root->accV);
                
                root->tuplesig = RSA::condenseTwo(root->left->tuplesig,root->right->tuplesig);
                stringstream buffer;
                buffer << root->tuplesig;
                digest += buffer.str();
                
                root->h = sha224(digest);
                
                string cdigest = (root->left->h + root->right->h);
                cdigest = sha224(cdigest);
                root->childh = cdigest;
                
                root->h = sha224(root->h+root->childh);
                st.pop();
            }
            
            else{
                int m = (root->rindex-root->lindex)/2 + root->lindex;
                APBTree * leftnode = new APBTree();
                leftnode->lindex = root->lindex;
                leftnode->rindex = m;
                
                APBTree * rightnode = new APBTree();
                rightnode->lindex = m+1;
                rightnode->rindex = root->rindex;
                
                
                root->left = leftnode;
                root->right = rightnode;
                
                st.push(leftnode);
                st.push(rightnode);
            }
        }
    }
    
    return parent;
    
    
}

APBTree * createAPB (vector <double> &dataset, vector<double> &prob, vector<mpz_class> &sigs){
    
    APBTree * root;
    int n = int (dataset.size());
    
    root = createRoot(dataset, prob, sigs, 0, n-1);
    return root;
}

double selectionQuery(APBTree * root, vector <double> & dataset, double lower, double upper){
    
    queue<APBTree *> apb_q;
    apb_q.push(root);
    
    int MTcount = 0 ;
    int CTcount = 0 ;
    int MFcount = 0 ;
    int MHcount = 0 ;
    int MF = 0;
    int MH = 0 ;
    string vo = "";
    
    while( !apb_q.empty()){
        
        APBTree * node = apb_q.front();
        apb_q.pop();
        
        if ( node == NULL) continue;
        
        
        if ( node->lindex == node->rindex &&  dataset[node->lindex] >= lower &&
            dataset[node->lindex] <= upper){
            MTcount++;
            
            string tmp = "";
            tmp += "("; tmp += to_string(dataset[node->lindex]); tmp += ",";
            tmp += to_string(dataset[node->lindex]); tmp += ",";
            tmp += to_string(node->accV); tmp += ",";
            tmp += ")";
            vo += tmp;
            //mutplesize += tmp;
            
            continue;
        }
        
        if ( node->lindex == node->rindex &&  (dataset[node->lindex] < lower ||
                                               dataset[node->lindex] > upper ) ){
            CTcount++;
            continue;
        }
        
        if ( dataset[node->lindex]  >= lower && dataset[node->rindex] <= upper ){
            MHcount += (node->rindex-node->lindex+1);
            MH++;
            //MHcount++;
            
            string tmp = "<";
            for ( int i = node->lindex ; i <= node->rindex; ++ i ){
                tmp += to_string(dataset[i]); tmp += ",";
                tmp += to_string(dataset[i]); tmp += ",";
            }
            tmp += ">";
            vo += tmp;
            //mutplesize += tmp;
            
            continue;
        }
        
        if ( dataset[node->lindex]  > upper || dataset[node->rindex] < lower ){
            MFcount += (node->rindex-node->lindex+1);
            MF++;
            //MFcount++;
            continue;
        }
        
        if ( dataset[node->lindex] <= upper && dataset[node->rindex] >= lower ){
            
            apb_q.push(node->left);
            apb_q.push(node->right);
        }
        
    }
    
    //cout << MTcount+MHcount << endl;
    //cout << MTcount << "\t" << CTcount << "\t" << MF << "\t" << MH << endl;
    return  double (MTcount+MHcount) / double (MTcount+CTcount+MFcount+MHcount);
    
    
    //cout << double(MTcount+MHcount) / double (dataset.size() ) << "\t";
    
}


string selectionQueryVO(APBTree * node, vector <double> &dataset, vector <double> &prob, double lower, double upper){
    
    //ofstream myfile;
    //myfile.open ("/Users/Morbrick/Desktop/VO");
    string vo = "";
    
    if ( node == NULL) return "";
    
    //Match tuple
    if ( node->lindex == node->rindex &&  dataset[node->lindex] >= lower &&
        dataset[node->rindex] <= upper){
        
        vo += "("; vo += to_string(dataset[node->lindex]); vo += ",";
        vo += to_string(prob[node->lindex]);
        vo+= "," ; vo += to_string(node->accV) ; vo += ")";
        
        seedproServer += node->accV;
        probServer.push_back(prob[node->lindex]);
        //cout << "server :\n" << node->accV << endl;
        
        //cout << "(" << to_string(dataset[node->lindex]) << "," << to_string(prob[node->lindex]) << ")";
        return vo;
    }
    
    //Candidate tuple
    if ( node->lindex == node->rindex &&  (dataset[node->lindex] < lower ||
                                           dataset[node->lindex] > upper ) ){
        //cout << "(" << to_string(dataset[node->lindex]) << "," << to_string(prob[node->lindex]) << ")";
        vo += "("; vo += to_string(dataset[node->lindex]); vo += ",";
        vo += to_string(prob[node->lindex]);
        vo+= "," ; vo += to_string(node->accV) ; vo += ")";
        return vo;
    }
    
    //Maximum hit tuple [min,max,tuplehash, accV, lefthash, righthash]
    if ( lower <= dataset[node->lindex]  && dataset[node->rindex] <= upper ){
        
        
        vo += "["; vo += to_string(dataset[node->lindex]); vo += ","; vo += to_string(dataset[node->rindex]);
        vo += ","; vo += to_string(node->accV); vo += ",<" ;
        
        //Those are the query result;
        for ( int i = node->lindex; i <= node->rindex; ++ i ){
            vo += to_string(dataset[i]); vo += ",";
            vo += to_string(prob[i]); vo += ",";
        }
        
        vo += ">,"; vo += node->childh;
        vo += "]";
        //cout << "[" << to_string(dataset[node->lindex]) << "," << to_string(dataset[node->rindex]) << "," <<
        //node->accV << "," << node->tupleh << "," << node->childh << "]";
        
        seedproServer += node->accV;
        for ( int i = node->lindex ; i <= node->rindex; ++ i ) probServer.push_back(prob[i]);
        //cout << "node->accV server:\n" << node->accV << endl;
       
        return vo;
    }
    
    //Maximum Negative hit
    if ( dataset[node->lindex]  > upper || dataset[node->rindex] < lower ){
        
        vo += "["; vo += to_string(dataset[node->lindex]); vo += ","; vo += to_string(dataset[node->rindex]);
        vo += ","; vo += to_string(node->accV); vo += ",";
        
        stringstream buffer;
        buffer << node->tuplesig;
        vo += buffer.str();
        
        vo += "," ; vo += node->childh;
        vo += "]";
        
        return vo;
    }
    
    //cout << "(" << node->lindex << "," << node->rindex << ")" << endl;
    
    //Candidate Node
    if ( (dataset[node->lindex] < lower && lower <= dataset[node->rindex])
           || (upper < dataset[node->rindex] && upper >= dataset[node->lindex])  ){
        vo += "{";
        vo += to_string(node->accV) ; vo += ",";
        vo += selectionQueryVO(node->left,dataset,prob,lower,upper);
        vo += selectionQueryVO(node->right,dataset,prob,lower,upper);
        vo += "}";
    }
    
    return vo;
}

void pushDigest(string digest, stack<char> &st){
    
    for ( int i = int (digest.size()-1); i >= 0 ; --i){
        st.push(digest[i]);
    }
}

int getTupleLength(stack<char> & st, int down, int up){
    
    string value = "";
    double prob = 0, v = 0 ;
    string accV = "" ;
    int count = 0 ;
    int x = 0 ;
    while( !st.empty() ){
        
        char c = st.top(); st.pop();
        x ++;
        if ( c != ')' && c!= '(' && c != ',') value += c;
        
        if ( c == ',' && count == 0){
            v = stod(value);
            value = "";
            count++;
        }
        
        else if ( c == ',' && count == 1){
            prob = stod(value);
            value = "";
            count++;
        }
        
        else if (c == ',' && count == 2){
            accV = value;
            value = "";
            count++;
        }
        
        if ( c == ')' ) break;
    }
    
    
    if ( accV == "" ) accV = value;
    
    if ( down <= v && v <= up){
        return x;
    }
    
    return 0 ;
}


void verifyTuple(stack<char> & st, int down, int up){
    
    string value = "";
    double prob = 0, v = 0 ;
    string accV = "" ;
    int count = 0 ;
    while( !st.empty() ){
        
        char c = st.top(); st.pop();
        if ( c != ')' && c!= '(' && c != ',') value += c;
        
        if ( c == ',' && count == 0){
            v = stod(value);
            value = "";
            count++;
        }
        
        else if ( c == ',' && count == 1){
            prob = stod(value);
            value = "";
            count++;
        }
        
        else if (c == ',' && count == 2){
            accV = value;
            value = "";
            count++;
        }

        if ( c == ')' ) break;
    }
    
    
    if ( accV == "" ) accV = value;
    
    if ( down <= v && v <= up){
        clock_t begin = clock();
        seedproLocal += stod(accV);
        clock_t end = clock();
        seedtime += double(end - begin) / CLOCKS_PER_SEC;
        //cout << "local : node->accV :" << stod(accV) << endl;
    }
    
    string digest = "";
    //min max value are the same
    digest += to_string(v);
    digest += to_string(v);
    //acc Value just one;
    
    
    digest += accV;
    //cout << "acc:\n" << accV << endl;
    
    string tuplehash = to_string(v)+to_string(prob);
    mpz_class tuplesig = RSA::sign(tuplehash);
    
    stringstream buffer;
    buffer << tuplesig;
    digest += buffer.str();
    digest = sha224(digest);
    
   
    string otherinfo  = "";
    otherinfo += to_string(v); otherinfo+= ",";
    otherinfo += to_string(v); otherinfo+= ",";
    otherinfo += buffer.str(); otherinfo += ",";
    otherinfo += digest; otherinfo += ",";
    
    pushDigest(otherinfo,st);
}

void parseInternalNode(string &segment, double &minv, double &maxv, double & accV,
                       string &tuplehash, string & childhash ){
    
    string value = "";
    int count = 0;
    
    for (auto c : segment ) {
        
        if ( c!= ',' && c!= '[' && c!= ']' ){
            value += c;
        }
        
        if ( c == ',' && count == 0 ){
            minv = stod(value);
            value = "";
            count ++;
        }
        else if (c == ',' && count == 1 ){
            maxv = stod(value);
            value = "";
            count ++;
        }
        else if (c == ',' && count == 2 ){
            accV = stod(value);
            value = "";
            count ++;
        }
        else if (c == ',' && count == 3 ){
            tuplehash = value;
            value = "";
            count ++;
        }
    }
    
    childhash = value;
}

void parseInternal(string &segment, double &minv, double &maxv, double & accV,
                   string &tuplehash, string & childhash ){
    
    string value = "";
    int count = 0;
    mpz_class childsig = 1;
    string lefth = "";
    string righth = "";
    
    
    
    for (auto c : segment ) {
        
        if ( c != '{' && c!= ',' && c!= '}' ){
            value += c;
        }
        
        if ( c == ',' && count == 0 ){
            accV = stod(value);
            value = "";
            count ++;
        }
        else if (c == ',' && count == 1 ){
            minv = stod(value);
            value = "";
            count ++;
        }
        else if (c == ',' && count == 2 ){
            maxv = stod(value);
            value = "";
            count ++;
        }
        else if (c == ',' && count == 3 ){
            childsig = value;
            value = "";
            count ++;
        }
        else if (c == ',' && count == 4 ){
            lefth = value;
            value = "";
            count ++;
        }
        else if (c == ',' && count == 5 ){
            minv = min(minv,stod(value));
            value = "";
            count ++;
        }
        else if (c == ',' && count == 6 ){
            maxv = max(maxv,stod(value));
            value = "";
            count ++;
        }
        else if (c == ',' && count == 7 ){
            
            mpz_class v = 1;
            v = value;
            childsig = (childsig * v) % RSA::N;
            value = "";
            count ++;
        }
        else if (c == ',' && count == 8 ){
            righth = value;
            value = "";
            count ++;
        }
    }
    
    
    stringstream buffer;
    buffer << childsig;
    tuplehash = buffer.str();
    
    childhash = lefth + righth;
}

string combinehash(double &minv, double &maxv, double &accV, string &tuplehash,
                   string &childhash){
    
    string digest = "";
    digest += to_string(minv); digest += to_string(maxv);
    digest += to_string(accV); digest += tuplehash;
    digest = sha224(digest);
    
    //cout << digest << endl;
    digest += childhash;
    
    digest = sha224(digest);
    //cout << digest << endl;
    
    return digest;
}

void verifyNode(stack<char> & st, int down, int up){
    
    //[min,max,acc, tuplehash, childhash]
    double minv = 0 , maxv = 0 , accV = 0;
    string segment = "";
    string tuplehash = "";
    string childhash = "";
    
    while( !st.empty() ){
        
        char c = st.top();
        st.pop();
        segment += c;
        if ( c == ']' ) break;
    }
    
    parseInternalNode(segment, minv, maxv, accV, tuplehash, childhash);
    if ( minv >= down && up >= maxv ){
        clock_t begin = clock();
        seedproLocal += accV;
        clock_t end = clock();
        seedtime += double(end - begin) / CLOCKS_PER_SEC;
        //cout << minv << "," << maxv << endl;
        //cout << "local:\n" << accV << endl;
    }
    
    if ( !((minv >= down && up >= maxv) || (up < minv || down > maxv)) ) {
        cout << minv  << "," << maxv << endl;
        cout << down  << "," << up << endl;
        cout << "wrong vo from maximum node" << endl;
    }
    
    string digest = combinehash(minv, maxv, accV, tuplehash,childhash);
    
    string otherinfo = "";
    otherinfo += to_string(minv); otherinfo  += ",";
    otherinfo += to_string(maxv); otherinfo  += ",";
    otherinfo += tuplehash; otherinfo  += ","; otherinfo += digest;
    otherinfo += ",";
    pushDigest(otherinfo,st);
    
}

void verifyInternal(stack<char> & st, int down, int up){
    
    double minv = 0 , maxv = 0 , accV = 0;
    string segment = "";
    string tuplehash = "";
    string childhash = "";
    
    while( !st.empty() ){
        
        char c = st.top();
        st.pop();
        segment += c;
        if ( c == '}' ) break;
    }
    
    parseInternal(segment, minv, maxv, accV, tuplehash, childhash);
    
    if ( (down <= minv && maxv <= up) || (up < minv || down > maxv) )
    {
        cout << minv << "," << maxv << endl;
        cout << down << "," << up << endl;
        cout << "wrong vo from candiate node" << endl;
        
    }
    
    childhash = sha224(childhash);
    string digest = combinehash(minv, maxv, accV, tuplehash, childhash);
    
    string otherinfo = "";
    otherinfo += to_string(minv); otherinfo  += ",";
    otherinfo += to_string(maxv); otherinfo  += ",";
    otherinfo += tuplehash; otherinfo  += ","; otherinfo += digest;
    otherinfo += ",";
    
    
    
    pushDigest(otherinfo,st);
    
}

void printstack( stack<char> st){
    
    while( !st.empty()){
        cout << st.top();
        st.pop();
    }
    cout << endl;
}

void verifyMaximumHitTuple(stack<char> &st){
    
    vector<string> value;
    vector<string> prob;
    
    string tmp = "";
    int count =  0;
    while ( !st.empty()){
        
        char c = st.top();
        st.pop();
        
        if ( c != '<' && c!= '>' && c != ','){
            tmp += c;
        }
        
        if ( c == ',' && ((count % 2) == 0)){
            value.push_back(tmp);
            tmp = "";
            count ++;
        }
        
        else if ( c == ',' && ((count % 2) == 1)){
            prob.push_back(tmp);
            tmp = "";
            count ++;
        }
        
        if ( c == '>'){
            break;
        }
    }
    
    vector< mpz_class> sigs;
    for ( int i = 0 ; i < value.size(); ++ i){
        sigs.push_back(RSA::sign(value[i]+prob[i]));
    }
    
    mpz_class childsig = RSA::condense(sigs);
    stringstream buffer;
    buffer << childsig;
    
    string tuplehash = "";
    tuplehash += buffer.str();
    
    //cout << "tuplehash\n:" << tuplehash << endl;
    
    pushDigest(tuplehash,st);
    
}

string verifyVP(string vo, int down, int up){
    
    stack<char> st;
    int n = int (vo.size());
    
    for ( int i = n-1 ; i >= 0 ; -- i){
        
        st.push(vo[i]);
        //Match tuple or candidate tuple
        if ( vo[i] == '(' ){
            verifyTuple(st,down,up);
        }
        
        if ( vo[i] == '['){
            verifyNode(st,down,up);
        }
        
        if ( vo[i] == '{' ){
            verifyInternal(st,down,up);
        }
        
        if ( vo[i] == '<'){
            verifyMaximumHitTuple(st);
        }
        
    }
    
    string roothash = "";
    while( !st.empty()){
        roothash += st.top();
        st.pop();
    }
    
    return roothash;
    
}

void createAccTuple(vector<double> &dataset,  vector<double> &prob, vector<mpz_class> &sigs){
    
    int n = int (dataset.size());
    string vp = "";
    for ( int i = 0 ; i < n ; ++ i){
        vp += to_string(dataset[i]);
        vp += to_string(prob[i]);
        
        //cout << "vp: " << vp << endl;
        sigs.push_back(RSA::sign(vp));
        vp = "";
    }
    
}

void printroot(string & roothash){
    
    stack<char> st;
    for ( int i = int (roothash.size()-2); i>=0; --i){
        
        if ( roothash[i] == ',') break;
        st.push(roothash[i]);
    }
    
    while( !st.empty() ){
        cout << st.top() ;
        st.pop();
    }
    cout << "?=\n";
    
}

void addData(vector<double> & dataset, int time){
    
    int n = int (dataset.size());
    vector<double> tmp_dataset = dataset;
    
    for ( int i = 0 ; i < time; ++ i){
        for ( int j = 0 ; j < n ; ++ j){
            dataset.push_back(tmp_dataset[j]);
        }
    }
}


int calculate_resultsize( string & vo, int down, int up){
    
    int count = 0;
    
    stack<char> st;
    
    for ( auto itm : vo){
        
        st.push(itm);
        
        if ( itm == '>'){
            while( ! st.empty() ){
                
                char c = st.top();
                st.pop();
                count ++;
                if ( c == '<' ) break;
            }
        }
        
        if ( itm == ')'){
            count += getTupleLength(st, down, up);
        }
        
    }
    
    return count;
}


int minKey(vector<int> key, vector<bool> mstSet)
{
    // Initialize min value
    int min = INT_MAX, min_index = -1;
    int V = int (mstSet.size());
    
    for (int v = 0; v < V; v++){
        if (mstSet[v] == false && key[v] < min){
            min = key[v] ; min_index = v;
        }
    }
    
    return min_index;
}

// A utility function to print the constructed MST stored in parent[]
void printMST(vector<int> parent, int n, vector< vector<int> > & graph )
{
    int V = int (graph.size());
    for (int i = 1; i < V; i++){
        cout << parent[i] << "-" << i << "   " << graph[i][parent[i]] << endl;
    }
}

// Function to construct and print MST for a graph represented using adjacency
// matrix representation

void primMST (vector< vector< int > > &graph, VPResult & test, vector<double> & primtime,
              vector<double> & findRepresentative, int k)
{
    clock_t begin = clock();
    
    int V = int (graph.size());
    vector<int> parent(V,0);; // Array to store constructed MST
    vector<int> key(V,INT_MAX);   // Key values used to pick minimum weight edge in cut
    vector<bool> mstSet(V, false);  // To represent set of vertices not yet included in MST
    //cout << "V" << V << endl;
    
    
    // Always include first 1st vertex in MST.
    key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
    parent[0] = -1; // First node is always root of MST
    
    // The MST will have V vertices
    for (int count = 0; count < V-1; count++)
    {
        // Pick the minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet);
        
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
        
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < V; v++){
            
            // graph[u][v] is non zero only for adjacent vertices of m
            // mstSet[v] is false for vertices not yet included in MST
            // Update the key only if graph[u][v] is smaller than key[v]
            if (graph[u][v] && mstSet[v] == false && graph[u][v] <  key[v]){
                parent[v]  = u; key[v] = graph[u][v];
            }
        }
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    primtime[k] += elapsed_secs;
    
    //for ( auto itm : parent ) cout << itm << " ";
    //cout << endl;
    
    // print the constructed MST
    
    vector< vector<bool> > topkresult;
    vector<double> prob;
    for ( auto itm : test.result ) {
        topkresult.push_back(itm.first);
        prob.push_back(itm.second);
    }
    
    int num = int ( topkresult[0].size() );
    //begin = clock();
    
    for (int i = 1 ; i < V; i++){
        
        int next = parent[i];
        vector<bool> connectkey = topkresult[next];
        
        //debug code
        vector<bool> cur_key = topkresult[i];
        
        for ( int j = 0 ; j < num; ++ j ){
            
            
            if ( cur_key[j] == 1 && connectkey[j] == 0 ){

                
                cur_key[j] = 0 ;
                
                auto it = test.result.find(cur_key);
                if ( it == test.result.end() ){
                    
                    begin = clock();
                    double newprob = prob[i]+ log10( (1-probServer[j])/probServer[j]);
                    test.insert_rp(cur_key, newprob);
                    test.represent.insert(make_pair(cur_key, newprob));
                    end = clock();
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    findRepresentative[k] += elapsed_secs;
                }
                
                
            }
            
            
            if ( cur_key[j] == 0 && connectkey[j] == 1 ){
                
                
                
                cur_key[j] = 1 ;
                auto it = test.result.find(cur_key);
                if ( it == test.result.end() ){
                    begin = clock();
                    double newprob = prob[i]+ log10( probServer[j] / (1-probServer[j]) );
                    
                    test.insert_rp(cur_key, newprob);
                    test.represent.insert(make_pair(cur_key, newprob));
                    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    findRepresentative[k] += elapsed_secs;
                }
            }
            
            /*
            if ( topkresult[i][j] == 1 && connectkey[j] == 0 ){
                vector<bool> tmp_key = topkresult[i];
                tmp_key[j] = 0 ;
                
                auto it = test.result.find(tmp_key);
                if ( it == test.result.end() ){
                    //begin = clock();
                    double newprob = prob[i]+ log10( (1-probServer[j])/probServer[j]);
                    //end = clock();
                    //elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    //findRepresentative[k] += elapsed_secs;
                    
                    test.insert_rp(tmp_key, newprob);
                    test.represent.insert(make_pair(tmp_key, newprob));
                    
                    
                    
                }
            }
            
            if ( topkresult[i][j] == 0 && connectkey[j] == 1 ){
                vector<bool> tmp_key = topkresult[i];
                tmp_key[j] = 1 ;
                
                auto it = test.result.find(tmp_key);
                if ( it == test.result.end() ){
                    //begin = clock();
                    double newprob = prob[i]+ log10( probServer[j] / (1-probServer[j]) );
                    //end = clock();
                    //elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                    //findRepresentative[k] += elapsed_secs;
                    
                    test.insert_rp(tmp_key, newprob);
                    test.represent.insert(make_pair(tmp_key, newprob));
                    
                }
            }
            */
            
        }
    }
    
    //end = clock();
    //elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //findRepresentative[k] += elapsed_secs;
    
    
    
    
}

int getDistance( vector<bool> a, vector<bool>  b){
    
    int m = int (a.size());
    //int count = 0 ;
    //int n1 = 0, n2 = 0 ;
    
    
    vector<int> location;
    int tmp_dis = 0 ;
    for ( int i = 0 ; i < m ; ++ i ){
        
        /*
        if ( a[i] == 1 ) n1 ++;
        if ( b[i] == 1 ) n2 ++;
        if ( a[i] == 1 && b[i] == 1 ) count++;
        */
        
        //Simulate the time to store the difference bit location
        if ( a[i] != b[i] ) {
            location.push_back(i);
            tmp_dis ++ ;
        }
        
    }
    
    
    //int dis = n1+n2-count-count;
    
    //if ( dis != tmp_dis ) cout << "cout 1" << endl;
    
    return tmp_dis;
}

void find_representative( int hitnum, VPResult & test,
    vector<double> &primtime, vector<double> &distancetime,
                           vector<double> & findRepresentative, int k ){
    
    clock_t begin = clock();
    map <vector<bool>, double>  result = test.result;
    int m = int (result.size());
    
    vector< vector< int > > graph(m, vector< int> (m,0) );
    map <vector<bool>, double> copy_result = result;
    
    int i = 0, j = 0 ;
    
    for (auto itm : result){
        j = 0 ;
        for ( auto jtm : copy_result ){
            
            if ( i != j ){
                graph[i][j] = getDistance(itm.first, jtm.first);
            }
            
            j++;
        }
        i++;
    }
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    distancetime[k] += elapsed_secs;
    
    primMST(graph,test, primtime, findRepresentative, k);
    
}

void assignRandom( vector<double> &ran){
    
    double AOD;
    srand((unsigned)time(NULL));
    for (int i = 0; i < 30000000 ; ++i){
        AOD = rand();
        ran.push_back(AOD);
    }
}

int main(int argc, const char * argv[]) {
    
    
    string inputfile = "/Users/Morbrick/Desktop/AuthPDB-exp/IIP2002-2015.csv";
    //string inputfile = "/Users/Morbrick/Desktop/AuthPDB-exp/uservisit";
    //cout << inputfile.size() << " " << inputfile.length() << endl;
    
    //read data from file
    vector<double> dataset;
    vector<double> prob;
    readData3(inputfile,dataset,prob); //1-254 0.004
    //readData(inputfile, dataset); //0-2359 0.00004
    
    
    /*
    //check the ratio for each unique value;
    map<double, int> mp;
    for ( auto itm : dataset ){
        if ( mp.find(itm) != mp.end() ){
            mp[itm] ++;
        }
        else {
            mp[itm] = 1;
        }
    }
    
    //vector<double> ratios;
    double s = 0 ;
    for ( auto itm : mp ){
        double tmp =  double(itm.second) / double(dataset.size()) ;
        s += tmp;
        cout << itm.first << "   " << s << "    " << tmp << endl;
    }
    */
    
    
    
    int time = 10;
    //addData(dataset, time);
    
    int base2 = 1400;
    
    //change the location of the hit set : 1-254
    for ( int i = 0 ; i <= 4200; i++){
        for ( int j = 1 ; j <= 30 ; ++ j ){
            dataset.push_back(2359+j);
        }
    }
    
    vector<double> dataset2;
    for ( int i = 0 ; i < dataset.size(); ++ i){
        if ( dataset[i] < 30+base2 || dataset[i] > 60+base2 ){
            dataset2.push_back(dataset[i]);
        }
    }
    for ( int i = 30+base2 ; i <= 60+base2 ; ++ i ){
        dataset2.push_back(i);
    }
    dataset.clear();
    dataset = dataset2;
    
    
    sort(dataset.begin(), dataset.end());
    //assignNormalProb(prob, int (dataset.size()), 0.5, 0.2);
    assignUniformProb(prob, int (dataset.size()) );
    
    
    RSA::initialize();
    vector<mpz_class> sigs;
    cout << "create tuple hash....." << endl;
    clock_t begin = clock();
    createAccTuple(dataset, prob, sigs);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << elapsed_secs << endl;
    
    cout << "create root hash....." << endl;
    begin = clock();
    APBTree * root = createAPB(dataset, prob, sigs);  //9.67383s
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << elapsed_secs << endl;
    
    cout << "selection ratio and query time\n";
    vector<double> ratio(300,0.0);
    vector<double> qtime(300,0.0);
    //vector<double> resultsize(100,0.0);
    
    time = 3;
    //int range = 60;
    //int interval = 2;
    //int range = 1500;
    //int interval = 50;
    //int range = 26;
    //int interval = 10;
    int k = 0 ;
    double base = 30+base2;
    
    
    for ( int j = 0 ; j < time; j ++){
        
        k = 0 ;
        
        //test for TP-results
        for ( int i = 46+base2 ; i <= 56+base2 ; i = i + 2 ){
            //mutplesize = "";
            begin = clock();
            ratio[k] = selectionQuery(root, dataset,
                                      base, i); //change base here
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            qtime[k] += elapsed_secs;
            //resultsize[k] += ( double ( mutplesize.size())/1024/1024 );
            k++;
        }
    }
    
    for ( int j = 0 ; j < k; j ++){
        qtime[j] = qtime[j] / time;
        //resultsize[j] = resultsize[j] / time;
    }
    
    
    /*
    for ( int j = 0 ; j < k; j ++){
        cout << ratio[j]  << "\t" << qtime[j] << endl;
    }
    */
   
    /*
    cout << "VO C Time, VO Size, Verify Time, V Seed Time " << endl;
    vector<double> vosize(100,0.0);
    vector<double> votime(100,0.0);
    vector<double> verifytime(100,0.0);
    vector<double> verifytime2(100,0.0);
    vector<double> seedTime(100,0.0);
    vector<double> resultsize(100,0.0);
    vector<double> bruteforceTime(100,0.0);
    //Brute force to get optimal answers
    for ( int j = 0 ; j < time ; ++ j){
    
        k = 0 ;
        seedtime = 0 ;
        
        //test for TP result
        
        for ( int i = 46+base2 ; i <= 56+base2 ; i = i + 2 ){
            
            probServer.clear();
            begin = clock();
            string vo = selectionQueryVO(root,dataset,prob,
                                         base, i); //change base here
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            votime[k] += elapsed_secs;
            vosize[k] += double(vo.size())/1024/1024;
            resultsize[k] += (double(calculate_resultsize(vo,base,i)) / 1024 / 1024); //change base here
        
     
            begin = clock();
            verifyVP(vo, base, i); //change base here
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            verifytime[k] += elapsed_secs ;
            seedTime[k] += seedtime;
        
            
            int m = int (probServer.size());
            //cout << m << endl;
            
            VPResult test(m);
            begin = clock();
            test.create_allrp(probServer);
            end = clock();
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            bruteforceTime[k] += elapsed_secs;
            
            
            begin = clock();
            test.checkallgroupresult2(probServer);
            end = clock();
            //cout << test.result.size() << endl;
            
            elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            verifytime[k] += elapsed_secs ;
            verifytime2[k] += elapsed_secs;
     
            k++;
        }
        
        
    }
    
    
    for ( int j = 0 ; j < k; j++ ){
        vosize[j] = vosize[j] / time;
        resultsize[j] = resultsize[j] / time;
        votime[j] = votime[j] / time;
        //verifytime[j] = verifytime[j] / time;
        //seedTime[j] = seedTime[j] / time;
        //bruteforceTime[j] = bruteforceTime[j] / time;
        //verifytime2[j] = verifytime2[j] / time;
    }
    
    
    //[0]    double    1.2051820000000002
    
    for ( int j = 0 ; j < k; j ++){
        cout << verifytime[j] <<  "\t" << verifytime2[j] << "\t" << bruteforceTime[j] << "\t"
        << seedTime[j] << "\t" << vosize[j] << "\t" << resultsize[j] << "\t" << qtime[j] << "\t" << votime[j] << endl;
        //cout << vosize[j] << "\t" << resultsize[j] << "\t" << qtime[j] << "\t" << votime[j]  << "\t" << ratio[j]  << endl;
    }
    
    
    
    //MC simulation results and this is included in the verification time
    */
    
    
    //vector<double> prob;
    assignNormalProb(prob, 5000, 0.5, 0.2);
    //assignUniformProb(prob, 5000);
    //prob[0] = 0.5; prob[1] = 0.6; prob[2] = 0.7; prob[3] = 0.8;
    
    
    vector<int> topk = {300};
    //vector<int> topk = {800};
    //vector<int> topk2 = {40,120,200,280,360};
    //vector<int> topk2 = {40};
    //int time = 20;
    
    //for ( auto itm : prob ) cout << itm << endl;
    //vector<double> ran;
    //assignRandom(ran);
    
    
    
    for ( int i = 14 ; i <= 26;  i = i + 2 ){
        
        vector<double> bruteforce_time(20,0.0);
        vector<double> verify_time(20,0.0);
        vector<double> verify_time2(20,0.0);
        vector<double> vo_time(20,0.0);
        
        vector<double> primtime(20,0.0);
        vector<double> distancetime(20,0.0);
        vector<double> findRepresentative(20,0.0);
        vector<double> vo_size(20,0.0);
        vector<double> bruteforce_less(20,0.0);
        
        probServer.clear();
        
        vector<bool> seed(i,1);
        double proball = 0;
        for ( int j = 0 ; j < i ; ++ j) {
            probServer.push_back(prob[j]);
            proball += log10(prob[j]);
        }
     
        
        
        
        /* // For the test of real ranking of the lower bound
        vector <double> tmp;
        brute_force2(probServer, tmp);
        sort(tmp.begin(),tmp.end());
        */
        
        //cout << "size of hit set :" << i << endl;
        /*
        clock_t begin = clock();
        vector< Answer > result  = Merge_TopK(probServer,0,i-1, topk[0]);
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        bruteforce_time[0] += elapsed_secs;
        */
    
        for ( int k = 0 ; k < topk.size() ; ++ k ){
  
            /*
            //cout << "top-" << topk[k] << ":" << endl;
            VPResult test(i);
            clock_t begin = clock();
            mc_simulationvector(probServer,test,0.99, topk[k]);
            mc_simulationvector(probServer,test,0.85, topk[k]);
            clock_t end = clock();
            double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
            simulate_time[k] += elapsed_secs;
            */
            
            
            
            for ( int j = 0 ; j < time; ++ j ){
              
                //Calculate the real ranking of the lower bound
                /*
                double plb = getLowerbound(probServer, topk[k], 0);
                for ( int index = tmp.size()-1; index >= 0; --index ){
                    if ( fabs(tmp[index]-plb) < 0.0000001 ){
                        cout << tmp.size() - index << endl;
                        break;
                    }
                }
                 
                VPResult test2(i);
                //Two Bruteforce Methods
                clock_t begin = clock();
                brute_force(probServer, test2, topk[k]);
                deterministic(probServer,test2,topk[k]);
                clock_t end = clock();
                double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                bruteforce_time[k] += elapsed_secs;
                */
            
                
                int k1 = 0.2 * topk[k] ;
                //int k1 = topk2[k];
                
               
                clock_t begin = clock();
                vector< Answer > result1 = Merge_TopK(probServer,0,i-1, k1);
                clock_t end = clock();
                double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                bruteforce_less[k] += elapsed_secs;
                
                
                begin = clock();
                vector< Answer > result  = Merge_TopK(probServer,0,i-1, topk[k]);
                end = clock();
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                bruteforce_time[k] += elapsed_secs;
                
                
                VPResult test(i);
                for ( int l = 0 ; l < result.size(); ++ l ){
                    test.insert_rp(result[l].key, result[l].value);
                }
                
                
                //map_compare (test2.result,test.result);
                
                /* All anwer verification
                test2.create_allrp(probServer);
                begin = clock();
                test2.checkallgroupresult(probServer);
                end = clock();
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                verify_time[k] += elapsed_secs;
                
            
                begin = clock();
                int x = estimateK(topk[k],i);
                //int x = estimateK2(topk[k],i);
                VPResult test(i);
                deterministic(probServer,test,x);
                end = clock();
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                primtime[k] += elapsed_secs;
                realK[k] = x;
                */
                
                //continue to find the representativfe
    
                
                test.insert_rp(seed, proball);
                begin = clock();
                find_representative(i, test, primtime,distancetime, findRepresentative, k);
                end = clock();
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                vo_time[k] += elapsed_secs;
                
                //check the authenticity of the answer
                begin = clock();
                test.checkgroupresult(probServer);
                end = clock();
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                verify_time[k] += elapsed_secs;
                vo_size[k] += ( double (test.represent.size()) /1024/1024 * (1+i+8) );
            
      
                //check the validness of the K-th answer
                /*
                begin = clock();
                int countFlip = 0 ;
                int indexForresult = 0 ;
                for ( auto itm : test.result ){
                    
                    if ( indexForresult == 0 ){
                        for ( int l = 0 ; l < i ; ++ l ){
                            if ( itm.first[l] == true ){
                                countFlip ++;
                            }
                        }
                        if ( countFlip > x ) cout << "";
                        indexForresult ++ ;
                    }
                }
                end = clock();
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                verify_time2[k] += elapsed_secs;
                verify_time[k] += elapsed_secs;
                */

                
                /* //get TOP-K probability summation
                vector<double> tmp;
                for( auto itm : test2.result ){
                    tmp.push_back(itm.second);
                }
                sort(tmp.begin(),tmp.end());
                double sp = 0 ;
                for ( int y = int (tmp.size()-1) ; y >= 0 ; --y ){
                    sp += pow(10,tmp[y]);
                    x--;
                    if ( x < 0 ) break;
                }
                cout << "--" << sp << endl;
                */
                
                /*
                begin = clock();
                //map_compare (test2.result,test.result);
                end = clock();
                
                elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
                verify_time2[k] += elapsed_secs;
                */
               
            }
            
            /*
            cout << bruteforce_time[k] / time << "\t" <<  bruteforce_less[k] / time << "\t" <<
            (vo_time[k]) / time << "\t" << vo_size[k] / time << "\t" << verify_time[k] / time << endl;
            */
        }
        
        
        
        for ( int k = 0 ; k < topk.size(); ++ k){
            cout << bruteforce_time[k] / time << "\t" <<  bruteforce_less[k] / time << "\t" <<
            (vo_time[k]) / time << "\t" << vo_size[k] / time << "\t" << verify_time[k] / time << endl;
         
            //cout << (primtime[k] + distancetime[k] + findRepresentative[k] ) / time  << "\t" << vo_size[k] / time << "\t" << findRepresentative[k]/ time  << "\t " <<
            //(primtime[k] + distancetime[k] ) / time << endl;
         
        }
        
        
        
    }
    
    return 0;
}

