#include <iostream>
#include "networks.h"
#include <math.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <string>
#include <stdexcept>
#include <queue>
#include <ctime>
#include <algorithm>
#include <random>
#include <sstream>
#include <fstream>
#include <ostream>
#include <time.h>
#include <functional>
#include <cassert>
//#include <omp.h>

#define  P(A) cout << #A << ": " <<(A) << endl;

using namespace std;
typedef std::mersenne_twister_engine<uint32_t, 32, 351,
        175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253> mt11213b;

enum OS {
    LINUX, WIN, MAC
};

template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<T>());
    return result;
}

//global variables:
float k_total = 4;              //k_inter + k_intra = k_total
int const M = 1000;
int const COMM_SIZE = 1000;
int const LATTICE_LENGTH = 100;
int const LATTICE_SIZE = LATTICE_LENGTH * LATTICE_LENGTH;
long int N = COMM_SIZE * M;
double p_min = 0.02;           // alpha parameter, alpha = k_intra/k_inter
double p_max = 1;
double dp = 0.02;
int num_of__sim = 50;

//functions:
vector<double> build_alpha_vector(double min, double max, double da);

vector<double> build_p_vector(double min, double max, double dp);

vector<float> calculateKinterKintra(long double a);

template<class T>
void printVector(vector<T> &v);

template<class T>
T getVectorAvg(vector<T> &v);

template<class T>
vector<T> getVectorOfVectorsAverage(vector<vector<T> > &v);

void save_adj_list(OS os, vector< vector<int> >& adj);

int main() {
    vector<double> p_v;
    vector<int> avg__gcc_nodes;
    vector<int> avg_secondGCC_nodes;
    vector<int> avg_gcc_modules;
    int number_of_inter_links = M * (M - 1);

    double p;

    /*****************************************/
    cout << "-------starting--------" << endl;
    clock_t tStart = clock();

    P(k_total);
    P(M);
    P(N);
    p_v = build_p_vector(p_min, p_max, dp);
    //printVector(p_v);
    for (int p_idx = 0; p_idx < p_v.size(); p_idx++) {

        vector<int> gcc_nodes_storage;
        vector<int> secGcc_nodes_storage;
        vector<int> gcc_modules_storage;
        p = p_v[p_idx];

//        omp_set_num_threads(6);
#pragma omp parallel for
        for (int s = 0; s < num_of__sim; s++) {
            cout << "running simualtion " << s << " out of " << num_of__sim << endl;

            network myNetwork(ER);
            //vector<float> k_s = calculateKinterKintra(k_s[0]);
            myNetwork.build_er_network_with_k(4);  //i need to change this function
            myNetwork.assign_inter_links(number_of_inter_links);
            cout << "finished building communities and assigning inter-links, starting attack with" << endl;
            //save_adj_list(WIN, myNetwork.adj_list);
            //attack
            P(p);
            myNetwork.attack_random(p);
            //myNetwork.attack_highest_degree_nodes(p);
            //myNetwork.attack_highest_degree_modules(p);
            vector<int> results;
            results.resize(3);
            cout<<"testing connectivity"<<endl;
            results = myNetwork.testConnectivity();
            gcc_nodes_storage.push_back(results[0]);
            gcc_modules_storage.push_back(results[1]);
            secGcc_nodes_storage.push_back(results[2]);
            printVector(results);
        }
        //averaging for all simulations of current alpha
        cout << "storing average values" << endl;
        avg__gcc_nodes.push_back(getVectorAvg(gcc_nodes_storage));
        avg_gcc_modules.push_back(getVectorAvg(gcc_modules_storage));
        avg_secondGCC_nodes.push_back((getVectorAvg(secGcc_nodes_storage)));
    }
    cout << "saving data" << endl;
    ostringstream oss;
    oss << M << "_" << COMM_SIZE <<"_"<<number_of_inter_links<<"_num_inter_links"
            "_ER_modules_GCC_SGCC_RANDOM_ATTACK_"<<num_of__sim<<"_simulations.txt";
    string str = oss.str();
    ofstream data;
    string path = "G:\\dropbox\\Dropbox\\M.Sc\\results\\" + str;
    //mac:        "/Users/eitanasher/Desktop/" + str;
    //linux:       "/home/eitan/Dropbox/M.Sc/results/" + str;
    //windows:    "G:\\dropbox\\Dropbox\\M.Sc\\results" + str;
    cout << "path = " << path << endl;
    const char *c_path = path.c_str();
    data.open(c_path, ofstream::out);
    data << "p      nodes in giant      nodes in second    " << endl;

    for (int j = 0; j < p_v.size(); j++) {
        data << p_v[j] << "        " << avg__gcc_nodes[j] << "          " << avg_secondGCC_nodes[j] << endl;
    }
    data.close();
    cout << "data saved successfully" << endl;

    printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

    return 0;
}//end of main


vector<double> build_alpha_vector(double min, double max, double da) {
    double step = da;
    vector<double> v;
    v.push_back(min);
    double n = (max - min) / da;
    for (int i = 1; i < n; i++) {
        v.push_back(v[0] + da);
        da += step;
    }
    cout << "finished building alpha vector" << endl;
    return v;
}

vector<double> build_p_vector(double min, double max, double dp) {
    double step = dp;
    vector<double> v;
    v.push_back(min);
    int n = int((max - min) / dp);
    for (int i = 0; i < n; i++) {
        v.push_back(min + dp);
        dp += step;
    }
    return v;
}

vector<float> calculateKinterKintra(long double a) {
    //gets alpha as input and calculates k_inter and k_intra, assuming k_total is 4
    vector<float> res;
    float intra = float(k_total * a / (a + M - 1));
    float inter = float(k_total - intra);
    res.push_back(intra);
    res.push_back(inter);
    return res;
}

template<class T>
void printVector(vector<T> &v) {
    cout << "results:";
    long int n = v.size();
    for (int i = 0; i < n; i++) {
        cout << v[i] << " ";
        if (i % 100 == 0 && i>0)
            cout << endl;

    }
    cout << endl;
}

//template  <class T>
//T getVectorAvg(vector<T>& v){
//    T average = accumulate(v.begin(), v.end(), 0)/v.size();
//    return  average;
//}
//template <class T>
//T getVectorOfVectorsAverage(vector<vector<T> >& v) {
//    int n = v[0].size;
//    vector<int> sum(n, 0);
//    for (int i = 0; i < n; i++) {
//        sum += v[i];
//    }
//
//}

template<class T>
vector<float> getVectorOfVectorsAverage(vector<vector<T> > &v) {
    long int n = v.size();
    vector<float> sum(n, 0);
    vector<float> res(n, 0);
    for (int i = 0; i < n; i++)
        sum += v[i];

    for (int i = 0; i < n; i++)
        res[i] = sum[i] / n;

    return res;
}

template<class T>
T getVectorAvg(vector<T> &v) {
    T average = accumulate(v.begin(), v.end(), 0) / v.size();
    return average;
}
//set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -fopenmp")


void save_adj_list(OS os, vector< vector<int> >& adj) {
    //this function saves adjancency list
    cout << "saving adjancency list data..." << endl;
    ostringstream oss;
    oss << "N" << N << "_M" << M << "_adjList.txt";
    string str = oss.str();

    ofstream data;
    string path;
    switch (os) {
        case (WIN):
            path = "G:\\dropbox\\Dropbox\\M.Sc\\Results\\" + str;
            break;
        case (LINUX):
            path = "/home/eitan/Dropbox/M.Sc/results/" + str;
            break;
        case (MAC):
            path = "/Users/eitanasher/Desktop/" + str;
            break;
    }
    const char * c_path  = path.c_str();
    data.open(c_path, ofstream::out);
    for (int i = 0; i < N ; i++ ) {
            data << i ;
            for (int j=0 ; j<adj[i].size() ; j++) {
                data << "	" << adj[i][j];
            }
            data << endl;
        }
        data.close();
        cout<<" path = " << path<<endl;
        cout<< "saved succesfully"<<endl;

    }
//            ****** adjancency list******
//        cout<<"saving data..." << endl;
//
//        ostringstream oss;
//        oss <<"N" <<N<<"_M"<<M<<"_kt4_kinter"<<k_inter<<"_p"<<"0.001"<<"_permitted_distances_adjList.txt";
//        string str = oss.str();
//
//        ofstream data;
//        string path ="F:\\master's\\shlomo\\communities\\results\\" + str;
//        const char * c_path  = path.c_str();
//        data.open(c_path, ofstream::out);
//        for (int i = 0; i < N ; i++ ) {
//            data << i ;
//            for (int j=0 ; j<modules[1]->adj_list[i].size() ; j++) {
//                data << "	" << modules[1]->adj_list[i][j];
//            }
//            data << endl;
//        }
//        data.close();
//        cout<<" path = " << path<<endl;
//        cout<< "saved succesfully"<<endl;
