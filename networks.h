// Created by eitan and la on 27/05/2017.
#ifndef NETWORKSFIX_NETWORKS_H
#define NETWORKSFIX_NETWORKS_H
#define  P(A) cout << #A << ": " <<(A) << endl;

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
#include <assert.h>

using namespace std;

typedef std::mersenne_twister_engine<uint32_t, 32, 351,
        175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 1812433253> mt11213b;

enum networkType {
    ER, LATTICE
};
extern const int M, COMM_SIZE;
extern const int LATTICE_LENGTH, LATTICE_SIZE;
extern long int N;
extern float k_total;
mt11213b rnd;

template<class T>
vector<float> getVectorOfVectorsAverage(vector<vector<T> > &v);


class network {
public:
    //variables:
    networkType type_;
    vector<vector<int> > adj_list;
    vector<long int> distances;
    vector<long int> components;                   //for rvery node in the network, holds its module number
    vector<int> connected;
    vector<int> active;
    vector<int> degrees;
    vector<int> modules_inter_links_number;
    double k_inter, k_intra;

    //functions:m
    network(networkType type, double a);

    network(networkType type);

    void build_er_network();

    void build_er_network_with_k(float k);

    void build_lattices();

    void assign_inter_links();

    void assignKinterKintra(double a);

    void addEdge(int u, int v);

    void removeEdge(int u, int v);

    template<class T>
    T vector_sum(vector<T> &v);

    template<class T>
    T print_vector(vector<T> &v);

    bool isConnected(int u, int v);

    int getModule(int node_id);

    vector<long int> bfs(int s);

    long double findMeanPathLength();                       //calculate the networks' mean path length
    void setNetworkType(networkType t);

    networkType getNetworkType();

    void assign_inter_links(vector<vector<int> > &adj_list);

    void assign_inter_links(int num);

    int getDegree(vector<vector<int> > &adjlist, int u);

    vector<int> testConnectivity();

    vector<int> getSecondGiant(int max_clstr);

    vector<int> createRandomVector(int max);

    int getLongestPath();

    int getDegree(int u);

    float getAverageDegree();

    template<class T>
    T getVectorAvg(vector<T> &v);

    void deleteAll();

    vector<int> getGeometricalPosition(int node_id);

    float getGeometricalDistance(vector<int> x, vector<int> y);

    int getfirstLatticeNodeId(int community_id);

    vector<int> getGeometricalVsChemicalDist(float percent_of_nodes_to_check);

    void attack_highest_degree_nodes(double p);

    void attack_random(double p);

    void attack_highest_degree_modules(double p);

    void get_all_degrees();

    template<typename T>
    vector<size_t> sort_indexes(const vector<T> &v);

    template<typename T>
    vector<size_t> ordered(std::vector<T> const &values);

    ~network();
};

network::network(networkType t, double alpha) {
    setNetworkType(t);
    assignKinterKintra(alpha);
    switch (t) {
        case (ER):
            adj_list.resize(N);
            distances.resize(N);
            components.resize(N);
            connected.resize(N);
            active.resize(N, 1);
            degrees.resize(N);
            modules_inter_links_number.resize(M,0);
            break;
        case (LATTICE):
            adj_list.resize(N);
            distances.resize(N);
            components.resize(N);
            connected.resize(N);
            active.resize(N, 1);
            degrees.resize(N);
            modules_inter_links_number.resize(M,0);
            break;
    }
}

network::network(networkType t) {
    setNetworkType(t);
    switch (t) {
        case (ER):
            adj_list.resize(N);
            distances.resize(N);
            components.resize(N);
            connected.resize(N);
            active.resize(N, 1);
            degrees.resize(N);
            modules_inter_links_number.resize(M,1);
            break;

        case (LATTICE):
            adj_list.resize(N);
            distances.resize(N);
            components.resize(N);
            connected.resize(N);
            active.resize(N, 1);
            degrees.resize(N);
            modules_inter_links_number.resize(M,1);
            break;
    }
}

void network::build_lattices() {
    int n = LATTICE_LENGTH;
    P(LATTICE_SIZE);
    for (int CN = 0; CN < M; CN++) {
        adj_list[(CN * n * n)].push_back((CN * n * n) + 1);
        adj_list[(CN * n * n) + 1].push_back((CN * n * n));
        adj_list[(CN * n * n)].push_back((CN * n * n) + n);
        adj_list[(CN * n * n) + n].push_back((CN * n * n));

        //upper right
        adj_list[(CN * n * n) + n - 1].push_back((CN * n * n) + n - 2);
        adj_list[(CN * n * n) + n - 2].push_back((CN * n * n) + n - 1);
        adj_list[(CN * n * n) + n - 1].push_back((CN * n * n) + 2 * n - 1);
        adj_list[(CN * n * n) + 2 * n - 1].push_back((CN * n * n) + n - 1);

        //lower left - n*n-1
        adj_list[(CN * n * n) + n * (n - 1)].push_back((CN * n * n) + n * (n - 1) + 1);
        adj_list[(CN * n * n) + n * (n - 1) + 1].push_back((CN * n * n) + n * (n - 1));
        adj_list[(CN * n * n) + n * (n - 1)].push_back((CN * n * n) + n * (n - 2));
        adj_list[(CN * n * n) + n * (n - 2)].push_back((CN * n * n) + n * (n - 1));

        //lower right
        adj_list[(CN * n * n) + n * n - 1].push_back((CN * n * n) + n * n - 2);
        adj_list[(CN * n * n) + n * n - 2].push_back((CN * n * n) + n * n - 1);
        adj_list[(CN * n * n) + n * n - 1].push_back((CN * n * n) + n * n - n - 1);
        adj_list[(CN * n * n) + n * n - n - 1].push_back((CN * n * n) + n * n - 1);
        //upper and lower row
        for (int i = 1; i < n - 1; i++) {
            adj_list[(CN * n * n) + i].push_back((CN * n * n) + n + i);
            adj_list[(CN * n * n) + n + i].push_back((CN * n * n) + i);
            adj_list[(CN * n * n) + i].push_back((CN * n * n) + i + 1);
            adj_list[(CN * n * n) + i + 1].push_back((CN * n * n) + i);
            adj_list[(CN * n * n) + i].push_back((CN * n * n) + i - 1);
            adj_list[(CN * n * n) + i - 1].push_back((CN * n * n) + i);
            adj_list[(CN * n * n) + n * (n - 1) + i].push_back((CN * n * n) + n * (n - 2) + i);
            adj_list[(CN * n * n) + n * (n - 2) + i].push_back((CN * n * n) + n * (n - 1) + i);
            adj_list[(CN * n * n) + n * (n - 1) + i].push_back((CN * n) + n * (n - 1) + i - 1);
            adj_list[(CN * n * n) + n * (n - 1) + i - 1].push_back((CN * n) + n * (n - 1) + i);
            adj_list[(CN * n * n) + n * (n - 1) + i].push_back((CN * n * n) + n * (n - 1) + i + 1);
            adj_list[(CN * n * n) + n * (n - 1) + i + 1].push_back((CN * n * n) + n * (n - 1) + i);
        }
        //left and right columns
        for (int i = 1; i < n - 1; i++) {
            adj_list[(CN * n * n) + n * i].push_back((CN * n * n) + n * i + 1);
            adj_list[(CN * n * n) + n * i + 1].push_back((CN * n * n) + n * i);
            adj_list[(CN * n * n) + n * i].push_back((CN * n * n) + n * (i - 1));
            adj_list[(CN * n * n) + n * (i - 1)].push_back((CN * n * n) + n * i);
            adj_list[(CN * n * n) + n * i].push_back((CN * n * n) + n * (i + 1));
            adj_list[(CN * n * n) + n * (i + 1)].push_back((CN * n * n) + n * i);

            adj_list[(CN * n * n) + n * (i + 1) - 1].push_back((CN * n * n) + n * (i + 1) - 1 - 1);
            adj_list[(CN * n * n) + n * (i + 1) - 1 - 1].push_back((CN * n * n) + n * (i + 1) - 1);
            adj_list[(CN * n * n) + n * (i + 1) - 1].push_back((CN * n * n) + n * i - 1);
            adj_list[(CN * n * n) + n * i - 1].push_back((CN * n * n) + n * (i + 1) - 1);
            adj_list[(CN * n * n) + n * (i + 1) - 1].push_back((CN * n * n) + n * (i + 2) - 1);
            adj_list[(CN * n * n) + n * (i + 2) - 1].push_back((CN * n * n) + n * (i + 1) - 1);

        }
        //all internal vertices
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < n - 1; j++) {
                adj_list[(CN * n * n) + i * n + j].push_back(((CN * n * n) + i * n + j + 1));
                adj_list[(CN * n * n) + i * n + j + 1].push_back((CN * n * n) + i * n + j);
                adj_list[(CN * n * n) + i * n + j].push_back((CN * n * n) + i * n + j - 1);
                adj_list[(CN * n * n) + i * n + j - 1].push_back((CN * n * n) + i * n + j);

                adj_list[(CN * n * n) + i * n + j].push_back((CN * n) + (i - 1) * n + j);
                adj_list[(CN * n * n) + (i - 1) * n + j].push_back((CN * n * n) + i * n + j);
                adj_list[(CN * n * n) + i * n + j].push_back((CN * n * n) + (i + 1) * n + j);
                adj_list[(CN * n * n) + (i + 1) * n + j].push_back((CN * n * n) + i * n + j);
            }
        }
    }
}

void network::setNetworkType(networkType t) {
    type_ = t;
    //cout<<"network typre is "<<network.type<<endl;
}

networkType network::getNetworkType() {
    return type_;
}

void network::build_er_network() {
    rnd.seed(time(0));
    for (int m = 0; m < M; m++) {
        int num_of_edges;
        num_of_edges = int(round(COMM_SIZE * k_intra / 2));
        for (int i = 0; i < num_of_edges; i++) {
            int r1 = (m * COMM_SIZE) + (rnd() % COMM_SIZE);
            int r2 = (m * COMM_SIZE) + (rnd() % COMM_SIZE);
            //cout<<"r1 = "<<r1<<" r2 = "<<r2<<endl;
            if ((r1 != r2) && !isConnected(r1, r2)) {
                adj_list[r1].push_back(r2);
                adj_list[r2].push_back(r1);
                //cout<<"edge assigned, i = "<<i<<endl;
                i++;
            }
            i--;
        }
    }
    assign_inter_links();
}

void network::assign_inter_links(int num) {

    networkType t = getNetworkType();
    switch (t) {
        case (ER):
            cout << "number of inter-links in the network = " << num << endl;
            //double inter = round(k_inter * N / 2);
            for (int l = 0; l < num; l++) {
                int r1 = rnd() % N;
                int r2 = rnd() % N;
                //cout<<"r1 , r2 = "<<r1<<" & "<<r2<<endl;
                int m1 = getModule(r1);
                int m2 = getModule(r2);
                if (m1 !=m2 && !isConnected(r1, r2)) {
                    adj_list[r1].push_back(r2);
                    adj_list[r2].push_back(r1);
                    //keeping list of modules inter-connected links
                    modules_inter_links_number[m1]++;
                    modules_inter_links_number[m2]++;
                }
            }
            break;

        case (LATTICE):
            cout << "LATTICE" << endl;
            cout << "number of inter-links = " << round(k_inter * N / 2) << endl;
            for (int l = 0; l < round(k_inter * N / 2); l++) {
                int r1 = rnd() % N;
                int r2 = rnd() % N;
                //cout<<"r1 , r2 = "<<r1<<" & "<<r2<<endl;
                if (getModule(r1) != getModule(r2) && !isConnected(r1, r2)) {
                    adj_list[r1].push_back(r2);
                    adj_list[r2].push_back(r1);
                    //here need to add network of networks list
                }
            }
            cout << "finished assigning inter-links" << endl;
            break;
    }
}

void network::assign_inter_links() {

    networkType t = getNetworkType();
    switch (t) {
        case (ER):
            cout << "number of inter-links = " << round(k_inter * N / 2) << endl;
            //double inter = round(k_inter * N / 2);
            for (int l = 0; l < round(k_inter * N / 2); l++) {
                int r1 = rnd() % N;
                int r2 = rnd() % N;
                //cout<<"r1 , r2 = "<<r1<<" & "<<r2<<endl;
                if (getModule(r1) != getModule(r2) && !isConnected(r1, r2)) {
                    adj_list[r1].push_back(r2);
                    adj_list[r2].push_back(r1);
                    //here need to add network of networks list
                }
            }
            cout << "finished assigning inter-links" << endl;
            break;

        case (LATTICE):
            cout << "LATTICE" << endl;
            cout << "number of inter-links = " << round(k_inter * N / 2) << endl;
            for (int l = 0; l < round(k_inter * N / 2); l++) {
                int r1 = rnd() % N;
                int r2 = rnd() % N;
                //cout<<"r1 , r2 = "<<r1<<" & "<<r2<<endl;
                if (getModule(r1) != getModule(r2) && !isConnected(r1, r2)) {
                    adj_list[r1].push_back(r2);
                    adj_list[r2].push_back(r1);
                    //here need to add network of networks list
                }
            }
            cout << "finished assigning inter-links" << endl;
            break;
    }
}

int network::getModule(int u) {
    if (getNetworkType() == ER)
        return u / COMM_SIZE;
    else if (getNetworkType() == LATTICE)
        return u / LATTICE_SIZE;
    return 1;
}

vector<long int> network::bfs(int s) {
    vector<long int> results;
    bool visited[N];
    int longest_distance, num_p_degree, num_n_degree;
    long int degree = 1;
    long int sum = 0;
    long int counter = 0;
    //visited vector initialization:
    for (int k = 0; k < N; k++) visited[k] = 0;
    visited[s] = 1;

    queue<int> bfsq;
    bfsq.push(s);
    num_p_degree = 1;
    num_n_degree = 0;

    while (!bfsq.empty()) {
        int u = bfsq.front();
        bfsq.pop();
        ////ask bnaya whats n and p
        if (num_p_degree == 0) {
            degree++;
            num_p_degree = num_n_degree;
            num_n_degree = 0;
        }
        //looking at all neighbors of u
        for (int i = 0; i < adj_list[u].size(); i++) {
            int v = adj_list[u][i];

            if (!visited[v]) {
                visited[v] = 1;
                bfsq.push(v);
                num_n_degree++;
                sum += degree;
                counter += 1;
                distances[degree]++;
            }
        }
        num_p_degree--;
    }
    results.push_back(sum);
    results.push_back(counter);
    results.push_back(degree);
    return results;
}

vector<int> network::testConnectivity() {
    //this function uses BFS and returns 3 components: 1.number of modules in gcc, 2.number of nodes in gcc
    // 3.second giant component size
    vector<int> results;
    queue<int> Q;
    int size, expCount = 0, cluster_id = 0, giant_size = 0, giant_id;
    //int explored[N];         // 0-white, 1-grey, 2-black
    vector<int> explored(N, 0);
    //int *explored = new int[N];// 0-white, 1-grey, 2-black
    vector<int> community_gcc(M, 0);             //this vector holds the communities that have nodes that belong to the gcc

    for (long int i = 0; i < N; i++) {
        if (active[i] == 0 || getDegree(i) == 0) {
            explored[i] = 2;
            expCount++;
            components[i] = cluster_id++;
        } else
            explored[i] = 0;
    }
    int node_id = 0;
    while (expCount < N) {
        while (explored[node_id] != 0 || active[node_id]==0) {     //find a source to begin BFS
            node_id++;
        }
        cluster_id++;
        Q.push(node_id);
        components[node_id] = cluster_id;
        size = 1;
        explored[node_id] = 1;
        expCount++;
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            explored[u] = 2;
            for (int i = 0; i < adj_list[u].size(); i++) {
                if ( (explored[adj_list[u][i]]==0) && (active[adj_list[u][i]] == 1) ){
                    //cout<<"explored[list[u][i] = "<< explored[list[u][i]]<<endl;
                    size++;
                    Q.push(adj_list[u][i]);
                    components[adj_list[u][i]] = cluster_id;
                    explored[adj_list[u][i]] = 1;
                    expCount++;
                }
            }
        }
        if (size > giant_size) {
            giant_size = size;
            giant_id = cluster_id;
        }
    }
    //marking nodes and communities that belong to GCC
    for (int k = 0; k < N; k++) {
            if (components[k] == giant_id) {
            connected[k] = 1;
            community_gcc[getModule(k)] = 1;
        } else {
            connected[k] = 0;
        }
    }

    vector<int> second_gcc_results = getSecondGiant(cluster_id);
    int second_id = second_gcc_results[0];
    int second_count = second_gcc_results[1];

    cout<<"giant size = "<<giant_size<<" , second size = "<<second_count<<endl;
    results.push_back(giant_size);
    results.push_back(vector_sum(community_gcc));
    results.push_back(second_count);

    second_gcc_results.clear();
    explored.clear();
    community_gcc.clear();
    return results;
}

bool network::isConnected(int u, int v) {
    //l stands for list
    for (int k = 0; k < adj_list[u].size(); k++) {
        if (adj_list[u][k] == v) return true;
    }
    return false;
}

void network::assignKinterKintra(double a) {
    //gets alpha as input and sets k_inter and k_intra, assuming k_total is 4
    k_intra = k_total * a / (a + M - 1);
    k_inter = k_total - k_intra;
}

vector<int> network::getSecondGiant(int max_cluster) {
    //cout<<"calculating second giant"<<endl;
    vector<int> results;
    vector<int> count(max_cluster + 1, 0);
    for (int i = 0; i < N; i++) {
        count[components[i]]++;
    }
    // Traverse through the count[] and find second highest element.
    int first = 0, second_id = 0;
    for (int i = 0; i < max_cluster; i++) {
        /* If current element is higher than first then update both
          first and second */
        if (count[i] > count[first]) {
            second_id = first;
            first = i;
        }
            /* If count[i] is in between first and second then update second  */
        else if (count[i] > count[second_id] &&
                 count[i] != count[first])
            second_id = i;
    }
    results.push_back(second_id);    //second is the second gcc id
    results.push_back(count[second_id]);

    count.clear();
    //cout<<"finished calculating second"<<endl;
    return results;  //this is the number of nodes in second gcc

}

int network::getDegree(int u) {
    return adj_list[u].size();
}

template<class T>
T network::vector_sum(vector<T> &v) {
    T sum = 0;
    int n = v.size();
    for (int i = 0; i < n; i++) {
        sum += v[i];
    }
    return sum;
}

long double network::findMeanPathLength() {
    cout << "calculating mean path length.." << endl;
    long double sum = 0, counter = 0;
    vector<long int> mpl_v;
    for (int i = 0; i < 0.05 * N; i++) {          //calculating BFS only for 0.05 nodes
        int r = rnd() % N;
        mpl_v = bfs(r);                         //BFS(rand()%N);
        sum += mpl_v[0];
        counter += mpl_v[1];
        if (i % 1000 == 0) { cout << "BFS for node " << i << endl; }
    }
    return (long double) sum / counter;
}

network::~network() {
    adj_list.clear();
    distances.clear();
    components.clear();
    connected.clear();
}

int network::getLongestPath() {
    int longest = 0, dist;
    for (int i = 0; i < N; ++i) {
        if (i % 5000 == 0) { cout << "traversing node " << i << endl; }
        dist = network::bfs(i)[2];
        //cout<<"dist = "<<dist<<endl;
        if (dist > longest) {
            longest = dist;
            cout << "longest path found is = " << longest << endl;
        }
    }
    return longest;
}

float network::getAverageDegree() {
    int sum = 0;
    for (int i = 0; i < N; i++)
        sum += adj_list[i].size();
    return sum / N;
}

int network::getfirstLatticeNodeId(int community_id) {
    //this function returns the upper left node id in the lattice
    return community_id * LATTICE_SIZE;
}

vector<int> network::getGeometricalPosition(int node_id) {
    int i_in, i_out, j_in, j_out, s, i_f, j_f, community_number;
    vector<int> res;
    community_number = getModule(node_id);
    //relative position of nodes inside specific lattice
    s = node_id - getfirstLatticeNodeId(community_number);
    j_in = s % LATTICE_SIZE;
    i_in = (s - j_in) / LATTICE_SIZE;
    //position of the the specific lattice community relative to the lattice of lattices
    //need to
    int large_lattice_length = sqrt(M);
    //if ( float(sqrt(M)) % 1 != 0) cout<<"warning!  number of lattices is not a squred integer"<<endl;
    j_out = community_number % large_lattice_length;
    i_out = s / large_lattice_length;

    i_f = i_out * LATTICE_LENGTH + i_in;
    j_f = j_out * LATTICE_LENGTH + j_in;
    res.push_back(i_f);
    res.push_back(j_f);
    return res;
}

float network::getGeometricalDistance(vector<int> u, vector<int> v) {
    int u_x = u[0];
    int u_y = u[1];
    int v_x = v[0];
    int v_y = v[1];
    //cout<<"inside function (u_x-u_y)^2 + (v_x - v_y)^2 = " << (u_x-u_y)^2 + (v_x - v_y)^2<<endl;
    return sqrt(pow(u_x - v_x, 2) + pow(u_y - v_y, 2));
}

vector<int> network::getGeometricalVsChemicalDist(float percent_of_nodes_to_check) {
    //change the bfs algorithm in a way that suites our purpose.
    //1. calculate number of nodes to traverse from p
    //2. for each node run BFS
    //3. at each degree jump in bfs save all the nodes at this level
    //4. for each node calculate geometrical distance
    //5. store in array
    cout << "calculating geometrical vs. chemical" << endl;
    int number_of_nodes_to_check = percent_of_nodes_to_check * N;
    cout << "number_of_nodes_to_check = " << number_of_nodes_to_check << endl;
    vector<int> degree(N, 0);
    vector<int> geometricalDistance;
    geometricalDistance.push_back(
            1);   //storing 1 in the zeroth element because only one node is at this distance(self)

    vector<float> distances_storage;
    for (int idx = 0; idx < number_of_nodes_to_check; idx++) {
        if (idx % 5000 == 0) cout << "traversing node = " << idx << endl;
        int s = rnd() % N;  //source for bfs
        //cout<<"source node = "<<s<<endl;
        int num_p_degree = 1, num_n_degree = 0;
        degree[s] = 0;
        bool visited[N];
        long int degree = 1, counter = 0;
        queue<int> current_degree_nodes_storage;
        vector<int> geom_pos_1 = getGeometricalPosition(
                s);   // s is the source node. we make all our calculations relative to it, thus need to calculate its geometrical position

        for (int k = 0; k < N; k++) visited[k] = 0;
        visited[s] = 1;

        queue<int> bfsq;
        bfsq.push(s);
        num_p_degree = 1;
        num_n_degree = 0;

        while (!bfsq.empty()) {

            int u = bfsq.front();
            bfsq.pop();

            if (num_p_degree == 0) {
                //meaning we finished scaning a level, and now calculating distance from s to all nodes in current level
                degree++;
                num_p_degree = num_n_degree;
                num_n_degree = 0;

                for (int i = 0; current_degree_nodes_storage.size(); i++) {
                    int v = current_degree_nodes_storage.front();      // current node
                    vector<int> geom_pos_2 = getGeometricalPosition(v);
                    float distance = getGeometricalDistance(geom_pos_1, geom_pos_2);
                    distances_storage.push_back(distance);
                    current_degree_nodes_storage.pop();
                }
                //after calculating distances from every node in current degree,
                // we average and store it in geometrical distances at its right place
                assert(current_degree_nodes_storage.empty());
                geometricalDistance.push_back(getVectorAvg(distances_storage));
                distances_storage.clear();
                assert(distances_storage.empty());
            }

            //looking at all u neighbors
            for (int i = 0; i < adj_list[u].size(); i++) {
                int v = adj_list[u][i];

                if (!visited[v]) {
                    visited[v] = 1;
                    bfsq.push(v);
                    num_n_degree++;
                    current_degree_nodes_storage.push(v);
                }
            }
            num_p_degree--;
        }
    }

    vector<int>::const_iterator first = geometricalDistance.begin();
    vector<int>::const_iterator last = first + 100;
    vector<int> newVec(first, last);
    return newVec;
}

template<class T>
T network::getVectorAvg(vector<T> &v) {
    T average = accumulate(v.begin(), v.end(), 0) / v.size();
    return average;
}

void network::build_er_network_with_k(float k) {
    rnd.seed(time(0));
    for (int m = 0; m < M; m++) {
        int num_of_edges;
        num_of_edges = int(round(COMM_SIZE * k / 2));
        for (int i = 0; i < num_of_edges; i++) {
            int r1 = (m * COMM_SIZE) + (rnd() % COMM_SIZE);
            int r2 = (m * COMM_SIZE) + (rnd() % COMM_SIZE);
            //cout<<"r1 = "<<r1<<" r2 = "<<r2<<endl;
            if ((r1 != r2) && !isConnected(r1, r2)) {
                adj_list[r1].push_back(r2);
                adj_list[r2].push_back(r1);
                i++;
            }
            i--;
        }
    }
    for (int i=0; i<N; i++)
        active[i]=1;
    // assign_inter_links();
}

void network::attack_highest_degree_nodes(double p) {
    cout<<"attacking"<<endl;
    get_all_degrees();
    int num_of_nodes_to_kill = round(p * COMM_SIZE);
    cout << "number of nodes to kill per community= " << num_of_nodes_to_kill << endl;
    for (int m = 0; m < M; m++) {
        int begin_idx = m * COMM_SIZE;  // this holds the first index of the begining of current community
        int end_idx = (m + 1) * COMM_SIZE;
        //vector<int>::const_iterator begin = arr.begin();
        vector<int>::const_iterator first = degrees.begin() + begin_idx;
        vector<int>::const_iterator last = degrees.begin() + end_idx;

        vector<int> curr(first, last);

        //vector<int> curr (&adj_list[begin_idx],&adj_list[end_idx]); //current community adj list
        vector<size_t> indices = ordered(curr);

        for (int i = 0; i < num_of_nodes_to_kill; i++) {
            active[indices[i] + m*COMM_SIZE] = 0;
        }
    }
}

void network::attack_random(double p){
    int num_of_nodes_to_kill = round(p*N);
    cout<<"number of nodes to randomly kill in the network = "<<num_of_nodes_to_kill<<endl;
    for (int i=0; i<num_of_nodes_to_kill; i++){
        int r = rnd()%N;
        if (active[r]==1){
            active[r]=0;
            i++;
        }
        i--;
    }
}

void network::attack_highest_degree_modules(double p){

}
void network::get_all_degrees() {
    for (int i = 0; i < N; i++)
        degrees[i] = adj_list[i].size();
}

template<typename T>
vector<size_t> network::sort_indexes(const vector<T> &v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

    return idx;
}

template<typename T>
vector<size_t> network::ordered(std::vector<T> const &values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

template<class T>
T network::print_vector(vector<T> &v) {
    cout << "results:";
    long int n = v.size();
    for (int i = 0; i < n; i++) {
        cout << v[i] << " ";
        if (i % 100 == 0 && i > 0)
            cout << endl;

    }
    cout << endl;
}
#endif //NETWORKSFIX_NETWORKS_H
