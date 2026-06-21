// build and write to the disk all Turan EX(n,G) graphs for some bipartite graph G.
// Here we assume that G is C4, C6, C8, K23, K24, K25, K26, K33, K34, K35, K36, K44 or Q3.
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <thread>
#include <chrono>
#include <mutex>

#include "Graph.h"

#define VERBOSE 1

const std::vector<std::string> names = {"C4", "C6", "C8", "K22", "K23", "K24", "K25", "K26", "K33", "K34", "K35", "K36", "K44", "Q3"};

int num_threads = 8;

class Turan {
public:
    Turan(const std::string& graph_name) : graph_name(graph_name) {
        initH();

        mkdir("../data/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir("../data/TURAN/", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        address = "../data/TURAN/" + graph_name + "/";
        mkdir(address.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        std::string cr_address = address + "Critical/";
        std::string ex_address = address + "Extremal/";

        mkdir(cr_address.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(ex_address.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        ln.resize(Hn);
        un.resize(Hn);
        qx.resize(Hn);
        ex.resize(Hn);
        for (int n = 1; n < Hn; ++n) {
            ex[n] = n * (n - 1) / 2;
        }
 
        ln[Hn - 1] = ((Hn - 1) * (Hn - 2)) / 2;
        un[Hn - 1] = ((Hn - 1) * (Hn - 2)) / 2;

        GraphSet one(Hn - 1);
        one.insert(K(Hn - 1).certify());
        std::string filename = address + "Critical/Cr(" + std::to_string(Hn - 1) + ", " + graph_name + ").gr";
        one.write(filename);

        DG.resize(num_threads);
        cycles.resize(num_threads);
    }

    void initH() {
        if (graph_name == "C4") {
            H = C(4);
        } else if (graph_name == "C6") {
            H = C(6);
        } else if (graph_name == "C8") {
            H = C(8);
        } else if (graph_name == "K23") {
            H = K(2, 3);
        } else if (graph_name == "K24") {
            H = K(2, 4);
        } else if (graph_name == "K25") {
            H = K(2, 5);
        } else if (graph_name == "K26") {
            H = K(2, 6);
        } else if (graph_name == "K33") {
            H = K(3, 3);
        } else if (graph_name == "K34") {
            H = K(3, 4);
        } else if (graph_name == "K35") {
            H = K(3, 5);
        } else if (graph_name == "K36") {
            H = K(3, 6);
        } else if (graph_name == "K44") {
            H = K(4, 4);
        } else if (graph_name == "Q3") {
            H = Q(3);
        }

        dH = H.deg();
        Hn = H.size();
        He = H.edges();

        He1.assign(He, 0);
        He2.assign(He, 0);

        int k = 0;
        for (int i = 0; i < Hn; i++) {
            for (int j = i + 1; j < Hn; j++) {
                if (H.edge(i, j)) {
                    He1[k] = i;
                    He2[k] = j;
                    k++;
                }
            }
        }
    }

    int start() const {
        return Hn;
    }

    std::pair<int, int> compute(int n) {
        N = n;
        while (N >= ln.size()) {
            ln.push_back(0);
            un.push_back(0);
            qx.push_back(0);
            ex.push_back(0);
        }

        CR.resize(N);
        EX.resize(N);

        un[N] = un[N - 1];
        while ((un[N] + 1) - (2 * un[N] + 2) / N <= un[N - 1]) {
            un[N]++;
        }

        ln[N] = un[N];
        getBounds(N - 1);
        getGraphs();

        while (EX.empty()) {
            un[N]--;
            ln[N]--;
 
            getBounds(N - 1);
            getGraphs();
        }

        std::string path = address + "Extremal/EX(" + std::to_string(n) + ", " + graph_name + ").s6";
        qx[N] = EX.size();
        ex[N] = ln[N];
        EX.writeSparse6(N, path);
        EX.clear();

        return {ex[N], qx[N]};
    }

private:

    void getCycles(const Graph& G, int th) {
        size_t n = G.size();
        cycles[th].clear();
        std::vector<size_t> g(Hn);
        if (graph_name == "C4") { // k -- i -- (n - 1) -- j -- k
            for (size_t i = 0; i < n - 1; i++) if (G.edge(i, n - 1))
            for (size_t j = i + 1; j < n - 1; j++) if (G.edge(j, n - 1))
            for (size_t k = 0; k < n - 1; k++) if (G.edge(i, k) && G.edge(j, k)) {
                 g[0] = n - 1;
                 g[1] = i;
                 g[2] = k;
                 g[3] = j;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "C6") { //   k -- j1 -- i1 -- (n - 1) -- i2 -- j2 -- k
            for (size_t i1 = 0; i1 < n - 1; i1++) if (G.edge(i1, n - 1))
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++) if (G.edge(i2, n - 1))
            for (size_t j1 = 0; j1 < n - 1; j1++)      if (G.edge(j1, i1) && j1 != i2)
            for (size_t j2 = 0; j2 < n - 1; j2++)      if (G.edge(j2, i2) && j2 != i1 && j2 != j1)
            for (size_t k = 0; k < n - 1; k++)	if (G.edge(k, j1) && G.edge(k, j2) && k != i1 && k != i2) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = j1;
                 g[3] = k;
                 g[4] = j2;
                 g[5] = i2;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "C8") { //   k -- l1 -- j1 -- i1 -- (n - 1) -- i2 -- j2 -- l2 -- k
            for (size_t i1 = 0; i1 < n - 1; i1++) if (G.edge(i1, n - 1))
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++) if (G.edge(i2, n - 1))
            for (size_t j1 = 0; j1 < n - 1; j1++)      if (G.edge(j1, i1) && j1 != i2)
            for (size_t j2 = 0; j2 < n - 1; j2++)      if (G.edge(j2, i2) && j2 != i1 && j2 != j1)
            for (size_t l1 = 0; l1 < n - 1; l1++)      if (G.edge(l1, j1) && l1 != i1 && l1 != i2 && l1 != j2)
            for (size_t l2 = 0; l2 < n - 1; l2++)      if (G.edge(l2, j2) && l2 != i1 && l2 != i2 && l2 != j1 && l1 != l2)
            for (size_t k = 0; k < n - 1; k++)	if (G.edge(k, l1) && G.edge(k, l2) && k != i1 && k != i2 && k != j1 && k != j2) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = j1;
                 g[3] = l1;
                 g[4] = k;
                 g[5] = l2;
                 g[6] = j2;
                 g[7] = i2;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K23") { // (n - 1), i =-= j0, j1, j2
            for (size_t i = 0; i < n - 1; i++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i)) {
                 g[0] = n - 1;
                 g[1] = i;
                 g[2] = j0;
                 g[3] = j1;
                 g[4] = j2;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2 =-= j0, j1
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = n - 1;
                 g[3] = i1;
                 g[4] = i2;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K24") { // (n - 1), i =-= j0, j1, j2, j3
            for (size_t i = 0; i < n - 1; i++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i)) {
                 g[0] = n - 1;
                 g[1] = i;
                 g[2] = j0;
                 g[3] = j1;
                 g[4] = j2;
                 g[5] = j3;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2, i3 =-= j0, j1
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = n - 1;
                 g[3] = i1;
                 g[4] = i2;
                 g[5] = i3;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K25") { // (n - 1), i =-= j0, j1, j2, j3, j4
            for (size_t i = 0; i < n - 1; i++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i))
            for (size_t j4 = j3 + 1; j4 < n - 1; j4++)	if (G.edge(j4, n - 1) && G.edge(j4, i)) {
                 g[0] = n - 1;
                 g[1] = i;
                 g[2] = j0;
                 g[3] = j1;
                 g[4] = j2;
                 g[5] = j3;
                 g[6] = j4;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2, i3 =-= j0, j1, j2
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2) && G.edge(j2, i3)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = j2;
                 g[3] = n - 1;
                 g[4] = i1;
                 g[5] = i2;
                 g[6] = i3;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K26") { // (n - 1), i =-= j0, j1, j2, j3, j4, j5
            for (size_t i = 0; i < n - 1; i++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i))
            for (size_t j4 = j3 + 1; j4 < n - 1; j4++)	if (G.edge(j4, n - 1) && G.edge(j4, i))
            for (size_t j5 = j4 + 1; j5 < n - 1; j5++)	if (G.edge(j5, n - 1) && G.edge(j5, i)) {
                 g[0] = n - 1;
                 g[1] = i;
                 g[2] = j0;
                 g[3] = j1;
                 g[4] = j2;
                 g[5] = j3;
                 g[6] = j4;
                 g[7] = j5;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2, i3, i4, i5 =-= j0, j1
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t i4 = i3 + 1; i4 < n - 1; i4++)
            for (size_t i5 = i4 + 1; i5 < n - 1; i5++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3) && G.edge(j0, i4) && G.edge(j0, i5))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3) && G.edge(j1, i4) && G.edge(j1, i5)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = n - 1;
                 g[3] = i1;
                 g[4] = i2;
                 g[5] = i3;
                 g[6] = i4;
                 g[7] = i5;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K33") { // (n - 1), i1, i2 =-= j0, j1, j2
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2)) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = i2;
                 g[3] = j0;
                 g[4] = j1;
                 g[5] = j2;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K34") { // (n - 1), i1, i2 =-= j0, j1, j2, j3
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i1) && G.edge(j3, i2)) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = i2;
                 g[3] = j0;
                 g[4] = j1;
                 g[5] = j2;
                 g[6] = j3;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2, i3 =-= j0, j1, j2
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2) && G.edge(j2, i3)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = j2;
                 g[3] = n - 1;
                 g[4] = i1;
                 g[5] = i2;
                 g[6] = i3;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K35") { // (n - 1), i1, i2 =-= j0, j1, j2, j3, j4
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i1) && G.edge(j3, i2))
            for (size_t j4 = j3 + 1; j4 < n - 1; j4++)	if (G.edge(j4, n - 1) && G.edge(j4, i1) && G.edge(j4, i2)) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = i2;
                 g[3] = j0;
                 g[4] = j1;
                 g[5] = j2;
                 g[6] = j3;
                 g[7] = j4;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2, i3, i4 =-= j0, j1, j2
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t i4 = i3 + 1; i4 < n - 1; i4++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3) && G.edge(j0, i4))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3) && G.edge(j1, i4))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2) && G.edge(j2, i3) && G.edge(j2, i4)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = j2;
                 g[3] = n - 1;
                 g[4] = i1;
                 g[5] = i2;
                 g[6] = i3;
                 g[7] = i4;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K36") { // (n - 1), i1, i2 =-= j0, j1, j2, j3, j4, j5
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i1) && G.edge(j3, i2))
            for (size_t j4 = j3 + 1; j4 < n - 1; j4++)	if (G.edge(j4, n - 1) && G.edge(j4, i1) && G.edge(j4, i2))
            for (size_t j5 = j4 + 1; j5 < n - 1; j5++)	if (G.edge(j5, n - 1) && G.edge(j5, i1) && G.edge(j5, i2)) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = i2;
                 g[3] = j0;
                 g[4] = j1;
                 g[5] = j2;
                 g[6] = j3;
                 g[7] = j4;
                 g[8] = j5;
                 cycles[th].emplace_back(g);
            }                             // (n - 1), i1, i2, i3, i4, i5 =-= j0, j1, j2
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t i4 = i3 + 1; i4 < n - 1; i4++)
            for (size_t i5 = i4 + 1; i5 < n - 1; i5++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3) && G.edge(j0, i4) && G.edge(j0, i5))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3) && G.edge(j1, i4) && G.edge(j1, i5))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2) && G.edge(j2, i3) && G.edge(j2, i4) && G.edge(j2, i5)) {
                 g[0] = j0;
                 g[1] = j1;
                 g[2] = j2;
                 g[3] = n - 1;
                 g[4] = i1;
                 g[5] = i2;
                 g[6] = i3;
                 g[7] = i4;
                 g[8] = i5;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "K44") { // (n - 1), i1, i2, i3 =-= j0, j1, j2, j3
            for (size_t i1 = 0; i1 < n - 1; i1++)
            for (size_t i2 = i1 + 1; i2 < n - 1; i2++)
            for (size_t i3 = i2 + 1; i3 < n - 1; i3++)
            for (size_t j0 = 0; j0 < n - 1; j0++)	if (G.edge(j0, n - 1) && G.edge(j0, i1) && G.edge(j0, i2) && G.edge(j0, i3))
            for (size_t j1 = j0 + 1; j1 < n - 1; j1++)	if (G.edge(j1, n - 1) && G.edge(j1, i1) && G.edge(j1, i2) && G.edge(j1, i3))
            for (size_t j2 = j1 + 1; j2 < n - 1; j2++)	if (G.edge(j2, n - 1) && G.edge(j2, i1) && G.edge(j2, i2) && G.edge(j2, i3))
            for (size_t j3 = j2 + 1; j3 < n - 1; j3++)	if (G.edge(j3, n - 1) && G.edge(j3, i1) && G.edge(j3, i2) && G.edge(j3, i3)) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = i2;
                 g[3] = i3;
                 g[4] = j0;
                 g[5] = j1;
                 g[6] = j2;
                 g[7] = j3;
                 cycles[th].emplace_back(g);
            }
        } else if (graph_name == "Q3") {                                     
            for (size_t i1 = 0; i1 < n - 1; ++i1) if (G.edge(i1, n - 1))     
            for (size_t i2 = i1 + 1; i2 < n - 1; ++i2) if (G.edge(i2, n - 1))
            for (size_t i3 = i2 + 1; i3 < n - 1; ++i3) if (G.edge(i3, n - 1))
            for (size_t j1 = 0; j1 < n - 1; ++j1) if (G.edge(j1, i1) && G.edge(j1, i2) && j1 != i3)  
            for (size_t j2 = 0; j2 < n - 1; ++j2) if (G.edge(j2, i1) && G.edge(j2, i3) && j2 != j1 && j2 != i2)
            for (size_t j3 = 0; j3 < n - 1; ++j3) if (G.edge(j3, i2) && G.edge(j3, i3) && j3 != j1 && j3 != j2 && j3 != i1)  
            for (size_t k = 0; k < n - 1; ++k) if (G.edge(k, j1) && G.edge(k, j2) && G.edge(k, j3) && k != i1 && k != i2 && k != i3) {
                 g[0] = n - 1;
                 g[1] = i1;
                 g[2] = i2;
                 g[3] = j1;
                 g[4] = i3;
                 g[5] = j2;
                 g[6] = j3;
                 g[7] = k;
                 cycles[th].emplace_back(g);
            }
        }
    }

    void next(int level, int n) {
        if (level < perm.size()) {
            for (int i = perm[level - 1] + 1; i < n; i++) {
                perm[level] = i;
                next(level + 1, n);
            }
        } else {
            cliques.push_back(perm);
        }
    }

    void getCliques(int n, int p) {
        cliques.clear();
        if (n == Hn - 1) {
            // in this case there is only one critical graph and one clique of each size up to isomorphism
            std::vector<size_t> g(p);
            for (int i = 0; i < p; i++) {
                g[i] = i;
            }
            cliques.emplace_back(std::move(g));
        } else {
            perm.resize(p);
            for (int i = 0; i <= n - p; i++) {
                perm[0] = i;
                next(1, n);
            }
        }
    }

    bool critical(Graph& G) {
        size_t n = G.size();
        for (int ii = 0; ii < n; ii++) {
            for (int jj = ii + 1; jj < n; jj++) {
                if (G.edge(ii, jj)) {
                    continue;
                }

                G.addEdge(ii, jj);

                bool BB = false;
                if (graph_name == "C4") { // i -- ii -- jj -- j
                    for (int i = 0; i < n && !BB; i++) if (G.edge(i, ii))
                    for (int j = 0; j < n && !BB; j++) if (G.edge(j, jj) && G.edge(i, j)) {
                        BB = true;
                    }
                } else if (graph_name == "C6") { //  i2 -- i1 -- ii -- jj -- j1 -- j2 -- i2
                    for (int i1 = 0; i1 < n && !BB; i1++) if (G.edge(i1, ii) && i1 != jj)
                    for (int j1 = 0; j1 < n && !BB; j1++) if (G.edge(j1, jj) && j1 != ii && j1 != i1)
                    for (int i2 = 0; i2 < n && !BB; i2++) if (G.edge(i1, i2) && i2 != j1 && i2 != ii && i2 != jj)
                    for (int j2 = 0; j2 < n && !BB; j2++) if (G.edge(j1, j2) && G.edge(i2, j2) && j2 != i1 && j2 != ii && j2 != jj) {
                        BB = true;
                    }
                } else if (graph_name == "C8") { //  i3 -- i2 -- i1 -- ii -- jj -- j1 -- j2 -- i2 -- i3
                    for (int i1 = 0; i1 < n && !BB; i1++) if (G.edge(i1, ii) && i1 != jj)
                    for (int j1 = 0; j1 < n && !BB; j1++) if (G.edge(j1, jj) && j1 != ii && j1 != i1)
                    for (int i2 = 0; i2 < n && !BB; i2++) if (G.edge(i1, i2) && i2 != j1 && i2 != ii && i2 != jj)
                    for (int j2 = 0; j2 < n && !BB; j2++) if (G.edge(j1, j2) && j2 != i2 && j2 != i1 && j2 != ii && j2 != jj)
                    for (int i3 = 0; i3 < n && !BB; i3++) if (G.edge(i2, i3) && i3 != j2 && i3 != i2 && i3 != j1 && i3 != ii && i3 != jj)
                    for (int j3 = 0; j3 < n && !BB; j3++) if (G.edge(j2, j3) && G.edge(i3, j3) && j3 != i2 && j3 != j1 && j3 != ii && j3 != jj) {
                        BB = true;
                    }
                } else if (graph_name == "K23") { // i1, ii =-= jj, j1, j2
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2)) {
                        BB = true;
                    }                             // i1, i2, ii =-= jj, j1
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1)) {
                        BB = true;
                    }
                } else if (graph_name == "K24") { // i1, ii =-= jj, j1, j2, j3
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3)) {
                        BB = true;
                    }                             // i1, i2, i3, ii =-= jj, j1
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1)) {
                        BB = true;
                    }
                } else if (graph_name == "K25") { // i1, ii =-= jj, j1, j2, j3, j4
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3))
                    for (int j4 = j3 + 1; j4 < n && !BB; j4++) if (jj != j4 && G.edge(ii, j4) && G.edge(i1, j4)) {
                        BB = true;
                    }                             // i1, i2, i3, i4, ii =-= jj, j1
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int i4 = i3 + 1; i4 < n && !BB; i4++) if (i4 != ii && G.edge(i4, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1) && G.edge(i4, j1)) {
                        BB = true;
                    }
                } else if (graph_name == "K26") { // i1, ii =-= jj, j1, j2, j3, j4, j5
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3))
                    for (int j4 = j3 + 1; j4 < n && !BB; j4++) if (jj != j4 && G.edge(ii, j4) && G.edge(i1, j4))
                    for (int j5 = j4 + 1; j5 < n && !BB; j5++) if (jj != j5 && G.edge(ii, j5) && G.edge(i1, j5)) {
                        BB = true;
                    }                             // i1, i2, i3, i4, i5, ii =-= jj, j1
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int i4 = i3 + 1; i4 < n && !BB; i4++) if (i4 != ii && G.edge(i4, jj))
                    for (int i5 = i4 + 1; i5 < n && !BB; i5++) if (i5 != ii && G.edge(i5, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1) && G.edge(i4, j1) && G.edge(i5, j1)) {
                        BB = true;
                    }
                } else if (graph_name == "K33") { // i1, i2, ii =-= jj, j1, j2
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2)) {
                        BB = true;
                    }
                } else if (graph_name == "K34") { // i1, i2, ii =-= jj, j1, j2, j3
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3) && G.edge(i2, j3)) {
                        BB = true;
                    }                             // i1, i2, i3, ii =-= jj, j1, j2
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2) && G.edge(i3, j2)) {
                        BB = true;
                    }
                } else if (graph_name == "K35") { // i1, i2, ii =-= jj, j1, j2, j3, j4
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3) && G.edge(i2, j3))
                    for (int j4 = j3 + 1; j4 < n && !BB; j4++) if (jj != j4 && G.edge(ii, j4) && G.edge(i1, j4) && G.edge(i2, j4)) {
                        BB = true;
                    }                             // i1, i2, i3, i4, ii =-= jj, j1, j2
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int i4 = i3 + 1; i4 < n && !BB; i4++) if (i4 != ii && G.edge(i4, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1) && G.edge(i4, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2) && G.edge(i3, j2) && G.edge(i4, j2)) {
                        BB = true;
                    }
                } else if (graph_name == "K36") { // i1, i2, ii =-= jj, j1, j2, j3, j4, j5
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3) && G.edge(i2, j3))
                    for (int j4 = j3 + 1; j4 < n && !BB; j4++) if (jj != j4 && G.edge(ii, j4) && G.edge(i1, j4) && G.edge(i2, j4))
                    for (int j5 = j4 + 1; j5 < n && !BB; j5++) if (jj != j5 && G.edge(ii, j5) && G.edge(i1, j5) && G.edge(i2, j5)) {
                        BB = true;
                    }                             // i1, i2, i3, i4, i5, ii =-= jj, j1, j2
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int i4 = i3 + 1; i4 < n && !BB; i4++) if (i4 != ii && G.edge(i4, jj))
                    for (int i5 = i4 + 1; i5 < n && !BB; i5++) if (i5 != ii && G.edge(i5, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1) && G.edge(i4, j1) && G.edge(i5, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2) && G.edge(i3, j2) && G.edge(i4, j2) && G.edge(i5, j2)) {
                        BB = true;
                    }
                } else if (graph_name == "K44") { // i1, i2, i3, ii =-= jj, j1, j2, j3
                    for (int i1 = 0; i1 < n && !BB; i1++) if (i1 != ii && G.edge(i1, jj))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (i2 != ii && G.edge(i2, jj))
                    for (int i3 = i2 + 1; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, jj))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (jj != j1 && G.edge(ii, j1) && G.edge(i1, j1) && G.edge(i2, j1) && G.edge(i3, j1))
                    for (int j2 = j1 + 1; j2 < n && !BB; j2++) if (jj != j2 && G.edge(ii, j2) && G.edge(i1, j2) && G.edge(i2, j2) && G.edge(i3, j2))
                    for (int j3 = j2 + 1; j3 < n && !BB; j3++) if (jj != j3 && G.edge(ii, j3) && G.edge(i1, j3) && G.edge(i2, j3) && G.edge(i3, j3)) {
                        BB = true;
                    }
                } else if (graph_name == "Q3") { // i3 < i1, i2 > ii =-= jj < j1, j2 > j3
                    for (int i1 = 0; i1 < n && !BB; i1++) if (G.edge(i1, ii))
                    for (int i2 = i1 + 1; i2 < n && !BB; i2++) if (G.edge(i2, ii))
                    for (int i3 = 0; i3 < n && !BB; i3++) if (i3 != ii && G.edge(i3, i1) && G.edge(i3, i2))
                    for (int j1 = 0; j1 < n && !BB; j1++) if (j1 != ii && j1 != i2 && j1 != i3 && G.edge(j1, jj) && G.edge(j1, i1))
                    for (int j2 = 0; j2 < n && !BB; j2++) if (j2 != ii && j2 != i1 && j2 != i3 && j2 != j1 && G.edge(j2, jj) && G.edge(j2, i2))
                    for (int j3 = 0; j3 < n && !BB; j3++) if (j3 != ii && j3 != jj && j3 != i1 && j3 != i2 && G.edge(j3, j1) && G.edge(j3, i2) && G.edge(j3, i3)) {
                        BB = true;
                    }
                }
                G.killEdge(ii, jj);
                if (!BB) {
                    return false;
                }
            }
        }
        return true;
    }

    void nextCycle(Graph& G, int th, int level) {
        int n = G.size();
        if (G.edges() < ln[n] || G.edges() > un[n] + cycles[th].size() - level) {
            return;
        }
        const std::vector<size_t>& cycle = cycles[th][level];

        if (level < cycles[th].size()) {
            // do we still need to check this cycle???
            bool B = true;
            for (int ii = 0; ii < He; ii++) {
                if (!G.edge(cycle[He1[ii]], cycle[He2[ii]])) {
                    B = false;
                    break;
                }
            }

            if (!B) {
                nextCycle(G, th, level + 1);
            } else { // we need to clean this cycle off
                for (int ii = 0; ii < He; ii++) {
                    if (cycle[He1[ii]] != n - 1 && cycle[He2[ii]] != n - 1) {
                        if (DG[th][cycle[He1[ii]]] > d && DG[th][cycle[He2[ii]]] > d) {
                            G.killEdge(cycle[He1[ii]], cycle[He2[ii]]);
                            DG[th][cycle[He1[ii]]]--;
                            DG[th][cycle[He2[ii]]]--;

                            nextCycle(G, th, level + 1);

                            G.addEdge(cycle[He1[ii]], cycle[He2[ii]]);
                            DG[th][cycle[He1[ii]]]++;
                            DG[th][cycle[He2[ii]]]++;
                        }
                    }
                }
            }
        } else {
            if (G.deg() == d) {
                if (critical(G)) {
                    G.certify();
                    CR.insert(G);
                    if (n == N && G.edges() == ln[N]) {
                        EX.insert(G);
                    }
                }
            }
        }
    }

    void deleteCycles(Graph& G, int th) {
        DG[th] = G.getDegrees();
        nextCycle(G, th, 0);
    }

    void getGraphs() {
        for (int n = Hn; n <= N; n++) {
            if (ln[n] > un[n]) {
                continue;
            }

            for (d = n - 1; d >= dH - 1; --d) {   // adding new vertex of degree  d > 0 // vertex number [n-1]
                std::fstream stream;
                stream.open(address + "Critical/Cr(" + std::to_string(n - 1) + ", " + graph_name + ").gr", std::ios::in | std::ios::binary);
                getCliques(n - 1, d);
                std::mutex read_mutex;
                std::vector<std::thread> threads;
                for (int th = 0; th < num_threads; ++th) {
                    threads.emplace_back([this, n, th, &stream, &read_mutex] {
                    Graph G(n - 1);
                    for (;;) {
                        read_mutex.lock();
                        bool end = readGraph(stream, G);
                        read_mutex.unlock();
                        if (end == false) {
                            return;
                        }
		                if (G.deg() + 1 < d || G.edges() + d < ln[n]) { //if minimal degree is small enough && there are enough edges
		                    continue;
		                }
                        Graph F = G + 1;
                        for (const std::vector<size_t>& clique : cliques) {
                            for (int x : clique) { // adding  a vertex [n-1] of degree d to the clique [i]
                                F.addEdge(x, n - 1);
                            }
                            if (F.deg() >= d) {
                                getCycles(F, th);
                                deleteCycles(F, th); // deleting all cycles and adding new critical and extremal graphs
                            }
                            for (int x : clique) { // adding  a vertex [n-1] of degree d to the clique [i]
                                F.killEdge(x, n - 1);
                            } 
                        }
                    }
                    });
                }
                for (int th = 0; th < num_threads; ++th) {
                    threads[th].join();
                }
                
                //adding new graphs to lists
                std::string path = address + "Critical/Cr(" + std::to_string(n) + ", " + graph_name + ").gr";
                if (!CR.empty()) {
                    CR.write(path, true);
                }
                CR.clear();
            }
        }
    }

    void getBounds(int k) {
        for (int n = k; n >= Hn - 1; n--) {
            un[n] = ln[n] - 1;
        }

        while (k >= Hn) {
            if (ln[k] > ln[k + 1] - ((2 * ln[k + 1]) / (k + 1))) {
                ln[k] = ln[k + 1] - ((2 * ln[k + 1]) / (k + 1));
            }
            k--;
        }
    }

    Graph H;
    int N;
    int dH, Hn, He;
    int d;

    GraphSet CR;
    GraphSet EX;

    std::string graph_name;
    std::string address;

    std::vector<size_t> perm;
    std::vector<size_t> ln;
    std::vector<size_t> un;
    std::vector<size_t> ex;
    std::vector<size_t> qx;
    std::vector<size_t> He1;
    std::vector<size_t> He2;

    std::vector<std::vector<size_t>> DG;
    std::vector<std::vector<size_t>> cliques;
    std::vector<std::vector<std::vector<size_t>>> cycles;
};


int main(int argc, char** argv) {
    if (argc != 2 && argc != 3) {
        std::cout << "We expect graph name G to compute the set EX(n,G)" << std::endl;
        return 1;
    }

    std::string graph_name = argv[1];
    if (std::find(names.begin(), names.end(), graph_name) == names.end()) {
        std::cout << "Wrong graph name. We expect G to be C4, C6, C8, K23, K24, K25, K26, K33, K34, K35, K36, K44 or Q3" << std::endl;
        return 1;       
    }
    if (graph_name == "K22") {
        graph_name = "C4";
    }

    if (argc == 3) {
        num_threads = std::atoi(argv[2]);
    }

    // compute as much Turan numbers as possible
    Turan turan(graph_name);
    std::cout << "Values of Turan numbers for the graph " << graph_name << std::endl;
    std::cout << "--------------------------------------";
    for (char c : graph_name) {
        std::cout << "-";
    }
    std::cout << std::endl;

    for (int n = turan.start();; ++n) {
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        const auto [ex, e] = turan.compute(n);
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        std::string result = "ex(" + std::to_string(n) + ", " + graph_name + ") = " + std::to_string(ex) + "/" + std::to_string(e);
        while (result.size() < 32) {
            result.push_back(' ');
        }
        std::cout << result;
        if (VERBOSE) {
            std::chrono::duration<float> difference = end - start;
            std::cout << difference.count() << "s";
        }
        std::cout << std::endl;
    }
    return 0;
}
