#include "Graph.h"

#include <algorithm>
#include <string>
#include <vector>

Graph::Graph() : Structure(0), e(0) {
}

Graph::Graph(size_t n) : Structure(n), A(n * n, 0), e(0) {
}

Graph::Graph(size_t n, const Certificate& cert) : Structure(n, cert), A(n * n, 0), e(0) {
    size_t l = n * (n - 1) / 2;
    if (l % 8 == 0) {
        l >>= 3;
    } else {
        l >>= 3;
        l++;
    }
    // cert must be of length l
    size_t ii = 0;
    size_t jj = 1;
	
    for (size_t i = 0; i < l; i++) {
        byte b = static_cast<byte>(cert[i]);
        for (int j = 7; j >= 0; j--) {
            byte c = (b >> j) & 1;

            if (c) {                
                addEdge(ii, jj);
            }

            jj++;

            if (jj >= n) {
                ii++;
                jj = ii + 1;
            }
	    if (ii >= n - 1) {
		return;
            }
        }
    }
}

bool Graph::edge(size_t i, size_t j) const {
    return A[n * i + j] != 0;
}

void Graph::resize(size_t m) {
    n = m;
    A.assign(m * m, 0);
    e = 0;
}

std::vector<size_t> Graph::getDegrees() const {
    std::vector<size_t> degrees(n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (A[i * n + j]) {
                degrees[i]++;
            }
        }
    }
    return degrees;
}

void Graph::addEdge(size_t i, size_t j) {
    byte c = A[n * i + j];
    if (c == 0) {
        A[n * i + j] = 1;
        A[n * j + i] = 1;
        e++;
    }
}

void Graph::killEdge(size_t i, size_t j) {
    byte c = A[n * i + j];
    if (c == 1) {
        A[n * i + j] = 0;
        A[n * j + i] = 0;
        e--;
    }
}

size_t Graph::degsize() const {
    return 1;
}

int Graph::color(size_t i, size_t j) const {
    return static_cast<int>(A[n * i + j]);
}

size_t Graph::edges() const {
    return e;
}

int Graph::compareOrders(const Perm& F, const Perm& B, size_t p, size_t q) const {
    for (size_t i = p; i < q; i++) {
        for (size_t j = 0; j < i; j++) {
             if (A[n * F[i] + F[j]] > A[n * B[i] + B[j]]) {
                 return 1;
             } else {
                 if (A[n * F[i] + F[j]] < A[n * B[i] + B[j]]) {
                     return -1;
                 }
             }
        }
    }
    return 0;
}

Certificate Graph::getCertificate(const Perm& P) const {
    size_t l = certSize(n);
    Certificate C(l);

    size_t q = 0;
    byte b = 0;

    for (size_t i = 0; i + 1 < n; i++) {
        for (size_t j = i + 1; j < n; j++) {			
            b = b * 2 + A[n * P[i] + P[j]];
            q++;
            if (q % 8 == 0) {
                C[q / 8 - 1] = b;
                b = 0;
            }
        }
    }
    if (q < 8 * l) {
        while (q < 8 * l) {
            q++;
            b = b * 2;		
        }
        C[q / 8 - 1] = b;
    }
    return C;
}

size_t Graph::deg() const {
    if (n == 0) {
        return 0;
    }

    size_t d = n - 1;
    for (size_t i = 0; i < n; i++) {
        size_t dd = 0;
        for (size_t j = 0; j < n; j++) {
            if (edge(i, j)) {
                dd++;
                if (dd >= d) {
                    break;
                }
            }
        }
        if (dd < d) {
            d = dd;
        }
        if (dd == 0) {
           return 0;
        }
    }
    return d;
}

void Graph::clear() {
    e = 0;
    A.assign(n * n , 0);
}

bool Graph::subClique(size_t k) const {
    Perm Q(k);
    for (size_t i = 0; i + k <= n; i++) {
        Q[0] = i;
        if (nextS(1, Q)) {
            return true;
        }
    }
    return false;
}

bool Graph::nextS(int level, Perm& Q) const {
    if (level >= Q.size()) {
        return true;
    }

    if (level < Q.size()) {
        for (size_t i = Q[level - 1] + 1; i < n; i++) {
            bool B = true;
            for (size_t j = 0; j < level; j++) {
                if (edge(Q[j], i)) {
                    B = false;
                }
            }
            if (B) {
                Q[level] = i;
                if (nextS(level + 1, Q)) {
                    return true;
                }
            }           
        }
    } 
    return false;
}

bool readGraph(std::fstream& stream, Graph& G) {
    size_t n = G.size();
    size_t l = Graph::certSize(n);
    Certificate cert(l);
    Structure::readStruct(stream, cert);

    if (stream.eof()) {
        return false;
    }

    G = Graph(n, cert);
    return true;
}

bool readSparse6(std::fstream& stream, Graph& G) {
    std::string line;
    while (std::getline(stream, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.size() >= 11 && line.compare(0, 11, ">>sparse6<<") == 0) {
            line.erase(0, 11);
        }
        if (line.empty() || line[0] != ':') {
            return false;
        }

        std::vector<int> chars;
        chars.reserve(line.size());
        for (size_t i = 1; i < line.size(); ++i) {
            const unsigned char c = static_cast<unsigned char>(line[i]);
            if (c < 63 || c > 126) {
                return false;
            }
            chars.push_back(static_cast<int>(c - 63));
        }

        size_t pos = 0;
        size_t n = 0;
        if (pos >= chars.size()) {
            return false;
        }
        if (chars[pos] <= 62) {
            n = static_cast<size_t>(chars[pos++]);
        } else if (pos + 4 <= chars.size() && chars[pos + 1] <= 62) {
            n = (static_cast<size_t>(chars[pos + 1]) << 12)
                + (static_cast<size_t>(chars[pos + 2]) << 6)
                + static_cast<size_t>(chars[pos + 3]);
            pos += 4;
        } else if (pos + 8 <= chars.size()) {
            n = (static_cast<size_t>(chars[pos + 2]) << 30)
                + (static_cast<size_t>(chars[pos + 3]) << 24)
                + (static_cast<size_t>(chars[pos + 4]) << 18)
                + (static_cast<size_t>(chars[pos + 5]) << 12)
                + (static_cast<size_t>(chars[pos + 6]) << 6)
                + static_cast<size_t>(chars[pos + 7]);
            pos += 8;
        } else {
            return false;
        }

        G = Graph(n);
        if (n == 0) {
            return true;
        }

        size_t k = 1;
        while ((size_t{1} << k) < n) {
            ++k;
        }

        size_t data_pos = pos;
        int d = 0;
        int d_len = 0;
        size_t v = 0;
        while (true) {
            if (d_len < 1) {
                if (data_pos >= chars.size()) {
                    break;
                }
                d = chars[data_pos++];
                d_len = 6;
            }
            --d_len;
            int b = (d >> d_len) & 1;

            int x = d & ((1 << d_len) - 1);
            int x_len = d_len;
            while (x_len < static_cast<int>(k)) {
                if (data_pos >= chars.size()) {
                    return true;
                }
                d = chars[data_pos++];
                d_len = 6;
                x = (x << 6) + d;
                x_len += 6;
            }
            x = x >> (x_len - static_cast<int>(k));
            d_len = x_len - static_cast<int>(k);

            if (b == 1) {
                ++v;
            }
            if (static_cast<size_t>(x) >= n || v >= n) {
                break;
            }
            if (static_cast<size_t>(x) > v) {
                v = static_cast<size_t>(x);
            } else {
                G.addEdge(static_cast<size_t>(x), v);
            }
        }
        return true;
    }
    return false;
}

void writeSparse6(std::fstream& stream, const Graph& G, bool header) {
    if (header) {
        stream << ">>sparse6<<";
    }

    const size_t n = G.size();
    stream.put(':');

    std::vector<int> nd;
    if (n <= 62) {
        nd.push_back(static_cast<int>(n));
    } else if (n <= 258047) {
        nd = {
            63,
            static_cast<int>((n >> 12) & 0x3F),
            static_cast<int>((n >> 6) & 0x3F),
            static_cast<int>(n & 0x3F),
        };
    } else {
        nd = {
            63,
            63,
            static_cast<int>((n >> 30) & 0x3F),
            static_cast<int>((n >> 24) & 0x3F),
            static_cast<int>((n >> 18) & 0x3F),
            static_cast<int>((n >> 12) & 0x3F),
            static_cast<int>((n >> 6) & 0x3F),
            static_cast<int>(n & 0x3F),
        };
    }
    for (int d : nd) {
        stream.put(static_cast<char>(d + 63));
    }

    if (n == 0) {
        stream << '\n';
        return;
    }

    size_t k = 1;
    while ((size_t{1} << k) < n) {
        ++k;
    }

    std::vector<std::pair<size_t, size_t>> edges;
    edges.reserve(G.edges());
    for (size_t u = 0; u + 1 < n; ++u) {
        for (size_t v = u + 1; v < n; ++v) {
            if (G.edge(u, v)) {
                edges.emplace_back(v, u);
            }
        }
    }
    std::sort(edges.begin(), edges.end());

    std::vector<int> bits;
    size_t curv = 0;
    for (const auto& edge : edges) {
        const size_t ev = edge.first;
        const size_t eu = edge.second;
        if (ev == curv) {
            bits.push_back(0);
            for (size_t i = 0; i < k; ++i) {
                bits.push_back((eu >> (k - 1 - i)) & 1);
            }
        } else if (ev == curv + 1) {
            ++curv;
            bits.push_back(1);
            for (size_t i = 0; i < k; ++i) {
                bits.push_back((eu >> (k - 1 - i)) & 1);
            }
        } else {
            curv = ev;
            bits.push_back(1);
            for (size_t i = 0; i < k; ++i) {
                bits.push_back((ev >> (k - 1 - i)) & 1);
            }
            bits.push_back(0);
            for (size_t i = 0; i < k; ++i) {
                bits.push_back((eu >> (k - 1 - i)) & 1);
            }
        }
    }

    const size_t pad = bits.size() % 6 == 0 ? 0 : 6 - (bits.size() % 6);
    if (k < 6 && n == (size_t{1} << k) && pad >= k && curv < n - 1) {
        bits.push_back(0);
        bits.insert(bits.end(), pad, 1);
    } else {
        bits.insert(bits.end(), pad, 1);
    }

    for (size_t i = 0; i < bits.size(); i += 6) {
        int byte = 0;
        for (size_t j = 0; j < 6; ++j) {
            byte = (byte << 1) + bits[i + j];
        }
        stream.put(static_cast<char>(byte + 63));
    }
    stream << '\n';
}

size_t Graph::certSize(size_t n) {
    size_t l = n * (n - 1) / 2;
    if (l % 8 == 0) {
        l >>= 3;
    } else {
        l >>= 3;
        l++;
    }
    if (l == 0) {
        l = 1;
    }
    return l;
}

Graph operator+(const Graph& G, size_t m) {
    size_t n = G.n;
    Graph H(m + n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (G.edge(i, j)) {
                H.addEdge(i, j);
            }
        }
    }
    return H;
}

Graph operator+(size_t m, const Graph& G) {
    size_t n = G.n;
    Graph H(m + n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (G.edge(i, j)) {
                H.addEdge(i + m, j + m);
            }
        }
    }
    return H;
}

Graph operator+(const Graph& G, const Graph& H) {
    size_t n = G.n;
    size_t m = H.n;

    Graph F(n + m);

    for (size_t i = 0; i + 1 < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (G.edge(i, j)) {
                F.addEdge(i, j);
            }
        }
    }

    for (size_t i = 0; i + 1 < m; i++) {
        for (size_t j = i + 1; j < m; j++) {
            if (H.edge(i, j)) {
                F.addEdge(n + i, n + j);
            }
        }
    }
    return F;
}

Graph K(size_t n) {
    if (n <= 1) {
        return Graph(1);
    }
    Graph G(n);
    for (size_t i = 0; i + 1 < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            G.addEdge(i, j);
        }
    }
    return G;
}

Graph K(size_t n, size_t m) {
    Graph G(n + m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = n; j < n + m; j++) {
            G.addEdge(i, j);
        }
    }
    return G;
}

Graph C(size_t n) {
    Graph G(n);
    for (size_t i = 0; i + 1 < n; i++) {
        G.addEdge(i, i + 1);
    }
    G.addEdge(0, n - 1);
    return G;
}

Graph P(size_t n) {
    Graph G(n);
    for (size_t i = 0; i + 1 < n; i++) {
        G.addEdge(i, i + 1);
    }
    return G;
}

Graph Q(size_t n) {
    if (n == 0) {
        return Graph(1);
    }
    Graph G = Q(n - 1) + Q(n - 1);
    size_t m = G.size() >> 1;
    for (size_t i = 0; i < m; i++) {
        G.addEdge(i, i + m);
    }
    return G;
}
