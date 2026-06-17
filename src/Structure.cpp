#include <algorithm>
#include <iostream>

#include "Structure.h"

Cell::Cell() : Perm(), counted(false), discrete(false) {
};

Cell::Cell(size_t m): Perm(m), counted(false), discrete(false) {
};

Cell::Cell(const Cell& W, int s, int m): counted(false), discrete(false), Perm(m) {
    for (size_t i = 0; i < m; ++i) {
        data_[i] = W[s + i];
    }
};

void Cell::sort(const std::function<int(int,int)>& comp) {
    std::sort(data_, data_ + size_, comp);
}

Structure::Structure(size_t n): n(n), auto_group(nullptr) {
};

Structure::Structure(size_t n, const Certificate& cert): n(n), 
    cert(cert), auto_group(nullptr) {
};

size_t Structure::size() const {
    return n;
}

void Structure::writeStruct(std::fstream& stream, const Certificate& cert) {
    stream.write((char*)cert.data(), cert.size());
}

void Structure::readStruct(std::fstream& stream, Certificate& cert) {
    stream.read((char*)cert.data(), cert.size());
}

bool isomorphic(const Structure& s, const Structure& t) {
    if (s.n != t.n) {
        return false;
    }

    return 0 == compareCertificates(s.cert, t.cert);
}

Certifier::Certifier(const Structure* S) : S(S), degrees(S->n) {
    size_t n = S->n;
    Top = new SearchNode(n, this);
    B.id(n);	
    F.id(n);
    Bexists = false;
    BasisOK = 0;
    AutoFound = false;
    LastBaseChange = Top;

    Top->G = std::make_shared<Group>(n);

    size_t s = S->degsize();
    for (size_t i = 0; i < n; i++) {
        degrees[i].assign(s, 0);
    }
	
    Cell C(n);

    Top->P.push_back(C);
    Top->NFixed = 0;
    Top->Depth = 0;
    Top->CellOrbits.id(n);
    Top->stabilise();
}

Certifier::~Certifier() {
    delete Top;
}

Structure& Structure::certify() {
    Certifier certifier(this);
    cert = getCertificate(certifier.B);
    auto_group = std::move(certifier.Top->G);
    return *this;
}

Group Structure::aut() {
    if (!auto_group) {
        certify();
    }
    return *auto_group;
}

SearchNode::SearchNode(size_t n, Certifier* crt) : G(nullptr), Next(nullptr), OnBestPath(false), CellOrbits(n), crt(crt) {
};

SearchNode::~SearchNode() {
    delete Next;
};

// comparing for sort. Must return true if node x goes before y
bool Certifier::Compare(int x, int y) {
    size_t s = S->degsize();
    if (degrees[y].size() < s) {
        return false;
    }
    if (degrees[y].size() > s) {
        return true;
    }

    size_t i = 0;
    while (i < s && degrees[x][i] == degrees[y][i]) {
        i++;
    }

    if (i >= s) {
        return false;
    }

    if (degrees[x][i] < degrees[y][i]) {
        return true;
    }
    return false;
}

int SearchNode::orbitRep(size_t v) {
    if (CellOrbits[v] < 0) {
        return v;
    }
    int w = orbitRep(CellOrbits[v]);
    CellOrbits[v] = w;
    return w;
}

void SearchNode::merge(size_t u, size_t v) {
    int u_size = -CellOrbits[u];
    int v_size = -CellOrbits[v];

    int w;
    if (u_size < v_size) {
        CellOrbits[u] = v;
        w = v;
    } else {
        CellOrbits[v] = u;
        w = u;
    }
    CellOrbits[w] = -u_size - v_size;
}

void SearchNode::updateOrbits(const Perm& Q) {
    if (P.empty()) {
        return;
    }

    const Cell& C = P.front();
    for (size_t i = 0; i < C.size(); i++) {
        int u = C[i];
        int v = Q[u];
    	
        int uRep = orbitRep(u);
        int vRep = orbitRep(v);		

        if (uRep != vRep) {
            merge(uRep, vRep);
        }
    }
}

void SearchNode::addGen(const Perm& P) {
    if (!G->Gu) {
        G->Gu = std::make_shared<Group>(G->n);
    }
    if (Next != nullptr) {
        Next->G = G->Gu;
    }

    if (G->u == -1) {
        if (FixedPoint < G->n && FixedPoint >= 0) {
            G->u = FixedPoint;
        } else {
            G->u = 0;
            while (P[G->u] == G->u) {
                G->u++;
            }
        }
		
        G->Orbit[0] = G->u;
        G->Cosets[G->u].id(G->n);	
    }

    G->Generators.push_back(P);

    if (!CellOrbits.isConst(0))	{
        updateOrbits(P);
    }

    size_t M = G->NPoints;
    size_t k = 0;
    while (k < M) {
        int v = G->Orbit[k];
        int w = P[v];

        if (G->Cosets[w].empty()) {			
            G->Orbit[G->NPoints] = w;
            G->NPoints++;

            G->Cosets[w] = P * G->Cosets[v];
            G->Inverses[w] = !G->Cosets[w];
        } else {
            if (G->Inverses[w].empty())	{
                G->Inverses[w] = !G->Cosets[w];
            }

            Perm Q = Mult(G->Inverses[w], P, G->Cosets[v]);
            if (!G->Gu->contains(Q)) {
                if (Next) {
                    Next->addGen(Q);
                } else {
                    G->Gu->addGen(Q);
                }
            }
        }
	k++;
    }
    // apply new generators to all points
    while (k < G->NPoints) {
        int v = G->Orbit[k];		
        for (const Perm& Gen : G->Generators) {
            int w = Gen[v];
            if (G->Cosets[w].empty()) {
                G->Orbit[G->NPoints] = w;
                G->NPoints++;
                G->Cosets[w] = Gen * G->Cosets[v];
                G->Inverses[w] = !G->Cosets[w];
            } else {
                if (G->Inverses[w].empty()) {
                    G->Inverses[w] = !G->Cosets[w];
                }
                Perm Q = Mult(G->Inverses[w], Gen, G->Cosets[v]);
                if (!G->Gu->contains(Q)) {
                    if (Next) {
                        Next->addGen(Q);
                    } else {
                        G->Gu->addGen(Q);
                    }
                }
            }			
        }
        k++;
    }
}

void SearchNode::changeBase(int d) {
    SearchNode* node = crt->LastBaseChange;
    std::shared_ptr<Group> G = node->G;

    if (!G) {
        return;
    }
    if (!G->Gu) {
        return;
    }

    while (node) {
        if (node->Depth <= d) {
            CellOrbits = 0;
            const Cell& C = node->P.front();			
            for (size_t i = 0; i < C.size(); i++) {
                node->CellOrbits[C[i]] = -1;
            }
        }		
        node = node->Next;
    }
    node = crt->LastBaseChange;
    PermList Geners = G->Generators;

    size_t n = G->n;
		
    // clearing the stabilizer tower from LastBaseChange
    // it is important not to delete any groups for avoiding memory leak

    std::shared_ptr<Group> GG = G;
    while (GG) {
        for (size_t i = 0; i < n; i++) {
            GG->Cosets[i].clear();
            GG->Inverses[i].clear();
        }

        GG->Generators.clear();

        GG->u = -1;			
        GG->NPoints = 1;
        GG = GG->Gu;
    }	
		
    // backwards
    for (const Perm& P : Geners) {
        crt->LastBaseChange->addGen(P);
    }
}

void SearchNode::refine() {
    // chosing a first non-discrete cell to refine
    for (auto it = P.begin(); it != P.end(); ++it) {
        it->counted = false;
        it->discrete = false;
    }

    // if this partition is discrete, it is refined already	

    // if not discrete
    bool Stab = false;
    std::vector<Degree>& degrees = crt->degrees;

    size_t s = crt->S->degsize();
    size_t n = crt->S->n;
		
    for (size_t i = 0; i < n; i++) {
        degrees[i].assign(s, 0);
    }
		
    auto c = P.begin();	

    do { // repeat until we get a stable partition
        const Cell& C = *c;		
        for (size_t i = 0; i < C.size(); i++) {			
            for (size_t j = 0; j < n; j++) {
                int col = crt->S->color(C[i], j);
                if (col) {
                    degrees[j][col - 1]++;
                }
            }
        }

        c->counted = true;
        Part PP;

        for (Cell CC : P) {
	    // first we treat discrete cells by adding them to F array
	    if (CC.size() == 1 || CC.discrete) { //
                for (size_t i = 0; i < CC.size(); i++) {
                    crt->F[NFixed++] = CC[i];
                }
                continue;
            }

            // if this is a non-discrete cell					
            CC.sort(std::bind(&Certifier::Compare, crt, std::placeholders::_1, std::placeholders::_2));
            CC.discrete = true;	
            if (CC.size() > 1) {
                for (size_t i = 0; i + 1 < CC.size(); i++) {
                    if (degrees[CC[i]] == degrees[CC[i + 1]]) {
                        CC.discrete = false;
                        break;
                    }		
                }
            }
            // if this cell splits into one-cells, add them to F
            if (CC.discrete) {
                for (size_t i = 0; i < CC.size(); i++) {
                    crt->F[NFixed++] = CC[i];
                }
                continue;
            }

            // if there is still one cell
            if (degrees[CC[0]] == degrees[CC[CC.size() - 1]]) {
                PP.push_back(std::move(CC));
                continue;
            }
		
            // splitting C into new cells
            int l = 0;
            for (size_t i = 1; i < CC.size(); i++) {
                if (degrees[CC[i]] != degrees[CC[i - 1]]) {
                    // add a new cell					
                    PP.emplace_back(CC, l, i - l);
                    l = i;
                }
            }
            PP.emplace_back(CC, l, CC.size() - l);
        }

        P = PP;
        // evaluating wether P is a stable partition
        Stab = true;
        crt->IsDiscrete = true;

        c = P.begin();
        while (c != P.end()) {
            if (!c->counted) {
                Stab = false;				
                break;
            } else {
                ++c;
            }
        }

        crt->IsDiscrete = true;
        for (auto it = P.begin(); it != P.end(); ++it) {
            if (it->size() > 1 && it->discrete == false) {
                crt->IsDiscrete = false;
                break;
            }
        }
        if (crt->IsDiscrete) {
            Stab = true;
        }
    } while(!Stab);

    crt->IsDiscrete = true;
    for (auto it = P.begin(); it != P.end(); ++it) {
        if (!it->discrete) {
            crt->IsDiscrete = false;
            break;
        }
    }
}

void SearchNode::stabilise() {
    size_t m = NFixed;
    OnBestPath = false;

    refine();
    int res = 1;
    if (crt->Bexists) {
        res = crt->S->compareOrders(crt->F, crt->B, m, NFixed);
    }

    size_t n = crt->S->n;

    if (crt->IsDiscrete) {
        if (crt->Bexists) {
            if (res == 0) { // this means that we have afound an automorphism
                Perm Q(n);
                for (size_t i = 0; i < n; i++) {
                    Q[crt->F[i]] = crt->B[i];
                }
                if (!crt->Top->G->contains(Q)) {
                    crt->Top->addGen(Q);
                    crt->AutoFound = true;
                }
            } else if (res == 1) { // if this ordering is better
                crt->B = crt->F;
                SearchNode* Node = crt->Top;
                while (Node != nullptr) {
                    Node->OnBestPath = true;
                    Node = Node->Next;
                }
            }
        } else {
            crt->B = crt->F;			
            SearchNode* Node = crt->Top;
            while (Node != nullptr) {
                Node->OnBestPath = true;
                Node = Node->Next;
            }
            crt->Bexists = true;
        }
        goto Finish;
    } else {
        if (res == 1) {
            crt->Bexists = false;
        } else if (res == -1) {
            goto Finish;
        }
        // we get here only if result is 0
        if (Next == nullptr) {
            Next = new SearchNode(n, crt);
            Next->Depth = Depth + 1;
            if (G) {
                Next->G = G->Gu;
            }
        }

        CellOrbits = 0;		
        SearchNode* Su = Next;

        //CellOrbits = new Perm(n);
        const Cell& C = P.front();
        for (size_t i = 0; i < C.size(); i++) {
            CellOrbits[C[i]] = -1;
        }
        int u;
        size_t jj = 0;
				
        while (jj < C.size()) {			
            u = C[jj];
            FixedPoint = u;
            // splitting the first cell into {u}{****}
            Part Pu;
            Cell C1(1);
            C1[0] = u;

            Pu.push_back(std::move(C1));
            Cell C2(C.size() - 1);
            for (size_t j = 0; j < jj; j++) {
                C2[j] = C[j];
            }
            for (size_t j = jj + 1; j < C.size(); j++) {
                C2[j - 1] = C[j];
            }
            Pu.push_back(std::move(C2));
					
            //for (int j = 1; j < P.size(); j++)
            auto it = P.begin();
            ++it;
            while (it != P.end()) {
                Pu.push_back(*it);
                ++it;
            }

            if (Depth > crt->BasisOK) {
                changeBase(Depth);
            }

            Su->P = std::move(Pu);
            Su->NFixed = NFixed;

            crt->BasisOK = Depth;
            crt->LastBaseChange = this;

            Next->stabilise();

            CellOrbits[orbitRep(u)] -= n;

            if (crt->AutoFound) {
                if (!OnBestPath) {
                    goto Finish;
                }
                crt->AutoFound = false;
            }				
            while (jj < C.size() && CellOrbits[orbitRep(C[jj])] < -n) {
                jj++;
            }
        }		
    }

Finish:
    // zero the degrees array 
    size_t s = crt->S->degsize();
    for (size_t i = m; i < NFixed; i++) {
        crt->degrees[i].assign(s, 0);
        FixedPoint = -1;
        P.clear();
    }
}


void StructSet::insert(const Structure& s) {
    std::lock_guard<std::mutex> lock(mut_);
    data_.insert(s.cert);
}

void StructSet::write(const std::string& path, bool append) const {
    std::lock_guard<std::mutex> lock(mut_);
    std::fstream stream;
    if (append) {
        stream.open(path, std::ios::app | std::ios::out | std::ios::binary);
    } else {
        stream.open(path, std::ios::out | std::ios::binary);
    }
    for (const Certificate& cert : data_) {
        Structure::writeStruct(stream, cert);
    }
    stream.close();
}

size_t StructSet::size() const {
    std::lock_guard<std::mutex> lock(mut_);
    return data_.size();
}

bool StructSet::empty() const {
    std::lock_guard<std::mutex> lock(mut_);
    return data_.empty();
}

void StructSet::clear() {
    std::lock_guard<std::mutex> lock(mut_);
    data_.clear();
}

bool StructSet::contains(const Structure& s) const {
    std::lock_guard<std::mutex> lock(mut_);
    return data_.count(s.cert);
}


StructHashSet::StructHashSet() : data_(bsize), empty_(true) {
}

void StructHashSet::extend(const Structure& s) {
    empty_ = false;
    int b = std::hash<Certificate>()(s.cert) % bsize;
    std::mutex& m = mut_[b];
    //std::lock_guard<std::mutex> lock(mut_[b]);
    m.lock();
    for (const Certificate& cert : data_[b]) {
        if (cert == s.cert) {
            m.unlock();
            return;
        }
    }   
    data_[b].emplace_back(s.cert);
    m.unlock();
}

bool StructHashSet::empty() const {
    return empty_;
}

void StructHashSet::write(const std::string& path, bool append) const {
    std::fstream stream;
    if (append) {
        stream.open(path, std::ios::app | std::ios::out | std::ios::binary);
    } else {
        stream.open(path, std::ios::out | std::ios::binary);
    }
    for (const auto& bucket : data_) {
		for (const Certificate& cert : bucket) {
		    Structure::writeStruct(stream, cert);
		}
    }
    stream.close();
}

void StructHashSet::clear() {
    data_.assign(bsize, std::vector<Certificate>());
    empty_ = true;
}

bool StructHashSet::contains(const Structure& s) const {
    int b = std::hash<Certificate>()(s.cert) % bsize;
    std::lock_guard<std::mutex> lock(mut_[b]);  
    for (const Certificate& cert : data_[b]) {
        if (cert == s.cert) {
            return true;
        }
    }   
    return false;
}
