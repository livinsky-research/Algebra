#pragma once

#include <fstream>
#include <unordered_set>
#include <functional>
#include <vector>
#include <string>
#include <memory>
#include <list>
#include <mutex>
#include <atomic>

#include "Group.h"
#include "Certificate.h"

typedef uint8_t byte;
typedef std::vector<int> Degree;

struct Cell : public Perm {
    Cell();
    explicit Cell(size_t m);
    Cell(const Cell& W, int s, int m);
    void sort(const std::function<int(int,int)>& comp);

    bool counted;
    bool discrete;
};

typedef std::list<Cell> Part;

class SearchNode;
// abstract class containing all algebraic structures such as
// graphs, digraphs, hypergraphs, semigroups, posets, lattices
class Structure {
friend class SearchNode;
friend class StructSet;
friend class StructHashSet;
friend class Certifier;
public:
    explicit Structure(size_t n);
    Structure(size_t n, const Certificate& cert);
    size_t size() const;
    Structure& certify();
    Group aut();

    static void writeStruct(std::fstream&, const Certificate& cert);
    static void readStruct(std::fstream&, Certificate& cert);

protected:
    virtual int compareOrders(const Perm& P, const Perm& Q, size_t p, size_t q) const = 0;
    virtual size_t degsize() const = 0;
    virtual int color(size_t i, size_t j) const = 0;
    virtual Certificate getCertificate(const Perm& P) const = 0;

    friend bool isomorphic(const Structure& s, const Structure& t);

    size_t n;
    Certificate cert;
    std::shared_ptr<Group> auto_group;
};

class Certifier {
friend class Structure;
friend class SearchNode;
public:
    explicit Certifier(const Structure* S);
    ~Certifier();

private:	
    Perm F;
    Perm B;
    std::vector<Degree> degrees;
    const Structure* S;
    bool AutoFound;
    bool Bexists;
    int BasisOK;
    bool IsDiscrete;
    SearchNode* Top;
    SearchNode* LastBaseChange;

    bool Compare(int x, int y);
};

class SearchNode {
friend class Certifier;
friend class Structure;
public:
    SearchNode(size_t n, Certifier* crt);
    ~SearchNode();
	
    int orbitRep(size_t v);
    void merge(size_t u, size_t v);
    void updateOrbits(const Perm& Q);
    void addGen(const Perm& Q);
    void refine();
    void stabilise();
    void changeBase(int d);

private:
    Part P;
    int FixedPoint;
    std::shared_ptr<Group> G;
    Perm CellOrbits;
    size_t Depth;
    size_t NFixed;
    bool OnBestPath;
    SearchNode* Next;
    Certifier* crt;
};

// class modelling an unordered collection of non-isomorphic structures. 
class StructSet {
public:
    void insert(const Structure& s);
    size_t size() const;
    bool empty() const;
    void write(const std::string& path, bool append = false) const;
    void clear();
    bool contains(const Structure& s) const;

protected:
    std::unordered_set<Certificate> data_;
    mutable std::mutex mut_;
};

// class modelling an unordered collection of non-isomorphic structures. 
class StructHashSet {
public:
    StructHashSet();
    void extend(const Structure& s);
    bool empty() const;
    void write(const std::string& path, bool append = false) const;
    void clear();
    bool contains(const Structure& s) const;
    
protected:
    std::vector<std::vector<Certificate>> data_;
    //mutable std::vector<std::mutex> mut_;
    mutable std::mutex mut_[10000];
    std::atomic<bool> empty_;
    const size_t bsize = 10000;
};
