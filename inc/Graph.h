#pragma once

#include "Structure.h"
#include <fstream>
#include <iostream>

// a class representing simple graph
class Graph : public Structure {
public:
    Graph();
    explicit Graph(size_t n);
    Graph(size_t n, const Certificate& cert);

    size_t edges() const;
    bool edge(size_t i, size_t j) const;
    size_t deg() const;
    void clear();
    bool subClique(size_t k) const;
    void resize(size_t m);
    std::vector<size_t> getDegrees() const;
    void addEdge(size_t i, size_t j);
    void killEdge(size_t i, size_t j);

    friend Graph operator+(const Graph& G, const Graph& H);
    friend Graph operator+(const Graph& G, size_t m);
    friend Graph operator+(size_t m, const Graph& G);
    friend bool readGraph(Graph& G);
    static size_t certSize(size_t n);

protected:
    size_t e;
    std::vector<byte> A;
    bool nextS(int level, Perm& Q) const;

    virtual size_t degsize() const override;
    virtual int color(size_t i, size_t j) const override;
    virtual int compareOrders(const Perm& F, const Perm& B, size_t p, size_t q) const override;
    virtual Certificate getCertificate(const Perm& P) const override;
};

bool readGraph(std::fstream&, Graph& G);
bool readSparse6(std::fstream& stream, Graph& G);
void writeSparse6(std::fstream& stream, const Graph& G, bool header = false);

Graph operator+(const Graph& G, size_t m);
Graph operator+(size_t m, const Graph& G);
Graph operator+(const Graph& G, const Graph& H);

// some special types of graphs
Graph K(size_t n);
Graph K(size_t n, size_t m);
Graph C(size_t n);
Graph P(size_t n);
Graph Q(size_t n);

// a class for a collection of pairwise non-isomorphic
// simple graphs of same size
class GraphSet : public StructSet {
public:
    GraphSet() = default;
    GraphSet(size_t n) : n(n) {
    }

    void resize(int m) {
        n = m;
    }

    std::vector<Certificate> getList() const {
        return std::vector<Certificate>(data_.begin(), data_.end());
    }
    
    void writeSparse6(int n, const std::string& path) const;

private:
    size_t n;
};

class GraphHashSet : public StructHashSet {
public:
    GraphHashSet() = default;
    GraphHashSet(size_t n) : n(n) {
    }

    void resize(int m) {
        n = m;
    }

    //std::vector<Certificate> getList() const {
    //    return std::vector<Certificate>(data_.begin(), data_.end());
    //}

private:
    size_t n;
};
