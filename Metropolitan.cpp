#include "json.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>


using json = nlohmann::json;


class Metropolitan {
    struct Vertex {
        double x, y;
    };

    struct IndeciesEdge {
        int p1, p2;
    };

    struct RefP {
        double orient; // x<0.5 => horizontal; x>0.5 => verticle;
        double Ori;
        double Min;
        double Max;
        double extraCoordinate;
    };

    size_t numRefPs;
    std::vector<Vertex> vertices;
    std::vector<IndeciesEdge> edges;
    std::vector<RefP> refPoints;
    std::vector< std::vector<int> > indGroups;
    std::vector<IndeciesEdge> filledObs;

    void _parseEdges();
    void _parseRefPoints();
    void _parseRefGroups();

    public:
    json met;
    Metropolitan(std::string filename);
    void parseJSON();
    void resetVertices();
    template<typename Func> void lineApply(Func f);
};



Metropolitan::Metropolitan(std::string filename){
    std::ifstream file(filename);
    if(!file.is_open()) {
        throw std::runtime_error("not a file");
    } else {
        file >> met;
    }
    parseJSON();
}

void Metropolitan::parseJSON() {
    resetVertices();
    _parseEdges();
    _parseRefPoints();
    _parseRefGroups();
}

void Metropolitan::resetVertices() {
    vertices.clear();
    const json& jVerts = met["vertexes"];
    for (auto it=jVerts.begin();
            it < jVerts.end();
            it += 2) {
        vertices.push_back(Vertex{(double) *it, (double)*(it+1)});
    }
}

void Metropolitan::_parseEdges() {
    const json& jEdges = met["edges"];
    for (auto it=jEdges.begin();
            it < jEdges.end();
            it += 2) {
        edges.push_back(IndeciesEdge{
                (int) *it,
                (int) *(it+1)
                });
    }
}

void Metropolitan::_parseRefPoints() {
    const json& jRefPs = met["refPoints"];
    for (auto it=jRefPs.begin();
            it < jRefPs.end();
            it += 5) {
        refPoints.push_back(RefP{
                (double) *it,
                (double) *(it+1),
                (double) *(it+2),
                (double) *(it+3),
                (double) *(it+4)
                });
    }
}

void Metropolitan::_parseRefGroups() {
    const json& jRefGs = met["selectedVertices"];
    for (auto& group: jRefGs){
        std::vector<int> refG;
        for (auto& vert: group) {
            refG.push_back((int) vert);
        }
    }
}

int main()
{
    Metropolitan m("j_metropolitan.json");
    return 0;
}
