#include "json.hpp"

#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#define PI 3.14159265359

using json = nlohmann::json;


class Metropolitan {
    typedef int VertIndex;
    typedef std::vector<VertIndex> IndexGroup;
    struct Vertex {
        double x, y;
    };

    struct IndeciesEdge {
        VertIndex p1, p2;
    };

    struct RefPoint {
        double orient; // x<0.5 => horizontal; x>0.5 => verticle;
        double ref;
        double min;
        double max;
        double extraCoordinate;
    };

    json metropolitan;
    size_t numRefPoints;
    std::vector<Vertex> vertices;
    std::vector<IndeciesEdge> edges;
    std::vector<RefPoint> refPoints;
    std::vector<IndexGroup> indGroups;
    std::vector<IndeciesEdge> filledObs;

    void parseJSON();
    void resetVertices();
    void _parseEdges();
    void _parseRefPointoints();
    void _parseRefGroups();
    void _parseFilledObstacles();

    public:
    std::vector<double> parameters;
    Metropolitan(std::string filename);
    void randomParameter();
    void applyParameter();
    // Func f(p0.x, p0.y, p1.x, p1.y);
    template<typename Func> void lineApply(Func f);
    // adapting SteerSuite
    // Func f(centerX, centerZ, lengthX, theta)
    template<typename Func> void orientedObstacleApply(Func f);
    // Func f(xmin, xmax, ymin, ymax);
    template<typename Func> void obstacleBoxApply(Func f);
};



Metropolitan::Metropolitan(std::string filename){
    std::ifstream file(filename);
    if(!file.is_open()) {
        throw std::runtime_error("not a file");
    } else {
        file >> metropolitan;
    }
    parseJSON();
}

void Metropolitan::parseJSON() {
    resetVertices();
    _parseEdges();
    _parseRefPointoints();
    _parseRefGroups();
    _parseFilledObstacles();
}

void Metropolitan::resetVertices() {
    vertices.clear();
    const json& jVerts = metropolitan["vertexes"];
    for (auto it=jVerts.begin();
            it < jVerts.end();
            it += 2) {
        vertices.push_back(Vertex{(double) *it, (double)*(it+1)});
    }
}

void Metropolitan::_parseEdges() {
    const json& jEdges = metropolitan["edges"];
    for (auto it=jEdges.begin();
            it < jEdges.end();
            it += 2) {
        edges.push_back(IndeciesEdge{
                (int) *it,
                (int) *(it+1)
                });
    }
}

void Metropolitan::_parseRefPointoints() {
    const json& jRefPoints = metropolitan["refPoints"];
    for (auto it=jRefPoints.begin();
            it < jRefPoints.end();
            it += 5) {
        refPoints.push_back(RefPoint{
                (double) *it,
                (double) *(it+1),
                (double) *(it+2),
                (double) *(it+3),
                (double) *(it+4)
                });
    }
    numRefPoints = refPoints.size();
}

void Metropolitan::_parseRefGroups() {
    const json& jRefGs = metropolitan["selectedVertices"];
    for (auto& group: jRefGs){
        std::vector<int> refG;
        for (auto& vert: group) {
            refG.push_back((int) vert);
        }
        indGroups.push_back(refG);
    }
    assert(indGroups.size() == numRefPoints);
}

void Metropolitan::_parseFilledObstacles() {
    const json& jFilledObs = metropolitan["filledObs"];
    for (auto it=jFilledObs.begin();
            it < jFilledObs.end();
            it += 2) {
        filledObs.push_back(IndeciesEdge{
                (int)*it,
                (int)*(it+1)
                });
    }
}

void Metropolitan::randomParameter() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> rand(0., 1.);

    for (size_t i=0; i<numRefPoints; i++) {
        RefPoint& refP = refPoints[i];
        double param = refP.min + rand(gen) * (refP.max - refP.min);
        parameters.push_back(param);
    }
}

void Metropolitan::applyParameter() {
    resetVertices();
    for (size_t i=0; i<numRefPoints; i++) {
        RefPoint& refP = refPoints[i];
        double diff = parameters[i] - refP.ref;
        if (refP.orient < 0.5) {
            for (auto vIdx: indGroups[i]) {
                vertices[vIdx].x += diff;
            }
        } else {
            for (auto vIdx: indGroups[i]) {
                vertices[vIdx].y += diff;
            }
        }
    }
}

template<typename Func>
void Metropolitan::lineApply(Func f) {
    for (auto& e: edges) {
        const Vertex &v1=vertices[e.p1];
        const Vertex &v2=vertices[e.p2];
        f(v1.x, v1.y, v2.x, v2.y);
    }
}

template<typename Func>
void Metropolitan::obstacleBoxApply(Func f) {
    for (auto& e: edges) {
        const Vertex &v1=vertices[e.p1];
        const Vertex &v2=vertices[e.p2];
        double xmin = std::min(v1.x, v2.x);
        double xmax = std::max(v1.x, v2.x);
        double ymin = std::min(v1.y, v2.y);
        double ymax = std::max(v1.y, v2.y);
        f(xmin, xmax, ymin, ymax);
    }
}

template<typename Func>
void Metropolitan::orientedObstacleApply(Func f) {
    auto warp = [f](double x0, double y0, double x1, double y1){
        double centerX = (x0 + x1)/2;
        double centerZ = (y0 + y1)/2;
        double diffX = (x0 - x1);
        double diffY = (y0 - y1);
        double lengthX = std::sqrt(diffX*diffX + diffY*diffY);
        double theta = std::atan2(diffX, diffY)*180/PI;
        f(centerX, centerZ, lengthX, theta);
    };
    lineApply(warp);
}

int main()
{
    Metropolitan m("j_metropolitan.json");
    m.randomParameter();
    m.applyParameter();
    m.lineApply([](double a, double b, double c, double d){
            std::cout << "vert" << a << " " << b << " " << c << " " << d << std::endl;
            });
    return 0;
}
