// Super-cell Grouping Module (k = 4 default)
// Includes logic level computation without needing precomputed topological sort
#include "grouping.hpp"

#include <iostream>
#include <limits>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <queue>
#include <cmath>
#include <algorithm>
#include <fstream>

#include <CGAL/Min_circle_2.h>
#include <CGAL/convex_hull_2.h>

#include "netlist.hpp"
#include "partitioning.hpp"
#include "supercells.hpp"
#include "util.hpp"

using namespace std;

const int GROUP_SIZE_K = 2; // default grouping size

static const int MAX_SEARCH_LEVEL = 2;
static const int NORMALIZATION_FACTOR = MAX_SEARCH_LEVEL;
// TODO: what use for Beta?
static const int DISTANCE_NORMALIZATION_FACTOR = 2;

// Map from level to all gates at that level
map<int, vector<int>> levelToNodes;
unordered_map<int, int> nodeLevelMap;

unordered_map<int, Subgraph*> subgraphs;

/**
 * Makes a subgraph from the main Netlist. Only adds nodes
 * of the specified logic level to the map
 */
static Subgraph* makeSubgraph(Netlist &netlist, int logicLevel) {
    Subgraph *graph = nullptr;

    if (levelToNodes.find(logicLevel) != levelToNodes.end()) {
        graph = new Subgraph(logicLevel);

        // Add verticies from all nodes at this logic level
        for (const int gateID : levelToNodes.at(logicLevel)) {
            SubgraphVertex vert;
            vert.id = gateID;

            // Only work with moveable cells in the subgraph
            if (!netlist.at(gateID).isPrimaryInput && !netlist.at(gateID).isPrimaryOutput) {
                graph->addVertex(vert);
            }
        }
    }

    return graph;
}

static double euclidean_distance(const Point& a, const Point& b)
{
    double dx = a.x() - b.x();
    double dy = a.y() - b.y();
    return std::sqrt(dx * dx + dy * dy);
}

Subgraph::~Subgraph() {
    for (auto &edge : this->edges) {
        delete edge;
    }
}

void Subgraph::calcMinMaxCellDistances(Netlist &netlist) {
    // Find minimum distance between all points
    // TODO: This is brute force method. Update this to a faster
    // version if necessary
    this->minCellDistance = std::numeric_limits<double>::max();

    for (const auto &[id1, v1] : this->graph) {
        if (netlist.at(id1).placement.isPlaced) {
            for (const auto &[id2, v2] : this->graph) {
                NetlistNode &node1 = netlist.at(id1);
                NetlistNode &node2 = netlist.at(id2);

                if (id2 != id1 && node2.placement.isPlaced) {
                    double dist = sqrt(pow(node1.placement.p.x() - node2.placement.p.x(), 2) +
                                       pow(node1.placement.p.y() - node2.placement.p.y(), 2));

                    this->minCellDistance = std::min(this->minCellDistance, dist);
                }
            }
        }
    }

    this->findMaxDistance(netlist);

    // Workaround to division by zero in distance graph processing
    // We need these two values to be different
    if (this->minCellDistance == this->maxCellDistance) {
        // Make this distance slightly bigger
        this->maxCellDistance += 0.01;
    }
}

void Subgraph::findMaxDistance(Netlist &netlist) {
    std::vector<Point> hull;
    std::vector<Point> points;

    // Gather all cells with placement info
    for (const auto &[id, v] : this->graph) {
        if (netlist.at(id).placement.isPlaced) {
            points.push_back(netlist.at(id).placement.p);
        }
    }

    // Finds maximum distance between all points using Convex Hull from CGAL
    // (ChatGPT generated code)

    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(hull));

    double max_dist = 0.0;
    std::pair<Point, Point> farthest_pair;

    int n = hull.size();
    if (n < 2) {
        this->maxCellDistance = 0; 
    } else { 
        int j = 1;
        for (int i = 0; i < n; ++i) {
            while (euclidean_distance(hull[i], hull[(j + 1) % n]) > euclidean_distance(hull[i], hull[j]))
                j = (j + 1) % n;

            double d = euclidean_distance(hull[i], hull[j]);
            if (d > max_dist) {
                max_dist = d;
                farthest_pair = {hull[i], hull[j]};
            }
        }

        this->maxCellDistance = max_dist;

#if 1
        std::cout << "[" << this->level << "]: "; 
        std::cout << "The maximum distance is " << this->maxCellDistance << " between: ";
        std::cout << farthest_pair.first.x() << " " << farthest_pair.first.y() << " and ";
        std::cout << farthest_pair.second.x() << " " << farthest_pair.second.y() << std::endl;
#endif
    }
}

void Subgraph::addVertex(const SubgraphVertex &vertex) {
    this->graph[vertex.id] = vertex;
}

SubgraphVertex& Subgraph::getVertex(int id) {
    return this->graph.at(id);
}

void Subgraph::addEdge(int origin, int target, double weight) {
    SubgraphEdge *edge;

    // Ensure both vertices of this edge are in the graph
    // (they might not be if one is an I/O pad)
    if (this->graph.find(target) != this->graph.end() 
            && this->graph.find(origin) != this->graph.end()) {
        SubgraphVertex &originVert = this->graph.at(origin);
        SubgraphVertex &targetVert = this->graph.at(target);

        // Check if an edge already exists in the graph
        if (originVert.connections.find(target) != originVert.connections.end()) {
            // Edge exists origin -> target. Update it
            edge = originVert.connections.at(target);

            // Update weight by adding
            edge->weight += weight;
        } else if (targetVert.connections.find(origin) != targetVert.connections.end()) {
            // Edge exists target -> origin.
            // We don't care about direction, so this is the same edge. update it
            edge = targetVert.connections.at(origin);

            edge->weight += weight;
        } else {
            edge = new SubgraphEdge;
            edge->originNode = origin;
            edge->targetNode = target;
            edge->weight = weight;

            // Add edge to graph
            this->edges.insert(edge);
            originVert.connections[target] = edge;
            targetVert.connections[origin] = edge;
        }
    }
}

void Subgraph::dump(std::ostream &out) const {
    std::cout << "Maximum/Minimum Cell Distance: " << this->maxCellDistance << " ";
    std::cout << this->minCellDistance << std::endl;

    for (const auto &[id, vert] : this->graph) {
        std::cout << id << " connections and weights:" << std::endl;
        for (const auto &[target, con] : vert.connections) {
            std::cout << "  -> " << target << " (weight: " << con->weight << ")" << std::endl;
        }
    }
}

std::ostream& operator<<(std::ostream &out, const Subgraph &subgraph) {
    subgraph.dump(out);
    return out;
}

std::unordered_set<NetlistNode*> findNeighbors(Netlist &netlist, int baseNode, int maxSearchLevel) {
    std::unordered_set<NetlistNode*> neighbors;
    std::queue<NetlistNode*> Q;

    Q.push(&netlist.at(baseNode));

    while(!Q.empty()) {
        int levelDiff;
        NetlistNode *front = Q.front();
        Q.pop();

        levelDiff = abs(nodeLevelMap[front->id] - nodeLevelMap[baseNode]);
        if (levelDiff <= maxSearchLevel) {
            // Base node is not a neighbor of itself
            if (front->id != baseNode) {
                neighbors.insert(front);
            }

            // Foreach child node
            for (const int child : front->fanOutList) {
                if (abs(nodeLevelMap.at(child) - nodeLevelMap.at(baseNode)) <= maxSearchLevel
                        && neighbors.find(&netlist.at(child)) == neighbors.end()) {
                    Q.push(&netlist.at(child));
                }
            }
            
            // Foreach parent node
            for (const int parent : front->fanInList) {
                if (abs(nodeLevelMap.at(parent) - nodeLevelMap.at(baseNode)) <= maxSearchLevel
                        && neighbors.find(&netlist.at(parent)) == neighbors.end()) {
                    Q.push(&netlist.at(parent));
                }
            }
        }
    }

    return neighbors;
}

void computeLogicLevels(Netlist &netlist) {
    unordered_map<int, int> inDegree;
    for (const auto &[id, gate] : netlist) {
        inDegree[id] = gate.fanInList.size();
        nodeLevelMap[id] = 0;
    }

    queue<int> q;
    for (const auto &[id, deg] : inDegree)
        if (deg == 0) q.push(id);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        int uLevel = nodeLevelMap[u];
        for (int v : netlist[u].fanOutList) {
            nodeLevelMap[v] = max(nodeLevelMap[v], uLevel + 1);
            if (--inDegree[v] == 0) q.push(v);
        }
    }

    for (const auto &[id, level] : nodeLevelMap) {
        std::cout << "Node: " << id << ", level: " << level << std::endl;
        levelToNodes[level].push_back(id);
    }
}

// Super-cell ID assignment: cellID -> groupID
unordered_map<int, int> superCellMap;

void connectivityGraphProcessing(Netlist &netlist) {
    for (const auto &[level, nodes] : levelToNodes) {
        for (const int u : nodes) {
            for (const int v : nodes) {
                if (u != v) {
                    std::unordered_set<NetlistNode*> uNeighbors;
                    std::unordered_set<NetlistNode*> vNeighbors;
                    std::unordered_set<NetlistNode*> commonNeighbors;

                    NetlistNode *uNode = &netlist.at(u);
                    NetlistNode *vNode = &netlist.at(v);

                    // TODO: cache neighbors for each node ahead of time
                    uNeighbors = findNeighbors(netlist, uNode->id, MAX_SEARCH_LEVEL);
                    vNeighbors = findNeighbors(netlist, vNode->id, MAX_SEARCH_LEVEL);

                    commonNeighbors = uNeighbors;
                    commonNeighbors.merge(vNeighbors);

                    for (const auto &neighbor : commonNeighbors) {
                        if (u != neighbor->id && v != neighbor->id) {
                            int levelDiff = abs(nodeLevelMap[u] - nodeLevelMap[neighbor->id]);
                            double edgeWeight = NORMALIZATION_FACTOR / (levelDiff * 1.0);

                            subgraphs.at(level)->addEdge(u, v, edgeWeight);
#if 0
                            std::cout << "Edge Weight added between " << u << " (" << nodeLevelMap[u] << ")";
                            std::cout << " and " << neighbor->id << " (" << nodeLevelMap[neighbor->id] << ")";
                            std::cout << " is " << edgeWeight;
                            std::cout << " (v is " << v << ")" << std::endl;
#endif
                        }
                    }
                }
            }
        }
    }
}

void distanceGraphProcessing(Netlist &netlist) {
    for (const auto &[level, subgraph] : subgraphs) {
        // TODO: is Ximin and Ximax only X/Y dimension or euclidean distance?

        for (const auto &[id, vertex] : subgraph->getVertices()) {
            for (const auto &[id2, vertex2] : subgraph->getVertices()) {
                if (id != id2 && netlist.at(id).placement.isPlaced 
                        && netlist.at(id2).placement.isPlaced) {
                    // Calculate distance between cells
                    double Xu = netlist.at(id).placement.p.x();
                    double Xv = netlist.at(id2).placement.p.y();
                    double Yu = netlist.at(id).placement.p.x();
                    double Yv = netlist.at(id2).placement.p.y();

                    double Wx = DISTANCE_NORMALIZATION_FACTOR * (
                            1 - ((abs(Xu - Xv) - subgraph->getMinCellDistance())
                                / (subgraph->getMaxCellDistance() - subgraph->getMinCellDistance()))
                            );
                    double Wy = DISTANCE_NORMALIZATION_FACTOR * (
                            1 - ((abs(Yu - Yv) - subgraph->getMinCellDistance())
                                / (subgraph->getMaxCellDistance() - subgraph->getMinCellDistance()))
                            );

#if 0
                    std::cout << "Xu: " << Xu << ", " << Xv;
                    std::cout << " Min Max: " << subgraph->getMinCellDistance() << " ";
                    std::cout << subgraph->getMaxCellDistance() << std::endl;
#endif

                    double W = floor(sqrt(pow(Wx, 2) + pow(Wy, 2)));

                    subgraph->addEdge(id, id2, W);
                }
            }
        }
    }
}

void groupCells(Netlist &netlist) {
    int superCellID = 0;
    for (const auto &[level, nodes] : levelToNodes) {
        for (size_t i = 0; i < nodes.size(); i += GROUP_SIZE_K) {
            for (size_t j = i; j < min(i + GROUP_SIZE_K, nodes.size()); ++j) {
                superCellMap[nodes[j]] = superCellID;
            }
            superCellID++;
        }
    }
}

void writeGroupingToFile(const string &filename) {
    ofstream out(filename);
    for (const auto &[node, group] : superCellMap) {
        out << "n" << node << " -> SuperCell " << group << "\n";
    }
    out.close();
    cout << "Super-cell mapping written to " << filename << endl;
}

void doGrouping(Netlist &netlist) {
    // Dummy parser: insert gate parsing code here or link to your existing parser
    // Example: parsingCircuitFile("b15_1.isc", netlist);

    std::cout << "Computing Logic levels..." << std::endl;
    computeLogicLevels(netlist);


    // Create the subgraphs for each logic level
    for (const auto &[level, nodes] : levelToNodes) {
        std::cout << "Creating subgraph for level " << level << std::endl;
        subgraphs[level] = makeSubgraph(netlist, level);
    }

    std::cout << "Running connectivity-based graph processing" << std::endl;

    // Run connectivity-based graph processing step
    connectivityGraphProcessing(netlist);

    // Calculate min/max cell distances
    for (const auto &[level, subgraph] : subgraphs) {
        std::cout << "Calculating min/max cell distances for subgraph " << level << std::endl;
        subgraph->calcMinMaxCellDistances(netlist);
    }

#if 0
    std::cout << "Subgraph 3 (Connectivity Processed)" << std::endl << *(subgraphs.at(3)) << std::endl;
#endif

    // Run Distance-based Graph Processing step
    std::cout << "Running distance-based graph processing" << std::endl;
    distanceGraphProcessing(netlist);

#if 0
    std::cout << "Subgraph 3 (Distance Processed)" << std::endl << *(subgraphs.at(3)) << std::endl;
#endif

    SupercellsPlacer supercells(&netlist, &subgraphs);
    supercells.process();
    supercells.displaySupercells(std::cout);
}
