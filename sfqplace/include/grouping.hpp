#ifndef SFQPLACE_GROUPING_HPP
#define SFQPLACE_GROUPING_HPP

#include "netlist.hpp"
#include <ostream>

struct SubgraphEdge {
    int originNode;
    int targetNode;
    double weight;
};

struct SubgraphVertex {
    int id;
    std::unordered_map<int, SubgraphEdge*> connections;
};

class Subgraph {
    public:
        Subgraph(int logicLevel) {
            this->level = logicLevel;
        }

        ~Subgraph();

        /**
         * Calculates the minimum and maximum cell distances in both
         * X and Y dimensions from the original netlist placement information
         * for all cells in this subgraph (and by extension this logic level)
         */
        void calcMinMaxCellDistances(Netlist &netlist);

        void addVertex(const SubgraphVertex &vert);
        void addEdge(int origin, int target, double weight);

        SubgraphVertex& getVertex(int id);

        std::unordered_map<int, SubgraphVertex>& getVertices(void) {
            return this->graph;
        };

        std::unordered_set<SubgraphEdge*> getEdges() { return this->edges; };

        void dump(std::ostream &out) const;

        int getLogicLevel(void) const { return level; };
        double getMinCellDistance(void) const { return minCellDistance; };
        double getMaxCellDistance(void) const { return maxCellDistance; };
    private:
        int level;

        double minCellDistance;
        double maxCellDistance;

        std::unordered_map<int, SubgraphVertex> graph;
        std::unordered_set<SubgraphEdge*> edges;

        void findMaxDistance(Netlist &netlist);
};

std::ostream& operator<<(std::ostream &out, const Subgraph &subgraph);

void doGrouping(Netlist &netlist);

#endif //SFQPLACE_GROUPING_HPP
