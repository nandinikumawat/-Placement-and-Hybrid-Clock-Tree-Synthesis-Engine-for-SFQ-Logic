#ifndef SFQPLACE_SUPERCELLS_HPP
#define SFQPLACE_SUPERCELLS_HPP

#include "netlist.hpp"
#include "grouping.hpp"

#include <ostream>
#include <unordered_map>

class SupercellsPlacer {
public:
    SupercellsPlacer(Netlist *originalNetlist, std::unordered_map<int, Subgraph*> *subgraphs);

    void process();

    void displaySupercells(std::ostream &out);
private:
    Netlist *originalNetlist;
    Netlist supercellNetlist;
    std::unordered_map<int, Subgraph*> *subgraphs;

    std::unordered_map<int, std::unordered_set<int>> supercells;
    // Master mapping of vertex IDs to supercell IDs
    std::unordered_map<int, int> verticesToSupercells;

    /**
     * Create a single supercell for ALL the vertices in a subgraph.
     * This is used when the subgraph is too small to partition, or partitioning failed.
     */
    void createSupercell(int id, Subgraph &subgraph);
    /**
     * Create supercells from a mapping of vertices to grouping IDs.
     * The individual grouping IDs are NOT supercell IDs, they are relative to the provided
     * map only.
     *
     * As this may create multiple supercells, it accepts a starting ID and will increment accordingly
     * until all supercells for the group mappings are created.
     *
     * It then returns the number of supercells created
     */
    int createSupercells(int startingId, std::unordered_map<int, int> verticesToSupercells);

    /**
     * Generates a netlist of supercells, adding edges between them corresponding to the 
     * edges between the individual cells inside.
     */
    void createSupercellNetlist();
};

#endif // SFQPLACE_SUPERCELLS_HPP
