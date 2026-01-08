#ifndef SFQPLACE_PARTITIONING_HPP
#define SFQPLACE_PARTITIONING_HPP

#include "grouping.hpp"

// HMETIS C symbols
extern "C" {
    void HMETIS_PartRecursive(int nvtxs, int nhedges, int *vwgts, int *eptr, int *eind, int *hewgts, int nparts, int ubfactor, int *options, int *part, int *edgecut);

    void HMETIS_PartKway(int nvtxs, int nhedges, int *vwgts, int *eptr, int *eind, int *hewgts, int nparts, int ubfactor, int *options, int *part, int *edgecut);
}

class PWayPartitioner {
public:
    /**
     * Creates a new P-way partitioner for the given subgraph,
     * preparing to partition it into P groups
     */
    PWayPartitioner(Subgraph *subgraph, int groups);

    ~PWayPartitioner();

    /**
     * Converts the subgraph internally to an input format readable by HMETIS,
     * then invokes HMETIS to perform the partitioning.
     * Will save the results internally. To obtain which vertexes belong to which partition,
     * call getPartitions() after this.
     *
     * Returns the total number of partitions created.
     */
    int doPartition(void);

    /**
     * Returns a mapping from each vertex ID to the partition it belongs to.
     * doPartition() must have been called before this.
     */
    std::unordered_map<int, int> getPartitions();

private:
    Subgraph *subgraph;
    int desiredPartitionCount;

    std::string hmetisInputFilename;

    int nvtxs;
    int nhedges;

    // Hyperedge weights array. hewgts[i] = weight of hyperedge i
    int *hewgts;
    // Map of hyperedge ID to the beginning of it's vertex list in "eind"
    int *eptr;
    // List of verticies in each hyperedge
    int *eind;

    // Array that maps the index to which partition it belongs to
    int *partitionedData;

    // Mapping of each vertex to its corresponding partition
    std::unordered_map<int, int> verticesToSupercells;

    // Mapping from the subgraph vertex IDs to IDs for HMETIS
    // IDs in the subgraph are same as the netlist, so they might not start at zero
    // HMETIS needs consecutive IDs starting at 0
    std::unordered_map<int, int> hmetisIdsMap;
    // Mapping from HMETIS IDs to subgraph vertex IDs
    std::unordered_map<int, int> sgraphIdsMap;

    void writeHMETISInput(const std::unordered_map<int, std::unordered_set<int>> &hedges,
                          const std::vector<int> &weights);
    int invokeHMETIS(void);

    // Frees all structures allocated by/for HMETIS.
    void freeHMETISStructures(void);
};

#endif // SFQPLACE_PARTITIONING_HPP
