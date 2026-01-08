#ifndef SFQPLACE_NETLIST_HPP
#define SFQPLACE_NETLIST_HPP

#include <functional>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

static const std::string ISCAS85_NODE_TYPE_INPUT = "inpt";
static const std::string ISCAS85_NODE_TYPE_FANOUT_BRANCH = "from";
static const std::string NODE_TYPE_OUTPUT = "otpt";


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;

struct CellPlacementData {
    // True if placement data is available, false if not
    bool isPlaced;

    Point p;
};

struct NetlistNode {
    // identifier from the ISCAS file
    int id;
    // identifier for hypergraph format, used by suraj_parser and Fastplace
    int hypergraphId;

    std::string name;
    // Type of the node, see ISCAS_85_NODE_TYPE_*
    std::string nodeType;
    std::unordered_set<int> fanInList;
    std::unordered_set<int> fanOutList;

    CellPlacementData placement;

    bool isPrimaryOutput = false;
    bool isPrimaryInput = false;
};

class Netlist : public std::unordered_map<int, NetlistNode> {
public:
    /**
     * Reads a given netlist in the ISCAS '85 format.
     * If the file cannot be read this function will return false.
     * Otherwise, a successful load will return true.
     */
    bool loadFromDisk(const std::string &filename);

    /**
     * Load placement information for this netlist from FastPlace .kiaPad output.
     *
     * Do not include the file extension, that will be added automatically.
     */
    bool loadPlacementKiaPad(const std::string &filePrefix);

    /**
     * Saves the netlist in the hypergraph format used by the FastPlace
     * implementation (.net files). If successful this function returns true,
     * otherwise a return value of false indicates an error occurred.
     *
     * @param genPadFile Generate a .kiaPad file with auto-spaced pad locations
     *                   for all input and output cells in the netlist.
     *
     * Do not include the file extension, that will be added automatically.
     */
    bool saveHypergraphFile(const std::string &outFilePrefix, bool genPadFile);

    /**
     * Calculates the distance in terms of logic levels between two nodes.
     *
     * Returns the distance, or -1 if one or both of the node IDs are invalid.
     */
    int levelsBetween(int startId, int endId);
private:
    int nextId;

    // Mapping of hypergraph Ids (from fastplace) to the main NetlistNode ids
    // Does not include I/O pads!
    std::unordered_map<int, int> hyperIdMappings;

    void hypergraphWriteNode(const NetlistNode &node, std::ostream &areaFile, std::ostream &graphFile);

    /**
     * Generates a .kiaPad file with auto-spaced pad locations for the given set of
     * cell IDs.
     */
    bool generatePadFile(const std::string &outFilePrefix, 
                         const std::set<int, std::greater<int>> &inPads,
                         const std::set<int, std::greater<int>> &outPads);

    /**
     * Iterates through the netlist removing "fanout branches",
     * a special type of node defined by the ISCAS85 format, but useless.
     */
    void eliminateFanoutBranches(void);

    /**
     * Assigns a "hypergraphId" to each node that does not skip numbers,
     * used to create the hypergraph netlist format.
     */
    void consolidateIds(void);

    /**
     * Adds output pad nodes to any cell with zero fan out.
     */
    void addOutputs(void);
};

std::ostream& operator<<(std::ostream &out, const Netlist &netlist);

#endif // SFQPLACE_NETLIST_HPP
