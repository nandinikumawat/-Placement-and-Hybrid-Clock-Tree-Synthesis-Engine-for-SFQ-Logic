#include "supercells.hpp"
#include "partitioning.hpp"

#include <unordered_map>

static const int GROUP_SIZE_K = 4;

SupercellsPlacer::SupercellsPlacer(Netlist *ogNet, std::unordered_map<int, Subgraph*> *subgraphs) {
    this->originalNetlist = ogNet;
    this->subgraphs = subgraphs;
}

void SupercellsPlacer::process() {
    this->supercells.clear();

    std::cout << "Begin subgraph partitioning" << std::endl;

    // Perform partitioning
    int nextSupercellId = 0;
    for (const auto &[level, subgraph] : *(this->subgraphs)) {
        int p = std::ceil(subgraph->getVertices().size() / (1.0 * GROUP_SIZE_K));

        if (p <= 1) {
            // Too small, place all cells into a single supercell
            std::cout << "Logic level" << subgraph->getLogicLevel();
            std::cout << " is too small to partition. All cells in single supercell";
            std::cout << std::endl;

            this->createSupercell(nextSupercellId++, *subgraph);
        } else {
            PWayPartitioner partitioner(subgraph, p);

            std::cout << "Partitioning logic level " << subgraph->getLogicLevel();
            std::cout << " into " << p << " parts" << std::endl;
            if (partitioner.doPartition() < 0) {
                // Partitioning failed. Place all these cells into a single supercell
                std::cerr << "WARNING: couldn't partition logic level ";
                std::cerr << subgraph->getLogicLevel() << std::endl;

                this->createSupercell(nextSupercellId++, *subgraph);
            } else {
                // Partitioning succeeded. Create supercells using the group mappings
                // from the partitioner

                int numSCells = this->createSupercells(nextSupercellId, partitioner.getPartitions());
                nextSupercellId += numSCells;
            }
        }
    }

    std::cout << "Creating supercell netlist" << std::endl;
    this->createSupercellNetlist();

    // Write FastPlace-format input files for this supercell netlist
    this->supercellNetlist.saveHypergraphFile("supercells", true);

    // Invoke FastPlace 
    // run khmetis via system call feeding it options and the temporary input file
    // then read the output file partition data and return number of partitions
    std::string cmd = "./PA3 supercells";
    std::cout << "Running FastPlace with command: " << cmd << std::endl;

    int ret = system(cmd.c_str());
    std::cout << "FastPlace finished with code: " << ret << std::endl;
}

void SupercellsPlacer::createSupercell(int id, Subgraph &subgraph) {
    for (const auto &[logicLevel, vertex] : subgraph.getVertices()) {
        this->supercells[id].insert(vertex.id);
        this->verticesToSupercells[vertex.id] = id;
    }
}

int SupercellsPlacer::createSupercells(int id, std::unordered_map<int, int> groupings) {
    int largestGroupId = 0;

    for (const auto &[vertex, group] : groupings) {
        if (group > largestGroupId) {
            largestGroupId = group;
        }

        this->supercells[id + group].insert(vertex);
        this->verticesToSupercells[vertex] = id + group;
    }

    // + 1 supercells created, because group IDs start at 0
    return largestGroupId + 1;
}

void SupercellsPlacer::displaySupercells(std::ostream &out) {
    out << "Supercell Information:" << std::endl;
    for (const auto &[supercell, members] : this->supercells) {
        out << "Supercell ID: " << supercell << ", members:" << std::endl;
        out << "\t";
        for (const int vertex : members) {
            out << vertex << ' ';
        }
        out << std::endl;
    }
}

void SupercellsPlacer::createSupercellNetlist() {
    // First copy I/O pads from the original netlist, as they are not present
    // in the supercells

    for (auto [id, node] : *(this->originalNetlist)) {
        if (node.isPrimaryInput || node.isPrimaryOutput) {
            // Clear the fan in/out lists because IDs change (only supercell IDs now)
            node.fanOutList.clear();
            node.fanInList.clear();
            this->supercellNetlist[id] = node;
        }
    }

    for (const auto &[supercell, members] : this->supercells) {
        NetlistNode supercellNode;

        supercellNode.id = supercell;
        supercellNode.isPrimaryOutput = false;
        supercellNode.isPrimaryInput = false;

        // Add connections between this supercell and others
        // based on the connections of the vertices inside it
        for (const int member : members) {
            for (const int dependency : this->originalNetlist->at(member).fanOutList) {
                NetlistNode originalNode = this->originalNetlist->at(dependency);

                // Check if this dependency is an I/O pad or normal cell
                if (originalNode.isPrimaryInput || originalNode.isPrimaryOutput) {
                    supercellNode.fanOutList.insert(dependency);
                    
                    // This supercell fans out to this I/O pad, so add ourself to it's fanin
                    this->supercellNetlist[dependency].fanInList.insert(supercell);
                } else {
                    if (this->verticesToSupercells.find(dependency) != this->verticesToSupercells.end()) {
                        supercellNode.fanOutList.insert(this->verticesToSupercells.at(dependency));
                    }
                }

                //std::cout << "Supercell " << supercell << " fans out to " << dependency << std::endl;
            }

            // Go through the Fan-In list too, but only look at I/O pads
            // We need to make sure we copy the input pads and connect them to the supercells
            // TODO: this approach is extremely unoptimized
            for (const int dependency : this->originalNetlist->at(member).fanInList) {
                NetlistNode originalNode = this->originalNetlist->at(dependency);

                // Check if this dependency is an I/O pad or normal cell
                if (originalNode.isPrimaryInput || originalNode.isPrimaryOutput) {
                    // This I/O pad fans out to ourself (supercell)
                    this->supercellNetlist[dependency].fanOutList.insert(supercell);
                    //std::cout << "I/O pad " << dependency << " fans out to supercell " << supercell << std::endl;
                }
            }
        }

        // Add supercell to the netlist
        this->supercellNetlist[supercell] = supercellNode;
    }

}
