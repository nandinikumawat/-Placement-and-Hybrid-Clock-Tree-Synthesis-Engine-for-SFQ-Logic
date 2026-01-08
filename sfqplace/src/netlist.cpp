#include "netlist.hpp"
#include "suraj_parser.h"

#include <climits>
#include <functional>
#include <iostream>
#include <fstream>
#include <ostream>
#include <unordered_set>
#include <queue>
#include <sstream>
#include <string>

static const int DEFAULT_CHIP_WIDTH = 6;
static const int DEFAULT_CHIP_HEIGHT = 6;
static const int DEFAULT_CELL_AREA = 4;

bool Netlist::loadFromDisk(const std::string &filename) {
    // Modified ChatGPT generated code
    
    bool status = true;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "ISCAS85 Parser: Error opening file: " << filename << std::endl;
        status = false;
    } else {
        std::string line;
        std::unordered_map<std::string, int> nodeNamesToIds;

        this->nextId = 0;

        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '*') { // Skip comments
                std::istringstream iss(line);
                NetlistNode currentNode;
                int numFanOut, numFanIn;

                // Read in node ID
                iss >> currentNode.id;

                if (currentNode.id >= this->nextId) {
                    this->nextId = currentNode.id + 1;
                }

                // Check if this node is already in the netlist
                if (this->find(currentNode.id) != this->end()) {
                    // Obtain existing data (fanout list?)
                    currentNode = this->at(currentNode.id);
                }

                iss >> currentNode.name >> currentNode.nodeType;

                // Associate the name to the ID 
                nodeNamesToIds[currentNode.name] = currentNode.id;

                if (ISCAS85_NODE_TYPE_FANOUT_BRANCH == currentNode.nodeType) {
                    // fanout branch type is special
                    std::string fanInNodeName;
                    iss >> fanInNodeName;

                    std::string fault;
                    while (iss >> fault) {
                        // Ignore faults
                    }

                    // Link the fanout node to this fanout branch
                    this->at(nodeNamesToIds.at(fanInNodeName)).fanOutList.insert(currentNode.id);
                    // Add the originating fanout node to our fanin list
                    currentNode.fanInList.insert(nodeNamesToIds.at(fanInNodeName));
                } else {
                    iss >> numFanOut >> numFanIn;

                    if (ISCAS85_NODE_TYPE_INPUT == currentNode.nodeType) {
                        currentNode.isPrimaryInput = true;
                    }
                    
                    std::string fault;
                    while (iss >> fault) {
                        // Ignore faults
                    }
                    
                    if (numFanIn > 0) {
                        std::getline(file, line);

                        std::istringstream fanin_stream(line);
                        int fanin_node;

                        while (fanin_stream >> fanin_node) {
                            currentNode.fanInList.insert(fanin_node);

                            // Look up the fanin node and add us to it's fanout list
                            if (this->find(fanin_node) != this->end()) {
                                // Fanin Node already in the netlist, add us to it
                                this->at(fanin_node).fanOutList.insert(currentNode.id);
                            } else {
                                // Fanin Node not already in the netlist
                                NetlistNode fanInNode;
                                fanInNode.id = fanin_node;
                                fanInNode.fanOutList.insert(currentNode.id);
                            }
                        }
                    }
                }

                (*this)[currentNode.id] = currentNode;
            }
        }

        // Add output "pads" for any cell that has a zero fanout
        this->addOutputs();

        file.close();
    }

    return status;
}

bool Netlist::loadPlacementKiaPad(const std::string &filePrefix) {
    bool status = true;
    std::ifstream file(filePrefix + ".kiaPad");
    std::string line;

    if (!file.is_open()) {
        std::cerr << "ISCAS85 Parser: Error opening file: " << filePrefix << std::endl;
        status = false;
    } else {
        while(std::getline(file, line)) {
            std::istringstream lineStream(line);
            std::string idStr;
            int id;
            bool isPad;
            double x;
            double y;
            
            // Get ID
            lineStream >> idStr;
            isPad = (idStr.at(0) == 'p');

            if (isPad) {
                // Pad found
                id = std::stoi(idStr.substr(1)); // Skip pad prefix
            } else {
                id = std::stoi(idStr);
            }

            // Get coordinates
            lineStream >> x >> y;

            if (!isPad) {
                // Ignore I/O pads, they don't move
                
                if (this->hyperIdMappings.find(id) == this->hyperIdMappings.end()) {
                    std::cerr << "Failed to find node with ID: " << id;
                    std::cerr << ", while loading placement data" << std::endl;
                } else {
                    int mappedId = this->hyperIdMappings.at(id);
                    this->at(mappedId).placement.isPlaced = true;
                    this->at(mappedId).placement.p = Point(x, y);
                }
            }
        }

        file.close();
    }

    return status;
}

bool Netlist::saveHypergraphFile(const std::string &outputFilename, bool genPadFile) {
    // Modified ChatGPT generated code
    
    bool status = true;
    std::ofstream outFile(outputFilename + ".net");
    std::ofstream areaFile(outputFilename + ".are");

    // Assign "hyperId" to each node
    // Removes gaps in the id space, parser may not work without it
    this->consolidateIds();

    if (!outFile.is_open() || !areaFile.is_open()) {
        std::cerr << "Error opening output file: " << outputFilename << std::endl;
        status = false;
    } else {
        int totalEndpoints = 0;
        int totalCells = this->size();
        std::unordered_set<int> uniqueHyperedges;
        std::unordered_set<int> gateCells;
        std::set<int, std::greater<int>> inPadCells;
        std::set<int, std::greater<int>> outPadCells;
        
        for (const auto &pair : *this) {
            std::cout << "Node " << pair.first << "type: " << pair.second.nodeType << std::endl;

            if (!pair.second.isPrimaryInput && !pair.second.isPrimaryOutput) {
                gateCells.insert(pair.first);
            } else if (pair.second.isPrimaryInput) {
                inPadCells.insert(pair.first);
            } else {
                outPadCells.insert(pair.first);
            }

            if (!pair.second.fanOutList.empty()) {
                totalEndpoints += pair.second.fanOutList.size() + 1;
                uniqueHyperedges.insert(pair.first);
            }
        }
        
        outFile << 0 << std::endl;
        outFile << totalEndpoints << "\n";
        outFile << uniqueHyperedges.size() << "\n";
        outFile << totalCells << "\n";
        outFile << gateCells.size() - 1 << "\n"; // suraj_Parser adds 1 for some reason

        for (int gateId : gateCells) {
            this->hypergraphWriteNode(this->at(gateId), areaFile, outFile);
        }

        for (int padId : inPadCells) {
            this->hypergraphWriteNode(this->at(padId), areaFile, outFile);
        }

        for (int padId : outPadCells) {
            this->hypergraphWriteNode(this->at(padId), areaFile, outFile);
        }
       
        outFile.close();
        areaFile.close();

        if (genPadFile) {
            this->generatePadFile(outputFilename, inPadCells, outPadCells);
        }
    }

    return status;
}

int Netlist::levelsBetween(int startId, int endId) {
    int levels = -1;

    // Calculate the number of logic levels between two nodes:
    // Djikstra's algorithm

    std::unordered_map<int, int> distances;
    // min-priority queue
    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        std::greater<>> djikstraQueue;
    std::unordered_map<int, int> prev;

    // Set all nodes to unvisited, no previous pointers
    for (const auto &pair : *this) {
        distances[pair.first] = INT_MAX;
        prev[pair.first] = INT_MIN;
    }

    // Insert source itself in queue and set distance to zero
    djikstraQueue.push({0, startId});
    distances.at(startId) = 0;

    // Continue while queue is non-empty and the next node isn't the end node
    while (!djikstraQueue.empty() && djikstraQueue.top().second != endId){
        int distance = djikstraQueue.top().first;
        int currentNodeId = djikstraQueue.top().second;

        djikstraQueue.pop();

        // Only continue if we have a longer best-distance stored
        if (distance <= distances.at(currentNodeId)) {
            for (const int v : this->at(currentNodeId).fanOutList) {
                if (distances.at(v) > distances.at(currentNodeId) + 1) {
                    distances.at(v) = distances.at(currentNodeId) + 1;
                    djikstraQueue.push({distances.at(v), v});
                }
            }
        }
    }

    if (distances.at(endId) != INT_MAX) {
        levels = distances.at(endId);
    }

    return levels;
}

void Netlist::hypergraphWriteNode(const NetlistNode &node, std::ostream &areaFile, std::ostream &graphFile) {
    bool cellIsPad = node.isPrimaryOutput || node.isPrimaryInput;
    std::string prefix = (cellIsPad) ? "p" : "a";

    // Write cell to area file
    // Just give every cell the same area
    areaFile << prefix << node.hypergraphId;
    areaFile << " " << (cellIsPad ? 0 : DEFAULT_CELL_AREA) << std::endl;
    
    if (node.fanOutList.size() > 0) {
        graphFile << prefix << node.hypergraphId << " s 1\n";

        for (int fanoutNode : node.fanOutList) {
            std::string fanoutPrefix = (this->at(fanoutNode).isPrimaryInput) ? "p" : "a";
            graphFile << fanoutPrefix << this->at(fanoutNode).hypergraphId << " l\n";
        }
    }
}

bool Netlist::generatePadFile(const std::string &outFilePrefix, 
                              const std::set<int, std::greater<int>> &inPads,
                              const std::set<int, std::greater<int>> &outPads) {
    bool success = true;
    std::ofstream out(outFilePrefix + ".kiaPad");

    if (!out.is_open()) {
        std::cerr << "Failed to open output file for generating .kiaPad" << std::endl;
        success = false;
    } else {
        int x = 0;
        int y = 1;

        x = 0; 
        y = DEFAULT_CHIP_HEIGHT;

        for (const int pad : outPads) {
            out << 'p' << this->at(pad).hypergraphId << ' ' << x << ' ' << y << std::endl;
            x += 1;
        }

        x = 0;
        y = 1;

        for (const int pad : inPads) {
            out << 'p' << this->at(pad).hypergraphId << ' ' << x << ' ' << y << std::endl;
            x += 1;
        }


        out.close();
    }

    return success;
}

void Netlist::eliminateFanoutBranches(void) {
    unordered_set<int> fanoutBranches;

    for (const auto &pair : *this) {
        const NetlistNode &node = pair.second;

        if (ISCAS85_NODE_TYPE_FANOUT_BRANCH == node.nodeType) {
            // Eliminate this node 
            fanoutBranches.insert(pair.first);

            // Fanout branches have a fanin and fanout of 1, so set is only 1 element
            int fanInNode = *(node.fanInList.begin());
            int fanOutNode = *(node.fanOutList.begin());

            // Manipulate fanin, fanout lists of these two nodes to link them together

            // Disconnect us from fanInNode's fanOut list
            this->at(fanInNode).fanOutList.erase(node.id);
            this->at(fanInNode).fanOutList.insert(fanOutNode);

            // Disconnect us from fanOutNode's fanIn list
            this->at(fanOutNode).fanInList.erase(node.id);
            this->at(fanOutNode).fanInList.insert(fanInNode);
        }
    }

    for (const auto &it : fanoutBranches) {
        this->erase(it);
    }
}

void Netlist::addOutputs(void) {
    stringstream outputNodeName;

    // First eliminate fanout branches
    this->eliminateFanoutBranches();

    // Iterate and find nodes with zero fanOut
    for (auto &pair : *this) {
        if (pair.second.fanOutList.empty() && pair.second.nodeType != NODE_TYPE_OUTPUT) {
            NetlistNode outputNode;

            outputNode.id = this->nextId++;

            outputNodeName.clear();
            outputNodeName << outputNode.id << "OUT";

            // Add a new output node (represents an I/O pad)
            // and connect the cell to it
            outputNode.isPrimaryOutput = true;
            outputNode.fanInList.insert(pair.first);
            outputNode.name = outputNodeName.str();
            outputNode.nodeType = NODE_TYPE_OUTPUT;

            pair.second.fanOutList.insert(outputNode.id);

            (*this)[outputNode.id] = outputNode;
        }
    }
}

void Netlist::consolidateIds(void) {
    int moveableCellCounter = 0;
    int padCounter = 1;

    for (auto &pair : *this) {
        if (pair.second.isPrimaryInput || pair.second.isPrimaryOutput) {
            pair.second.hypergraphId = padCounter++;
        } else {
            pair.second.hypergraphId = moveableCellCounter++;

            // Only map moveable cells
            this->hyperIdMappings[pair.second.hypergraphId] = pair.second.id;
        }
    }
}


std::ostream& operator<<(std::ostream &out, const Netlist &netlist) {
    // Modified ChatGPT generated code

    for (const auto &pair : netlist) {
        const NetlistNode &node = pair.second;

        out << "Node " << node.id << "[" << node.nodeType << "] ";
        out << "(" << node.name << ", hypergraph ID: " << node.hypergraphId << "): "
                  << " Fanin: " << node.fanInList.size() << " Fanout: " 
                  << node.fanOutList.size() << "\n";
        out << "  Placed(yes/no: " << node.placement.isPlaced << "): ";
        out << node.placement.p.x() << " " << node.placement.p.y() << "\n";

        if (!node.fanInList.empty()) {
            out << "  Fanin nodes: ";
            for (int fn : node.fanInList) {
                out << fn << " ";
            }
            out << "\n";
        }
    }

    return out;
}

