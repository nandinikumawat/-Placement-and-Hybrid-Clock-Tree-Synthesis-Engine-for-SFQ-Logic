// Based on fastplace main.cpp
#include <iostream>
#include <cstring>
#include <ostream>

#include "suraj_parser.h"
#include "placer.hpp"
#include "netlist.hpp"
#include "grouping.hpp"


int main(int argv, char *argc[])
{
    char inareFileName[100];
    char innetFileName[100];
    char iniscasFileName[100];
    char inPadLocationFileName[100];

    if (argv!=2) {
        std::cout << "Please provide a circuit file name with no extension." << std::endl;
        return 1;
    }

    strcpy (inareFileName, argc[1]);
    strcat(inareFileName, ".are");
    strcpy(innetFileName,argc[1]);
    strcat(innetFileName,".net");
    strcpy(iniscasFileName, argc[1]);
    strcat(iniscasFileName, ".isc");
    strcpy(inPadLocationFileName,argc[1]);
    strcat(inPadLocationFileName,".kiaPad");

    std::cout << "Reading ISCAS circuit file " << iniscasFileName << std::endl;

    Netlist *netlist = new Netlist();
    netlist->loadFromDisk(iniscasFileName);
    
    std::cout << "Converting to hypergraph format." << std::endl;
    netlist->saveHypergraphFile(argc[1], true);

    std::cout << "Feeding hypergraph to FastPlace" << std::endl;
    int success = parseIbmFile(inareFileName, innetFileName, inPadLocationFileName);
    if (success == -1) {
        cout << "Error reading input file(s)" << endl;
        return 0;
    }

    printf("\nNumber of vertices,hyper = %d %d\n",numCellsAndPads,numhyper);

    PA3Placement::AnalyticPlacer placer;

    std::cout << "Invoking FastPlace" << std::endl;
    placer.doPlacement(argc[1]);

    std::cout << "Loading initial Placement data" << std::endl;

    // Placement is saved in "xx_spread.kiaPad"
    netlist->loadPlacementKiaPad(std::string(argc[1]) + "_spread");

    std::cout << "NETLIST OBJECT DUMP:" << std::endl;
    std::cout << *netlist << std::endl;

    doGrouping(*netlist);

    delete netlist;

    free(pinLocations);
    free(hEdge_idxToFirstEntryInPinArray);
    free(cellPinArray);
    free(hyperwts);
}
