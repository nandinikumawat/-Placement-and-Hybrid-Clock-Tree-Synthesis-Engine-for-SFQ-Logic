// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)


#include<iostream>
#include<stdio.h>
#include<vector>
#include <string.h>
#include "suraj_parser.h"

#include "placer.hpp"

using namespace std;


int main(int argv, char *argc[])
{
    char inareFileName[100];
    char innetFileName[100];
    char inPadLocationFileName[100];

    if (argv!=2) {
        cout << "Please provide a circuit file name with no extension." << endl;
        return 1;
    }

    cout << "Reading circuit file " << argc[1] << endl;

    strcpy (inareFileName, argc[1]);
    strcat(inareFileName, ".are");
    strcpy(innetFileName,argc[1]);
    strcat(innetFileName,".net");
    strcpy(inPadLocationFileName,argc[1]);
    strcat(inPadLocationFileName,".kiaPad");

    int success = parseIbmFile(inareFileName, innetFileName, inPadLocationFileName);
    if (success == -1) {
        cout << "Error reading input file(s)" << endl;
        return 0;
    }

    printf("\nNumber of vertices,hyper = %d %d\n",numCellsAndPads,numhyper);


    // call function(s) dealing with creating the Q matrix, placement, etc.

    PA3Placement::AnalyticPlacer placer;
    placer.doPlacement(argc[1]);

    free(pinLocations);
    free(hEdge_idxToFirstEntryInPinArray);
    free(cellPinArray);
    free(hyperwts);
}
