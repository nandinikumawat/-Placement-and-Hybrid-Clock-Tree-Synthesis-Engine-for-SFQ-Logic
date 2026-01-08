//
// Created by Alejandro Zeise on 12/11/23.
//

#include "placer.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <cassert>
#include <cmath>

namespace PA3Placement
{
    static const unsigned int NUM_BINS = 5;

    static const double SPREADING_ALPHA = 0.8;
    static const double SPREADING_SIGMA = 1.5;

    static const double CONJ_GRADIENT_TOLERANCE = 1e-6;
    static const int CONJ_GRADIENT_ITERATIONS = 1000;

    double AnalyticPlacer::calculateTotalWirelength(std::vector<std::pair<double, double>> &cellLocations) const
    {
        double sum = 0.0;

        if (nullptr != this->matrixQ)
        {
            std::unordered_set<std::pair<int, int>, pairHashInteger> visitedPairs;

            for (unsigned int currentCell = 0; currentCell < this->matrixQ->getCellConnectionsList()->size(); currentCell++)
            {
                // Only looking at non I/O pads as our "current cell" (exclude star nodes)
                if (currentCell < numCells_noPads) {
                    auto currentCellCoords = cellLocations.at(currentCell);
                    auto connectedCells = this->matrixQ->getCellConnectionsList()->at(currentCell);

                    // First element in pair is the cell number, second element is the connection weight
                    for (auto &connectedCell: connectedCells) {
                        std::pair<double, double> connectedCellCoords;
                        double connectionWeight = connectedCell.second;

                        // Only visit a pair once
                        if (visitedPairs.find({currentCell, connectedCell.first}) == visitedPairs.end()
                            && visitedPairs.find({connectedCell.first, currentCell}) == visitedPairs.end())
                        {

                            // Check if the cell is an I/O pad or not (include Star nodes)
                            if (connectedCell.first < numCells_noPads + this->matrixQ->getStarNodeCount())
                            {
                                // Not an I/O pad, and also not a star node
                                connectedCellCoords = cellLocations.at(connectedCell.first);
                            }
                            else
                            {
#if 0
                                // Is an I/O pad, get its coordinates from the pinLocations list
                                connectedCellCoords.first = this->matrixDx->getIOPadLocation(connectedCell.first)->x;
                                connectedCellCoords.second = this->matrixDx->getIOPadLocation(connectedCell.first)->y;
#else
                                continue; // Only calculating wirelength between movable cells
#endif
                            }

                            // w_i * (x1-x2)^2 + (y1-y2)^2
                            double diffx = currentCellCoords.first - connectedCellCoords.first;
                            double diffy = currentCellCoords.second - connectedCellCoords.second;

                            double x1x2 = pow(diffx, 2);
                            double y1y2 = pow(diffy, 2);

                            sum += connectionWeight * (x1x2 + y1y2);

                            visitedPairs.insert({currentCell, connectedCell.first});
                            //std::cout << "Visited pair: " << currentCell << " " << connectedCell.first << std::endl;
                        }
                    }
                }
            }
        }

        return sum;
    }

    void AnalyticPlacer::saveCellLocationsToDisk(std::string filename)
    {
        std::ofstream fout(filename);

        // Only save moveable cells (no star nodes, no I/O pads)
        for (unsigned int i = 0; i < this->cellLocations.size() - this->matrixQ->getStarNodeCount(); i++)
        {
            fout << i << " " << this->cellLocations.at(i).first << " " << this->cellLocations.at(i).second << std::endl;
        }

        // Now add the I/O pads
        for (int i = 0; i < numCellsAndPads - numCells_noPads; i++)
        {
            fout << "p" << i << " " << pinLocations[i].x << " " << pinLocations[i].y << std::endl;
        }
    }

    void AnalyticPlacer::saveSpreadedCellsToDisk(std::string filename)
    {
        std::ofstream fout(filename);

        // Spreaded cells list only has cells, no I/O pads. Exclude star nodes
        for (unsigned int i = 0; i < this->spreadedCellLocations.size() - this->matrixQ->getStarNodeCount(); i++)
        {
            fout << i << " " << this->spreadedCellLocations.at(i).first << " " << this->spreadedCellLocations.at(i).second << std::endl;
        }

        // Now add the I/O pads
        for (int i = 0; i < numCellsAndPads - numCells_noPads; i++)
        {
            fout << "p" << i << " " << pinLocations[i].x << " " << pinLocations[i].y << std::endl;
        }
    }

    void AnalyticPlacer::calculateCellLocations()
    {
        // Solve the Matrices to get X/Y coordinates for each cell
#if 0
        std::cout << "Solving for X coordinates..." << std::endl;
        auto resultX = solveMatrixGradientDescent(0.01, 100, *this->matrixQ, this->matrixDx->convertTo2D());
        std::cout << "Solving for Y coordinates..." << std::endl;
        auto resultY = solveMatrixGradientDescent(0.01, 100, *this->matrixQ, this->matrixDy->convertTo2D());
#else
        pthread_t xSolverTid, ySolverTid;
        ColumnMatrix<double> *resultX = new ColumnMatrix<double>(1, 0);
        ColumnMatrix<double> *resultY = new ColumnMatrix<double>(1, 0);

        MatrixSolverParams *xParams = new MatrixSolverParams();
        MatrixSolverParams *yParams;
        xParams->id = 1;
        xParams->tolerance = CONJ_GRADIENT_TOLERANCE;
        xParams->maxIterations = CONJ_GRADIENT_ITERATIONS;
        xParams->Q = this->matrixQ;
        xParams->Dx = this->matrixDx;
        xParams->xAnswer = resultX;

        yParams = new MatrixSolverParams(*xParams);
        yParams->id = 2;
        yParams->Dx = this->matrixDy;
        yParams->xAnswer = resultY;

        std::cout << "Solving for X coordinates..." << std::endl;
        pthread_create(&xSolverTid, NULL, (void *(*)(void *)) (matrixSolverThread), (void*) xParams);
        std::cout << "Solving for Y coordinates..." << std::endl;
        pthread_create(&ySolverTid, NULL, (void *(*)(void *)) (matrixSolverThread), (void*) yParams);

        // Wait for X and Y solver threads to finish
        pthread_join(xSolverTid, NULL);
        pthread_join(ySolverTid, NULL);
#endif

        // Place the x/y coordinates into our vector
        for (int i = 0; i < matrixQ->getHeight(); i++)
        {
            this->cellLocations.push_back({resultX->get(0, i), resultY->get(0, i)});
        }

        saveCellLocationsToDisk("preSpread.kiaPad");
    }

    std::pair<double, double> AnalyticPlacer::calculateChipDimensions()
    {
        std::pair<double, double> dimensions = {0, 0};

        // Look at movable cell coordinates
        for (auto &pos : this->cellLocations)
        {
            if (pos.first > dimensions.first)
            {
                dimensions.first = pos.first;
            }

            if (pos.second > dimensions.second)
            {
                dimensions.second = pos.second;
            }
        }

        // Check I/O pad coordinates now
        for (int i = 0; i < numCellsAndPads - numCells_noPads; i++)
        {
            if (pinLocations[i].x > dimensions.first)
            {
                dimensions.first = pinLocations[i].x;
            }

            if (pinLocations[i].y > dimensions.second)
            {
                dimensions.second = pinLocations[i].y;
            }
        }

        return dimensions;
    }

    void AnalyticPlacer::doPlacement(std::string filePrefix)
    {
        std::cout << "Constructing Matrices..." << std::endl;
        // Create Q, Dx, Dy matrices
        this->matrixQ = new QMatrix(numCells_noPads, numCellsAndPads, numhyper, cellPinArray, hEdge_idxToFirstEntryInPinArray, hyperwts);
        this->matrixDx = new DMatrix(DMatrix::Dimension::X, pinLocations, numCells_noPads, this->matrixQ->getStarNodeCount(), this->matrixQ->getCellConnectionsList());
        this->matrixDy = new DMatrix(DMatrix::Dimension::Y, pinLocations, numCells_noPads, this->matrixQ->getStarNodeCount(), this->matrixQ->getCellConnectionsList());

        //std::cout << *this->matrixQ << std::endl << std::endl;
#if 0
        double diagonal;
        double nonDiag = 0;
        std::multiset<double> elementsInRow;
        for (int i =0; i < this->matrixQ->getHeight(); i++)
        {
            nonDiag = 0;
            elementsInRow.clear();
            diagonal = this->matrixQ->get(i, i);
            for (int j = 0; j < this->matrixQ->getWidth(); j++)
            {
                if (this->matrixQ->get(j, i) != 0)
                {
                    //std::cout << "Connection between " << i << " and " << j << " " << this->matrixQ->get(j, i) << std::endl;
                    elementsInRow.insert(this->matrixQ->get(j, i));
                }

                if (j != i) {
                    nonDiag += this->matrixQ->get(j, i);
                    //elementsInRow.insert(this->matrixQ->get(j, i));
                }
            }
            double dx = this->matrixDx->get(0, i);
            double dy = this->matrixDy->get(0, i);
            if (dx == 0 && dy == 0 && -nonDiag == diagonal)
            {

            }
            else if (dx == 0 && dy == 0)
            {
                auto connections = this->matrixQ->getCellConnectionsList()->at(i);
                std::cerr << "Problem at " << i << std::endl;
                std::cerr << diagonal << " " << nonDiag << "( Dx/Dy: " << this->matrixDx->get(0, i) << " " << this->matrixDy->get(0, i) << ")" << std::endl;

            }
            //else
              //  std::cout << diagonal << " " << nonDiag << "( Dx/Dy: " << this->matrixDx->get(0, i) << " " << this->matrixDy->get(0, i) << ")" << std::endl;
        }

        for (int i = 0; i < this->matrixDx->getHeight(); i++)
        {
            //std::cout << this->matrixDy->get(0, i) << std::endl;
        }

#else
        std::cout << "Solving Matrices..." << std::endl;
#if 1
        this->calculateCellLocations();

        double wirelength = this->calculateTotalWirelength(this->cellLocations);

        std::cout << "Total Wirelength: " << wirelength << std::endl;
        std::cout << "Sqrt of total Wirelength: " << sqrt(wirelength) << std::endl;

        std::cout << "Spreading..." << std::endl;
        this->doSpreading(filePrefix);

        std::cout << "Sqrt of total Wirelength (post-spreading): " << sqrt(this->calculateTotalWirelength(this->spreadedCellLocations)) << std::endl;
#else
        std::cout << *this->matrixQ * *this->matrixDx << std::endl;

        AcceleratedConjugateGradientSolver accelSolver(*this->matrixQ, *this->matrixDx);
        accelSolver.solve();
#endif
#endif
    }

    void AnalyticPlacer::createBins()
    {
        const auto chipDimensions = this->calculateChipDimensions();
        const double xBinSize = chipDimensions.first / NUM_BINS;
        const double yBinSize = chipDimensions.second / NUM_BINS;

        bins.clear();

        // Divide the chip into bins
        double x = xBinSize;
        double y = yBinSize;
        for (unsigned int i = 0; i < NUM_BINS; i++)
        {
            for (unsigned int j = 0; j < NUM_BINS; j++)
            {
                Bin bin = {x - xBinSize, x, y - yBinSize, y, 0};

                x += xBinSize;

                bins.push_back(bin);
            }

            x = xBinSize;
            y += yBinSize;
        }
    }

    void AnalyticPlacer::calculateSpreadedCellLocations()
    {
        const auto chipDimensions = this->calculateChipDimensions();

        std::vector<Bin> unequalBins;

        for (unsigned int i = 0; i < bins.size(); i ++)
        {
            auto bin = bins.at(i);

            double OBx = bin.upperX;
            double OBy = bin.upperY;
            double OBxMinus1 = bin.lowerX;
            double OByMinus1 = bin.lowerY;
            double OBxPlus1;
            double OByPlus1;
            double UxPlus1;
            double UyPlus1;

            double Ux = bin.utilization;
            double Uy = bin.utilization;

            // Assume we are the first bin in the row or column
            double NBxMinus1 = 0;
            double NByMinus1 = 0;
            double NBx;
            double NBy;

            if (bin.lowerX != 0)
            {
                // We are not the first bin in the row
                NBxMinus1 = unequalBins.at(i - 1).upperX;
            }

            if (bin.upperX == chipDimensions.first)
            {
                // We are the last bin in the row
                NBx = bin.upperX;
            }
            else
            {
                // In middle of row
                OBxPlus1 = bins.at(i + 1).upperX;
                UxPlus1 = bins.at(i + 1).utilization;

                double top = (OBxMinus1 * (UxPlus1 + SPREADING_SIGMA) + OBxPlus1 * (Ux + SPREADING_SIGMA));
                double bottom = (Ux + UxPlus1 + (2 * SPREADING_SIGMA));

                NBx = top / bottom;
            }

            if (bin.lowerY != 0)
            {
                // We are not the first bin in the column
                NByMinus1 = unequalBins.at(i - NUM_BINS).upperY;
            }

            if (bin.upperY == chipDimensions.second)
            {
                // We are the last bin in the column
                NBy = bin.upperY;
            }
            else
            {
                // In middle of column
                OByPlus1 = bins.at(i + NUM_BINS).upperY;
                UyPlus1 = bins.at(i + NUM_BINS).utilization;

                double top = (OByMinus1 * (UyPlus1 + SPREADING_SIGMA) + OByPlus1 * (Uy + SPREADING_SIGMA));
                double bottom = (Uy + UyPlus1 + (2 * SPREADING_SIGMA));

                NBy = top / bottom;
            }

            for (auto &cell : bin.memberCells)
            {
                auto cellPos = this->cellLocations.at(cell);

                double originalX = cellPos.first;
                double originalY = cellPos.second;

                double topX = (NBx * (originalX - OBxMinus1) + NBxMinus1 * (OBx - originalX));
                double topY = (NBy * (originalY - OByMinus1) + NByMinus1 * (OBy - originalY));
                double bottomX = (OBx - OBxMinus1);
                double bottomY = (OBy - OByMinus1);
                double newCoordX = topX / bottomX;
                double newCoordY = topY / bottomY;

                double distanceX = (SPREADING_ALPHA * abs(newCoordX - originalX));
                double distanceY = (SPREADING_ALPHA * abs(newCoordY - originalY));

                // Set direction to move based on where the new coordinate was supposed to be compared to the current
                if (newCoordX < originalX)
                {
                    distanceX = -distanceX;
                }

                if (newCoordY < originalY)
                {
                    distanceY = -distanceY;
                }

                spreadedCellLocations.at(cell) = {cellPos.first + distanceX,
                                                 cellPos.second + distanceY};
            }

            unequalBins.push_back({NBxMinus1, NBx, NByMinus1, NBy, 0});
            std::cout << "Bin " << i << "old: (X: " << bin.lowerX << " to " << bin.upperX << ", Y: " << bin.lowerY << " to " << bin.upperY;
            std::cout << ") new: (X: " << NBxMinus1 << " to " << NBx << ", Y: " << NByMinus1 << " to " << NBy << ")" << std::endl;
        }
    }

    int AnalyticPlacer::getBinIndex(std::pair<double, double> pos, std::vector<Bin> &binsList)
    {
        int foundBin = -1;

        for (unsigned int i = 0; i < binsList.size(); i++)
        {
            auto bin = binsList.at(i);

            if (pos.first >= bin.lowerX && pos.first <= bin.upperX
                && pos.second >= bin.lowerY && pos.second <= bin.upperY)
            {
                foundBin = i;
                break;
            }
        }

        return foundBin;
    }

    // Get the bin a specific point belongs in
    Bin* AnalyticPlacer::getBin(std::pair<double, double> pos, std::vector<Bin> &binsList)
    {
        int index = this->getBinIndex(pos, binsList);
        if (index < 0)
        {
            return nullptr;
        }
        else
        {
            return &binsList.at(index);
        }
    }

    void AnalyticPlacer::updateBinUtilizations()
    {
        const int cellArea = 1; // Assuming each cell is 1x1

        // Iterate through each cell, determine which bin it is in and update that bin's area
        for (unsigned int i = 0; i < cellLocations.size(); i++)
        {
            // Only care about movable cells, no I/O pads or star nodes
            if (i < numCells_noPads)
            {
                {
                    Bin *bin = getBin(cellLocations.at(i), bins);

                    if (bin != nullptr) {
                        // Add the cell to this bin
                        bin->memberCells.push_back(i);
                    } else {
                        std::cerr << "Failed to determine bin for cell at " << cellLocations.at(i).first << ", "
                                  << cellLocations.at(i).second;
                        std::cerr << " Cell is probably outside the chip area" << std::endl;
                    }
                }
            }
        }

        for (auto &bin : bins)
        {
            double binArea = (bin.upperX - bin.lowerX) * (bin.upperY - bin.lowerY);
            bin.utilization = bin.memberCells.size() * cellArea / binArea;
        }
    }

    void AnalyticPlacer::doSpreading(std::string filePrefix)
    {
        this->spreadedCellLocations.resize(cellLocations.size());
        this->spreadedCellLocations.clear();

        // Only copy movable cell locations to our spreaded cell list
        for (unsigned int i = 0; i < cellLocations.size(); i++)
        {
            // Only care about movable cells, no I/O pads
            if (i < numCells_noPads + this->matrixQ->getStarNodeCount())
            {
                this->spreadedCellLocations.push_back(cellLocations.at(i));
            }
        }

        this->createBins();
        this->updateBinUtilizations();

        // Creates the uneqal bins and computes the spreaded cell locations
        this->calculateSpreadedCellLocations();

        this->saveSpreadedCellsToDisk(filePrefix + "_spread.kiaPad");
    }
}
