//
// Created by Alejandro Zeise on 12/11/23.
//

#ifndef PA3ANALYTICPLACEMENT_PLACER_HPP
#define PA3ANALYTICPLACEMENT_PLACER_HPP

#include <vector>
#include "matrix.hpp"

namespace PA3Placement
{
    struct Bin
    {
        double lowerX;
        double upperX;
        double lowerY;
        double upperY;

        int utilization;

        std::vector<int> memberCells;
    };

    class AnalyticPlacer
    {
    private:
        QMatrix *matrixQ;
        DMatrix *matrixDx;
        DMatrix *matrixDy;

        std::vector<std::pair<double, double>> cellLocations;
        std::vector<std::pair<double, double>> spreadedCellLocations;

        std::vector<Bin> bins;

        /**
         * Obtain the sum of all wirelengths in the circuit.
         * Calculated by summing the weight of each hyperedge,
         * which is modified itself due to being either a clique or star.
         * @return The sum of the weights of all hyperedges
         */
        double calculateTotalWirelength(std::vector<std::pair<double, double>> &cellLocations) const;

        /**
         * Solves the equation for the X and Y coordinates of each cell.
         * The equation is in the form of Q * x = Dx for x coordinates
         * and Q * y = Dy for y coordinates, where Q and Dx are known, but x is unknown (same for y variant)
         */
        void calculateCellLocations();

        /**
         * Calculates the dimensions of the chip by finding the maximum I/O pad or cell X coordinate
         * and maximum I/O pad or cell Y coordinate
         * @return Dimensions of the chip as a pair, where X is the first element and Y is the second
         */
        std::pair<double, double> calculateChipDimensions();

        void saveCellLocationsToDisk(std::string filename);
        void saveSpreadedCellsToDisk(std::string filename);

        void doSpreading(std::string filePrefix);
        void calculateSpreadedCellLocations();
        void createBins();
        void updateBinUtilizations();

        int getBinIndex(std::pair<double, double> coordinates, std::vector<Bin> &binsList);
        Bin* getBin(std::pair<double, double> coordinates, std::vector<Bin> &binsList);
    public:
        void doPlacement(std::string filePrefix);
    };
}

#endif //PA3ANALYTICPLACEMENT_PLACER_HPP
