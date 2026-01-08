//
// Created by Alejandro Zeise on 12/11/23.
//
#include "matrix.hpp"

#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cassert>

#define DEBUG

namespace PA3Placement
{
    /**
     * The threshold of which a degree of a hyperedge must exceed in order
     * to be represented with a Star model instead of a clique.
     */
    static const int STAR_MODEL_THRESHOLD = 3;

    template <typename T>
    T Matrix<T>::norm() const
    {
        assert(this->getHeight() >= 1);
        assert(this->getWidth() >= 1);

        T sum = this->get(0, 0) * this->get(0, 0);
        for (int i = 0; i < this->getHeight(); i++)
        {
            for (int j = 1; j < this->getWidth(); j++)
            {
                sum += this->get(j, i) * this->get(j, i);
            }
        }

        return sqrt(static_cast<double>(sum));
    }

    template <typename T>
    T Matrix<T>::dot(const Matrix<T> &rhs) const
    {
        // Check that the dot product is defined for these matrices
        assert(!(this->getWidth() != 1 || rhs.getWidth() != 1 || this->getHeight() != rhs.getHeight()));

        T result = this->get(0, 0) * rhs.get(0, 0);

        // Calculate the dot product
        for (unsigned int i = 0; i < this->getHeight(); ++i) {
            result += this->get(0, i) * rhs.get(0, i);
        }

        return result;
    }

    template <typename T>
    T* Matrix<T>::squashToRowMajorArray() const
    {
        T *array = new T[this->getHeight() * this->getWidth()];

        for (int i = 0; i < this->getHeight(); i++)
        {
            for (int j = 0; j < this->getWidth(); j++)
            {
                array[(i * this->getWidth()) + j] = this->get(j, i);
            }
        }

        return array;
    }

    template <typename T>
    Matrix2D<T>::Matrix2D(long width, long height)
    {
        this->data = new std::vector<T>(height * width);
        this->height = height;
        this->width = width;
    }

    template <typename T>
    Matrix2D<T>::Matrix2D(const Matrix2D<T> &other)
    {
        // Make a deep copy of the data structure

        this->height = other.height;
        this->width = other.width;

        this->data = new std::vector<T>(height * width);
        this->data->clear();
        std::copy(other.data->begin(), other.data->end(), std::back_inserter(*this->data));
    }

    template <typename T>
    Matrix2D<T>::~Matrix2D()
    {
        delete this->data;
    }

    template <typename T>
    Matrix2D<T>& Matrix2D<T>::operator=(const Matrix2D<T> &rhs)
    {
        if (this != &rhs) {
            // Make a deep copy of the data structure
            this->~Matrix2D();

            this->height = rhs.height;
            this->width = rhs.width;

            this->data = new std::vector<T>(height * width);
            this->data->clear();
            std::copy(rhs.data->begin(), rhs.data->end(), std::back_inserter(*this->data));
        }

        return *this;
    }

    template <typename T>
    Matrix2D<T> Matrix2D<T>::operator*(const Matrix2D<T> &rhs) const
    {
        Matrix2D<T> result = Matrix2D<T>(1, rhs.getHeight());

        assert(this->getWidth() == rhs.getHeight());

        for(int i = 0; i < this->getHeight(); i++) // for every row in A
        {
            for (int j = 0; j < rhs.getWidth(); j++) // for every column in B
            {
                T sum = result.get(j, i);
                for (int k = 0; k < this->getWidth(); k++) // for every row in B
                {
                    sum += this->get(k, i) * rhs.get(j, k);
                }
                result.set(j, i, sum);
            }
        }

        return result;
    }

    template <typename T>
    ColumnMatrix<T> Matrix2D<T>::operator*(const ColumnMatrix<T> &rhs) const
    {
        // Check if the matrices can be multiplied
        assert(this->getHeight() == rhs.getHeight());

        // Get dimensions of the matrices
        //int numRows = matrix.size();
        //int numCols = columnMatrix.size();

        // Initialize the result matrix with zeros
        //std::vector<std::vector<int>> result(numRows, std::vector<int>(1, 0));
        ColumnMatrix<T> result(this->getHeight(), 0);

        // Perform matrix multiplication
        for (int i = 0; i < this->getHeight(); ++i) {
            for (int j = 0; j < rhs.getHeight(); ++j) {
                T existing = result.get(0, i);
                T multiplied = this->get(j, i) * rhs.get(0, j);
                result.set(0, i, existing + multiplied);
                //result[i][0] += matrix[i][j] * columnMatrix[j];
            }
        }

        return result;
    }

    template <typename T>
    Matrix2D<T> Matrix2D<T>::operator*(const T &rhs) const
    {
        Matrix2D<T> result(this->getWidth(), this->getHeight());

        // Multiply each entry of the Matrix by the given scalar
        for (int i = 0; i < this->getHeight(); i++)
        {
            for (int j = 0; j < this->getWidth(); j++)
            {
                result.set(j, i, this->get(j, i) * rhs);
            }
        }

        return result;
    }

    template <typename T>
    Matrix2D<T> Matrix2D<T>::operator+(const Matrix2D<T> &rhs) const
    {
        Matrix2D<T> result = Matrix2D<T>(this->getWidth(), this->getHeight());

        for (int i = 0; i < this->getHeight(); ++i) {
            for (int j = 0; j < this->getWidth(); ++j) {
                result.set(j, i, this->get(j, i) + rhs.get(j, i));
            }
        }
        return result;
    }

    template <typename T>
    Matrix2D<T> Matrix2D<T>::operator-(const Matrix2D<T> &rhs) const
    {
        Matrix2D<T> result = Matrix2D<T>(this->getWidth(), this->getHeight());

        for (int i = 0; i < this->getHeight(); ++i) {
            for (int j = 0; j < this->getWidth(); ++j) {
                result.set(j, i, this->get(j, i) - rhs.get(j, i));
            }
        }
        return result;
    }

    template <typename T>
    Matrix2D<T> Matrix2D<T>::transpose() const
    {
        // Get the number of rows and columns in the original matrix
        int rows = this->getHeight();
        int cols = this->getWidth();

        // Initialize a new matrix with swapped rows and columns
        Matrix2D<T> result(rows, cols);

        // Calculate the transpose
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.set(i, j, this->get(j, i));
            }
        }

        return result;
    }

    template <typename T>
    void Matrix2D<T>::addRow()
    {
        this->height++;

        for (unsigned int i = 0; i < this->width; i++)
        {
            this->data->push_back(0);
        }
    }

    template <typename T>
    void Matrix2D<T>::resizeWidth(long newWidth)
    {
        // TODO: actually resize here
        // Currently the only place this is used is when the matrix is empty, so just reallocate everything

        this->~Matrix2D();
        this->width = newWidth;

        this->data = new std::vector<T>(newWidth * this->height);
    }

    template <typename T>
    void Matrix2D<T>::set(long x, long y, T value)
    {
        assert(x >= 0);
        assert(y >= 0);

        this->data->at((y * this->width) + x) = value;
    }

    template <typename T>
    T Matrix2D<T>::get(long x, long y) const
    {
        assert(x >= 0);
        assert(y >= 0);

        return this->data->at((y * this->width) + x);
    }

    template <typename T>
    void Matrix2D<T>::display(std::ostream &out) const
    {
        for(int i = 0; i < this->getHeight(); i++)
        {
            for (int j = 0; j < this->getWidth(); j++)
            {
                out << std::setw(3) << std::setprecision(3) << std::fixed << this->get(j, i) << " ";
            }
            out << std::endl;
        }
    }

    template <typename T>
    long Matrix2D<T>::getWidth() const
    {
        return this->width;
    }

    template <typename T>
    long Matrix2D<T>::getHeight() const
    {
        return this->height;
    }

    QMatrix::QMatrix(int numCellsNoPads, int numCellsAndPads, int numHyperedges, int *cellPinArray,
                     int *hEdgesToFirstMemberCellArray, int *hEdgeWeights) : Matrix2D<double>(1, numCellsNoPads) // We will be resizing the Matrix anyway later
    {
        this->numCellsNoPads = numCellsNoPads;
        this->numCellsAndPads = numCellsAndPads;
        this->numHyperedges = numHyperedges;
        this->cellPinArray = cellPinArray;
        this->hEdgesToFirstMemberCellArray = hEdgesToFirstMemberCellArray;
        this->hEdgeWeights = hEdgeWeights;

        // Will be resized later to account for added star nodes
        cellConnectionsList = WeightedCellConnectionsList(numCellsAndPads);

        this->generate();
    }

    int QMatrix::getStarNodeCount() const
    {
        return this->numStars;
    }

    QMatrix::WeightedCellConnectionsList* QMatrix::getCellConnectionsList()
    {
        return &this->cellConnectionsList;
    }

    void QMatrix::calculateNumberOfStarNodes()
    {
#ifdef DEBUG
        std::cout << "Calculating number of star nodes to add." << std::endl;
#endif
        /*
         * First iterate through each hyperedge, and determine it's degree
         * If the degree is greater than 3 then we use a star model (which adds a node),
         * otherwise we use clique.
         */
        for (int i = 0; i < this->numHyperedges; i++)
        {
            // Degree of the hyperedge is the element at index i subtracted from element at i + 1
            // The last element in the hyperedge array is the total number of pins
            int degree = this->hEdgesToFirstMemberCellArray[i + 1] - this->hEdgesToFirstMemberCellArray[i];
            if (degree > STAR_MODEL_THRESHOLD)
            {
                // Add row to matrix, since we are adding a node
                this->addRow();
                this->numStars++;
            }
            // Else we use clique, doesn't affect matrix dimensions
        }

#ifdef DEBUG
        std::cout << "Resizing Matrix width" << std::endl;
#endif

        // Resize matrix (and cellConnectionsList) due to our added star nodes
        this->resizeWidth(this->numCellsNoPads + this->numStars);
        this->cellConnectionsList.resize(this->numCellsAndPads + this->numStars);
    }

    void QMatrix::generate()
    {
        // First star node will be after the last cell
        const int firstStarIndex = this->numCellsNoPads;

        // Calculate the number of star nodes we will need, and resize ourself to accomodate
        // This also creates the required row/cols in the matrix for the new star nodes
        this->calculateNumberOfStarNodes();

        // Iterate through each hyperedge, now populating the Q matrix
        for (int i = 0, currentStar = firstStarIndex; i < this->numHyperedges; i++)
        {
            const int indexFirstCellPinInHyperedge = this->hEdgesToFirstMemberCellArray[i];
            const int degree = this->hEdgesToFirstMemberCellArray[i + 1] - indexFirstCellPinInHyperedge;
            const double weight = this->hEdgeWeights[i];

            const bool isClique = (degree <= STAR_MODEL_THRESHOLD);

            // Iterate through each cell or I/O pad in the hyperedge
            for (int j = indexFirstCellPinInHyperedge; j < indexFirstCellPinInHyperedge + degree; j++)
            {
                // Action depends on if this hyperedge is a clique or not,
                if (isClique)
                {
                    this->processCliqueHyperedge(degree, weight, j, indexFirstCellPinInHyperedge);
                }
                else
                {
                    this->processStarHyperedge(degree, weight, j, currentStar);
                }
            }

            // If the current hyperedge was modeled as a star node, increment the currentStar index
            if (!isClique)
            {
                currentStar++;
            }
        }

        this->populateDiagonal();
    }

    void QMatrix::processCliqueHyperedge(int degree, double weight, int currentCellIndex, int firstCellInHEdgeIndex)
    {
        const int currentCell = this->cellPinArray[currentCellIndex];

        // The current cell is an I/O pin if it's number is larger than the last moveable cell number
        const bool currentCellIsIOPin = (currentCell >= this->numCellsNoPads);

        // This is a clique, so the weight is different now (1/(k-1)) * w
        const double modifiedWeight = (1.0 / (degree - 1.0)) * weight;

        // Every cell in this hyperedge needs to connect to the others in it.
        for (int k = currentCellIndex + 1; k < firstCellInHEdgeIndex + degree; k++)
        {
            int connectedCell = cellPinArray[k];

            // Check if connected cell/connected cell is an I/O pin or not
            // Account for potential offset if the cell is an I/O pin or not
            int currentCellWithOffset = currentCell;
            int connectedCellWithOffset = connectedCell;

            if (currentCellIsIOPin)
            {
                // cell is an I/O pin, add offset to it to account for added star nodes
                currentCellWithOffset += numStars;
            }

            if (connectedCell >= this->numCellsNoPads)
            {
                // cell is an I/O pin, add offset to it to account for added star nodes
                connectedCellWithOffset += numStars;
            }

            // Connect the cell or I/O pads together in our cell connections list
            auto &currentConnections = cellConnectionsList.at(currentCellWithOffset);
            auto &connectedConnections = cellConnectionsList.at(connectedCellWithOffset);

            // Set comparator only looks at first element in pair, so can leave the "weight" as zero for the query
            auto currentToCon = currentConnections.find({connectedCellWithOffset, 0});
            auto conToCurrent = connectedConnections.find({currentCellWithOffset, 0});

            // Check if there is an already existing connection
            if (currentToCon != currentConnections.end())
            {
                // "currentCell" already connected to "connectedCell". Replace the connection with an updated weight
                std::pair<double, double> newPair = {currentToCon->first, currentToCon->second + modifiedWeight};
                currentConnections.erase(currentToCon);
                currentConnections.insert(newPair);
            }
            else
            {
                currentConnections.insert({connectedCellWithOffset, modifiedWeight});
            }

            if (conToCurrent != connectedConnections.end())
            {
                // "connectedCell" already connected to "currentCell". Replace the connection with an updated weight
                std::pair<double, double> newPair = {conToCurrent->first, conToCurrent->second + modifiedWeight};
                connectedConnections.erase(conToCurrent);
                connectedConnections.insert(newPair);
            }
            else
            {
                connectedConnections.insert({currentCellWithOffset, modifiedWeight});
            }


            if (!currentCellIsIOPin && connectedCell < this->numCellsNoPads) // Avoid I/O pins
            {
                // Connect the cells together in the matrix (since they are not I/O pins)
                this->set(currentCell, connectedCell, this->get(currentCell, connectedCell) + -modifiedWeight);
                this->set(connectedCell, currentCell, this->get(connectedCell, currentCell) + -modifiedWeight);
            }
        }
    }

    void QMatrix::processStarHyperedge(int degree, double weight, int currentCellIndex, int currentStar)
    {
        const int currentCell = this->cellPinArray[currentCellIndex];

        // The current cell is an I/O pin if it's number is larger than the last moveable cell number
        const bool currentCellIsIOPin = (currentCell >= this->numCellsNoPads);

        // Hyperedge to be modeled as a Star
        const double modifiedWeight = (1.0 * degree / (degree - 1.0)) * weight;

        // Connect the current cell or I/O pad to the current star
        if (currentCellIsIOPin)
        {
            // Add offset to currentCell to account for added star nodes
            cellConnectionsList.at(currentStar).insert({currentCell + numStars, modifiedWeight});
            cellConnectionsList.at(currentCell + numStars).insert({currentStar, modifiedWeight});
        }
        else
        {
            cellConnectionsList.at(currentStar).insert({currentCell, modifiedWeight});
            cellConnectionsList.at(currentCell).insert({currentStar, modifiedWeight});

            // Current cell is NOT an I/O pad, so connect the cell to the star in the Matrix
            this->set(currentCell, currentStar, this->get(currentCell, currentStar) + -modifiedWeight);
            this->set(currentStar, currentCell, this->get(currentStar, currentCell) + -modifiedWeight);
        }
    }

    void QMatrix::populateDiagonal()
    {
        // Populate the diagonal of the Q matrix based on cell connections
        for (unsigned int i = 0; i < cellConnectionsList.size(); i++)
        {
            auto &ithCellConnections = cellConnectionsList.at(i);

            std::for_each(ithCellConnections.begin(), ithCellConnections.end(), [this, i] (auto connectedCellData)
                  {
                    // First element of pair is the connected cell number, second element is the weight of the connection

                      // Exclude I/O pins from Q matrix
                      // We include them as part of the weights (hence only checking i), but they are not actual rows/columns in the matrix
                      if (i < (this->numCellsNoPads + numStars))
                      {
                          this->set(i, i, this->get(i, i) + connectedCellData.second);
                      }
                  }
            );
        }
    }

    template <typename T>
    ColumnMatrix<T>::ColumnMatrix(std::vector<T> &data)
    {
        this->rows = new std::vector<T>(data.size());
        this->rows->clear();

        std::copy(data.begin(), data.end(), std::back_inserter(*this->rows));
    }

    template <typename T>
    ColumnMatrix<T>::ColumnMatrix(long rowCount, T defaultValue)
    {
        this->rows = new std::vector<T>(rowCount);
        this->defaultValue = defaultValue;
    }

    template <typename T>
    ColumnMatrix<T>::ColumnMatrix(const ColumnMatrix<T> &other)
    {
        // Make a deep copy of the data structure
        this->rows = new std::vector<T>(other.getHeight());
        this->rows->clear();

        std::copy(other.rows->begin(), other.rows->end(), std::back_inserter(*this->rows));
    }

    template <typename T>
    ColumnMatrix<T>::~ColumnMatrix()
    {
        delete this->rows;
    }

    template <typename T>
    ColumnMatrix<T>& ColumnMatrix<T>::operator=(const ColumnMatrix<T> &rhs)
    {
        if (this != &rhs) {
            // Make a deep copy of the data structure

            // Resize our rows to match the target amount
            this->rows->resize(rhs.getHeight());

            // Clear out old data
            this->rows->clear();

            // Copy the data from RHS
            std::copy(rhs.rows->begin(), rhs.rows->end(), std::back_inserter(*this->rows));
        }

        return *this;
    }

    template <typename T>
    ColumnMatrix<T> ColumnMatrix<T>::operator*(const ColumnMatrix<T> &rhs) const
    {
        assert(this->getHeight() == rhs.getHeight());

        // Initialize the result vector with zeros
        ColumnMatrix<T> result(this->getHeight(), 0);

        // Perform element-wise multiplication
        for (int i = 0; i < this->getHeight(); ++i) {
            result.set(0, i, this->get(0, i) * rhs.get(0, i));
        }

        return result;
    }

    template <typename T>
    ColumnMatrix<T> ColumnMatrix<T>::operator*(const T &rhs) const
    {
        // Initialize the result vector with zeros
        ColumnMatrix<T> result(this->getHeight(), 0);

        // Perform element-wise multiplication
        for (int i = 0; i < this->getHeight(); ++i) {
            result.set(0, i, this->get(0, i) * rhs);
        }

        return result;
    }

    template <typename T>
    ColumnMatrix<T> ColumnMatrix<T>::operator+(const ColumnMatrix<T> &rhs) const
    {
        // Initialize the result vector with zeros
        ColumnMatrix<T> result(this->getHeight(), 0);

        // Perform element-wise multiplication
        for (int i = 0; i < this->getHeight(); ++i) {
            result.set(0, i, this->get(0, i) + rhs.get(0, i));
        }

        return result;
    }

    template <typename T>
    ColumnMatrix<T> ColumnMatrix<T>::operator-(const ColumnMatrix<T> &rhs) const
    {
        // Initialize the result vector with zeros
        ColumnMatrix<T> result(this->getHeight(), 0);

        // Perform element-wise multiplication
        for (int i = 0; i < this->getHeight(); ++i) {
            result.set(0, i, this->get(0, i) - rhs.get(0, i));
        }

        return result;
    }

    template <typename T>
    void ColumnMatrix<T>::addRow()
    {
        this->rows->push_back(this->defaultValue);
    }

    template <typename T>
    void ColumnMatrix<T>::set(long x, long y, T value)
    {
        assert(x == 0);
        assert(y >= 0);

        this->rows->at(y) = value;
    }

    template <typename T>
    T ColumnMatrix<T>::get(long x, long y) const
    {
        assert(x == 0);
        assert(y >= 0);

        return this->rows->at(y);
    }

    template <typename T>
    void ColumnMatrix<T>::display(std::ostream &out) const
    {
        out << std::setprecision(3) << std::fixed;
        for (int i = 0; i < this->getHeight(); i++)
        {
            out << this->get(0, i) << " ";
        }
        out << std::endl;
    }

    template <typename T>
    long ColumnMatrix<T>::getHeight() const
    {
        return this->rows->size();
    }

    template <typename T>
    long ColumnMatrix<T>::getWidth() const
    {
        return 1;
    }

    DMatrix::DMatrix(Dimension dimension, const SPinLocation *pinLocations, int numCellsNoPads, int numStars, QMatrix::WeightedCellConnectionsList *cellConnectionsList)
        : ColumnMatrix<double>(numCellsNoPads + numStars, 0)
        // Number of rows in the D matrix is equal to the number of cells and star nodes
        // (everything except I/O pads, they are not present as rows or columns in any matrix)
    {
        this->dimension = dimension;
        this->pinLocations = pinLocations;
        this->numCellsNoPads = numCellsNoPads;
        this->numStars = numStars;
        this->cellConnectionsList = cellConnectionsList;

        this->generate();
    }

    const SPinLocation* DMatrix::getIOPadLocation(int ioPadNumberWithOffset) const
    {
        // According to suraj_parser.cpp: the pin location for an I/O pad is (padNumber - numCells_noPads)
        // but due to our offset we must subtract that as well (added star nodes)
        return this->pinLocations + (ioPadNumberWithOffset - this->numStars - this->numCellsNoPads);
    }

    void DMatrix::generate()
    {
        for (unsigned int currentCell = 0; currentCell < this->cellConnectionsList->size(); currentCell++)
        {
            // Only looking at non I/O pads as our "current cell"
            if (currentCell < this->numCellsNoPads + this->numStars) {
                const auto &connectedCells = this->cellConnectionsList->at(currentCell);

                for (auto &connectedCell: connectedCells) {
                    // Check if the connected cell is an I/O pad.
                    // The cell numbers in the cellConnectionsList have an offset based on how many star nodes have been added
                    // So, an I/O pad is an I/O pad if it's number is greater than the last movable cell plus the number of stars
                    // (since the next numbers after the last movable cell are the added star nodes)
                    if (connectedCell.first >= this->numCellsNoPads + this->numStars) {
                        // Get the x,y coordinates of the I/O pad
                        const SPinLocation *ioPadLoc = this->getIOPadLocation(connectedCell.first);
                        // Pick the corresponding coordinate depending on our dimension
                        int coord = (this->dimension == X ? ioPadLoc->x : ioPadLoc->y);

                        // Update the matrix row for our current cell with the I/O pad coordinate
                        // Summing here because multiple pads might be connected to a cell, and they should be summed in this case
                        // Multiply the coordinate value by the weight of the connection between the pad and the current cell
                        this->set(0, currentCell, this->get(0, currentCell) + (coord * connectedCell.second));
                    }
                }
            }
        }
    }

    // Explicitly Instantiate Templates
    template class Matrix2D<double>;
    template class Matrix2D<int>;
    template class ColumnMatrix<double>;
    template class ColumnMatrix<int>;

    // Based on ChatGPT code
    Matrix2D<double> solveMatrixGradientDescent(double learningRate, int iterations, const Matrix2D<double> &Q, const Matrix2D<double> &Dx)
    {
        Matrix2D<double> x(Dx.getWidth(), Q.getWidth()); // Initialize d with zeros

        for (int iter = 0; iter < iterations; ++iter) {
            Matrix2D<double> Qd = Q * x;  // Q * x
            Matrix2D<double> error = Qd - Dx;  // Error term

            // Update d using gradient descent
            for (int i = 0; i < Q.getWidth(); ++i) {
                for (int j = 0; j < Dx.getWidth(); ++j) {
                    x.set(j, i, x.get(j, i) - (learningRate * error.get(j, i)));
                }
            }
#ifdef DEBUG
            std::cout << "[Gradient Descent Solver]: Completed iteration: " << iter << std::endl;
#endif
        }

        return x;
    }

    // Based on ChatGPT code
    ColumnMatrix<double> solveMatrixConjugateGradient(double tolerance, int iterations, const Matrix2D<double> &Q, const ColumnMatrix<double> &Dx)
    {
#if 0
        // Initialize the solution vector x
        Matrix2D<double> x(Dx.getWidth(), Dx.getHeight());

        Matrix2D<double> r = Q * x - Dx;
        Matrix2D<double> p = r;

        for (int k = 0; k < iterations; ++k) {
            double alpha = r.dot(r) / p.dot(Q * p);
            x = x - p * alpha;
            Matrix2D<double> new_r = Q * x - Dx;

            if (new_r.norm() < tolerance) {
#ifdef DEBUG
                std::cout << "[Conjugate Gradient Solver]: Converged in " << k + 1 << " iterations." << std::endl;
#endif
                break;
            }

            double beta = new_r.dot(new_r) / r.dot(r);
            p = new_r + p * beta;
            r = new_r;

#ifdef DEBUG
            std::cout << "[Conjugate Gradient Solver]: Completed iteration: " << k << std::endl;
#endif
        }

        return x;
#else
        long n = Q.getHeight();
        ColumnMatrix<double> Ap(n , 0);
        ColumnMatrix<double> x(n, 0); // col

        ColumnMatrix<double> r = Dx - Q * x; // col
        ColumnMatrix<double> p = r; // col

        //double rs_old = inner_product(r, r);
        double rs_old = r.dot(r);

        for (int i = 0; i < iterations; i++) {
            Ap = Q * p; // col
            double alpha = rs_old / p.dot(Ap);
            x = x + (p * alpha);
            r = r - (Ap * alpha);
            double rs_new = r.dot(r);
            if (sqrt(rs_new) < tolerance) break;
            p = r + p * (rs_new / rs_old);
            rs_old = rs_new;
#ifdef DEBUG
            std::cout << "[Conjugate Gradient Solver]: Completed Iteration " << i << std::endl;
#endif
        }

        return x;
#endif
    }

    void matrixSolverThread(void *args)
    {
        MatrixSolverParams *params = (MatrixSolverParams*) args;

        std::cout << "[Matrix Solver Thread " << params->id << "]: Started." << std::endl;

#if 1
        *(params->xAnswer) = solveMatrixConjugateGradient(params->tolerance, params->maxIterations, *params->Q, *params->Dx);
        //*(params->xAnswer) = acceleratedSolveMatrixConjugateGradient(*params->Q, *params->Dx);
#else
        *(params->xAnswer) = solveMatrixGradientDescent(0.01, params->maxIterations, *params->Q, *params->Dx);
#endif

        std::cout << "[Matrix Solver Thread " << params->id << "]: Exit." << std::endl;
    }
}
