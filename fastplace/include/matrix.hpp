//
// Created by Alejandro Zeise on 12/11/23.
//

#ifndef PA3ANALYTICPLACEMENT_MATRIX_HPP
#define PA3ANALYTICPLACEMENT_MATRIX_HPP

#include "util.hpp"
#include "suraj_parser.h"

#include <vector>
#include <set>
#include <iostream>
#include <pthread.h>

namespace PA3Placement
{
    /**
     * Abstract base class for a matrix to hold values of type T.
     * @tparam T Type of values held in this matrix (usually int or double)
     */
    template <typename T>
    class Matrix
    {
    public:

        /**
         * Adds a new row to the matrix keeping the same width
         * as all other rows.
         */
        virtual void addRow() = 0;

        /**
         * Displays the matrix to the provided ostream
         * @param out output stream to display the matrix to
         */
        virtual void display(std::ostream &out) const = 0;

        /**
         * Sets a value in the matrix to the given value
         * @param x The row number (starting at 0) of the value
         * @param y The column number (starting at 0) of the value
         * @param value The value to be set at that location.
         */
        virtual void set(long x, long y, T value) = 0;

        /**
         * Gets a value in the matrix from a specific position.
         * @param x The row number (starting at 0) of the value.
         * @param y The column number (starting at 0) of the value.
         * @return The value
         */
        virtual T get(long x, long y) const = 0;

        /**
         * Get the height of the matrix (number of rows)
         * @return Number of rows in the matrix.
         */
        virtual long getHeight() const = 0;

        /**
         * Get the width of the matrix (number of columns)
         * @return Number of columns in the matrix.
         */
        virtual long getWidth() const = 0;

        /**
         * Calculates the Frobenius norm of this Matrix,
         * which is the square root of the sum of the squares of all elements
         * in the Matrix.
         *
         * May be different for a Column Matrix or Vector
         *
         * @return The Frobenius norm.
         */
        virtual T norm() const;

        /**
         * Compute the dot product of this Matrix with another
         * @param rhs The other Matrix to compute the dot product with
         * @return The dot product
         */
        virtual T dot(const Matrix<T> &rhs) const;

        /**
         * Creates a copy of this column matrix as a raw array of type T elements
         * on the heap. The array is 1 dimensional, stored in row-major order.
         *
         * Use getHeight and getWidth in order to calculate indices properly
         *
         * @return A copy of the matrix as a row-major ordered array.
         */
        virtual T* squashToRowMajorArray() const;
    };

    template <typename T>
    class ColumnMatrix : public Matrix<T>
    {
    private:
        T defaultValue;

        std::vector<T> *rows;
    public:
        ColumnMatrix(std::vector<T> &data);
        ColumnMatrix(long rowCount, T defaultValue);
        ColumnMatrix(const ColumnMatrix<T> &other);
        ~ColumnMatrix();

        ColumnMatrix<T>& operator=(const ColumnMatrix<T> &rhs);

        // Math operations
        ColumnMatrix<T> operator*(const ColumnMatrix<T> &rhs) const;
        ColumnMatrix<T> operator*(const T &rhs) const;
        ColumnMatrix<T> operator+(const ColumnMatrix<T> &rhs) const;
        ColumnMatrix<T> operator-(const ColumnMatrix<T> &rhs) const;

        // Implemented the following from Matrix
        virtual void addRow();
        virtual void display(std::ostream &out) const;
        virtual void set(long x, long y, T value);
        virtual T get(long x, long y) const;
        virtual long getHeight() const;
        virtual long getWidth() const;
    };


    template <typename T>
    class Matrix2D : public Matrix<T>
    {
    private:
        std::vector<T> *data;
        size_t height;
        size_t width;
    public:
        Matrix2D(long width, long height);
        Matrix2D(const Matrix2D<T> &other);
        ~Matrix2D();

        /**
         * Resizes all rows of the matrix to have a new width.
         * If the width is smaller than the  current width nothing will happen,
         * otherwise new columns will be added until the desired newWidth is reached.
         *
         * @param newWidth
         */
        virtual void resizeWidth(long newWidth);

        /**
         * Calculate the transpose of this Matrix.
         * @return The transpose of this matrix.
         */
        virtual Matrix2D<T> transpose() const;

        Matrix2D<T>& operator=(const Matrix2D<T> &rhs);

        // Math operations
        Matrix2D<T> operator*(const Matrix2D<T> &rhs) const;
        ColumnMatrix<T> operator*(const ColumnMatrix<T> &rhs) const;
        Matrix2D<T> operator*(const T &rhs) const;
        Matrix2D<T> operator+(const Matrix2D<T> &rhs) const;
        Matrix2D<T> operator-(const Matrix2D<T> &rhs) const;

        // Implemented the following from Matrix
        virtual void addRow();
        virtual void display(std::ostream &out) const;
        virtual void set(long x, long y, T value);
        virtual T get(long x, long y) const;
        virtual long getHeight() const;
        virtual long getWidth() const;
    };

    template <typename T>
    std::ostream& operator<<(std::ostream &out, const Matrix<T> &mat)
    {
        mat.display(out);
        return out;
    }

    class QMatrix : public Matrix2D<double>
    {
    public:
        typedef std::vector<std::set<std::pair<int, double>, firstElementOfPairComparator>> WeightedCellConnectionsList;

        QMatrix(int numCellsNoPads, int numCellsAndPads, int numHyperedges, int *cellPinArray, int *hEdgesToFirstMemberCellArray, int *hEdgeWeights);

        int getStarNodeCount() const;
        WeightedCellConnectionsList* getCellConnectionsList();
    private:
        int numCellsNoPads;
        int numCellsAndPads;
        int numHyperedges;

        int *cellPinArray;
        int *hEdgesToFirstMemberCellArray;
        int *hEdgeWeights;

        int numStars = 0;

        WeightedCellConnectionsList cellConnectionsList;

        void calculateNumberOfStarNodes();
        void generate();
        void populateDiagonal();
        void processCliqueHyperedge(int degree, double weight, int currentCellIndex, int firstCellInHEdgeIndex);
        void processStarHyperedge(int degree, double weight, int currentCellIndex, int currentStar);
    };

    class DMatrix : public ColumnMatrix<double>
    {
    public:
        enum Dimension
        {
            X,
            Y
        };

        /**
         * Dx/Dy matrix have only a single column
         *
         * The dimension (either X or Y) is denoted by the dimension parameter
         */
        DMatrix(Dimension dimension, const SPinLocation *pinLocations, int numCellsNoPads, int numStars, QMatrix::WeightedCellConnectionsList *cellConnectionsList);

        /**
         * Get the X/Y location of an I/O pad provided it's number, but
         * it's number with the offset due to added star nodes.
         *
         * This is used due to numbers from the cellConnectionsList.
         * @param ioPadNumberWithOffset
         * @return
         */
        const SPinLocation* getIOPadLocation(int ioPadNumberWithOffset) const;
    private:
        Dimension dimension;
        const SPinLocation *pinLocations;

        /**
         * List to keep track of which cells are connected to others.
         * The index is the cell number, value is a set of pairs of
         * cell numbers the cell is connected to (first) with its corresponding weight (second).
         *
         * Cell numbers will be slightly different than those in the hyperedge datastructures
         * All normal cells will be the same
         * Added star Nodes will be placed right after the last cell number, pushing I/O pad numbers
         * to start after the last star node number.
         *
         * Will incorporate both clique and star models.
         */
        QMatrix::WeightedCellConnectionsList *cellConnectionsList;

        int numCellsNoPads;
        int numStars;

        void generate();
    };

    struct MatrixSolverParams
    {
        int id;

        double tolerance;
        int maxIterations;

        Matrix2D<double> *Q;
        ColumnMatrix<double> *Dx;

        ColumnMatrix<double> *xAnswer;
    };

    /**
     * Solves the equation of the form Q*x = Dx
     * where Q, x, and Dx are all matrices. x is the matrix that is unknown.
     *
     * Uses gradient descent method.
     *
     * Credit to ChatGPT for writing this function, I modified it to work with my Matrix class(es).
     *
     * @param learningRate
     * @param iterations
     * @return x from the equation, the solved matrix
     */
    Matrix2D<double> solveMatrixGradientDescent(double learningRate, int iterations, const Matrix2D<double> &Q, const Matrix2D<double> &Dx);

    /**
     * Solves the equation of the form Q*x = Dx
     * where Q, x, and Dx are all matrices. x is the matrix that is unknown.
     *
     * Uses Conjugate Gradient method.
     *
     * Credit to ChatGPT for writing this function, I modified it to work with my Matrix class(es).
     *
     * @param tolerance
     * @param iterations
     * @return x from the equation, the solved matrix
     */
    ColumnMatrix<double> solveMatrixConjugateGradient(double tolerance, int iterations, const Matrix2D<double> &Q, const ColumnMatrix<double> &Dx);

    void matrixSolverThread(void *args);
}

#endif //PA3ANALYTICPLACEMENT_MATRIX_HPP
