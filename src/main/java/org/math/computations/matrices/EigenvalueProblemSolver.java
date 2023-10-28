package org.math.computations.matrices;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Class that provides different methods for solving partial and full eigenvalue problem
 */
public class EigenvalueProblemSolver {

    private EigenvalueProblemSolver() {}

    public static List<Matrix> solveFullEigenvalueProblem(Matrix matrixA, double precision) {

        double[][] dataA = matrixA.toArray();
        int[] matrixSize = matrixA.getSize();
        double[][] dataX = Utils.getIdentityMatrix(matrixSize[1]).toArray();
        double[][] prevA = matrixA.toArray();

        double maxUpperElement = 0;

        int iterationCounter = 0;
        while (maxUpperElement >= precision || iterationCounter == 0) {
            iterationCounter++;

            maxUpperElement = 0;
            int[] indices = new int[2];
            for (int i = 0; i < matrixSize[0]; i++) {
                for (int j = i + 1; j < matrixSize[1]; j++) {
                    if (Math.abs(dataA[i][j]) > maxUpperElement) {
                        maxUpperElement = Math.abs(dataA[i][j]);
                        indices = new int[]{i, j};
                    }
                }
            }

            double aii = dataA[indices[0]][indices[0]];
            double ajj = dataA[indices[1]][indices[1]];
            double aij = dataA[indices[0]][indices[1]];
            double d = Math.sqrt(Math.pow(aii - ajj, 2) + 4 * Math.pow(aij, 2));

            double absDifference = Math.abs(aii - ajj) / d;
            double cos = Math.sqrt((1 + absDifference) / 2);
            double sin = Math.signum(aij * (aii - ajj))*Math.sqrt((1 - absDifference) / 2);

            for (int i = 0; i < matrixSize[0]; i++) {
                if (i == indices[0] || i == indices[1])
                    continue;

                dataA[i][indices[0]] = cos * prevA[i][indices[0]] + sin * prevA[i][indices[1]];
                dataA[indices[0]][i] = dataA[i][indices[0]];

                dataA[i][indices[1]] = cos * prevA[i][indices[1]] - sin * prevA[i][indices[0]];
                dataA[indices[1]][i] = dataA[i][indices[1]];
            }

            dataA[indices[0]][indices[0]] = cos*cos * aii + 2 * cos * sin * aij + sin*sin * ajj;
            dataA[indices[1]][indices[1]] = sin*sin * aii - 2 * cos * sin * aij + cos*cos * ajj;

            dataA[indices[0]][indices[1]] = 0;
            dataA[indices[1]][indices[0]] = 0;

            for (int i = 0; i < matrixSize[0]; i++) {
                double xi = dataX[i][indices[0]];
                double xj = dataX[i][indices[1]];
                dataX[i][indices[0]] = cos * xi + sin * xj;
                dataX[i][indices[1]] = cos * xj - sin * xi;
            }

            for (int i = 0; i < matrixSize[0]; i++) {
                prevA[i] = Arrays.copyOf(dataA[i], matrixSize[0]);
            }
        }

        double[][] eigenvaluesData = new double[matrixSize[0]][1];
        for (int i = 0; i < matrixSize[0]; i++) {
            eigenvaluesData[i][0] = dataA[i][i];
        }

        Matrix eigenvaluesVector = new DenseMatrix(matrixSize[0], 1, eigenvaluesData);
        Matrix eigenvectorsMatrix = new DenseMatrix(matrixSize[0], matrixSize[1], dataX).transpose();
        List<Matrix> eigenvectors = new ArrayList<>();
        for (int i = 0; i < matrixSize[1]; i++) {
            double[][] vector = new double[matrixSize[0]][1];
            for (int j = 0; j < matrixSize[0]; j++) {
                vector[j][0] = eigenvectorsMatrix.getElement(j, i);
            }
            eigenvectors.add(new DenseMatrix(matrixSize[0], 1, vector));
        }

        eigenvectors.add(0, eigenvaluesVector);
        return eigenvectors;
    }

    public double getLargestEigenvalue() {
        return 0;
    }

}
