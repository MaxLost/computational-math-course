package org.math.computations.matrices;

public class Utils {

    public static Matrix getIdentityMatrix(int size) {

        double[][] matrixData = new double[size][size];
        for (int i = 0; i < size; i++) {
            matrixData[i][i] = 1;
        }

        return new DenseMatrix(size, size, matrixData);
    }
}