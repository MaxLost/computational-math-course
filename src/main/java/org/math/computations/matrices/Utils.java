package org.math.computations.matrices;

public class Utils {

    public static Matrix getIdentityMatrix(int size) {

        double[][] matrixData = new double[size][size];
        for (int i = 0; i < size; i++) {
            matrixData[i][i] = 1;
        }

        return new DenseMatrix(size, size, matrixData);
    }

    public static double euclideanMean(Matrix matrix) {

        int[] matrixSize = matrix.getSize();
        double mean = 0;

        for (int i = 0; i < matrixSize[0]; i++) {
            for (int j = 0; j < matrixSize[1]; j++) {
                mean += Math.pow(matrix.getElement(j, i), 2);
            }
        }
        mean = Math.sqrt(mean);

        return mean;
    }

    public static double vectorMaxMean(Matrix vector) {

        int[] matrixSize = vector.getSize();
        double[][] data = vector.toArray();
        double mean = 0;

        for (int i = 0; i < matrixSize[0]; i++) {
            mean = Math.max(mean, Math.abs(data[i][0]));
        }

        return mean;
    }
}