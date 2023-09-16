package org.math.computations.matrices;

public class Decomposer {

  public static Matrix[] getLuDecomposition(Matrix matrix) {

    int[] matrixSize = matrix.getSize();
    double[][] data = matrix.toArray();

    double[][] lowerMatrixData = new double[matrixSize[0]][matrixSize[1] - 1];
    double[][] upperMatrixData = new double[matrixSize[0]][matrixSize[1]];

    double sum;

    for (int i = 0; i < matrixSize[0]; i++) {
      for (int k = i; k < matrixSize[1]; k++) {
        sum = 0;
        for (int j = 0; j < i; j++) {
          sum += (lowerMatrixData[i][j] * upperMatrixData[j][k]);
        }

        upperMatrixData[i][k] = data[i][k] - sum;
      }

      for (int k = i; k < matrixSize[0]; k++) {
        if (i == k) {
          lowerMatrixData[i][i] = 1;
        } else {
          sum = 0;
          for (int j = 0; j < i; j++) {
            sum += (lowerMatrixData[k][j] * upperMatrixData[j][i]);
          }

          lowerMatrixData[k][i] = (data[k][i] - sum) / upperMatrixData[i][i];
        }
      }
    }

    return new Matrix[]{new DenseMatrix(matrixSize[0], matrixSize[0], lowerMatrixData),
        new DenseMatrix(matrixSize[0], matrixSize[1], upperMatrixData)};
  }

}
