package org.math.computations.matrices;

/**
 * Class that provides methods for various matrix decompositions.
 */
public class Decomposer {

  private Decomposer() {}

  /**
   * Method that performs LU decomposition of matrix.
   *
   * @param matrix matrix for which decomposition will be performed
   * @return Array of two matrices:
   *        - First matrix is lower-triangular matrix L, main diagonal elements equals 1
   *        - Second matrix is upper-triangular matrix U
   */
  public static Matrix[] getLuDecomposition(Matrix matrix) {

    int[] matrixSize = matrix.getSize();
    double[][] data = matrix.toArray();

    int lowerColCount = matrixSize[0] != matrixSize[1] ? matrixSize[1] - 1 : matrixSize[1];

    double[][] lowerMatrixData = new double[matrixSize[0]][lowerColCount];
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
