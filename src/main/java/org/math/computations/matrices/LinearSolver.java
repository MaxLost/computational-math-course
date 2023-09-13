package org.math.computations.matrices;

import java.util.Arrays;
import java.util.List;

/**
 * Class for finding unique solution X for linear system of equations: A @ X = B.
 */
public class LinearSolver {

  private final DenseMatrix matrixA;
  private final DenseMatrix matrixB;

  /**
   * Initialization of linear system solver for finding unique solution X: A @ X = B.
   * <br>
   * Where A - N x N matrix, X - N x 1 vector, B - N x 1 vector
   *
   * @param matrixA N x N matrix
   * @param matrixB N x 1 vector
   */
  public LinearSolver(DenseMatrix matrixA, DenseMatrix matrixB) {
    this.matrixA = matrixA;
    this.matrixB = (DenseMatrix) matrixB.transpose();
  }

  /**
   * Finds unique solution linear system: A @ X = B.
   * <br>
   * Where A - N x N matrix and det(A) != 0, X - N x 1 vector, B - N x 1 vector
   *
   * @return X vector as a list
   */
  public List<Double> solve() {

    double[] b = matrixB.toArray()[0];

    Double[] x = new Double[b.length];

    double det = matrixA.det();

    if (Math.abs(det) < 10e-32) {
      throw new RuntimeException(
          "Finding unique solution for this system is impossible because det(A) = 0");
    }

    double[][] matrix = matrixA.toArray();
    double sum = 0;

    for (int i = 0; i < matrixA.colCount - 1; i++) {
      for (int j = i + 1; j < matrixA.rowCount; j++) {
        for (int k = i + 1; k < matrixA.colCount; k++) {
          matrix[j][k] -= (matrix[i][k] * matrix[j][i]) / matrix[i][i];
        }
      }
    }

    x[0] = b[0];
    for (int i = 1; i < matrixA.colCount; i++) {
      sum = 0;
      for (int j = 0; j < i; j++) {
        sum += matrix[i][j] * x[j];
      }
      x[i] = b[i] - sum;
    }

    for (int i = matrixA.colCount - 1; i >= 0; i--) {
      sum = 0;
      for (int j = i + 1; j < matrixA.colCount; j++) {
        sum += matrix[i][j] * x[j];
      }
      x[i] = (x[i] - sum) / matrix[i][i];
    }

    return Arrays.asList(x);
  }

}
