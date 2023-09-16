package org.math.computations.matrices;

import java.util.Arrays;

/**
 * Class for finding unique solution X for linear system of equations: A @ X = B.
 */
public class LinearSystemSolver {

  private final Matrix augmentedMatrix;

  /**
   * Initialization of linear system solver for finding unique solution X: A @ X = B.
   * <br>
   * Where A - N x N matrix, X - N x 1 vector, B - N x 1 vector
   *
   * @param augmentedMatrix N x N+1 augmented matrix of linear system A @ X = B
   */
  public LinearSystemSolver(Matrix augmentedMatrix) {
    this.augmentedMatrix = augmentedMatrix;
  }

  /**
   * Method for getting augmented matrix of a system A @ X = B.
   * Where A - N x N, B - N x 1, X - N x 1
   *
   * @param matrixA matrix A from system denoted above
   * @param b vector B from system denoted above
   *
   * @return augmented matrix of a system A @ X = B
   */
  public static Matrix getAugmentedMatrix(Matrix matrixA, Matrix b) {

    int[] matrixSize = matrixA.getSize();
    double[][] data = new double[matrixSize[0]][matrixSize[1] + 1];

    for (int i = 0; i < matrixSize[0]; i++) {
      for (int j = 0; j < matrixSize[1]; j++) {
        data[i][j] = matrixA.getElement(j, i);
      }

      data[i][matrixSize[1]] = b.getElement(0, i);
    }

    return new DenseMatrix(matrixSize[0], matrixSize[1] + 1, data);
  }

  /**
   * Finds unique solution of linear system: A @ X = B.
   *
   * @param method method of solving system: <br>
   *               - Default method: LU decomposition <br>
   *               - DG: Gauss elimination <br>
   *               - IG: Gauss elimination with selection of leading element <br>
   *
   * @return solution as vector
   */
  public Matrix solve(String method) {

    return switch (method) {
      case "DG" -> solveDefaultGauss();
      case "IG" -> solveImprovedGauss();
      default -> solveLuDecomposition();
    };
  }

  private Matrix solveDefaultGauss() {

    int[] matrixSize = augmentedMatrix.getSize();
    double[][] data = augmentedMatrix.toArray();

    boolean isSmallLeadingElement = false;

    for (int i = 0; i < matrixSize[1] - 1; i++) {

      double leadingElement = data[i][i];
      if (Math.abs(leadingElement) < 10e-8) {
        isSmallLeadingElement = true;
      }

      for (int j = i + 1; j < matrixSize[0]; j++) {

        double c = data[j][i] / data[i][i];
        for (int k = i; k < matrixSize[1]; k++) {

          data[j][k] = data[j][k] - data[i][k] * c;
        }
      }
    }

    double[][] finalData = doGaussBackwardElimination(data);

    if (isSmallLeadingElement) {
      throw new RuntimeException("Result may be incorrect due to small leading elements");
    }

    return new DenseMatrix(matrixSize[0], 1, extractLastColumn(finalData));
  }

  private Matrix solveImprovedGauss() {

    int[] matrixSize = augmentedMatrix.getSize();
    double[][] data = augmentedMatrix.toArray();

    for (int i = 0; i < matrixSize[1] - 1; i++) {

      int leadingElementRow = i;
      for (int j = i; j < matrixSize[0]; j++) {

        if (Math.abs(data[j][i]) > Math.abs(data[leadingElementRow][i])) {
          leadingElementRow = j;
        }
      }

      double[] temp = data[i];
      data[i] = data[leadingElementRow];
      data[leadingElementRow] = temp;

      for (int j = i + 1; j < matrixSize[0]; j++) {

        double c = data[j][i] / data[i][i];
        for (int k = i; k < matrixSize[1]; k++) {

          data[j][k] = data[j][k] - data[i][k] * c;
        }
      }
    }

    double[][] finalData = doGaussBackwardElimination(data);

    return new DenseMatrix(matrixSize[0], 1, extractLastColumn(finalData));
  }

  private Matrix solveLuDecomposition() {

    Matrix[] decomposed = Decomposer.getLuDecomposition(augmentedMatrix);

    double[][] upperMatrix = decomposed[1].toArray();
    int[] matrixSize = new int[]{upperMatrix.length, upperMatrix[0].length};

    double[][] data = doGaussBackwardElimination(upperMatrix);

    return new DenseMatrix(matrixSize[0], 1, extractLastColumn(data));
  }

  private static double[][] extractLastColumn(double[][] data) {

    double[][] result = new double[data.length][1];
    for (int i = 0; i < data.length; i++) {
      result[i][0] = data[i][data[0].length - 1];
    }

    return result;
  }

  private static double[][] doGaussBackwardElimination(double[][] data) {

    int[] matrixSize = new int[]{data.length, data[0].length};

    for (int i = matrixSize[0] - 1; i >= 0; i--) {

      for (int j = matrixSize[0]; j >= 0; j--) {
        data[i][j] = data[i][j] / data[i][i];
      }

      for (int j = i - 1; j >= 0; j--) {

        double c = data[j][i] / data[i][i];
        for (int k = matrixSize[0]; k >= 0; k--) {
          data[j][k] = data[j][k] - data[i][k] * c;
        }
      }
    }

    return data;
  }

}
