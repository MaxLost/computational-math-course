package org.math.computations.matrices;

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
   * Finds unique solution of linear system: A @ X = B.
   *
   * @param method method of solving system
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

    if (isSmallLeadingElement) {
      throw new RuntimeException("Result may be incorrect due to small leading elements");
    }

    double[][] result = new double[matrixSize[0]][1];
    for (int i = 0; i < matrixSize[0]; i++) {
      result[i][0] = data[i][matrixSize[1] - 1];
    }

    return new DenseMatrix(matrixSize[0], 1, result);
  }

  private Matrix solveImprovedGauss() {
    return null;
  }

  private Matrix solveLuDecomposition() {
    return null;
  }

}
