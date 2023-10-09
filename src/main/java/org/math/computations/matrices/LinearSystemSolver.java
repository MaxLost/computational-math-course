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

  public Matrix solveIterative(String method, int approximation) {

    return switch (method) {
      case "SE" -> solveSeidelMethod(approximation);
      case "SET" -> solveSeidelMethodWithTransformation(approximation);
      case "RE" -> solveRelaxation(approximation);
      case "RET" -> solveRelaxationWithTransformation(approximation);
      case "SI" -> solveSimpleIteration(approximation);
      default -> solveSimpleIterationWithTransformation(approximation);
    };
  }

  private Matrix solveSimpleIteration(Matrix system, int approximation) {

    double[][][] data = deAugmentMatrix(system);

    Matrix matrixH = new DenseMatrix(data[0].length, data[0][0].length, data[0]);
    Matrix matrixG = new DenseMatrix(data[1].length, data[1][0].length, data[1]);

    double meanH = mean(matrixH);
    if (meanH > 1) {
      System.out.println("Simple iteration algorithm might not converge for this system.\n");
    }

    Matrix x = new DenseMatrix(data[1].length, 1, new double[data[1].length][1]);
    Matrix prevX = x;
    for (int i = 0; i < approximation; i++) {
      prevX = x;
      x = matrixH.mul(x).add(matrixG);
    }

    double postEstimation = meanH * mean(x.add(prevX.scalarMultiply(-1))) / (1 - meanH);
    System.out.println("Post estimation for || x^(7) - x || <= " + postEstimation + "\n");

    return x;
  }

  private Matrix solveSimpleIteration(int approximation) {

    return solveSimpleIteration(this.augmentedMatrix, approximation);
  }

  private Matrix solveSimpleIterationWithTransformation(int approximation) {

    Matrix[] transformedSystem = transformSystem(this.augmentedMatrix);
    Matrix system = getAugmentedMatrix(transformedSystem[0], transformedSystem[1]);

    return solveSimpleIteration(system, approximation);
  }

  private static Matrix solveSeidelMethod(Matrix system, int approximation) {

    int[] matrixSize = system.getSize();
    matrixSize[1]--;

    double[][] dataL = new double[matrixSize[0]][matrixSize[1]];
    double[][] dataR = new double[matrixSize[0]][matrixSize[1]];
    double[][] dataE = new double[matrixSize[0]][matrixSize[1]];

    for (int i = 0; i < matrixSize[0]; i++) {
      for (int j = 0; j < i; j++) {
        dataL[i][j] = system.getElement(j, i);
      }
      for (int j = i; j < matrixSize[1]; j++) {
        dataR[i][j] = system.getElement(j, i);
      }
      dataE[i][i] = 1;
    }

    Matrix matrixB = new DenseMatrix(matrixSize[0], 1, extractLastColumn(system.toArray()));
    Matrix matrixL = new DenseMatrix(matrixSize[0], matrixSize[1], dataL);
    Matrix matrixR = new DenseMatrix(matrixSize[0], matrixSize[1], dataR);
    Matrix matrixE = new DenseMatrix(matrixSize[0], matrixSize[1], dataE);
    Matrix matrixT = ((DenseMatrix) matrixE.add(matrixL.scalarMultiply(-1))).invert();
    Matrix matrixH = matrixT.mul(matrixR);
    Matrix matrixG = matrixT.mul(matrixB);

    double meanH = mean(matrixH);
    if (meanH > 1) {
      System.out.println("Seidel's algorithm might not converge for this system.\n");
    }

    Matrix x = new DenseMatrix(matrixSize[0], 1, new double[matrixSize[0]][1]);
    for (int i = 0; i < approximation; i++) {
      x = matrixH.mul(x).add(matrixG);
    }

    return x;
  }

  private Matrix solveSeidelMethod(int approximation) {

    return solveSeidelMethod(this.augmentedMatrix, approximation);
  }

  private Matrix solveSeidelMethodWithTransformation(int approximation) {

    Matrix[] transformedSystem = transformSystem(this.augmentedMatrix);
    Matrix system = getAugmentedMatrix(transformedSystem[0], transformedSystem[1]);

    return solveSeidelMethod(system, approximation);
  }

  private Matrix solveRelaxation(Matrix system, int approximation) {

    double[][][] data = deAugmentMatrix(system);
    int[] matrixSize = new int[] {data[0].length, data[0][0].length};

    double[][] dataH = data[0];
    double[][] dataG = data[1];
    double q = 1;

    double[][] x = new double[matrixSize[0]][1];
    for (int k = 0; k < approximation; k++) {
      for (int i = 0; i < matrixSize[0]; i++) {

        double sumCurrent = 0;
        for (int j = 0; j < i; j++) {
          sumCurrent += dataH[i][j] * x[j][0];
        }

        double sumPrevious = 0;
        for (int j = i + 1; j < matrixSize[1]; j++) {
          sumPrevious += dataH[i][j] * x[j][0];
        }

        x[i][0] = x[i][0] + q * (sumCurrent + sumPrevious - x[i][0] + dataG[i][0]);
      }
    }

    return new DenseMatrix(matrixSize[0], 1, x);
  }

  private Matrix solveRelaxation(int approximation) {
    return solveRelaxation(this.augmentedMatrix, approximation);
  }

  private Matrix solveRelaxationWithTransformation(int approximation) {

    Matrix[] transformedSystem = transformSystem(this.augmentedMatrix);
    Matrix system = getAugmentedMatrix(transformedSystem[0], transformedSystem[1]);

    return solveRelaxation(system, approximation);
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

  private static double[][][] deAugmentMatrix(Matrix m) {

    int[] augmentedSize = m.getSize();
    double[][] matrixA = new double[augmentedSize[0]][augmentedSize[1] - 1];
    double[][] b = new double[augmentedSize[0]][1];

    for (int i = 0; i < augmentedSize[0]; i++) {

      for (int j = 0; j < augmentedSize[1] - 1; j++) {
        matrixA[i][j] = m.getElement(j, i);
      }

      b[i][0] = m.getElement(augmentedSize[1] - 1, i);
    }

    return new double[][][]{matrixA, b};
  }

  private double[][][] deAugmentMatrix() {

    return deAugmentMatrix(this.augmentedMatrix);
  }

  public static Matrix[] transformSystem(Matrix system) {

    double[][][] data = deAugmentMatrix(system);
    Matrix matrixA = new DenseMatrix(data[0].length, data[0][0].length, data[0]);
    Matrix matrixB = new DenseMatrix(data[1].length, data[1][0].length, data[1]);

    int[] matrixSize = matrixA.getSize();

    double[][] dataE = new double[matrixSize[0]][matrixSize[1]];
    double[][] dataD = new double[matrixSize[0]][matrixSize[1]];
    for (int i = 0; i < matrixSize[0]; i++) {
      dataE[i][i] = 1;
      dataD[i][i] = matrixA.getElement(i, i);
    }

    Matrix inverseD = new DenseMatrix(dataD.length, dataD[0].length, dataD).invert();
    Matrix matrixH = new DenseMatrix(dataE.length, dataE[0].length, dataE).add(
        inverseD.mul(matrixA).scalarMultiply(-1));

    Matrix matrixG = inverseD.mul(matrixB);

    return new Matrix[] {matrixH, matrixG};
  }

  public static double mean(Matrix m) {

    double mean = 0;
    int[] matrixSize = m.getSize();
    for (int i = 0; i < matrixSize[0]; i++) {

      double sum = 0;
      for (int j = 0; j < matrixSize[1]; j++) {
        sum += Math.abs(m.getElement(j, i));
      }

      mean = Math.max(mean, sum);
    }

    return mean;
  }

}
