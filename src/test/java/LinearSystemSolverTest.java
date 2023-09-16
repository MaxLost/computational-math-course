import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.math.computations.matrices.*;

import org.junit.Test;

public class LinearSystemSolverTest {

  @Test
  public void defaultGaussTest() {
    DenseMatrix matrix = new DenseMatrix("linear_system1.txt");
    LinearSystemSolver solver = new LinearSystemSolver(matrix);
    Matrix result = solver.solve("DG");

    int rowCount = matrix.getSize()[0];
    int colCount = matrix.getSize()[1];

    double[][] data = new double[rowCount][colCount - 1];
    for (int i = 0; i < rowCount; i++) {
      for (int j = 0; j < colCount - 1; j++) {
        data[i][j] = matrix.getElement(j, i);
      }
    }
    DenseMatrix matrix1 = new DenseMatrix(rowCount, colCount - 1, data);
    DenseMatrix answer = (DenseMatrix) matrix1.mul(result);

    for (int i = 0; i < rowCount; i++) {
      assertTrue(
          Math.abs(answer.getElement(0, i) - matrix.getElement(colCount - 1, i)) < 10e-14);
    }
  }

  @Test
  public void improvedGaussTest() {
    DenseMatrix matrix = new DenseMatrix("linear_system1.txt");
    LinearSystemSolver solver = new LinearSystemSolver(matrix);
    Matrix result = solver.solve("IG");

    int rowCount = matrix.getSize()[0];
    int colCount = matrix.getSize()[1];

    double[][] data = new double[rowCount][colCount - 1];
    for (int i = 0; i < rowCount; i++) {
      for (int j = 0; j < colCount - 1; j++) {
        data[i][j] = matrix.getElement(j, i);
      }
    }
    DenseMatrix matrix1 = new DenseMatrix(rowCount, colCount - 1, data);
    DenseMatrix answer = (DenseMatrix) matrix1.mul(result);

    for (int i = 0; i < rowCount; i++) {
      assertTrue(
          Math.abs(answer.getElement(0, i) - matrix.getElement(colCount - 1, i)) < 10e-14);
    }
  }

  @Test
  public void solveLuDecompositionTest() {
    DenseMatrix matrix = new DenseMatrix("linear_system1.txt");
    LinearSystemSolver solver = new LinearSystemSolver(matrix);
    Matrix result = solver.solve("");

    int rowCount = matrix.getSize()[0];
    int colCount = matrix.getSize()[1];

    double[][] data = new double[rowCount][colCount - 1];
    for (int i = 0; i < rowCount; i++) {
      for (int j = 0; j < colCount - 1; j++) {
        data[i][j] = matrix.getElement(j, i);
      }
    }
    DenseMatrix matrix1 = new DenseMatrix(rowCount, colCount - 1, data);
    DenseMatrix answer = (DenseMatrix) matrix1.mul(result);

    for (int i = 0; i < rowCount; i++) {
      assertTrue(
          Math.abs(answer.getElement(0, i) - matrix.getElement(colCount - 1, i)) < 10e-14);
    }
  }

  @Test
  public void getAugmentedMatrixTest() {

    double[][] dataA = new double[][]{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    DenseMatrix matrixA = new DenseMatrix(3, 3, dataA);

    double[][] dataB = new double[][]{{-1.0}, {-2.0}, {-3.0}};
    DenseMatrix b = new DenseMatrix(3, 1, dataB);

    double[][] expectedData = new double[][]{{1.0, 2.0, 3.0, -1.0}, {4.0, 5.0, 6.0, -2.0},
        {7.0, 8.0, 9.0, -3.0}};
    DenseMatrix expected = new DenseMatrix(3, 4, expectedData);

    Matrix result = LinearSystemSolver.getAugmentedMatrix(matrixA, b);

    assertEquals(expected, result);
  }

}
