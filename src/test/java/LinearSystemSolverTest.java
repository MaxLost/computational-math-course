import static org.junit.Assert.assertTrue;

import org.junit.Assert;
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
          Math.abs(answer.getElement(0, i) - matrix.getElement(colCount - 1, i)) < 10e-12);
    }
  }

}
