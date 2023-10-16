package tasks.sem5;

import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;

public class Task1 {

  private static Matrix subtract(Matrix a, Matrix b) {

    int[] matrixSize = a.getSize();
    double[][] dataA = a.toArray();
    double[][] dataB = b.toArray();
    double[][] dataC = new double[matrixSize[0]][matrixSize[1]];

    for (int i = 0; i < matrixSize[0]; i++) {
      for (int j = 0; j < matrixSize[1]; j++) {
        dataC[i][j] = dataA[i][j] - dataB[i][j];
      }
    }

    return new DenseMatrix(matrixSize[0], matrixSize[1], dataC);
  }

  public static void main(String[] args) {

    DenseMatrix matrixA = new DenseMatrix("linear_system_solver/linear_system3_A.txt");
    DenseMatrix b = new DenseMatrix("linear_system_solver/linear_system3_b.txt");

    double[][] dataC = matrixA.toArray();
    dataC[0][0] *= 10e-8;
    DenseMatrix matrixC = new DenseMatrix(3, 3, dataC);

    Matrix firstSystem = LinearSystemSolver.getAugmentedMatrix(matrixA, b);
    LinearSystemSolver firstSolver = new LinearSystemSolver(firstSystem);
    //Matrix firstResultDG = firstSolver.solve("DG");
    Matrix firstResultIG = firstSolver.solve("IG");
    Matrix firstResultLU = firstSolver.solve("LU");

    //Matrix firstErrorDG = subtract(matrixA.mul(firstResultDG), b);
    Matrix firstErrorIG = subtract(matrixA.mul(firstResultIG), b);
    Matrix firstErrorLU = subtract(matrixA.mul(firstResultLU), b);

    Matrix secondSystem = LinearSystemSolver.getAugmentedMatrix(matrixC, b);
    LinearSystemSolver secondSolver = new LinearSystemSolver(secondSystem);
    //Matrix secondResultDG = secondSolver.solve("DG");
    Matrix secondResultIG = secondSolver.solve("IG");
    Matrix secondResultLU = secondSolver.solve("LU");

   // Matrix secondErrorDG = subtract(matrixC.mul(secondResultDG), b);
    Matrix secondErrorIG = subtract(matrixC.mul(secondResultIG), b);
    Matrix secondErrorLU = subtract(matrixC.mul(secondResultLU), b);

    //System.out.println("Расширенная матрица системы A@x = b:\n" + firstSystem);
    //System.out.println("Решение системы методом Гаусса:\n" + firstResultDG);
    //System.out.println("Вектор невязки:\n" + firstErrorDG + "\n");
    System.out.println("Решение системы методом Гаусса с выбором элемента:\n" + firstResultIG);
    System.out.println("Вектор невязки:\n" + firstErrorIG + "\n");
    System.out.println("Решение системы с помощью LU разложения:\n" + firstResultLU);
    System.out.println("Вектор невязки:\n" + firstErrorLU + "\n");

    //System.out.println("Расширенная матрица системы C@x = b:\n" + secondSystem);
    //System.out.println("Решение системы методом Гаусса:\n" + secondResultDG);
    //System.out.println("Вектор невязки:\n" + secondErrorDG + "\n");
    System.out.println("Решение системы методом Гаусса с выбором элемента:\n" + secondResultIG);
    System.out.println("Вектор невязки:\n" + secondErrorIG + "\n");
    System.out.println("Решение системы с помощью LU разложения:\n" + secondResultLU);
    System.out.println("Вектор невязки:\n" + secondErrorLU + "\n");

    Matrix inverseA = matrixA.invert();
    Matrix inverseC = matrixC.invert();
    DenseMatrix matrixE = new DenseMatrix(3, 3, new double[][]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});

    Matrix errorA = subtract(matrixA.mul(inverseA), matrixE);
    Matrix errorC = subtract(matrixC.mul(inverseC), matrixE);

    System.out.println("Матрица, обратная к матрице A:\n" + inverseA);
    System.out.println("Матрица невязки:\n" + errorA + "\n");
    System.out.println("Матрица, обратная к матрице C:\n" + inverseC);
    System.out.println("Матрица невязки:\n" + errorC + "\n");
  }

}
