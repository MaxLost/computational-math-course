package tasks.sem5;

import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;

// Вариант 7
public class Task2 {

  public static void main(String[] args) {

    DenseMatrix matrixA = new DenseMatrix("linear_system_solver/linear_system2_A.txt");
    DenseMatrix b = new DenseMatrix("linear_system_solver/linear_system2_b.txt");
    Matrix system = LinearSystemSolver.getAugmentedMatrix(matrixA, b);
    System.out.println("Matrix of a system:\n" + system);

    LinearSystemSolver solver = new LinearSystemSolver(system);
    Matrix expected = solver.solve("IG");
    System.out.println("Solution via Gauss method:\n" + expected);

    Matrix[] transformedSystem = LinearSystemSolver.transformSystem(system);
    Matrix matrixH = transformedSystem[0];
    Matrix matrixG = transformedSystem[1];

    double meanH = LinearSystemSolver.mean(matrixH);
    System.out.println("Mean of H = " + meanH);
    double meanG = LinearSystemSolver.mean(matrixG);
    double priorEstimation = Math.pow(meanH, 7) * meanG / (1 - meanH);
    System.out.println("Prior estimation for || x^(7) - x || <= " + priorEstimation);

    LinearSystemSolver transformedSolver = new LinearSystemSolver(
        LinearSystemSolver.getAugmentedMatrix(matrixH, matrixG));

    Matrix x = transformedSolver.solveIterative("SI", 7);
    System.out.println("\nSolution via simple iteration method:\n" + x + "\nDifference between solutions:\n"
            + expected.add(x.scalarMultiply(-1)));
    Matrix x1 = transformedSolver.solveIterative("SE", 7);
    System.out.println("\nSolution via Seidel's method:\n" + x1+ "\nDifference between solutions:\n"
        + expected.add(x1.scalarMultiply(-1)));
    Matrix x2 = transformedSolver.solveIterative("RE", 7);
    System.out.println("\nSolution via upper relaxation method:\n" + x2+ "\nDifference between solutions:\n"
        + expected.add(x2.scalarMultiply(-1)));
  }

}
