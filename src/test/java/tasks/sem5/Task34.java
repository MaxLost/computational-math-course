package tasks.sem5;

import java.util.List;
import java.util.Locale;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.EigenvalueProblemSolver;
import org.math.computations.matrices.Matrix;

public class Task34 {

    public static void main(String[] args) {

        DenseMatrix matrixA = new DenseMatrix("eigenvalue_test/matrix7.txt");
        List<Matrix> result = EigenvalueProblemSolver.solveFullEigenvalueProblem(matrixA, 1e-6);
        for (int i = 1; i < result.size(); i++) {
            System.out.printf(
                Locale.US, "Eigenvalue λ%d = %.6f\nEigenvector:\n", i, result.get(0).getElement(0, i - 1));
            for (int j = 0; j < matrixA.getSize()[0]; j++) {
                System.out.printf(Locale.US, "%.6f\n", result.get(i).getElement(0, j));
            }
            Matrix error = matrixA.mul(result.get(i))
                .add(result.get(i).scalarMultiply(-1 * result.get(0).getElement(0, i - 1)));
            System.out.print("\nEigenvector error:\n");
            System.out.println(error);
        }
    }

}
