package tasks.sem5;

import java.util.List;
import java.util.Locale;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.EigenvalueProblemSolver;
import org.math.computations.matrices.Matrix;
import org.math.computations.matrices.Utils;

public class Task34 {

    public static void main(String[] args) {

        DenseMatrix matrixA = new DenseMatrix("eigenvalue_test/matrix7.txt");
        List<Matrix> result = EigenvalueProblemSolver.solveFullEigenvalueProblem(matrixA, 1e-6);
        for (int i = 1; i < result.size(); i++) {
            System.out.printf(
                Locale.US, "Eigenvalue λ%d = %.6f\nEigenvector:\n", i, result.get(0).getElement(0, i - 1));
            System.out.println(result.get(i));
            Matrix error = matrixA.mul(result.get(i))
                .add(result.get(i).scalarMultiply(-1 * result.get(0).getElement(0, i - 1)));
            System.out.print("\nEigenvector error:\n");
            System.out.println(error);
        }

        List<Matrix> resultPowers = EigenvalueProblemSolver.getLargestEigenvalueAndEigenvector(
            matrixA, 1e-3);
        System.out.printf(Locale.US, "Eigenvalue λ1 = %.6f\n", resultPowers.get(0).getElement(0, 0));
        System.out.println(resultPowers.get(1));
        System.out.printf(Locale.US, "Error post estimation = %.6f", Utils.euclideanMean(
            matrixA.mul(resultPowers.get(1)).add(
                resultPowers.get(1).scalarMultiply(-1 * resultPowers.get(0).getElement(0, 0)))));
    }

}
