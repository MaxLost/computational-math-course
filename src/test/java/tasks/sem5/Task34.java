package tasks.sem5;

import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map.Entry;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.EigenvalueProblemSolver;
import org.math.computations.matrices.Matrix;

public class Task34 {

    public static void main(String[] args) {


        DenseMatrix matrixA = new DenseMatrix("eigenvalue_test/matrix7.txt");

        double[] eigenvalues = new double[]{-1.6901924778225, 0.82209423831600, -0.48300176049347};

        HashMap<Double, Matrix> eigenvectors = new HashMap<>();
        eigenvectors.put(-1.6901924778225,
            new DenseMatrix(3, 1, new double[][]{{-0.691729}, {-0.146656}, {0.707108}}));
        eigenvectors.put(0.82209423831600,
            new DenseMatrix(3, 1, new double[][]{{-0.207403}, {0.978256}, {1.45553e-6}}));
        eigenvectors.put(-0.48300176049347,
            new DenseMatrix(3, 1, new double[][]{{0.691733}, {0.146655}, {0.707105}}));

        System.out.println("\n--- WolframAlpha results ---\n");
        for (Entry<Double, Matrix> x : eigenvectors.entrySet()) {
            System.out.printf(
                Locale.US, "Eigenvalue λ = %.12f\nEigenvector:\n", x.getKey());
            System.out.println(x.getValue());
        }

        System.out.println("\n--- Jacobi's method ---\n");

        List<Matrix> result = EigenvalueProblemSolver.solveFullEigenvalueProblem(matrixA, 1e-6);
        for (int i = 1; i < result.size(); i++) {
            System.out.printf(
                Locale.US, "Eigenvalue λ%d = %.12f\nEigenvector error: %.5e\nEigenvector:\n",
                i,
                result.get(0).getElement(0, i - 1),
                Math.abs(eigenvalues[result.size() - i - 1]) - Math.abs(
                    result.get(0).getElement(0, i - 1)));
            System.out.println(result.get(i));
        }

        System.out.println("\n--- Method of power iterations ---\n");

        List<Matrix> resultPowers = EigenvalueProblemSolver.getLargestEigenvalueAndEigenvector(
            matrixA, 1e-3);
        System.out.printf(Locale.US, "Eigenvalue λ1 = %.12f\nEigenvalue error = %.5e\nEigenvector:\n",
            resultPowers.get(0).getElement(0, 0),
            Math.abs(eigenvalues[0]) - Math.abs(resultPowers.get(0).getElement(0, 0)));
        System.out.println(resultPowers.get(1));

        System.out.println("\n--- Method of scalar multiplications ---\n");

        double resultScalar = EigenvalueProblemSolver.getLargestEigenvalueScalarMultiplication(matrixA, 1e-3);
        System.out.printf(Locale.US, "Eigenvalue λ1 = %.12f\nEigenvalue error = %.5e\n",
            resultScalar,
            Math.abs(eigenvalues[0]) - Math.abs(resultScalar));

        System.out.println("\n--- Finding opposite bound of spectre ---\n");

        List<Matrix> resultOpposite = EigenvalueProblemSolver.getOppositeSpectreBound(matrixA,
            resultScalar, 1e-3);
        System.out.printf(Locale.US, "Eigenvalue λ2 = %.12f\nEigenvalue error = %.5e\nEigenvector:\n",
            resultOpposite.get(0).getElement(0, 0),
            Math.abs(eigenvalues[1]) - Math.abs(resultOpposite.get(0).getElement(0, 0)));
        System.out.println(resultOpposite.get(1));

        System.out.println("\n--- Wielandt's method ---\n");

        List<Matrix> resultWielandt = EigenvalueProblemSolver.getLargestEigenvalueAndEigenvectorWielandt(
            matrixA, 1e-3);
        System.out.printf(Locale.US, "Eigenvalue λ1 = %.12f\nEigenvalue error = %.5e\nEigenvector:\n",
            resultWielandt.get(0).getElement(0, 0),
            Math.abs(eigenvalues[2]) - Math.abs(resultWielandt.get(0).getElement(0, 0)));
        System.out.println(resultWielandt.get(1));

        System.out.println("\n--- Aitken's refinement method ---\n");

        double resultAitken = EigenvalueProblemSolver.getLargestEigenvalueWithAitkenRefinement(matrixA, 1e-3);
        System.out.printf(Locale.US, "Eigenvalue λ1 = %.12f\nEigenvalue error = %.5e\n",
            resultAitken,
            Math.abs(eigenvalues[0]) - Math.abs(resultAitken));
    }

}
