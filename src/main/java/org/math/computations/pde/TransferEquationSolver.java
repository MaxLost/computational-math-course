package org.math.computations.pde;

import static org.math.computations.functions.Utils.generateMesh;

import java.util.Arrays;
import java.util.List;
import org.math.computations.functions.FunctionalR2;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;

public class TransferEquationSolver {

    private final FunctionalR2 f;
    private final double C;
    private final double lowerBoundX;
    private final double upperBoundX;
    private final double lowerBoundT;
    private final double upperBoundT;


    public TransferEquationSolver(FunctionalR2 f, double c, double lowerBoundX, double upperBoundX,
        double lowerBoundT, double upperBoundT) {

        this.f = f;
        this.C = c;
        this.lowerBoundX = lowerBoundX;
        this.upperBoundX = upperBoundX;
        this.lowerBoundT = lowerBoundT;
        this.upperBoundT = upperBoundT;
    }

    public double[][] solveExplicitFull(double[] u, double stepX, double stepT) {

        int sizeX = (int) Math.ceil((upperBoundX + lowerBoundX) / stepX);
        int sizeT = (int) Math.ceil((upperBoundT + lowerBoundT) / stepT);

        double[][] grid = generateGrid(sizeX, sizeT);

        List<Double> x = generateMesh(lowerBoundX, upperBoundX, sizeX);
        List<Double> t = generateMesh(lowerBoundT, upperBoundT, sizeT);

        System.arraycopy(u, 0, grid[0], 0, u.length);

        for (int i = 1; i < sizeT; i++) {
            double[] uNew = u.clone();
            for (int j = 1; j < sizeX; j++) {
                uNew[j] = u[j] + stepT * (-1 * C * (u[j] - u[j-1]) / stepX + f.evaluate(x.get(j), t.get(i-1)));
            }
            u = uNew;
            System.arraycopy(u, 0, grid[i], 0, u.length);
        }

        return grid;
    }

    public double[][] solveFullyImplicit(double[] u, double stepX, double stepT) {

        int sizeX = (int) Math.ceil((upperBoundX + lowerBoundX) / stepX);
        int sizeT = (int) Math.ceil((upperBoundT + lowerBoundT) / stepT);

        double[][] grid = generateGrid(sizeX, sizeT);
        System.arraycopy(u, 0, grid[0], 0, u.length);

        List<Double> x = generateMesh(lowerBoundX, upperBoundX, sizeX);
        List<Double> t = generateMesh(lowerBoundT, upperBoundT, sizeT);

        double[] uNew = Arrays.copyOf(u, u.length);

        double[][] A = buildSchemeMatrix(stepX, stepT, sizeX);

        for (int i = 1; i < sizeT; i++) {
            double[] b = new double[sizeX - 1];

            for (int j = 1; j < sizeX; j++) {
                b[j - 1] = uNew[j] + stepT * f.evaluate(x.get(j), t.get(i - 1));
            }

            Matrix solution = solveLinearSystem(sizeX, A, b);

            for (int j = 1; j < sizeX; j++) {
                uNew[j] = solution.getElement(0, j - 1);
            }

            System.arraycopy(uNew, 0, grid[i], 0, uNew.length);
        }

        return grid;
    }

    public double[][] solveImplicit(double[] u, double stepX, double stepT) {

        int sizeX = (int) Math.ceil((upperBoundX + lowerBoundX) / stepX);
        int sizeT = (int) Math.ceil((upperBoundT + lowerBoundT) / stepT);

        double[][] grid = generateGrid(sizeX, sizeT);
        System.arraycopy(u, 0, grid[0], 0, u.length);

        List<Double> x = generateMesh(lowerBoundX, upperBoundX, sizeX);
        List<Double> t = generateMesh(lowerBoundT, upperBoundT, sizeT);

        double[] uNew = Arrays.copyOf(u, u.length);

        grid[0] = Arrays.copyOf(uNew, uNew.length);

        double[][] A = buildSchemeMatrix(stepX, stepT, sizeX);

        for (int i = 1; i < sizeT; i++) {
            double[] b = new double[sizeX - 1];
            for (int j = 1; j < sizeX; j++) {
                b[j-1] = uNew[j-1] + stepT * f.evaluate(x.get(j), t.get(i - 1));
            }

            Matrix solution = solveLinearSystem(sizeX, A, b);

            for (int j = 1; j < sizeX; j++) {
                uNew[j] = solution.getElement(0, j - 1);
            }

            System.arraycopy(uNew, 0, grid[i], 0, uNew.length);
        }

        return grid;
    }

    public double[][] solveSymmetric(double[] u, double stepX, double stepT) {

        int sizeX = (int) Math.ceil((upperBoundX + lowerBoundX) / stepX);
        int sizeT = (int) Math.ceil((upperBoundT + lowerBoundT) / stepT);

        double[][] grid = generateGrid(sizeX, sizeT);
        System.arraycopy(u, 0, grid[0], 0, u.length);

        List<Double> x = generateMesh(lowerBoundX, upperBoundX, sizeX);
        List<Double> t = generateMesh(lowerBoundT, upperBoundT, sizeT);

        double[] uNew = Arrays.copyOf(u, u.length);
        grid[0] = Arrays.copyOf(uNew, uNew.length);

        double[][] A = new double[sizeX - 1][sizeX - 1];

        for (int i = 1; i < sizeX - 2; i++) {

            A[i][i] = 2 + 2 * C * stepT / stepX;
            A[i][i - 1] = -C * stepT / stepX;
            A[i][i + 1] = -C * stepT / stepX;
        }

        A[0][0] = A[sizeX - 2][sizeX - 2] = 2 + 2 * C * stepT / stepX;
        A[0][1] = A[sizeX - 2][sizeX - 3] = -C * stepT / stepX;

        for (int i = 1; i < sizeT; i++) {

            double[] b = new double[sizeX - 1];
            for (int j = 1; j < sizeX; j++) {

                b[j - 1] = 2 * uNew[j - 1] + 2 * stepT * f.evaluate(x.get(j) + stepX / 2,
                    t.get(i) + stepT / 2);
            }

            Matrix solution = solveLinearSystem(sizeX, A, b);

            for (int j = 1; j < sizeX; j++) {
                uNew[j] = solution.getElement(0, j - 1);
            }

            grid[i] = Arrays.copyOf(uNew, uNew.length);
        }

        return grid;
    }

    private double[][] buildSchemeMatrix(double stepX, double stepT, int sizeX) {
        double[][] A = new double[sizeX - 1][sizeX - 1];
        for (int i = 0; i < sizeX -1; i++) {
            A[i][i] = 1 + C * stepT / stepX;
            if (i > 0) {
                A[i][i-1] = -C * stepT / stepX;
            }
        }
        return A;
    }

    private static Matrix solveLinearSystem(int sizeX, double[][] A, double[] b) {

        double[][] newB = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            newB[i][0] = b[i];
        }

        DenseMatrix matrixA = new DenseMatrix(sizeX - 1, sizeX - 1, A);
        DenseMatrix vectorB = new DenseMatrix(sizeX - 1, 1, newB);
        LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(matrixA, vectorB));

        return solver.solve("IG");
    }


    private double[][] generateGrid(double lowerBoundX, double upperBoundX, double lowerBoundT,
        double upperBoundT, double stepX, double stepT) {

        int sizeX = (int) Math.ceil((upperBoundX + lowerBoundX) / stepX);
        int sizeT = (int) Math.ceil((upperBoundT + lowerBoundT) / stepT);

        return new double[sizeX][sizeT];
    }

    private double[][] generateGrid(int sizeX, int sizeT) {

        return new double[sizeT][sizeX];
    }
}
