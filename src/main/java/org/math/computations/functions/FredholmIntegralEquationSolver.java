package org.math.computations.functions;

import static org.math.computations.functions.Utils.generateMesh;

import java.util.ArrayList;
import java.util.List;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;
import org.math.computations.polynomials.Polynomial;

public class FredholmIntegralEquationSolver {

    private final double lowerBound;
    private final double upperBound;
    private final Function f;
    private final FunctionalR2 ker;
    private final double lambda;

    private final FunctionalR2 delta = (x, y) -> Math.abs(Math.abs(x) - Math.abs(y)) < 10e-6 ? 1 : 0;

    public FredholmIntegralEquationSolver(FunctionalR2 ker, Function f, double lambda,
        double lowerBound, double upperBound) {

        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.f = f;
        this.ker = ker;
        this.lambda = lambda;
    }

    public Function solveSingularKernel(int expansionTermsAmount) {

        List<Function> An = new ArrayList<>();
        List<Function> Bn = new ArrayList<>();
        for (int i = 0; i < 2*expansionTermsAmount; i+=2) {

            int k = i;
            An.add(x -> Math.pow(x, k) / Utils.factorial(k));
            Bn.add(x -> Math.pow(x, k));
        }

        double[][] A = new double[expansionTermsAmount][expansionTermsAmount];
        double[][] b = new double[expansionTermsAmount][1];

        for (int i = 0; i < expansionTermsAmount; i++) {
            for (int j = 0; j < expansionTermsAmount; j++) {
                A[i][j] = delta.evaluate(i, j) - lambda * dotProductL2(Bn.get(i), An.get(j));
            }
            b[i][0] = dotProductL2(Bn.get(i), f);
        }

        DenseMatrix matrixA = new DenseMatrix(expansionTermsAmount, expansionTermsAmount, A);
        DenseMatrix vectorB = new DenseMatrix(expansionTermsAmount, 1, b);
        LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(matrixA, vectorB));

        Matrix solution = solver.solve("IG");

        Polynomial p = new Polynomial(List.of(1.0));
        Polynomial sum = new Polynomial(List.of(0.0));
        for (int i = 0; i < expansionTermsAmount; i++) {
            int k = 2*i;
            sum = sum.sum(p.scalarMultiply(solution.getElement(0, i) / Utils.factorial(k)));
            p = p.multiply(new Polynomial(List.of(0.0, 0.0, 1.0)));
        }

        Polynomial finalSum = sum;
        Function u = x -> f.evaluate(x) + lambda * finalSum.evaluate(x);

        return u;
    }

    public Function solveQuadrature(int nodesAmount) {

        double step = (upperBound - lowerBound) / (nodesAmount - 1);

        List<Double> mesh = generateMesh(lowerBound, upperBound, nodesAmount - 1);

        double[][] D = new double[mesh.size()][mesh.size()];
        double[][] b = new double[mesh.size()][1];

        for (int i = 0; i < mesh.size(); i++) {
            for (int j = 0; j < mesh.size(); j++) {
                D[i][j] = delta.evaluate(i, j) - lambda * ker.evaluate(mesh.get(j), mesh.get(i)) * step;
            }
            b[i][0] = f.evaluate(mesh.get(i));
        }

        DenseMatrix matrixD = new DenseMatrix(mesh.size(), mesh.size(), D);
        DenseMatrix vectorB = new DenseMatrix(mesh.size(), 1, b);

        LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(matrixD, vectorB));
        Matrix solution = solver.solveIterative("SIT", 20);
        // Matrix solution = solver.solve("IG");

        return x -> {
            double value = 0.0;
            for (int i = 0; i < solution.getSize()[0]; i++) {
                value += lambda * ker.evaluate(x, mesh.get(i)) * solution.getElement(0, i) * step;
            }
            return value + f.evaluate(x);
        };
    }

    private double dotProductL2(Function f, Function g) {

        double step = (upperBound - lowerBound) / 1000;
        double x = lowerBound + step / 2;

        double value = 0;
        while (x < upperBound) {
            value += f.evaluate(x) * g.evaluate(x);
            x += step;
        }

        return value * step;
    }

}
