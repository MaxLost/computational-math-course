package org.math.computations.ode;

import static org.math.computations.ode.Utils.generateMesh;
import static org.math.computations.ode.Utils.normalizeEuclidean;
import static org.math.computations.ode.Utils.vectorDotProduct;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.math.computations.functions.Function;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;
import org.math.computations.matrices.Utils;
import org.math.computations.polynomials.JacobiPolynomialGenerator;
import org.math.computations.polynomials.Polynomial;

public class SturmLiouvilleEigenvaluesFinder {

    /*
        L[u] = -(p * u')' + q * u

        Find a: L[u] = a * u
     */

    private final Function p;
    private final Function q;
    private final double[] boundaryConstraints; // [a_1, a_2, a, b_1, b_2, b]

    private final double lowerBound;
    private final double upperBound;

    public SturmLiouvilleEigenvaluesFinder(Function p, Function q, double lowerBound, double upperBound, double[] boundaryConstraints) {

        this.p = p;
        this.q = q;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.boundaryConstraints = boundaryConstraints;
    }

    /*
        Boundaries:
        a_1 * y(a) - a_2 * dy(a)/dx = a
        b_1 * y(b) + b_2 * dy(b)/dx = b
     */

    public Map<Double, Function> find(int coordinateFunctionsAmount, int eigenvaluesAmount) {

        // Calculate Г and Г_L matrices and find eigenvalues of Г^-1*Г_L
        List<Polynomial> coordinateFunctions = JacobiPolynomialGenerator.getJacobiPolynomials(2, 2, coordinateFunctionsAmount - 1);
        Polynomial factor = new Polynomial(List.of(1.0, 0.0, -1.0));
        for (int i = 0; i < coordinateFunctions.size(); i++) {

            double coeff = Math.sqrt((double) ((i + 3) * (i + 4) * (2 * i + 5)) / (2*(i+1)*(i+2))) / 4;
            coordinateFunctions.set(i, coordinateFunctions.get(i).multiply(factor).scalarMultiply(coeff));
        }

        double[][] gamma = new double[coordinateFunctionsAmount][coordinateFunctionsAmount];
        double[][] gammaL = new double[coordinateFunctionsAmount][coordinateFunctionsAmount];

        for (int i = 0; i < coordinateFunctionsAmount; i++) {
            for (int j = 0; j < coordinateFunctionsAmount; j++) {
                gamma[i][j] = dotProductL2(coordinateFunctions.get(i), coordinateFunctions.get(j));
                gammaL[i][j] = dotProductFuncSpace(coordinateFunctions.get(i), coordinateFunctions.get(j));
            }
        }

        List<Function> eigenfunctionApprox = new ArrayList<>();
        for (int i = 1; i <= eigenvaluesAmount; i++) {
            eigenfunctionApprox.add(getEigenfunction(i));
        }

        DenseMatrix G = new DenseMatrix(coordinateFunctionsAmount, coordinateFunctionsAmount, gamma);
        DenseMatrix Gl = new DenseMatrix(coordinateFunctionsAmount, coordinateFunctionsAmount, gammaL);
        Matrix inverseG = G.invert();
        Matrix A = inverseG.mul(Gl);

        List<Function> eigenfunctions = new ArrayList<>();

        for (int i = 0; i < eigenfunctionApprox.size(); i++) {
            Matrix decomposition = getEigenfunctionCoefficients(A,
                eigenfunctionApprox.get(i),
                coordinateFunctions);
            // Matrix decomposition = generateRandomVector(coordinateFunctionsAmount);

            Polynomial u = new Polynomial(List.of(0.0));
            for (int j = 0; j < coordinateFunctionsAmount; j++) {
                u = u.sum(coordinateFunctions.get(j).scalarMultiply(decomposition.getElement(0, j)));
            }

            eigenfunctions.add(u);
        }

        HashMap<Double, Function> result = new HashMap<>();
        for (int i = 0; i < eigenfunctions.size(); i++) {
            Function u = eigenfunctions.get(i);
            Double eigenvalue = dotProductFuncSpace(u, u) / dotProductL2(u, u);
            result.put(eigenvalue, u);
        }

        return result;
    }

    private Matrix getEigenfunctionCoefficients(Matrix matrixA, Function u, List<? extends Function> base) {

        Matrix vectorY = decomposeFunction(u, base);
        double eigenvalue = dotProductFuncSpace(u, u) / dotProductL2(u, u);
        Matrix matrixB = matrixA.add(
            Utils.getIdentityMatrix(matrixA.getSize()[0]).scalarMultiply(-1 * eigenvalue));

        for (int i = 0; i < 30; i++) {

            LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(matrixB, vectorY));
            Matrix newVectorY = solver.solve("IG");

            eigenvalue += vectorDotProduct(newVectorY, vectorY) / vectorDotProduct(newVectorY, newVectorY);
            matrixB = matrixA.add(
                Utils.getIdentityMatrix(matrixA.getSize()[0]).scalarMultiply(-1 * eigenvalue));
        }

        return normalizeEuclidean(vectorY);
    }

    private Matrix decomposeFunction(Function f, List<? extends Function> base) {

        List<Double> mesh = new ArrayList<>();
        for (int i = 0; i < base.size(); i++) {
            mesh.add(Math.cos((2 * (i + 1) - 1) * Math.PI / (2 * base.size())));
        }

        double[][] matrix = new double[mesh.size()][mesh.size()];
        double[][] fValues = new double[mesh.size()][1];
        for (int i = 0; i < mesh.size(); i++) {
            for (int j = 0; j < mesh.size(); j++) {
                matrix[i][j] = base.get(j).evaluate(mesh.get(i));
            }
            fValues[i][0] = f.evaluate(mesh.get(i));
        }

        Matrix A = new DenseMatrix(mesh.size(), mesh.size(), matrix);
        Matrix b = new DenseMatrix(mesh.size(), 1, fValues);

        LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(A, b));

        return solver.solve("IG");
    }

    private Function getEigenfunction(int k) {

        Function f = x -> Math.cos((Math.PI * k - Math.PI / 2) * x);
        double l2Mean = Math.sqrt(dotProductL2(f, f));

        return f;
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

    private double dotProductFuncSpace(Function f, Function g) {

        double step = (upperBound - lowerBound) / 1000;
        double x = lowerBound + step / 2;

        double value = 0;
        while (x < upperBound) {

            double df = (f.evaluate(x + step) - f.evaluate(x - step)) / (2*step);
            double dg = (g.evaluate(x + step) - g.evaluate(x - step)) / (2*step);
            value += (p.evaluate(x) * df * dg + q.evaluate(x) * f.evaluate(x) * g.evaluate(x));

            x += step;
        }

        return value * step;
    }
}
