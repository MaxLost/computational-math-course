package org.math.computations.ode;

import static org.math.computations.ode.Utils.generateMesh;

import java.util.ArrayList;
import java.util.List;
import org.math.computations.PlanePoint;
import org.math.computations.polynomials.Polynomial;
import org.math.computations.functions.Function;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;
import org.math.computations.polynomials.JacobiPolynomialGenerator;

public class SecondOrderBVPSolver {

    private final Function p;
    private final Function q;
    private final Function r;
    private final Function f;
    private final double[] boundaryConstraints; // [a_1, a_2, a, b_1, b_2, b]

    public SecondOrderBVPSolver(Function p, Function q, Function r, Function f, double[] boundaryConstraints) {

        this.p = p;
        this.q = q;
        this.r = r;
        this.f = f;
        this.boundaryConstraints = boundaryConstraints;
    }

    /*
        Boundaries:
        a_1 * y(a) - a_2 * dy(a)/dx = a
        b_1 * y(b) + b_2 * dy(b)/dx = b

        7th variant

        a = -1, b = 1

        -1 * dy(a)/dx = 0
        2 * y(b) + dy(b)/dx = 0

        a_1 = 0, a_2 = 1
        b_1 = 2, b_2 = 1

        - O(h) approximation:
            -1 * (y_1 - y_0) / h = 0

     */

    public List<PlanePoint> solveBvpDiff(double start, double end, int cellAmount,
        int approximationOrder, boolean useRichardsonExtrapolation) {

        List<PlanePoint> mesh = generateMesh(start, end, cellAmount);
        if (approximationOrder == 2) {
            if (useRichardsonExtrapolation) {
                mesh = solveBvpRichardsonSecondApproximation(mesh);
            }else {
                mesh = solveBvpSecondApproximation(mesh);
            }
        } else {
            if (useRichardsonExtrapolation) {
                mesh = solveBvpRichardsonFirstApproximation(mesh);
            } else {
                mesh = solveBvpFirstApproximation(mesh);
            }
        }

        return mesh;
    }

    // Projective methods
    public Function solveCollocation(int coordinateFunctionsAmount) {

        int nodesAmount = coordinateFunctionsAmount;

        Polynomial w1 = new Polynomial(List.of((double) -13/3, -2.0, 1.0)); // x^2 - 2x - 13/2
        Polynomial w2 = new Polynomial(List.of(-2.0, -3.0, 0.0, 1.0)); // x^3 - 3x - 2
        List<Polynomial> jacobiPolynomials = JacobiPolynomialGenerator.getJacobiPolynomials(2, 2, coordinateFunctionsAmount);
        List<Polynomial> coordinateSystem = new ArrayList<>(List.of(w1, w2));

        Polynomial factor = new Polynomial(List.of(1.0, 0.0, -1.0));
        factor = factor.multiply(factor);

        for (int i = 2; i < coordinateFunctionsAmount; i++) {
            coordinateSystem.add(jacobiPolynomials.get(i - 2).multiply(factor));
        }

        double[][] A = new double[coordinateFunctionsAmount][coordinateFunctionsAmount];
        double[][] b = new double[coordinateFunctionsAmount][1];

        List<Double> nodes = new ArrayList<>();
        for (int i = 0; i < nodesAmount; i++) {
            nodes.add(Math.cos((2 * (i + 1) - 1) * Math.PI / (2*nodesAmount)));
        }

        for (int i = 0; i < nodesAmount; i++) {
            for (int j = 0; j < coordinateFunctionsAmount; j++) {
                A[i][j] = LWValue(coordinateSystem.get(j), nodes.get(i));
            }
            b[i][0] = f.evaluate(nodes.get(i));
        }

        Matrix matrixA = new DenseMatrix(nodesAmount, coordinateFunctionsAmount, A);
        Matrix vectorB = new DenseMatrix(nodesAmount, 1, b);

        LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(matrixA, vectorB));
        Matrix c = solver.solve("");

        Polynomial solution = new Polynomial(List.of(0.0));

        for (int i = 0; i < c.getSize()[0]; i++) {
            solution = solution.sum(coordinateSystem.get(i).scalarMultiply(c.getElement(0, i)));
        }

        return solution;
    }

    private double LWValue(Function w, double x) {

        double step = 0.01;

        double du = (w.evaluate(x + step) - w.evaluate(x - step)) / (2*step);
        double ddu = (w.evaluate(x - step) - 2*w.evaluate(x) + w.evaluate(x + step)) / (step * step);

        return -1 * p.evaluate(x)*ddu + q.evaluate(x)*du + r.evaluate(x)*w.evaluate(x);
    }

    public Function solveMinSquareMethod(int coordinateFunctionsAmount) {

        Polynomial w1 = new Polynomial(List.of((double) -13/3, -2.0, 1.0)); // x^2 - 2x - 13/2
        Polynomial w2 = new Polynomial(List.of(-2.0, -3.0, 0.0, 1.0)); // x^3 - 3x - 2
//        Polynomial w1 = new Polynomial(List.of(1.0, 2.0, 1.0)); // x^2 + 2*x + 1
//        Polynomial w2 = new Polynomial(List.of(-2.0, -3.0, 0.0, 1.0)); // x^3 - 3*x - 2

        List<Polynomial> jacobiPolynomials = JacobiPolynomialGenerator.getJacobiPolynomials(2, 2, coordinateFunctionsAmount);
        List<Polynomial> coordinateSystem = new ArrayList<>(List.of(w1, w2));

        Polynomial factor = new Polynomial(List.of(1.0, 0.0, -1.0));
        factor = factor.multiply(factor);

        for (int i = 2; i < coordinateFunctionsAmount; i++) {
            coordinateSystem.add(jacobiPolynomials.get(i - 2).multiply(factor));
        }

        int N = coordinateFunctionsAmount;

        double[][] A = new double[N][N];
        double[][] b = new double[N][1];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {

                A[i][j] = dotProductCoordinate(coordinateSystem.get(i), coordinateSystem.get(j));
            }
            b[i][0] = dotProductCoordinateAndF(f, coordinateSystem.get(i));
        }

        Matrix matrixA = new DenseMatrix(N, N, A);
        Matrix vectorB = new DenseMatrix(N, 1, b);

        LinearSystemSolver solver = new LinearSystemSolver(LinearSystemSolver.getAugmentedMatrix(matrixA, vectorB));
        Matrix c = solver.solve("IG");

        Polynomial solution = new Polynomial(List.of(0.0));

        for (int i = 0; i < c.getSize()[0]; i++) {
            solution = solution.sum(coordinateSystem.get(i).scalarMultiply(c.getElement(0, i)));
        }

        return solution;
    }

    private double dotProductCoordinate(Function a, Function b) {

        double step = 0.05;
        List<Double> La = LWValues(a, -1 + step / 2, 1 - step / 2, step);
        List<Double> Lb = LWValues(b, -1 + step / 2, 1 - step / 2, step);

        double value = 0;
        for (int i = 0; i < La.size(); i++) {
            value += La.get(i) * Lb.get(i);
        }

        return value * step;
    }

    private double dotProductCoordinateAndF(Function f, Function w) {

        double step = 0.05;
        List<Double> Lw = LWValues(w, -1 + step / 2, 1 - step / 2, step);

        double value = 0;
        double x = -1 + step;
        for (int i = 0; i < Lw.size(); i++) {
            value += Lw.get(i) * f.evaluate(x + step * i);
        }

        return value * step;
    }

    private List<Double> LWValues(Function w, double lowerBound, double upperBound, double step) {

        List<Double> result = new ArrayList<>();
        double x = lowerBound;

        while(x < upperBound) {

            double du = (w.evaluate(x + step) - w.evaluate(x - step)) / (2*step);
            double ddu = (w.evaluate(x - step) - 2*w.evaluate(x) + w.evaluate(x + step)) / (step * step);

            result.add(-1 * p.evaluate(x)*ddu + q.evaluate(x)*du + r.evaluate(x)*w.evaluate(x));

            x += step;
        }
        return result;
    }

    // Difference methods
    private List<PlanePoint> solveBvpFirstApproximation(List<PlanePoint> mesh) {

        double step = mesh.get(1).getX() - mesh.get(0).getX();
        double[][] coefficients = new double[mesh.size() + 1][4];

        coefficients[0][0] = 0;
        coefficients[0][1] = step * boundaryConstraints[0] + boundaryConstraints[1];
        coefficients[0][2] = boundaryConstraints[1];
        coefficients[0][3] = -1 * step * boundaryConstraints[2];

        coefficients[mesh.size()][0] = boundaryConstraints[4];
        coefficients[mesh.size()][1] = step * boundaryConstraints[3] + boundaryConstraints[4];
        coefficients[mesh.size()][2] = 0;
        coefficients[mesh.size()][3] = -1 * step * boundaryConstraints[5];

        for (int i = 1; i < mesh.size(); i++) {
            double x = mesh.get(i).getX();
            coefficients[i][0] = -1 * p.evaluate(x) - q.evaluate(x) * step / 2;
            coefficients[i][2] = -1 * p.evaluate(x) + q.evaluate(x) * step / 2;
            coefficients[i][1] = coefficients[i][0] + coefficients[i][2] - Math.pow(step, 2) * r.evaluate(x);
            coefficients[i][3] = Math.pow(step, 2) * f.evaluate(x);
        }

        return recurrentRun(mesh, coefficients);
    }

    private List<PlanePoint> solveBvpSecondApproximation(List<PlanePoint> mesh) {

        double step = mesh.get(1).getX() - mesh.get(0).getX();
        List<PlanePoint> tempMesh = generateMesh(mesh.get(0).getX() - step / 2,
            mesh.get(mesh.size() - 1).getX() + step / 2, step);

        double[][] coefficients = new double[tempMesh.size() + 1][4];

        coefficients[0][0] = 0;
        coefficients[0][1] = step * boundaryConstraints[0] + 2 * boundaryConstraints[1];
        coefficients[0][2] = 2 * boundaryConstraints[1] - step * boundaryConstraints[0];
        coefficients[0][3] = -2 * step * boundaryConstraints[2];

        coefficients[tempMesh.size()][0] = 2 * boundaryConstraints[4] - step * boundaryConstraints[3];
        coefficients[tempMesh.size()][1] = step * boundaryConstraints[3] + 2 * boundaryConstraints[4];
        coefficients[tempMesh.size()][2] = 0;
        coefficients[tempMesh.size()][3] = -2 * step * boundaryConstraints[5];

        for (int i = 1; i < tempMesh.size(); i++) {
            double x = tempMesh.get(i).getX();
            coefficients[i][0] = -1 * p.evaluate(x - step / 2) - q.evaluate(x) * step / 2;
            coefficients[i][2] = -1 * p.evaluate(x + step / 2) + q.evaluate(x) * step / 2;
            coefficients[i][1] = coefficients[i][0] + coefficients[i][2] - Math.pow(step, 2) * r.evaluate(x);
            coefficients[i][3] = Math.pow(step, 2) * f.evaluate(x);
        }

        List<PlanePoint> tempResult = recurrentRun(tempMesh, coefficients);
        List<PlanePoint> result = new ArrayList<>();
        for (int i = 0; i < mesh.size(); i++) {
            result.add(new PlanePoint(mesh.get(i).getX(), (tempResult.get(i).getY() + tempResult.get(i + 1).getY()) / 2));
        }
        return result;
    }

    private List<PlanePoint> solveBvpRichardsonFirstApproximation(List<PlanePoint> mesh) {

        double step = mesh.get(1).getX() - mesh.get(0).getX();
        double start = mesh.get(0).getX();
        double end = mesh.get(mesh.size() - 1).getX();

        List<PlanePoint> y = solveBvpFirstApproximation(mesh);
        List<PlanePoint> y2 = solveBvpFirstApproximation(generateMesh(start, end, step / 2));

        double[] errors = new double[mesh.size()];

        for (int i = 0; i < mesh.size(); i++) {
            errors[i] = y2.get(i * 2).getY() - y.get(i).getY();
        }

        return calculateRichardsonExtrapolation(mesh, y2, errors);
    }

    private List<PlanePoint> solveBvpRichardsonSecondApproximation(List<PlanePoint> mesh) {

        double step = mesh.get(1).getX() - mesh.get(0).getX();
        double start = mesh.get(0).getX();
        double end = mesh.get(mesh.size() - 1).getX();

        List<PlanePoint> y = solveBvpSecondApproximation(mesh);
        List<PlanePoint> y2 = solveBvpSecondApproximation(generateMesh(start, end, step / 2));

        double[] errors = new double[mesh.size()];

        for (int i = 0; i < mesh.size(); i++) {
            errors[i] = (y2.get(i * 2).getY() - y.get(i).getY()) / 3;
        }

        return calculateRichardsonExtrapolation(mesh, y2, errors);
    }

    private List<PlanePoint> calculateRichardsonExtrapolation(List<PlanePoint> mesh,
        List<PlanePoint> y2, double[] errors) {

        List<PlanePoint> result = new ArrayList<>();
        for (int i = 0; i < mesh.size(); i++) {
            result.add(new PlanePoint(mesh.get(i).getX(), y2.get(i * 2).getY() + errors[i]));
        }

        return result;
    }

    private List<PlanePoint> recurrentRun(List<PlanePoint> mesh, double[][] coefficients) {

        double[] s = new double[mesh.size() + 1];
        double[] t = new double[mesh.size() + 1];

        s[0] = coefficients[0][2] / coefficients[0][1];
        t[0] = -1 * coefficients[0][3] / coefficients[0][1];

        for (int i = 1; i < mesh.size() + 1; i++) {
            s[i] = coefficients[i][2] / (coefficients[i][1] - coefficients[i][0] * s[i-1]);
            t[i] = (coefficients[i][0] * t[i-1] - coefficients[i][3])
                / (coefficients[i][1] - coefficients[i][0] * s[i-1]);
        }

        PlanePoint y = new PlanePoint(mesh.get(mesh.size() - 1).getX(), t[mesh.size()]);
        mesh.set(mesh.size() - 1, y);

        for (int i = mesh.size() - 2; i >= 0; i--) {
            y = new PlanePoint(mesh.get(i).getX(),s[i] * mesh.get(i + 1).getY() + t[i]);
            mesh.set(i, y);
        }

        return mesh;
    }

}
