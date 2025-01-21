package org.math.computations.pde;

import static org.math.computations.pde.Utils.getValuesOnMesh;

import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.math.computations.PlanePoint;
import org.math.computations.functions.FunctionalR2;

public class EllipticalDirichletSolver {

    private final int MAX_ITERATIONS = 20;
    private final double[] boundaries;
    private final FunctionalR2 p;
    private final FunctionalR2 q;
    private final FunctionalR2 f;
    private final FunctionalR2 mu;
    private double[] constants;
    private double[][] fValues;
    private double absErrorInitialApprox;
    private FunctionalR2 uExact;
    private double[] chebyshevParameters = null;
    private double[] alternatingTriangleParameters = null;

    public EllipticalDirichletSolver(FunctionalR2 p, FunctionalR2 q, FunctionalR2 f, FunctionalR2 mu,
        double[] boundaries) {

        this.p = p;
        this.q = q;
        this.f = f;
        this.mu = mu;
        this.boundaries = boundaries;
    }

    public EllipticalDirichletSolver(FunctionalR2 p, FunctionalR2 q, FunctionalR2 f, FunctionalR2 mu,
        double[] boundaries, FunctionalR2 uExact) {

        this.p = p;
        this.q = q;
        this.f = f;
        this.mu = mu;
        this.boundaries = boundaries;
        this.uExact = uExact;
    }

    private long approximateIterationsAmount(String method, double epsilon, double xi) {
        return approximateIterationsAmount(method, epsilon, xi, 1);
    }

    private long approximateIterationsAmount(String method, double epsilon, int stepAmount) {
        return approximateIterationsAmount(method, epsilon, 1, stepAmount);
    }

    private long approximateIterationsAmount(String method, double epsilon, double xi, int stepAmount) {
        long maxIterations;
        switch (method) {
            case "SI", "IOP":
                maxIterations = Math.round(Math.ceil(Math.log(1 / epsilon) / (2*xi)));
                break;
            case "SE":
                maxIterations = Math.round(Math.ceil(Math.log(1 / epsilon) / (4*xi)));
                break;
            case "UR":
                maxIterations = Math.round(Math.ceil(Math.log(1 / epsilon) / Math.sqrt(xi)));
                break;
            case "ICP":
                maxIterations = Math.round(Math.ceil(Math.log(2 / epsilon) / (2*Math.sqrt(xi))));
                maxIterations = (int) Math.pow(2, Math.floor(Math.log(maxIterations)/Math.log(2)));
                break;
            case "AT":
                maxIterations = Math.round(Math.ceil(Math.log(1 / epsilon) / Math.log((1 + xi) / (1 - xi))));
                break;
            case "ATCP":
                maxIterations = Math.round(Math.ceil(Math.log(2 / epsilon) / (2 * Math.sqrt(2) * Math.pow(xi, 0.25))));
                maxIterations = (int) (Math.pow(2, Math.floor(Math.log(maxIterations)/Math.log(2))));
                break;
            case "AD":
                maxIterations = Math.round(Math.ceil(stepAmount * Math.log(1 / epsilon) / (2 * Math.PI)));
                break;
            default:
                maxIterations = MAX_ITERATIONS;
                break;
        }

        return maxIterations;
    }

    public double[][] solve(String method, int[] stepAmounts, double epsilon) {

        double stepX = (boundaries[1] - boundaries[0]) / stepAmounts[0];
        double stepY = (boundaries[3] - boundaries[2]) / stepAmounts[1];

        PlanePoint[][] mesh = Utils.generate2dMesh(this.boundaries, stepAmounts);
        this.fValues = getValuesOnMesh(this.f, mesh);
        double[][] negativeFValues = getValuesOnMesh((x, y) -> -1 * this.f.evaluate(x, y), mesh);

        double[] minMaxP = getMaxAndMinOfFunctional(this.p, mesh);
        double[] minMaxQ = getMaxAndMinOfFunctional(this.q, mesh);
        this.constants = new double[] {minMaxP[0], minMaxP[1], minMaxQ[0], minMaxQ[1]};

        double minEigenvalue = getMinEigenvalue(stepX, stepY);
        double maxEigenvalue = getMaxEigenvalue(stepX, stepY);
        double spectralRadius = (maxEigenvalue - minEigenvalue) / (maxEigenvalue + minEigenvalue);

        double[][] u = setBoundaryConditions(mesh);

        double[][] Lu = getLuValues(u, mesh, stepX, stepY);
        double initialMeanDifference = getMean(Lu, negativeFValues);

        System.out.printf(Locale.US, "Mean of difference for initial approximation = %3.5e%n",
            initialMeanDifference);
        System.out.printf(Locale.US, "Spectral radius = %3.5e%n", spectralRadius);

        double xi = minEigenvalue / maxEigenvalue;
        long maxIterations = method.equals("AD")
            ? approximateIterationsAmount(method, epsilon, stepAmounts[0])
            : approximateIterationsAmount(method, epsilon, xi);
        System.out.printf("Approximate amount of iterations = %d%n", maxIterations);

        if (this.uExact == null) {
            System.out.printf(
                " K   | Mean diff. | Relative diff | Approx. diff | Est. error | Sp. rad. approx |%n");
        } else {
            double[][] exactValues = getLuValues(getValuesOnMesh(this.uExact, mesh), mesh, stepX, stepY);
            double meanApproximationExact = getMean(exactValues, negativeFValues);
            System.out.printf(Locale.US, "Mean of approximation on exact solution = %3.5e%n",
                meanApproximationExact);
            this.absErrorInitialApprox = getMean(u, getValuesOnMesh(uExact, mesh));
            System.out.printf(
                " K   | Mean diff. | Relative diff | Abs. error | Rel. error | Approx. diff | Est. error | Sp. rad. approx |%n");
        }

        if (method.equals("ICP")) {
            this.chebyshevParameters = getChebyshevParameters(maxIterations, minEigenvalue, maxEigenvalue);
        }

        if (method.equals("AT") || method.equals("ATCP")) {
            double omega = 2 / (Math.sqrt(maxEigenvalue * minEigenvalue));
            double gamma1 = minEigenvalue / (2 + 2 * Math.sqrt(minEigenvalue / maxEigenvalue));
            double gamma2 = minEigenvalue / (4 * Math.sqrt(minEigenvalue / maxEigenvalue));

            double kappa1 = omega / (stepX * stepX);
            double kappa2 = omega / (stepY * stepY);
            this.alternatingTriangleParameters = new double[] {omega, gamma1, gamma2, kappa1, kappa2};

            this.chebyshevParameters = method.equals("ATCP") ? getChebyshevParametersAlternatingTriangles(maxIterations, gamma1, gamma2) : null;
        }
        // Only useful when method is "AT" or "ATCP"
        double[][] upperMatrix = new double[mesh.length][mesh[0].length];
        double[][] lowerMatrix = new double[mesh.length][mesh[0].length];
        double[][] fMatrix = new double[mesh.length][mesh[0].length];

        int counter = 0;
        double errorEstimate = 100.0;
        double lastApproxDiff = 1;

        while (errorEstimate > epsilon && counter < maxIterations) {

            double[][] pu = copyApproximation(u);

            switch (method) {
                case "SI":
                    u = simpleIterationStep(mesh, pu, stepX, stepY);
                    break;
                case "IOP":
                    u = iterationWithOptimalParameterStep(mesh, pu, stepX, stepY,
                        2 / (minEigenvalue + maxEigenvalue));
                    break;
                case "SE":
                    u = seidelMethodStep(mesh, pu, stepX, stepY);
                    break;
                case "UR":
                    u = upperRelaxationStep(mesh, pu, stepX, stepY,
                        2 / (1 + Math.sqrt(1 - Math.pow(spectralRadius, 2))));
                    break;
                case "ICP":
                    u = iterationWithChebyshevParameters(mesh, pu, stepX, stepY, this.chebyshevParameters[counter]);
                    break;
                case "AT":
                    fMatrix = alternatingTrianglesUpdateF(mesh, Lu);
                    lowerMatrix = alternatingTrianglesUpdateLower(mesh, fMatrix, stepX, stepY);
                    upperMatrix = alternatingTrianglesUpdateUpper(mesh, lowerMatrix, stepX, stepY);
                    u = alternatingTrianglesUpdateApproximation(mesh, pu, upperMatrix,
                        2 / (this.alternatingTriangleParameters[1] + this.alternatingTriangleParameters[2]));
                    break;
                case "ATCP":
                    fMatrix = alternatingTrianglesUpdateF(mesh, Lu);
                    lowerMatrix = alternatingTrianglesUpdateLower(mesh, fMatrix, stepX, stepY);
                    upperMatrix = alternatingTrianglesUpdateUpper(mesh, lowerMatrix, stepX, stepY);
                    u = alternatingTrianglesUpdateApproximation(mesh, pu, upperMatrix, this.chebyshevParameters[counter]);
                    break;
                case "AD":
                    u = alternatingDirections(mesh, pu, stepX, stepY, Math.pow(stepX, 2) / Math.sin(Math.PI * stepX));
                    break;
                default: return null;
            }

            Lu = getLuValues(u, mesh, stepX, stepY);
            double meanDifference = getMean(Lu, negativeFValues);
            double relativeMeanDifference = meanDifference / initialMeanDifference;

            double approxDiff = getMean(u, pu);
            double approxSpectralRadius = approxDiff / lastApproxDiff;
            errorEstimate = spectralRadius * approxDiff / (1 - spectralRadius);

            if (this.uExact == null) {
                System.out.printf(Locale.US, "%d    | %3.4e  | %3.4e  | %3.4e  | %3.4e  | %3.4e |%n",
                    counter + 1, meanDifference, relativeMeanDifference, approxDiff, errorEstimate, approxSpectralRadius);
            } else {
                double absError = getMean(u, getValuesOnMesh(uExact, mesh));
                double relError = absError / this.absErrorInitialApprox;
                System.out.printf(Locale.US,
                    "%d    | %3.4e  | %3.4e  | %3.4e  | %3.4e | %3.4e  | %3.4e  | %3.4e |%n",
                    counter + 1, meanDifference, relativeMeanDifference, absError, relError, approxDiff,
                    errorEstimate, approxSpectralRadius);
            }

            lastApproxDiff = approxDiff;
            counter++;
        }

        return u;
    }

    private double getMean(double[][] a, double[][] b) {

        if (a.length != b.length || a[0].length != b[0].length) {
            throw new RuntimeException("Dimensions doesn't match.");
        }

        double mean = 0;
        for (int i = 1; i < a.length - 1; i++) {
            for (int j = 1; j < a[0].length - 1; j++) {
                mean = Math.max(mean, Math.abs(a[i][j] - b[i][j]));
            }
        }

        return mean;
    }

    private double[] getMaxAndMinOfFunctional(FunctionalR2 f, PlanePoint[][] mesh) {

        double[] minmax = new double[2];
        minmax[0] = f.evaluate(mesh[0][0].getX(), mesh[0][0].getY());
        minmax[1] = f.evaluate(mesh[0][0].getX(), mesh[0][0].getY());

        for (int i = 0; i < mesh.length; i++) {
            for (int j = 0; j < mesh[0].length; j++) {
                double value = f.evaluate(mesh[i][j].getX(), mesh[i][j].getY());
                minmax[0] = Math.min(minmax[0], value);
                minmax[1] = Math.max(minmax[1], value);
            }
        }

        return minmax;
    }

    private double getMinEigenvalue(double stepX, double stepY) {
        double lenX = this.boundaries[1] - this.boundaries[0];
        double lenY = this.boundaries[3] - this.boundaries[2];
        return this.constants[0] * 4 * Math.pow(Math.sin((Math.PI * stepX) / (2 * lenX)), 2) / Math.pow(stepX, 2)
            + this.constants[2] * 4 * Math.pow(Math.sin((Math.PI * stepY) / (2 * lenY)), 2) / Math.pow(stepY, 2);
    }

    private double getMaxEigenvalue(double stepX, double stepY) {
        double lenX = this.boundaries[1] - this.boundaries[0];
        double lenY = this.boundaries[3] - this.boundaries[2];
        return this.constants[1] * 4 * Math.pow(Math.cos((Math.PI * stepX) / (2 * lenX)), 2) / Math.pow(stepX, 2)
            + this.constants[3] * 4 * Math.pow(Math.cos((Math.PI * stepY) / (2 * lenY)), 2) / Math.pow(stepY, 2);
    }

    private double[][] getLuValues(double[][] approximation, PlanePoint[][] mesh, double stepX, double stepY) {

        double[][] result = new double[mesh.length][mesh[0].length];

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                result[i][j] = (coeffsRight[0] * (approximation[i + 1][j] - approximation[i][j]) -
                    coeffsLeft[0] * (approximation[i][j] - approximation[i - 1][j])) / Math.pow(stepX, 2) +
                    (coeffsRight[1] * (approximation[i][j + 1] - approximation[i][j]) -
                    coeffsLeft[1] * (approximation[i][j] - approximation[i][j - 1])) / Math.pow(stepY, 2);
            }
        }

        return result;
    }

    private double[][] simpleIterationStep(PlanePoint[][] mesh, double[][] pu, double stepX, double stepY) {

        double[][] u = setBoundaryConditions(mesh);

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                double numerator = this.fValues[i][j] + (coeffsLeft[0] * pu[i - 1][j] + coeffsRight[0] * pu[i + 1][j]) / (stepX * stepX)
                        + (coeffsLeft[1] * pu[i][j - 1] + coeffsRight[1] * pu[i][j + 1]) / (stepY * stepY);

                double denominator = (coeffsLeft[0] + coeffsRight[0]) / (stepX * stepX)
                        + (coeffsLeft[1] + coeffsRight[1]) / (stepY * stepY);

                u[i][j] = numerator / denominator;
            }
        }

        return u;
    }

    private double[][] iterationWithOptimalParameterStep(PlanePoint[][] mesh, double[][] pu, double stepX, double stepY, double parameter) {

        double[][] u = setBoundaryConditions(mesh);

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                double value = this.fValues[i][j] +
                    (coeffsRight[0] * (pu[i + 1][j] - pu[i][j]) - coeffsLeft[0] * (pu[i][j] - pu[i - 1][j])) / (stepX * stepX)
                    + (coeffsRight[1] * (pu[i][j + 1] - pu[i][j]) - coeffsLeft[1] * (pu[i][j] - pu[i][j - 1])) / (stepY * stepY);

                u[i][j] = pu[i][j] + parameter * value;
            }
        }

        return u;
    }

    private double[][] seidelMethodStep(PlanePoint[][] mesh, double[][] pu, double stepX, double stepY) {

        double[][] u = setBoundaryConditions(mesh);

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                double numerator = this.fValues[i][j] + (coeffsLeft[0] * u[i - 1][j] + coeffsRight[0] * pu[i + 1][j]) / (stepX * stepX)
                    + (coeffsLeft[1] * u[i][j - 1] + coeffsRight[1] * pu[i][j + 1]) / (stepY * stepY);

                double denominator = (coeffsLeft[0] + coeffsRight[0]) / (stepX * stepX)
                    + (coeffsLeft[1] + coeffsRight[1]) / (stepY * stepY);

                u[i][j] = numerator / denominator;
            }
        }

        return u;
    }

    private double[][] upperRelaxationStep(PlanePoint[][] mesh, double[][] pu, double stepX, double stepY, double parameter) {

        double[][] u = setBoundaryConditions(mesh);

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                double numerator = this.fValues[i][j] +
                    (coeffsRight[0] * (pu[i + 1][j] - pu[i][j]) - coeffsLeft[0] * (pu[i][j] - u[i - 1][j])) / (stepX * stepX)
                    + (coeffsRight[1] * (pu[i][j + 1] - pu[i][j]) - coeffsLeft[1] * (pu[i][j] - u[i][j - 1])) / (stepY * stepY);
                double denominator = (coeffsLeft[0] + coeffsRight[0]) / (stepX * stepX)
                    + (coeffsLeft[1] + coeffsRight[1]) / (stepY * stepY);

                u[i][j] = pu[i][j] + parameter * numerator / denominator;
            }
        }

        return u;
    }

    private double[][] iterationWithChebyshevParameters(PlanePoint[][] mesh, double[][] pu, double stepX, double stepY, double parameter) {

        double[][] u = setBoundaryConditions(mesh);

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                double value = this.fValues[i][j] +
                    (coeffsRight[0] * (pu[i + 1][j] - pu[i][j]) - coeffsLeft[0] * (pu[i][j] - pu[i - 1][j])) / (stepX * stepX)
                    + (coeffsRight[1] * (pu[i][j + 1] - pu[i][j]) - coeffsLeft[1] * (pu[i][j] - pu[i][j - 1])) / (stepY * stepY);

                u[i][j] = pu[i][j] + parameter * value;
            }
        }

        return u;
    }

    private double[][] alternatingTrianglesUpdateF(PlanePoint[][] mesh, double[][] Lu) {

        double[][] newF = new double[mesh.length][mesh[0].length];

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {
                newF[i][j] = Lu[i][j] + f.evaluate(mesh[i][j].getX(), mesh[i][j].getY());
            }
        }

        return newF;
    }

    private double[][] alternatingTrianglesUpdateLower(PlanePoint[][] mesh, double[][] fMatrix, double stepX, double stepY) {

        double[][] newLower = new double[mesh.length][mesh[0].length];

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();
                double kappa1 = this.alternatingTriangleParameters[3];
                double kappa2 = this.alternatingTriangleParameters[4];

                double pValue = p.evaluate(x - stepX / 2, y);
                double qValue = q.evaluate(x, y - stepY / 2);

                newLower[i][j] = (kappa1 * pValue * newLower[i - 1][j]
                    + kappa2 * qValue * newLower[i][j - 1] + fMatrix[i][j]) / (1
                    + kappa1 * pValue + kappa2 * qValue);
            }
        }

        return newLower;
    }

    private double[][] alternatingTrianglesUpdateUpper(PlanePoint[][] mesh, double[][] lowerMatrix, double stepX, double stepY) {

        double[][] newUpper = new double[mesh.length][mesh[0].length];

        for (int i = mesh.length - 2; i >= 1; i--) {
            for (int j = mesh[0].length - 2; j >= 1; j--) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();
                double kappa1 = this.alternatingTriangleParameters[3];
                double kappa2 = this.alternatingTriangleParameters[4];

                double pValue = p.evaluate(x + stepX / 2, y);
                double qValue = q.evaluate(x, y + stepY / 2);

                newUpper[i][j] = (kappa1 * pValue * newUpper[i + 1][j]
                    + kappa2 * qValue * newUpper[i][j + 1] + lowerMatrix[i][j]) / (1
                    + kappa1 * pValue + kappa2 * qValue);
            }
        }

        return newUpper;
    }

    private double[][] alternatingTrianglesUpdateApproximation(PlanePoint[][] mesh, double[][] pu, double[][] upperMatrix, double parameter) {

        double[][] u = setBoundaryConditions(mesh);

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {
                u[i][j] = pu[i][j] + parameter * upperMatrix[i][j];
            }
        }

        return u;
    }

    private double[][] alternatingDirections(PlanePoint[][] mesh, double[][] pu, double stepX, double stepY, double tau) {

        double[][] intermediateU = setBoundaryConditions(mesh);
        double[][] lambda2Values = getLambda2Values(mesh, pu, stepX, stepY);

        // First step, u^k -> u^(k+0.5)
        for (int j = 1; j < mesh[0].length - 1; j++) {

            double[][] coefficientsY = new double[mesh.length][4];

            coefficientsY[0][0] = 0;
            coefficientsY[0][1] = -1;
            coefficientsY[0][2] = 0;
            coefficientsY[0][3] = mu.evaluate(mesh[0][j].getX(), mesh[0][j].getY());
            coefficientsY[mesh.length - 1][0] = 0;
            coefficientsY[mesh.length - 1][1] = -1;
            coefficientsY[mesh.length - 1][2] = 0;
            coefficientsY[mesh.length - 1][3] = mu.evaluate(mesh[mesh.length - 1][j].getX(), mesh[mesh.length - 1][j].getY());

            for (int i = 1; i < mesh.length - 1; i++) {
                coefficientsY[i][0] = tau / (2 * Math.pow(stepX, 2));
                coefficientsY[i][1] = tau / Math.pow(stepX, 2) + 1;
                coefficientsY[i][2] = tau / (2 * Math.pow(stepX, 2));
                coefficientsY[i][3] = -pu[i][j] - tau * (lambda2Values[i][j] + fValues[i][j]) / 2;
            }

            double[][] partialResult = recurrentRunAlongX(pu, coefficientsY, j);
            for (int i = 1; i < mesh.length - 1; i++) {
                intermediateU[i][j] = partialResult[i][j];
            }
        }

        double[][] u = setBoundaryConditions(mesh);
        double[][] lambda1Values = getLambda1Values(mesh, intermediateU, stepX, stepY);

        // Second step, u^(k+0.5) -> u^(k+1)
        for (int i = 1; i < mesh.length - 1; i++) {

            double[][] coefficientsX = new double[mesh[0].length][4];

            coefficientsX[0][0] = 0;
            coefficientsX[0][1] = -1;
            coefficientsX[0][2] = 0;
            coefficientsX[0][3] = mu.evaluate(mesh[i][0].getX(), mesh[i][0].getY());
            coefficientsX[mesh[0].length - 1][0] = 0;
            coefficientsX[mesh[0].length - 1][1] = -1;
            coefficientsX[mesh[0].length - 1][2] = 0;
            coefficientsX[mesh[0].length - 1][3] = mu.evaluate(mesh[i][mesh[0].length - 1].getX(), mesh[i][mesh[0].length - 1].getY());

            for (int j = 1; j < mesh[0].length - 1; j++) {
                coefficientsX[j][0] = tau / (2 * Math.pow(stepY, 2));
                coefficientsX[j][1] = tau / Math.pow(stepY, 2) + 1;
                coefficientsX[j][2] = tau / (2 * Math.pow(stepY, 2));
                coefficientsX[j][3] = -intermediateU[i][j] - tau * (lambda1Values[i][j] + fValues[i][j]) / 2;
            }

            double[][] partialResult = recurrentRunAlongY(intermediateU, coefficientsX, i);
            for (int j = 1; j < mesh[0].length - 1; j++) {
                u[i][j] = partialResult[i][j];
            }
        }

        return u;
    }

    private double[][] recurrentRunAlongX(double[][] approximation, double[][] coefficients, int j) {

        double[][] result = new double[approximation.length][approximation[0].length];

        double[] s = new double[approximation.length + 1];
        double[] t = new double[approximation.length + 1];

        s[0] = coefficients[0][2] / coefficients[0][1];
        t[0] = -1 * coefficients[0][3] / coefficients[0][1];

        for (int i = 1; i < approximation.length; i++) {
            s[i] = coefficients[i][2] / (coefficients[i][1] - coefficients[i][0] * s[i - 1]);
            t[i] = (coefficients[i][0] * t[i - 1] - coefficients[i][3])
                / (coefficients[i][1] - coefficients[i][0] * s[i - 1]);
        }

        result[approximation.length - 1][j] = t[approximation.length - 1];

        for (int i = approximation.length - 2; i >= 0; i--) {
            result[i][j] = s[i] * result[i + 1][j] + t[i];
        }

        return result;
    }

    private double[][] recurrentRunAlongY(double[][] approximation, double[][] coefficients, int i) {

        double[][] result = new double[approximation.length][approximation[0].length];

        double[] s = new double[approximation[0].length + 1];
        double[] t = new double[approximation[0].length + 1];

        s[0] = coefficients[0][2] / coefficients[0][1];
        t[0] = -1 * coefficients[0][3] / coefficients[0][1];

        for (int k = 1; k < approximation[0].length; k++) {
            s[k] = coefficients[k][2] / (coefficients[k][1] - coefficients[k][0] * s[k - 1]);
            t[k] = (coefficients[k][0] * t[k - 1] - coefficients[k][3])
                / (coefficients[k][1] - coefficients[k][0] * s[k - 1]);
        }

        result[i][approximation[0].length - 1] = t[approximation[0].length - 1];

        for (int j = approximation[0].length - 2; j >= 0; j--) {
            result[i][j] = s[j] * result[i][j + 1] + t[j];
        }

        return result;
    }

    private double[][] getLambda1Values(PlanePoint[][] mesh, double[][] u, double stepX, double stepY) {

        double[][] result = new double[mesh.length][mesh[0].length];

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                result[i][j] = (coeffsRight[0] * (u[i + 1][j] - u[i][j]) -
                    coeffsLeft[0] * (u[i][j] - u[i - 1][j])) / Math.pow(stepX, 2);
            }
        }

        return result;
    }

    private double[][] getLambda2Values(PlanePoint[][] mesh, double[][] u, double stepX, double stepY) {

        double[][] result = new double[mesh.length][mesh[0].length];

        for (int i = 1; i < mesh.length - 1; i++) {
            for (int j = 1; j < mesh[0].length - 1; j++) {

                double x = mesh[i][j].getX();
                double y = mesh[i][j].getY();

                double[] coeffsLeft = {p.evaluate(x - stepX / 2, y), q.evaluate(x, y - stepY / 2)};
                double[] coeffsRight = {p.evaluate(x + stepX / 2, y), q.evaluate(x, y + stepY / 2)};

                result[i][j] = (coeffsRight[1] * (u[i][j + 1] - u[i][j]) -
                    coeffsLeft[1] * (u[i][j] - u[i][j - 1])) / Math.pow(stepY, 2);
            }
        }

        return result;
    }

    private double[][] setBoundaryConditions(PlanePoint[][] mesh) {

        double[][] u = new double[mesh.length][mesh[0].length];

        for (int i = 0; i < mesh.length; i++) {
            for (int j = 0; j < mesh[0].length; j++) {
                if (i == 0 || j == 0 || i == mesh.length - 1 || j == mesh[0].length - 1){
                    u[i][j] = this.mu.evaluate(mesh[i][j].getX(), mesh[i][j].getY());
                }
            }
        }

        return u;
    }

    private double[][] copyApproximation(double[][] approximation) {

        double[][] copy = new double[approximation.length][approximation[0].length];
        for (int i = 0; i < approximation.length; i++) {
            copy[i] = Arrays.copyOf(approximation[i], approximation[i].length);
        }

        return copy;
    }

    private double[] getChebyshevParameters(long n, double minEigenvalue, double maxEigenvalue) {

        int[] prevIndices = {1};

        for (int k = 1; k < n; k *= 2) {
            int[] indices = new int[2*k];
            for (int i = 0; i < k; i++) {
                indices[2 * i] = prevIndices[i];
                indices[2 * i + 1] = 4 * k - indices[2 * i];
            }
            prevIndices = Arrays.copyOf(indices, 2*k);
        }

        double[] parameters = new double[prevIndices.length];

        for (int i = 0; i < prevIndices.length; i++) {
            parameters[i] = 2 / (maxEigenvalue + minEigenvalue +
                (maxEigenvalue - minEigenvalue) * Math.cos((2 * i - 1) * Math.PI / (2 * prevIndices.length))
            );
        }

        return parameters;
    }

    private double[] getChebyshevParametersAlternatingTriangles(long n, double gamma1, double gamma2) {

        int[] prevIndices = {1};

        for (int k = 1; k < n; k *= 2) {
            int[] indices = new int[2*k];
            for (int i = 0; i < k; i++) {
                indices[2 * i] = prevIndices[i];
                indices[2 * i + 1] = 4 * k - indices[2 * i];
            }
            prevIndices = Arrays.copyOf(indices, 2*k);
        }

        double[] parameters = new double[prevIndices.length];

        for (int i = 0; i < prevIndices.length; i++) {
            parameters[i] = 2 / (gamma1 + gamma2 +
                (gamma2 - gamma1) * Math.cos((2 * i - 1) * Math.PI / (2 * prevIndices.length))
            );
        }

        return parameters;
    }

}
