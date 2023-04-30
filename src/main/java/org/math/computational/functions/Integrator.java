package org.math.computational.functions;

import org.math.computational.PlanePoint;
import org.math.computational.matrices.DenseMatrix;
import org.math.computational.matrices.LinearSystemSolver;

import java.util.ArrayList;
import java.util.List;

public class Integrator {

	private final Function f;
	private final double lowerBound;
	private final double upperBound;

	public Integrator(Function f, double lowerBound, double upperBound) {
		this.f = f;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}


	/**
	 * Return a value of definite integral of f(x) from A to B.
	 * Where f(x), A and B stored in Integrator class instance.
	 *
	 * @param method method of integration, by default - 3/8 method, "LR" - left rectangles method,
	 * "RR" - right rectangles method, "CR" - center rectangles method, "TR" - trapezium method,
	 * "SI" - Simpson's method
	 */
	public double integrate(String method) {

		return switch (method) {
			case "LR" -> integrateLeftRectangles();
			case "RR" -> integrateRightRectangles();
			case "CR" -> integrateCenterRectangles();
			case "TR" -> integrateTrapezium();
			case "SI" -> integrateSimpson();
			default -> integrate38();
		};
	}

	/**
	 * Return a value of definite integral of f(x) from A to B.
	 * Where f(x), A and B stored in Integrator class instance.
	 * Integral value will be calculated using composite quadrature formulas.
	 *
	 * @param method method of integration, by default - Simpson's method, "LR" - left rectangles method,
	 * "RR" - right rectangles method, "CR" - center rectangles method, "TR" - trapezium method
	 * @param subsegments number of subsegments on which integral values will be calculated
	 */
	public double integrate(String method, int subsegments) {

		double step = (upperBound - lowerBound) / subsegments;

		return switch (method) {
			case "LR" -> integrateLeftRectangles(step);
			case "RR" -> integrateRightRectangles(step);
			case "CR" -> integrateCenterRectangles(step);
			case "TR" -> integrateTrapezium(step);
			default -> integrateSimpson(step);
		};

	}

	private double integrateLeftRectangles() {
		return f.evaluate(lowerBound) * (upperBound - lowerBound);
	}

	private double integrateLeftRectangles(double step) {

		double value = 0;
		double x = lowerBound;

		while (Math.abs(x - upperBound) > 10e-6 && x < upperBound) {
			value += f.evaluate(x);
			x += step;
		}

		return step * value;
	}

	private double integrateRightRectangles() {
		return f.evaluate(upperBound) * (upperBound - lowerBound);
	}

	private double integrateRightRectangles(double step) {

		double value = 0;
		double x = lowerBound + step;

		while (Math.abs(x - upperBound) > 10e-6 && x < upperBound) {
			value += f.evaluate(x);
			x += step;
		}

		value += f.evaluate(upperBound);

		return step * value;
	}

	private double integrateCenterRectangles() {
		return f.evaluate((lowerBound + upperBound) / 2) * (upperBound - lowerBound);
	}

	private double integrateCenterRectangles(double step) {

		double value = 0;
		double x = lowerBound + (step / 2.0);

		while (x < upperBound) {
			value += f.evaluate(x);
			x += step;
		}

		return step * value;
	}

	private double integrateTrapezium() {
		return (f.evaluate(lowerBound) + f.evaluate(upperBound)) * (upperBound - lowerBound) / 2;
	}

	private double integrateTrapezium(double step) {

		double z = f.evaluate(lowerBound) + f.evaluate(upperBound);
		double w = 0;
		double x = lowerBound + step;

		while (Math.abs(x - upperBound) > 10e-10 && x < upperBound) {
			w += f.evaluate(x);
			x += step;
		}

		return (z + 2 * w) * step / 2;
	}

	private double integrateSimpson() {
		double center = (lowerBound + upperBound) / 2;
		return (upperBound - lowerBound) * (f.evaluate(lowerBound) + 4*f.evaluate(center) + f.evaluate(upperBound)) / 6;
	}

	private double integrateSimpson(double step) {

		double z = f.evaluate(lowerBound) + f.evaluate(upperBound);
		double w = 0;
		double q = 0;
		double h = step / 2.0;
		double x = lowerBound + h;

		while (x < upperBound) {
			q += f.evaluate(x);
			w += f.evaluate(x + h);
			x += step;
		}
		w -= f.evaluate(upperBound);

		return (z + (2 * w) + (4 * q)) * step / 6;
	}

	private double integrate38() {
		double oneThird = lowerBound + (upperBound - lowerBound) / 3;
		double twoThird = lowerBound + 2 * (upperBound - lowerBound) / 3;
		return (upperBound - lowerBound) * (f.evaluate(lowerBound) + 3 * f.evaluate(oneThird) + 3 * f.evaluate(twoThird)
				+ f.evaluate(upperBound)) / 8;
	}

	private List<PlanePoint> buildInterpolationQuadrature (int N, List<Double> nodes, List<Double> momentums) {

		double[][] m = new double[N][N];
		double[][] momentum_array = new double[N][1];

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				m[i][j] = Math.pow(nodes.get(j), i);
			}
			momentum_array[i][0] = momentums.get(i);
		}

		DenseMatrix M = new DenseMatrix(N, N, m);
		DenseMatrix mu = new DenseMatrix(N, 1, momentum_array);
		LinearSystemSolver linearSolver = new LinearSystemSolver(M, mu);
		List<Double> A = linearSolver.solve();

		List<PlanePoint> result = new ArrayList<>();

		for (int i = 0; i < N; i++) {
			result.add(new PlanePoint(nodes.get(i), A.get(i)));
		}

		return result;
	}

	public List<PlanePoint> buildInterpolationQuadrature(Function rho, int N, List<Double> nodes) {

		List<Double> momentums = new ArrayList<>();

		for (int i = 0; i < N; i++) {
			int k = i;
			Function g = t -> rho.evaluate(t) * Math.pow(t, k);
			Integrator integrator = new Integrator(g, lowerBound, upperBound);
			momentums.add(integrator.integrate("CR", 1000));
		}

		return buildInterpolationQuadrature(N, nodes, momentums);
	}

	public double integrateInterpolationQuadrature (Function weight, int subsegments, List<Double> nodes) {

		List<PlanePoint> A = buildInterpolationQuadrature(weight, subsegments, nodes);
		double value = 0;
		for (int i = 0; i < subsegments; i ++) {
			value += A.get(i).getY() * f.evaluate(nodes.get(i));
		}

		return value;
	}

}
