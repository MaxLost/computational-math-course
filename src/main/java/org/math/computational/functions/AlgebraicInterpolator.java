package org.math.computational.functions;

import org.math.computational.PlanePoint;
import org.math.computational.Polynomial;
import org.math.computational.matrices.DenseMatrix;
import org.math.computational.matrices.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

public class AlgebraicInterpolator {

	private final List<PlanePoint> nodes;

	public AlgebraicInterpolator(List<PlanePoint> nodes) {
		this.nodes = nodes;
	}

	public double interpolate(int N, double x) {

		nodes.sort((a, b) -> (Math.abs(a.getX() - x) - Math.abs(b.getX() - x)) < 10e-8
							? 0
							: (Math.abs(a.getX() - x) > Math.abs(b.getX() - x)) ? 1 : -1);


		return createInterpolationPolynomialNewton(N, x).evaluate(x);
	}

	public List<Double> interpolateLogged(int N, double x){

		nodes.sort((a, b) -> Math.abs(Math.abs(a.getX() - x) - Math.abs(b.getX() - x)) < 10e-8
							? 0
							: ((Math.abs(a.getX() - x) > Math.abs(b.getX() - x)) ? 1 : -1));

		System.out.println("Таблица узлов для построения интерполяционного многочлена степени не выше " + N
							+ " функции f: \n\tx\t | \tf(x)\t");

		for (int i = 0; i < N + 1; i++){
			System.out.printf(Locale.US, "%.6f | %.6f\n", nodes.get(i).getX(), nodes.get(i).getY());
		}

		double valueLagrange = createInterpolationPolynomialLagrange(N, x);
		System.out.println("\nЗначение интерполяционного многочлена в форме Лагранжа в точке X: " + valueLagrange);
		double valueNewton = createInterpolationPolynomialNewton(N, x).evaluate(x);
		System.out.println("Значение интерполяционного многочлена в форме Ньютона в точке X: " + valueNewton);

		List<Double> result = new ArrayList<>();
		result.add(valueLagrange);
		result.add(valueNewton);

		return result;
	}

	private double createInterpolationPolynomialLagrange(int N, double x){

		double value = 0;
		for (int k = 0; k <= N; k++){

			double numerator = 1;
			double nominator = 1;
			for (int i = 0; i <= N; i++) {

				if (i != k) {

					numerator *= x - nodes.get(i).getX();
					nominator *= nodes.get(k).getX() - nodes.get(i).getX();
				}
			}

			value += (numerator / nominator) * nodes.get(k).getY();
		}

		return value;
	}

	private Polynomial createInterpolationPolynomialNewton(int N, double x){

		double[][] dividedDifferences = new double [N][N + 1];

		for (int i = 0; i < N + 1; i++) {
			dividedDifferences[0][i] = nodes.get(i).getY();
		}
		for (int i = 1; i < N + 1; i++) {
			for (int j = 0; j < N - i; j++){
				dividedDifferences[i][j] = (dividedDifferences[i - 1][j + 1] - dividedDifferences[i - 1][j])
											/ (nodes.get(j + i).getX() - nodes.get(j).getX());
			}
		}

		Polynomial p = new Polynomial(List.of(dividedDifferences[0][0]));
		Polynomial g = new Polynomial(List.of(1.0));

		for (int i = 1; i < N; i++){

			g = g.multiply(new Polynomial(List.of(-1*nodes.get(i - 1).getX(), 1.0)));
			Polynomial t = g.scalarMultiply(dividedDifferences[i][0]);

			p = p.sum(t);
		}

		return p;
	}

	public Polynomial getInterpolationPolynomial(int N, double x) {
		return createInterpolationPolynomialNewton(N, x);
	}

	public Polynomial getInterpolationPolynomialParker(int N, double x) {

		nodes.sort((a, b) -> (Math.abs(Math.abs(a.getX() - x) - Math.abs(b.getX() - x)) < 10e-8
				? 0
				: (Math.abs(a.getX() - x) > Math.abs(b.getX() - x)) ? 1 : -1));

		for (int i = 0; i < N + 1; i++){
			System.out.printf(Locale.US, "%.6f | %.6f\n", nodes.get(i).getX(), nodes.get(i).getY());
		}

		Polynomial P = new Polynomial(List.of(1.0));

		for (int i = 0; i < N; i++) {
			P = P.multiply(new Polynomial(List.of(-1*nodes.get(i).getX(), 1.0)));
		}

		double[][] q = new double[N + 1][N + 1];
		double[][] matrix = new double[N + 1][N + 1];

		for (int i = 0; i < N + 1; i++) {

			q[0][i] = 1;
			double dP = 1;

			for (int j = 1; j < N + 1; j++) {

				q[j][i] = nodes.get(i).getX() * q[j - 1][i] + P.getCoefficients().get(N - j);
				if (j != i) {
					dP *= (nodes.get(i).getX() - nodes.get(j).getX());
				}
			}

			for (int j = 0; j < N + 1; j++) {

				matrix[j][i] = q[j][i] / dP;
			}
		}

		double[][] y = new double[N + 1][1];
		for (int i = 0; i < N + 1; i++) {
			y[i][0] = nodes.get(i).getY();
		}

		DenseMatrix V = new DenseMatrix(N + 1, N + 1, matrix);
		DenseMatrix Y = new DenseMatrix(N + 1, 1, y);

		Matrix A = V.mul(Y);

		List<Double> coefficients = new ArrayList<>();
		for (int i = N; i >= 0; i--) {
			coefficients.add(A.getElement(0, i));
		}

		return new Polynomial(coefficients);
	}
}
