package org.math.computational.functions;

import org.math.computational.PlanePoint;

import java.util.ArrayList;
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


		return createInterpolationPolynomialNewton(N, x);
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
		double valueNewton = createInterpolationPolynomialNewton(N, x);
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

	private double createInterpolationPolynomialNewton(int N, double x){

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

		double quotient = 1;
		double value = 0;
		for (int i = 0; i < N; i++){
			value += dividedDifferences[i][0] * quotient;
			quotient *= (x - nodes.get(i).getX());
		}

		return value;
	}
}
