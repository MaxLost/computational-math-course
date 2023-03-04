package org.math.computational.functions;

import org.math.computational.PlanePoint;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class AlgebraicInterpolator {

	private List<PlanePoint> nodes;

	public AlgebraicInterpolator(List<PlanePoint> nodes) {
		this.nodes = nodes;
	}

	public double interpolate(int N, double x) {

		nodes.sort((a, b) -> (Math.abs(a.getX() - x) - Math.abs(b.getX() - x)) < 10e-8
							? 0
							: (Math.abs(a.getX() - x) > Math.abs(b.getX() - x)) ? 1 : -1);


		Function polynomialLagrange = createInterpolationPolynomialLagrange(N, x);
		Function polynomialNewton = createInterpolationPolynomialNewton(N, x);
		return 0;
	}

	public double interpolateLogged(int N, double x){

		nodes.sort((a, b) -> Math.abs(Math.abs(a.getX() - x) - Math.abs(b.getX() - x)) < 10e-8
							? 0
							: ((Math.abs(a.getX() - x) > Math.abs(b.getX() - x)) ? 1 : -1));

		System.out.println("Таблица узлов для построения интерполяционного многочлена степени " + N
							+ " функции f: \n\tx\t | \tf(x)\t");

		for (int i = 0; i < N; i++){
			System.out.printf(Locale.US, "%.6f | %.6f\n", nodes.get(i).getX(), nodes.get(i).getY());
		}

		Function polynomialLagrange = createInterpolationPolynomialLagrange(N, x);
		Function polynomialNewton = createInterpolationPolynomialNewton(N, x);
		return 0;
	}

	private Function createInterpolationPolynomialLagrange(int N, double x){


		return null;
	}

	private Function createInterpolationPolynomialNewton(int N, double x){


		return null;
	}
}
