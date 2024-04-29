package org.math.computations.functions;

import java.util.ArrayList;
import java.util.List;

public class Utils {

	public static double computeDerivative(double point, Function function, double h) {
		return (function.evaluate(point + h) - function.evaluate(point - h)) / (2*h);
	}

	public static double computeDerivative(double point, Function function) {
		double h = 1e-8;
		return computeDerivative(point, function, h);
	};

	public static double computeSecondDerivative(double point, Function function) {
		double h = 1e-8;
		return (function.evaluate(point + h) - 2* function.evaluate(point) + function.evaluate(point - h)) / (h * h);
	}

	public static double computeThirdDerivative(double point, Function function) {
		double h = 1e-8;
		return (function.evaluate(point + 2*h) - 2*function.evaluate(point + h) +
				2*function.evaluate(point - h) - function.evaluate(point - 2*h)) / (2 * Math.pow(h, 3));
	}

	public static double factorial(int n) {

		double result = 1;

		for (int i = 1; i <= n; i++) {
			result *= i;
		}

		return result;
	}

	public static List<Double> generateMesh(double lowerBound, double upperBound, int nodesAmount) {

		double step = (upperBound - lowerBound) / (nodesAmount - 1);
		double x = lowerBound;

		List<Double> mesh = new ArrayList<>();

		for (int i = 0; i < nodesAmount; i++) {
			mesh.add(x);
			x += step;
		}

		return mesh;
	}

}
