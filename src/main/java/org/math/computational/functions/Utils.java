package org.math.computational.functions;

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

}
