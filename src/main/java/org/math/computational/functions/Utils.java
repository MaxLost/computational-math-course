package org.math.computational.functions;

public class Utils {

	public static double computeDerivative(double point, Function function, double h) {
		return (function.evaluate(point + h) - function.evaluate(point - h)) / h;
	}

	public static double computeDerivative(double point, Function function) {
		double h = 1e-8;
		return computeDerivative(point, function, h);
	};

}
