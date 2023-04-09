package org.math.computational.functions;

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
	 * Return a value of definite integral of f(x) form A to B.
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

	private double integrateLeftRectangles() {
		return f.evaluate(lowerBound) * (upperBound - lowerBound);
	}

	private double integrateRightRectangles() {
		return f.evaluate(upperBound) * (upperBound - lowerBound);
	}

	private double integrateCenterRectangles() {
		return f.evaluate((lowerBound + upperBound) / 2) * (upperBound - lowerBound);
	}

	private double integrateTrapezium() {
		return (f.evaluate(lowerBound) + f.evaluate(upperBound)) * (upperBound - lowerBound) / 2;
	}

	private double integrateSimpson() {
		double center = (lowerBound + upperBound) / 2;
		return (upperBound - lowerBound) * (f.evaluate(lowerBound) + 4*f.evaluate(center) + f.evaluate(upperBound)) / 6;
	}

	private double integrate38() {
		double oneThird = lowerBound + (upperBound - lowerBound) / 3;
		double twoThird = lowerBound + 2 * (upperBound - lowerBound) / 3;
		return (upperBound - lowerBound) * (f.evaluate(lowerBound) + 3 * f.evaluate(oneThird) + 3 * f.evaluate(twoThird)
				+ f.evaluate(upperBound)) / 8;
	}

}
