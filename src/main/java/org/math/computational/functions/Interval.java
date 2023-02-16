package org.math.computational.functions;

public class Interval {

	private final double left_bound;
	private final double right_bound;

	public Interval(double a, double b){

		left_bound = a;
		right_bound = b;
	}

	public double getLeftBound(){
		return left_bound;
	}

	public double getRightBound(){
		return right_bound;
	}

}
