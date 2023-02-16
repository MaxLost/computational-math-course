package org.math.computational.functions;

import java.util.Random;

public class Segment {

	private final double lowerBound;
	private final double upperBound;

	public Segment(double a, double b){

		lowerBound = a;
		upperBound = b;
	}

	public double getLeftBound(){
		return lowerBound;
	}

	public double getRightBound(){
		return upperBound;
	}

	public double getRandomPoint() {

		return this.lowerBound + (new Random().nextDouble()) *
				(this.upperBound - this.lowerBound);

	}

}
