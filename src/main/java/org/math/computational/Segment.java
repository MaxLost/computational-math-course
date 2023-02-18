package org.math.computational;

import java.util.Random;

public class Segment {

	private final double lowerBound;
	private final double upperBound;

	public Segment(double a, double b){

		lowerBound = a;
		upperBound = b;
	}

	public double getLowerBound(){
		return lowerBound;
	}

	public double getUpperBound(){
		return upperBound;
	}

	public double getRandomPoint() {

		return this.lowerBound + (new Random().nextDouble()) *
				(this.upperBound - this.lowerBound);

	}

	public String toString(){
		return "[" + lowerBound + ";" + upperBound + "]";
	}

}
