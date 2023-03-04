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

	@Override
	public String toString(){
		return "[" + lowerBound + ";" + upperBound + "]";
	}

	@Override
	public boolean equals(Object o) {

		if (o instanceof Segment){

			if (o == this) {
				return true;
			}
			return Math.abs(this.lowerBound - ((Segment) o).getLowerBound()) < 10e-8 &&
					Math.abs(this.upperBound - ((Segment) o).getUpperBound()) < 10e-8;
		}
		return false;
	}

	@Override
	public int hashCode(){
		return (int) (Math.round(Math.ceil(upperBound)) ^ Math.round(Math.floor(lowerBound)));
	}

}
