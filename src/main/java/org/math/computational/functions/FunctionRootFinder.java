package org.math.computational.functions;

import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.Scanner;

public class FunctionRootFinder {

	private final Function f;
	private final double lowerBound;
	private final double upperBound;
	private final double epsilon;

	public FunctionRootFinder(Function f, double lowerBound, double upperBound, double epsilon) {
		this.f = f;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.epsilon = epsilon;
	}

	public ArrayList<Double> findRoots(double step_size){

		ArrayList<Segment> segments = getRootSegments(step_size);

		for (Segment segment : segments) {
			// Call methods for root approximation
		}

		return null;
	}

	public ArrayList<Double> findRoots(){

		double step_size = getStepSize();

		return findRoots(step_size);
	}

	public ArrayList<Segment> getRootSegments(double stepSize) {

		ArrayList<Segment> rootSegments = new ArrayList<>();

		double start = this.lowerBound;
		double end = start + stepSize;
		while (end <= this.upperBound) {

			double start_value = this.f.evaluate(start);
			double end_value = this.f.evaluate(end);

			if (start_value * end_value < 0) {
				rootSegments.add(new Segment(start, end));
			}

			start = end;
			end += stepSize;

		}

		return rootSegments;
	}

	private double getStepSize() {

		try (Scanner input = new Scanner(System.in)) {
			System.out.println("Enter number of intervals into which original interval will be divided for " +
					"roots search: ");
			int N = input.nextInt();

			if (N < 2) {
				System.out.println("There are not enough intervals for correct roots search, try bigger amount");
				return getStepSize();
			}

			double stepSize = (this.upperBound - this.lowerBound) / N;

			if (stepSize > 0.01) {
				System.out.println("Warning: intervals may be too big for correct roots search.");
				System.out.println("Would you like to choose bigger amount of intervals? [yes/no]: ");
				try {
					String clearance = input.next("[a-z]");
					if (clearance.equals("yes")) {
						return getStepSize();
					} else if (clearance.equals("no")) {
						return stepSize;
					}
				} catch (NoSuchElementException e) {
					System.out.println("Unexpected answer, try again");
					return getStepSize();
				}
			}
		}
		return 1e-4;
	}

	public double approximateRootBisect(Segment segment) {

		double left_bound = segment.getLeftBound();
		double right_bound = segment.getRightBound();
		double left_value = this.f.evaluate(left_bound);

		while (right_bound - left_bound > this.epsilon) {
			double center = (right_bound - left_bound) / 2;
			double center_value = this.f.evaluate(center);

			if (center_value * left_value <= 0) {

				right_bound = center;
			}
			else {

				left_bound = center;
				left_value = center_value;
			}
		}

		return this.f.evaluate((right_bound - left_bound) / 2);
	}

	public double approximateRootNewtonManual(Segment segment, Function firstDerivative, Function secondDerivative) {

		return 0;

	}

}
