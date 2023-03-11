package org.math.computational.functions;

import java.util.ArrayList;
import java.util.NoSuchElementException;
import java.util.Scanner;

import org.math.computational.Segment;

public class FunctionRootFinder {

	private final Function f;
	private final Function df;
	private final Function ddf;
	private final double lowerBound;
	private final double upperBound;
	private final double epsilon;

	public FunctionRootFinder(Function f, double lowerBound, double upperBound, double epsilon) {
		this.f = f;
		this.df = null;
		this.ddf = null;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.epsilon = epsilon;
	}

	public FunctionRootFinder(Function f, Function df, double lowerBound, double upperBound, double epsilon) {
		this.f = f;
		this.df = df;
		this.ddf = null;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.epsilon = epsilon;
	}

	public FunctionRootFinder(Function f, Function df, Function ddf,
	                          double lowerBound, double upperBound, double epsilon) {
		this.f = f;
		this.df = df;
		this.ddf = ddf;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.epsilon = epsilon;
	}

	public ArrayList<Double> findRootsLogged(double step_size){

		ArrayList<Segment> segments = getRootSegments(step_size);
		ArrayList<Double> roots = new ArrayList<>();

		System.out.println("После отделения корней получено " + segments.size() + " отрезков с возможными корнями.");
		for (int i = 0; i < segments.size(); i++){
			System.out.println(segments.get(i));
		}

		for (Segment segment : segments) {
			System.out.println("\n############################################\n");
			System.out.println(segment);
			System.out.println("\nМетод бисекции:");
			System.out.println("\tПриближённый корень x = " + approximateRootBisect(segment));
			System.out.println("\nМетод Ньютона:");
			System.out.println("\tПриближённый корень x = " + approximateRootNewton(segment));
			System.out.println("\nУлучшенный метод Ньютона:");
			double root = approximateRootNewtonImproved(segment);
			roots.add(root);
			System.out.println("\tПриближённый корень x = " + root);
			System.out.println("\nМетод секущих:");
			System.out.println("\tПриближённый корень x = " + approximateRootSecantMethod(segment));
		}

		return roots;
	}

	public ArrayList<Double> findRootsLogged(){

		double step_size = getStepSize();

		return findRootsLogged(step_size);
	}

	public ArrayList<Double> findRoots(double step_size){

		ArrayList<Segment> segments = getRootSegments(step_size);
		ArrayList<Double> roots = new ArrayList<>();

		for (Segment segment : segments) {
			roots.add(approximateRootNewtonImproved(segment));
		}

		return roots;
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

			double startValue = this.f.evaluate(start);
			double endValue = this.f.evaluate(end);

			if (startValue * endValue <= 0) {

				if (startValue == 0 &&
						this.f.evaluate(start - this.epsilon) * this.f.evaluate(start + this.epsilon) < 0) {
					rootSegments.add(new Segment(start - this.epsilon, start + this.epsilon));
					end = start + this.epsilon;
				} else if (endValue == 0 &&
						this.f.evaluate(end - this.epsilon) * this.f.evaluate(end + this.epsilon) < 0) {
					rootSegments.add(new Segment(end - this.epsilon, end + this.epsilon));
					end = end + this.epsilon;
				} else {
					rootSegments.add(new Segment(start, end));
				}
			}

			start = end;
			end += stepSize;

		}

		if (this.f.evaluate(start) * this.f.evaluate(this.upperBound) < 0) {
			rootSegments.add(new Segment(start, upperBound));
		}

		return rootSegments;
	}

	private double getStepSize() {

		try (Scanner input = new Scanner(System.in)) {
			System.out.println("Введите количество отрезков на которые будет разделён изначальный отрезок для поиска " +
					"корней: ");
			int N = input.nextInt();

			if (N < 2) {
				System.out.println("Недостаточно отрезков для корректного поиска корней, попробуйте ввести большее " +
						"число");
				return getStepSize();
			}

			double stepSize = (this.upperBound - this.lowerBound) / N;

			if (stepSize > 0.01) {
				System.out.println("Внимание! Отрезки могут оказаться слишком большими для корректного поиска корней");
				System.out.println("Хотите делить на большее число отрезков? [yes/no]: ");
				try {
					String clearance = input.next();
					if (clearance.equals("yes")) {
						return getStepSize();
					} else if (clearance.equals("no")) {
						return stepSize;
					} else {
						throw new NoSuchElementException();
					}
				} catch (NoSuchElementException e) {
					System.out.println("Неизвестный ответ, попробуйте снова");
					return getStepSize();
				}
			}
		}
		return 1e-4;
	}

	private double getStepSizeEN() {

		try (Scanner input = new Scanner(System.in)) {
			System.out.println("Enter number of segments into which original interval will be divided for " +
					"roots search: ");
			int N = input.nextInt();

			if (N < 2) {
				System.out.println("There are not enough segments for correct roots search, try bigger amount");
				return getStepSize();
			}

			double stepSize = (this.upperBound - this.lowerBound) / N;

			if (stepSize > 0.01) {
				System.out.println("Warning: segments may be too big for correct roots search.");
				System.out.println("Would you like to choose bigger amount of segments? [yes/no]: ");
				try {
					String clearance = input.next();
					if (clearance.equals("yes")) {
						return getStepSize();
					} else if (clearance.equals("no")) {
						return stepSize;
					} else {
						throw new NoSuchElementException();
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

		double start = segment.getLowerBound();
		double end = segment.getUpperBound();
		double startValue = this.f.evaluate(start);
		int iterationCounter = 0;

		System.out.println("\tНачальное приближение: " + (start + (end - start) / 2));

		while (end - start > this.epsilon) {

			double center = start + (end - start) / 2;
			double center_value = this.f.evaluate(center);

			if (center_value * startValue <= 0) {

				end = center;
			} else {

				start = center;
				startValue = center_value;
			}

			iterationCounter++;
		}
		System.out.println("\tКоличество итераций: " + iterationCounter);
		System.out.println("\tРазница между последними приближениями: " + (end - start));
		System.out.println("\tАбсолютное значение невязки: " + Math.abs(this.f.evaluate(start + (end - start) / 2)));
		return start + (end - start) / 2;
	}

	private boolean isNewtonAlgorithmConverges(double point) {

		double value = this.f.evaluate(point);
		double dfValue;
		double ddfValue;

		if (this.df == null) {
			dfValue = Utils.computeDerivative(point, this.f);
			ddfValue = Utils.computeDerivative(point, t -> Utils.computeDerivative(t, this.f));
		}
		else if (this.ddf == null) {
			dfValue = this.df.evaluate(point);
			ddfValue = Utils.computeDerivative(point, this.df);
		}
		else {
			dfValue = this.df.evaluate(point);
			ddfValue = this.ddf.evaluate(point);
		}

		return (Math.abs(dfValue) > this.epsilon) && (value * ddfValue >= -this.epsilon);
	}

	public double approximateRootNewton(Segment segment) {

		double x0 = segment.getRandomPoint();
		while(!isNewtonAlgorithmConverges(x0)) {
			x0 = segment.getRandomPoint();
		}

		System.out.println("\tНачальное приближение: " + x0);

		int p = 1;

		while (p < 10) {

			int iterationCounter = 1;

			double value = this.f.evaluate(x0);
			double dfValue;
			if (this.df == null) {
				dfValue = Utils.computeDerivative(x0, this.f);
			} else {
				dfValue = this.df.evaluate(x0);
			}
			double x = x0 - p * value / dfValue;

			while(Math.abs(x - x0) > this.epsilon && iterationCounter < 30) {

				x0 = x;
				if (this.df == null) {
					dfValue = Utils.computeDerivative(x0, this.f);
				} else {
					dfValue = this.df.evaluate(x0);
				}
				value = this.f.evaluate(x0);

				x = x0 - p * value / dfValue;
				iterationCounter++;
			}
			if (Math.abs(x - x0) < this.epsilon) {
				System.out.println("\tКоличество итераций: " + iterationCounter);
				System.out.println("\tРазница между последними приближениями: " + Math.abs(x - x0));
				System.out.println("\tАбсолютное значение невязки: " + Math.abs(this.f.evaluate(x)));
				return x;
			}
			p += 2;
		}

		return Double.NaN;
	}

	public double approximateRootNewtonImproved(Segment segment) {

		double x0 = segment.getRandomPoint();
		while(!isNewtonAlgorithmConverges(x0)) {
			x0 = segment.getRandomPoint();
		}

		System.out.println("\tНачальное приближение: " + x0);

		int p = 1;

		while (p < 10) {

			int iterationCounter = 1;
			double dfValue;

			if (this.df == null) {
				dfValue = Utils.computeDerivative(x0, this.f);
			} else {
				dfValue = this.df.evaluate(x0);
			}

			double value = this.f.evaluate(x0);
			double x = x0 - p * value / dfValue;

			while(Math.abs(x - x0) > this.epsilon && iterationCounter < 30) {
				x0 = x;
				value = this.f.evaluate(x0);

				x = x0 - p * value / dfValue;
				iterationCounter++;
			}
			if (Math.abs(x - x0) < this.epsilon) {
				System.out.println("\tКоличество итераций: " + iterationCounter);
				System.out.println("\tРазница между последними приближениями: " + Math.abs(x - x0));
				System.out.println("\tАбсолютное значение невязки: " + Math.abs(this.f.evaluate(x)));
				return x;
			}
			p += 2;
		}

		return Double.NaN;
	}

	public double approximateRootSecantMethod(Segment segment) {

		int p = 1;

		while (p < 10) {

			int iterationCounter = 0;

			double x0 = segment.getLowerBound();
			double x1 = segment.getUpperBound();
			double x1Value = this.f.evaluate(x1);

			System.out.println("\tНачальные приближения: x0 = " + x0 + ", x1 = " + x1);

			double x2 = x1 - p * x1Value * (x1 - x0) / (x1Value - this.f.evaluate(x0));

			while(Math.abs(x2 - x1) > this.epsilon && iterationCounter < 30) {

				x0 = x1;
				x1 = x2;
				x1Value = this.f.evaluate(x1);

				x2 = x1 - p * x1Value * (x1 - x0) / (x1Value - this.f.evaluate(x0));

				iterationCounter++;
			}
			if (Math.abs(x2 - x1) < this.epsilon) {
				System.out.println("\tКоличество итераций: " + iterationCounter);
				System.out.println("\tРазница между последними приближениями: " + Math.abs(x2 - x1));
				System.out.println("\tАбсолютное значение невязки: " + Math.abs(this.f.evaluate(x2)));
				return x2;
			}
			p += 2;
		}

		return Double.NaN;
	}

}
