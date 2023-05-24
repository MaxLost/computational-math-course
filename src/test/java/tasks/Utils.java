package tasks;

import org.math.computations.PlanePoint;
import org.math.computations.Segment;
import org.math.computations.functions.Function;
import org.math.computations.functions.NumericalDifferentiator;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;

public class Utils {

	public static List<PlanePoint> divideSegmentRandom(Function f, double left, double right, int N) {

		List<PlanePoint> result = new ArrayList<>();

		double step = (right - left) / N;
		double start = left;
		double end = start + step;

		for (int i = 0; i < N; i++) {

			double x = (new Segment(start, end)).getRandomPoint();
			result.add(new PlanePoint(x, f.evaluate(x)));

			start = end;
			end += step;
		}

		return result;
	}

	public static List<PlanePoint> divideSegmentEqual(Function f, double left, double right, int N) {

		List<PlanePoint> result = new ArrayList<>();

		double step = (right - left) / N;
		double start = left;
		double end = start + step;
		result.add(new PlanePoint(start, f.evaluate(start)));

		for (int i = 0; i < N; i++) {
			result.add(new PlanePoint(end, f.evaluate(end)));
			start = end;
			end += step;
		}

		return result;
	}

	public static double dfMax(Function f, double A, double B) {

		List<PlanePoint> nodes = divideSegmentEqual(f, A, B, 10000);
		NumericalDifferentiator task = new NumericalDifferentiator(nodes);
		List<PlanePoint> dfValues = task.computeDerivatives();
		Optional<PlanePoint> value = dfValues.stream().max(Comparator.comparingDouble(t -> Math.abs(t.getY())));
		if (value.isPresent()) {
			return Math.abs(value.get().getY());
		} else {
			return 0;
		}
	}

	public static double ddfMax(Function f, double A, double B) {

		List<PlanePoint> nodes = divideSegmentEqual(f, A, B, 10000);
		NumericalDifferentiator task = new NumericalDifferentiator(nodes);
		List<PlanePoint> ddfValues = task.computeSecondDerivatives();
		Optional<PlanePoint> value = ddfValues.stream().max(Comparator.comparingDouble(t -> Math.abs(t.getY())));
		if (value.isPresent()) {
			return Math.abs(value.get().getY());
		} else {
			return 0;
		}
	}

}
