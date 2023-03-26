package tasks;

import org.math.computational.PlanePoint;
import org.math.computational.Segment;
import org.math.computational.functions.Function;

import java.util.ArrayList;
import java.util.List;

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

}
