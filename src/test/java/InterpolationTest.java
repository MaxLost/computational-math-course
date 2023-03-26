import org.junit.Test;
import org.math.computational.PlanePoint;
import org.math.computational.Polynomial;
import org.math.computational.functions.AlgebraicInterpolator;
import org.math.computational.functions.Function;

import java.util.List;

import static org.junit.Assert.assertTrue;
import static tasks.Utils.divideSegmentEqual;

public class InterpolationTest {

	@Test
	public void getInterpolationPolynomialTest() {

		Function f = t -> Math.sin(t);

		List<PlanePoint> table = divideSegmentEqual(f, -1, 1, 10);

		AlgebraicInterpolator task = new AlgebraicInterpolator(table);

		Polynomial p = task.getInterpolationPolynomial(9, 0.5);
		System.out.println(p.evaluate(0.5));
		assertTrue(Math.abs(p.evaluate(0.5) - 0.4794255386042) < 10e-8);
	}



}
