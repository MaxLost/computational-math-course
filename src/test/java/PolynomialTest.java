import org.junit.Test;
import org.math.computational.Polynomial;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

public class PolynomialTest {

	@Test
	public void multiplicationTest() {
		Polynomial a = new Polynomial(Arrays.asList(-1.0, 1.0));
		Polynomial b = new Polynomial(Arrays.asList(1.0, 1.0));
		Polynomial result = a.multiply(b);
		Polynomial expected = new Polynomial(Arrays.asList(-1.0, 0.0, 1.0));
		assertEquals(expected, result);
	}

	@Test
	public void scalarMultiplicationTest() {
		Polynomial a = new Polynomial(Arrays.asList(-1.0, 1.0));
		Polynomial result = a.scalarMultiply(-1);
		Polynomial expected = new Polynomial(Arrays.asList(1.0, -1.0));
		assertEquals(expected, result);
	}

}
