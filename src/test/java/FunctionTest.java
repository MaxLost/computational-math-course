import org.math.computational.functions.*;

import org.junit.Test;

public class FunctionTest {

	@Test
	public void variant11(){

		Function f = t -> Math.tan(t)+Math.exp(t);
		Function df = t -> 1/Math.pow(Math.cos(t), 2) + Math.exp(t);
		Function ddf = t -> 2*Math.tan(t)/Math.pow(Math.cos(t), 2) + Math.exp(t);
		double epsilon = 1e-6;

		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -50, 0, epsilon);
	}

	@Test
	public void test(){

		double x = 2;

		Function eq = t -> Math.pow(t, 2) + t;
		System.out.println(eq.evaluate(x));
		System.out.println(1e-3);
	}

}
