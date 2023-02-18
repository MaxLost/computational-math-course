import org.math.computational.functions.*;

import org.junit.Test;

public class FunctionTest {

	@Test
	public void test1A(){

		System.out.println("Non-linear equation solving\nA = -8, B = 10\nf(x) = 2^x-2*cos(x)\nPrecision = 10^-6");
		Function f = t -> Math.pow(2, t) + 2*Math.cos(t);
		Function df = t -> Math.pow(2, t)*Math.log(2) + 2*Math.sin(t);
		Function ddf = t -> Math.pow(2, t)*Math.pow(Math.log(2), 2) + 2*Math.cos(t);

		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -8, 10, 10e-6);
		double stepSize = (10.0 + 8.0) / 10000.0; // Let's divide segment into 10000 subsegments
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test(){

		double x = 2;

		Function eq = t -> Math.pow(t, 2) + t;
		System.out.println(eq.evaluate(x));
		System.out.println(10e-6);
	}

}
