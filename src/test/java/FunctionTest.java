import org.math.computational.functions.*;

import org.junit.Test;

public class FunctionTest {

	@Test
	public void test1A(){

		System.out.println("Non-linear equation solving\nA = -8, B = 10\nf(x) = 2^x-2*cos(x)\nPrecision = 10^-6");
		Function f = t -> Math.pow(2, t) - 2*Math.cos(t);
		Function df = t -> Math.pow(2, t)*Math.log(2) + 2*Math.sin(t);
		Function ddf = t -> Math.pow(2, t)*Math.pow(Math.log(2), 2) + 2*Math.cos(t);

		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -8, 10, 10e-6);
		double stepSize = (10.0 + 8.0) / 1000.0; // Let's divide segment into 1000 subsegments
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test2A(){

		System.out.println("Non-linear equation solving\nA = -5, B = 3\nf(x) = x-10*sin(x)\nPrecision = 10^-6");
		Function f = t -> t - 10*Math.sin(t);
		Function df = t -> 1 - 10*Math.cos(t);
		Function ddf = t -> 10*Math.sin(t);
		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -5, 3, 10e-6);
		double stepSize = (3 + 5) / 1000.0;
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test2B(){

		System.out.println("Non-linear equation solving\nA = -5, B = 3\nf(x) = x-10*sin(x)\nPrecision = 10^-6");
		Function f = t -> t - 10*Math.sin(t);
		FunctionRootFinder task = new FunctionRootFinder(f, -5, 3, 10e-6);
		double stepSize = (3 + 5) / 1000.0;
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
