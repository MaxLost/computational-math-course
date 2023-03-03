import org.math.computational.functions.*;

import org.junit.Test;

public class FunctionTest {

	@Test
	public void test1A(){

		System.out.println("""
				Non-linear equation solving
				Find all roots of equation: f(x) = 0 with odd order
				A = -8, B = 10
				f(x) = 2^x-2*cos(x)
				Required precision = 10^-6""");
		Function f = t -> Math.pow(2, t) - 2*Math.cos(t);
		Function df = t -> Math.pow(2, t)*Math.log(2) + 2*Math.sin(t);
		Function ddf = t -> Math.pow(2, t)*Math.pow(Math.log(2), 2) + 2*Math.cos(t);

		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -8, 10, 10e-6);
		double stepSize = (10.0 + 8.0) / 1000.0; // Let's divide segment into 1000 subsegments
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test2A(){

		System.out.println("""
				Non-linear equation solving
				Find all roots of equation: f(x) = 0 with odd order
				A = -5, B = 3
				f(x) = x-10*sin(x)
				Required precision = 10^-6""");
		Function f = t -> t - 10*Math.sin(t);
		Function df = t -> 1 - 10*Math.cos(t);
		Function ddf = t -> 10*Math.sin(t);
		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -5, 3, 10e-6);
		double stepSize = (3 + 5) / 1000.0;
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test2B(){

		System.out.println("""
				Non-linear equation solving
				Find all roots of equation: f(x) = 0 with odd order
				A = -5, B = 3
				f(x) = x-10*sin(x)
				Required precision = 10^-6""");
		Function f = t -> t - 10*Math.sin(t);
		FunctionRootFinder task = new FunctionRootFinder(f, -5, 3, 10e-6);
		double stepSize = (3 + 5) / 1000.0;
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test3B(){

		System.out.println("""
				Non-linear equation solving
				Find all roots of equation: f(x) = 0 with odd order
				A = -5, B = 5
				f(x) = (1,2)*x^4+2*x^3‒13*x^2‒14,2*x‒24,1
				Required precision = 10^-6""");
		Function f = t -> 1.2*Math.pow(t, 4) + 2*Math.pow(t, 3) - 13*Math.pow(t, 2) - 14.2*t - 24.1;
		FunctionRootFinder task = new FunctionRootFinder(f, -5, 5, 10e-6);
		double stepSize = (5 + 5) / 1000.0;
		task.findRootsLogged(stepSize);
	}

	@Test
	public void test4B(){

		System.out.println("""
				Non-linear equation solving
				Find all roots of equation: f(x) = 0 with odd order
				A = -3, B = 3
				f(x) = x^2
				Required precision = 10^-8""");
		Function f = t -> Math.pow(t, 2);
		FunctionRootFinder task = new FunctionRootFinder(f, -3, 3, 10e-8);
		double stepSize = (3 + 3) / 1000.0;
		task.findRootsLogged(stepSize);
	}
}
