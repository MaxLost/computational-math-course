package tasks;

import org.math.computational.functions.Function;
import org.math.computational.functions.FunctionRootFinder;

public class Task1 {

	public static void main(String[] args) {

		System.out.println("""
				Non-linear equation solving
				Find all roots of equation: f(x) = 0 with odd order
				A = -1, B = 3
				f(x) = (x-1)^2 - exp(-x)
				Required precision = 10^-8""");
		Function f = t -> Math.pow(t - 1, 2) - Math.exp(-t);
		Function df = t -> 2*(t - 1) + Math.exp(-t);
		Function ddf = t -> 2 - Math.exp(-t);
		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -3, 3, 10e-8);
		task.findRootsLogged();
	}
}
