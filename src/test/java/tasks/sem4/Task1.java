package tasks.sem4;

import org.math.computations.functions.Function;
import org.math.computations.functions.FunctionRootFinder;

public class Task1 {

	public static void main(String[] args) {

		System.out.println("""
				Решение нелинейного уравнения
				Найти все корни нечётной кратности уравнения: f(x) = 0
				
				Вариант 15
				
				f(x) = (x-1)^2 - exp(-x)
				A = -1, B = 3
				
				Точность = 10^-8""");
		Function f = t -> Math.pow(t - 1, 2) - Math.exp(-t);
		Function df = t -> 2*(t - 1) + Math.exp(-t);
		Function ddf = t -> 2 - Math.exp(-t);
		FunctionRootFinder task = new FunctionRootFinder(f, df, ddf, -1, 3, 10e-8);
		task.findRootsLogged();
	}
}
