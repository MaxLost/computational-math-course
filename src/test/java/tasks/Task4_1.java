package tasks;

import org.math.computations.functions.Function;
import org.math.computations.functions.Integrator;

import java.util.Locale;
import java.util.Scanner;

public class Task4_1 {

	public static void integrationTest(Function f, Function F, double A, double B) {

		double exactValue = F.evaluate(B) - F.evaluate(A);
		System.out.print("\nТочное значение интеграла: " + exactValue);

		Integrator task = new Integrator(f, A, B);

		double value = task.integrate("LR");
		System.out.printf(Locale.US, """
						\n\nФормула левых прямоугольников:
						\tЗначение интеграла: %10.3e
						\tАбсолютная погрешность: %10.3e""",
							value, Math.abs(value - exactValue));

		value = task.integrate("RR");
		System.out.printf(Locale.US, """
						\n\nФормула правых прямоугольников:
						\tЗначение интеграла: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, Math.abs(value - exactValue));

		value = task.integrate("СR");
		System.out.printf(Locale.US, """
						\n\nФормула серединных прямоугольников:
						\tЗначение интеграла: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, Math.abs(value - exactValue));

		value = task.integrate("TR");
		System.out.printf(Locale.US, """
						\n\nФормула трапеций:
						\tЗначение интеграла: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, Math.abs(value - exactValue));

		value = task.integrate("SI");
		System.out.printf(Locale.US, """
						\n\nФормула Симпсона:
						\tЗначение интеграла: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, Math.abs(value - exactValue));

		value = task.integrate("");
		System.out.printf(Locale.US, """
						\n\nФормула 3/8:
						\tЗначение интеграла: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, Math.abs(value - exactValue));

	}

	public static void main(String[] args) {

		Function F = t -> Math.sin(t * t) + Math.pow(t, 4) / 4;
		Function f = t -> 2 * t * Math.cos(t * t) + Math.pow(t, 3);

		Function c = t -> 1;
		Function C = t -> t;

		Function x = t -> t;
		Function X = t -> Math.pow(t, 2) / 2;

		Function x3 = t -> Math.pow(t, 3);
		Function X3 = t -> Math.pow(t, 4) / 4;

		System.out.println("""
				Задача приближённого вычисления интеграла с помощью интерполяционных квадратурных формул
				                
				f(x) = 2t * cos(t^2) + t^3
				Найти интеграл от f(x) по [a; b]
				                
				Введите A - левую границу промежутка интегрирования:""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)) {

			double A = input.nextDouble();

			System.out.println("Введите B - правую границу промежутка интегрирования:");
			double B = input.nextDouble();

			integrationTest(f, F, A, B);
			//integrationTest(c, C, A, B);
			//integrationTest(x, X, A, B);
			//integrationTest(x3, X3, A, B);
		}

	}
}
