package tasks;

import org.math.computational.functions.Function;
import org.math.computational.functions.Integrator;

import static tasks.Utils.dfMax;
import static tasks.Utils.ddfMax;

import java.util.*;

public class Task4_2 {

	public static void integrationTest(Function f, Function F, Function ddf, double A, double B, int m) {

		double exactValue = F.evaluate(B) - F.evaluate(A);
		System.out.print("\nТочное значение интеграла: " + exactValue);

		Integrator task = new Integrator(f, A, B);

		double dfMax = dfMax(f, A, B);
		double ddfMax = ddfMax(f, A, B);
		double ddddfMax = ddfMax(ddf, A, B);

		System.out.printf(Locale.US, "\nЗначение h: %5.3e", (B - A) / m);

		double value = task.integrate("LR", m);
		double error = Math.pow(B - A, 2) * dfMax / (2 * m);
		System.out.printf(Locale.US, """
						\n\nФормула левых прямоугольников:
						\tЗначение интеграла: %14.8f
						\tТеоретическая погрешность: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, error, Math.abs(value - exactValue));

		value = task.integrate("RR", m);
		System.out.printf(Locale.US, """
						\n\nФормула правых прямоугольников:
						\tЗначение интеграла: %14.8f
						\tТеоретическая погрешность: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, error, Math.abs(value - exactValue));

		value = task.integrate("СR", m);
		error = Math.pow(B - A, 3) * ddfMax / (24 * Math.pow(m, 2));
		System.out.printf(Locale.US, """
						\n\nФормула серединных прямоугольников:
						\tЗначение интеграла: %14.8f
						\tТеоретическая погрешность: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, error, Math.abs(value - exactValue));

		value = task.integrate("TR", m);
		error = (Math.pow(B - A, 3) * ddfMax) / (12 * Math.pow(m, 2));
		System.out.printf(Locale.US, """
						\n\nФормула трапеций:
						\tЗначение интеграла: %14.8f
						\tТеоретическая погрешность: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, error, Math.abs(value - exactValue));

		value = task.integrate("SI", m);
		error = (Math.pow(B - A, 5) * ddddfMax) / (2880 * Math.pow(m, 4));
		System.out.printf(Locale.US, """
						\n\nФормула Симпсона:
						\tЗначение интеграла: %14.8f
						\tТеоретическая погрешность: %10.3e
						\tАбсолютная погрешность: %10.3e""",
				value, error, Math.abs(value - exactValue));
	}

	public static void main(String[] args) {

		Function F = t -> Math.sin(t * t) + Math.pow(t, 4) / 4;
		Function f = t -> 2 * t * Math.cos(t * t) + Math.pow(t, 3);
		Function ddf = t -> -12*t * Math.sin(t * t) - 8*Math.pow(t, 3) * Math.cos(t * t) + 6*t;

		Function ddc = t -> 0;
		Function c = t -> 1;
		Function C = t -> t;

		Function ddx = t -> 0;
		Function x = t -> t;
		Function X = t -> Math.pow(t, 2) / 2;

		Function ddx3 = t -> 6*t;
		Function x3 = t -> Math.pow(t, 3);
		Function X3 = t -> Math.pow(t, 4) / 4;

		System.out.println("""
				Задача приближённого вычисления интеграла с помощью составных интерполяционных квадратурных формул
				                
				f(x) = 2t * cos(t^2) + t^3
				Найти интеграл от f(x) по [a; b]
				                
				Введите A - левую границу промежутка интегрирования:""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)) {

			double A = input.nextDouble();

			System.out.println("Введите B - правую границу промежутка интегрирования:");
			double B = input.nextDouble();

			int m = 0;
			while (m < 1) {
				System.out.println("Введите M >= 1 - количество отрезков на которые хотите разделить изначальный");
				m = input.nextInt();
				if (m < 1) {
					System.out.println("M должно быть натуральным числом, пропробуйте ещё раз.");
				}
			}
			integrationTest(f, F, ddf, A, B, m);
			//integrationTest(c, C, ddc, A, B, m);
			//integrationTest(x, X, ddx, A, B, m);
			//integrationTest(x3, X3, ddx3, A, B, m);
		}

	}
}
