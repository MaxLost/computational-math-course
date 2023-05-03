package tasks;

import org.math.computational.functions.Function;
import org.math.computational.functions.Integrator;

import java.util.Locale;
import java.util.Scanner;

public class Task5_3 {

	public static void main(String[] args) {

		System.out.println("""
                Вычисление значения интеграла от f(x) с помощью формул составной формулы Гаусса

                15 вариант
                
                Границы интегрирования: [0, 1]
                
                f(x) = sin(x)
                ρ(x) = sqrt(x / (1 - x))""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)) {

			int N = -1;
			while (N < 1) {
				System.out.println("\n\nВведите N - число узлов для построения формулы Гаусса: ");
				N = input.nextInt();
				if (N < 1) {
					System.out.println("Слишком малое N, попробуйте ввести большее N");
				}
			}

			System.out.println("Введите A - левую границу промежутка");
			double A = input.nextDouble();
			System.out.println("Введите B - правую границу промежутка");
			double B = input.nextDouble();

			int M = -1;
			while (M < 2) {
				System.out.printf(Locale.US,
						"Введите M - число разбиений для изначального отрезка [%.6f, %.6f]:\n", A, B);
				M = input.nextInt();
				if (M < 2) {
					System.out.println("Слишком малое M, попробуйте ввести большее M");
				}
			}

			Function f = t -> Math.sin(t) / Math.sqrt(t / (1 - t));
			double exactValueF = 0.37277488169851963;

			Integrator integrator = new Integrator(f, A, B);
			double fValue = integrator.integrateCompositeGaussianQuadrature(N, M);

			System.out.printf(Locale.US,
					"\nЗначение интеграла: %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
					fValue, Math.abs(fValue - exactValueF),
					Math.abs(fValue - exactValueF) / fValue * 100);

		}
	}

}
