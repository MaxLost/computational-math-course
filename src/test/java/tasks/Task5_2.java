package tasks;

import org.math.computational.PlanePoint;
import org.math.computational.functions.Function;
import org.math.computational.functions.Integrator;

import java.util.List;
import java.util.Locale;
import java.util.Scanner;

public class Task5_2 {

	public static void main(String[] args) {
		System.out.println("""
                Вычисление значения интеграла от f(x) с помощью формул Гаусса и Мёлера

                1 вариант
                
                Границы интегрирования: [0, 2]
                Значения N: 2, 5, 6, 8
                
                Для формулы Гаусса:
                f(x) = sin(x) / x
                
                Для формулы Мёлера:
                f(x) = cos(x)""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)) {

			int N = -1;
			while (N < 1) {
				System.out.println("Введите N - число узлов для построения формул Гаусса и Мёлера: ");
				N = input.nextInt();
				if (N < 1) {
					System.out.println("Слишком малое N, попробуйте ввести большее N");
				}
			}

			Function f = t -> Math.sin(t) / t;
			double exactValueF = 1.605412976802694; // Exact value of integral of f(x) from 0 to 2

			String repeat = "yes";

			while (!repeat.equals("no")) {
				System.out.println("Введите A - левую границу промежутка");
				double A = input.nextDouble();
				System.out.println("Введите B - правую границу промежутка");
				double B = input.nextDouble();

				for (int i = 1; i <= N; i++) {
					int finalN = i;
					Function p = t -> Math.pow(t, 2 * finalN - 1);
					Function P = t -> Math.pow(t, 2 * finalN) / (2 * finalN);
					double exactValueP = P.evaluate(B) - P.evaluate(A);

					Integrator integratorF = new Integrator(f, A, B);
					Integrator integratorX = new Integrator(p, A, B);

					List<PlanePoint> Af = integratorF.buildGaussQuadrature(i);
					List<PlanePoint> Ax = integratorX.buildGaussQuadrature(i);

					System.out.println("\nУзлы и коэффициенты КФ Гаусса при N = " + i + "\n k |     x_k     |     A_k    ");
					for (int j = 0; j < i; j++) {
						System.out.printf(Locale.US, " %d | %11.6f | %11.6f\n", j, Af.get(j).getX(), Af.get(j).getY());
					}

					double pValue = 0;
					double fValue = 0;

					for (int j = 0; j < i; j++) {
						pValue += Ax.get(j).getY() * p.evaluate(Ax.get(j).getX());
						fValue += Af.get(j).getY() * f.evaluate(Af.get(j).getX());
					}

					System.out.printf(Locale.US,
							"\nx^%d : %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
							2*i - 1, pValue, Math.abs(pValue - exactValueP),
							Math.abs(pValue - exactValueP) / pValue * 100);
					System.out.printf(Locale.US,
							"\nf(x) : %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
							fValue, Math.abs(fValue - exactValueF),
							Math.abs(fValue - exactValueF) / fValue * 100);
				}

				input.nextLine();
				System.out.println("Хотите ли вы выбрать другие границы интегрирования [yes/no]?");
				repeat = input.nextLine();
			}
		}
	}

}