package tasks;

import org.math.computations.PlanePoint;
import org.math.computations.functions.Function;
import org.math.computations.functions.Integrator;

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

			System.out.println("\n\nКФ Гаусса\n\n");

			Function f = t -> Math.sin(t) / t;
			double exactValueF = 1.605412976802694; // Exact value of integral of f(x) from 0 to 2

			String repeat = "yes";

			while (!repeat.equals("no")) {

				int N = -1;
				while (N < 1) {
					System.out.println("Введите N - число узлов для построения формулы Гаусса: ");
					N = input.nextInt();
					if (N < 1) {
						System.out.println("Слишком малое N, попробуйте ввести большее N");
					}
				}

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

					List<PlanePoint> Af = integratorF.buildGaussianQuadrature(i);
					List<PlanePoint> d = integratorF.buildDefaultGaussianQuadrature(i);
					List<PlanePoint> Ax = integratorX.buildGaussianQuadrature(i);

					System.out.println("\nУзлы и коэффициенты изначальной КФ Гаусса\n k |     x_k     |     A_k    ");
					for (int j = 0; j < i; j++) {
						System.out.printf(Locale.US, " %d | %11.6f | %11.6f\n", j, d.get(j).getX(), d.get(j).getY());
					}

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
							"\nЗначение интеграла при f(x) = x^%d: %.13f\nАбсолютная погрешность: %.5e\n" +
									"Относительная погрешность: %.4f%%\n",
							2*i - 1, pValue, Math.abs(pValue - exactValueP),
							Math.abs(pValue - exactValueP) / pValue * 100);
					System.out.printf(Locale.US,
							"\nЗначение интеграла при f(x) = sin(x) / x : %.13f\nАбсолютная погрешность: %.5e\n" +
									"Относительная погрешность: %.4f%%\n",
							fValue, Math.abs(fValue - exactValueF),
							Math.abs(fValue - exactValueF) / fValue * 100);
				}

				input.nextLine();
				System.out.println("Хотите ли вы выбрать другие границы интегрирования [yes/no]?");
				repeat = input.nextLine();
			}

			System.out.println("\n\nКФ Мёлера\n\n");

			int[] N1 = {-1, -1, -1};
			for (int i = 0; i < 3; i++) {
				while (N1[i] < 2) {
					System.out.println("Введите N" + (i + 1) + " - число узлов для построения формулы Мёлера: ");
					N1[i] = input.nextInt();
					if (N1[i] < 1) {
						System.out.println("Слишком малое N, попробуйте ввести большее N");
					}
				}
			}

			Function g = t -> Math.cos(t);
			Function h = t -> g.evaluate(t) / Math.sqrt(1 - Math.pow(t, 2));
			double exactValueH = 2.403939430634; // Exact value of integral of g(x) / sqrt(1-x^2) from -1 to 1
			repeat = "yes";
			while(!repeat.equals("no")) {

				System.out.println("Введите A - левую границу промежутка");
				double A = input.nextDouble();
				System.out.println("Введите B - правую границу промежутка");
				double B = input.nextDouble();

				for (int i = 0; i < 3; i++) {
					int finalN1 = N1[i];
					Function x = t -> Math.pow(t, 2 * finalN1 - 1) / Math.sqrt(1 - Math.pow(t, 2));

					Integrator integratorH = new Integrator(h, A, B);
					Integrator integratorX = new Integrator(x, A, B);

					List<PlanePoint> Cg = integratorH.buildMohlerQuadrature(finalN1);
					List<PlanePoint> Cx = integratorX.buildMohlerQuadrature(finalN1);

					System.out.println("\nУзлы и коэффициенты КФ Мёлера при N = " + finalN1 + "\n k |     x_k     |     A_k    ");
					for (int j = 0; j < finalN1; j++) {
						System.out.printf(Locale.US, " %d | %11.6f | %11.6f\n", j, Cg.get(j).getX(), Cg.get(j).getY());
					}

					double xValue = 0;
					double hValue = 0;

					for (int j = 0; j < finalN1; j++) {
						xValue += Cx.get(j).getY() * x.evaluate(Cx.get(j).getX());
						hValue += Cg.get(j).getY() * g.evaluate(Cg.get(j).getX());
					}

					System.out.printf(Locale.US,
							"\nЗначение интеграла при f(x) = x^%d: %.13f\nАбсолютная погрешность: %.8f\n" +
									"Относительная погрешность: %.4f%%\n",
							2*finalN1 - 1, xValue, Math.abs(xValue),
							Math.abs(xValue) / xValue * 100);
					System.out.printf(Locale.US,
							"\nЗначение интеграла при f(x) = cos(x): %.13f\nАбсолютная погрешность: %.8f\n" +
									"Относительная погрешность: %.4f%%\n",
							hValue, Math.abs(hValue - exactValueH),
							Math.abs(hValue - exactValueH) / hValue * 100);
				}

				input.nextLine();
				System.out.println("Хотите ли вы выбрать другие границы интегрирования [yes/no]?");
				repeat = input.nextLine();

			}
		}
	}

}
