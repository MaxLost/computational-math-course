package tasks.sem4;

import org.math.computations.PlanePoint;
import org.math.computations.Segment;
import org.math.computations.functions.Function;
import org.math.computations.functions.Integrator;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Scanner;

public class Task5_1 {

	public static void main(String[] args) {
		System.out.println("""
                Вычисление значения интеграла по ИКФ и КФНАСТ
                
                15 вариант
                f(x) = sin(x)
                ρ(x) = sqrt(x / (1 - x))
                
                Введите A - левую границу промежутка:""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)) {

			Function f = t -> Math.sin(t);
			Function rho = t -> Math.sqrt(t / (1 - t));
			//Function rho = t -> Math.sqrt(t);

			//double exactValueF = 0.36422193203213236407385175;
			double exactValueF = 1.04071063342533; // Value of integral of f(x)*ρ(x) from 0 to 1
			double exactValueX = 1.17809724509617; // Value of integral of x*ρ(x) from 0 to 1

			double A = input.nextDouble();
			System.out.println("Введите B - правую границу промежутка");
			double B = input.nextDouble();

			int N = -1;
			while (N < 1) {
				System.out.println("Введите N - число узлов для построения ИКФ и КФНАСТ: ");
				N = input.nextInt();
				if (N < 1) {
					System.out.println("Слишком малое N, попробуйте ввести большее N");
				}
			}

			int finalN = N;
			Function x = t -> Math.pow(t, finalN - 1);
			Function x3 = t -> Math.pow(t, 2*finalN - 1);

			List<Double> nodes = new ArrayList<>();
			String repeat = "yes";

			input.nextLine();

			Integrator integratorX = new Integrator(x, A, B);
			Integrator integratorX3 = new Integrator(x3, A, B);
			Integrator integratorF = new Integrator(f, A, B);
			List<Double> momentums = integratorF.calculateMomentums(rho, 2*N);

			while (!repeat.equals("no")) {
				String status = "";
				System.out.println("Хотите ли определить узлы интерполирования для ИКФ автоматически [yes/no]?");
				while (!status.equals("yes") && !status.equals("no")) {
					if (!status.equals("")) {
						System.out.println("Хотите ли определить узлы интерполирования для ИКФ автоматически [yes/no]?");
					}
					status = input.nextLine();
					if (status.equals("no")) {
						System.out.println("Введите " + N + " узлов по-одному:");
						for (int i = 0; i < N; i++) {
							nodes.add(input.nextDouble());
						}
					} else if (status.equals("yes")) {
						nodes = (new Segment(A, B)).getUniformlySpacedPoints(N - 1);
					}
				}

				System.out.println("\n\nИКФ:\nМоменты:");
				for(int i = 0; i < N; i++) {
					System.out.printf(Locale.US, "µ_%d = %.6f\n", i, momentums.get(i));
				}

				List<PlanePoint> Ax = integratorX.buildInterpolationQuadrature(N, nodes, momentums);

				List<PlanePoint> Af = integratorF.buildInterpolationQuadrature(N, nodes, momentums);

				System.out.println("\nf(x) = sin(x)\n k |     x_k     |     A_k    ");
				for (int i = 0; i < N; i++) {
					System.out.printf(Locale.US," %d | %11.6f | %11.6f\n", i, Af.get(i).getX(), Af.get(i).getY());
				}

				double xValue = 0;
				double fValue = 0;

				for (int i = 0; i < N; i++) {
					xValue += Ax.get(i).getY() * x.evaluate(Ax.get(i).getX());
					fValue += Af.get(i).getY() * f.evaluate(Af.get(i).getX());
				}
				System.out.printf(Locale.US,
						"\nf(x) = x^%d : %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
						N - 1, xValue, Math.abs(xValue - momentums.get(N - 1)),
						Math.abs(xValue - momentums.get(N - 1)) / xValue * 100);
				System.out.printf(Locale.US,
							"\nf(x) = sin(x) : %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
						fValue, Math.abs(fValue - exactValueF),
						Math.abs(fValue - exactValueF) / fValue * 100);

				System.out.println("\n\nКФНАСТ:\nМоменты:");
				for(int i = 0; i < 2*N; i++) {
					System.out.printf(Locale.US, "µ_%d = %.6f\n", i, momentums.get(i));
				}

				List<PlanePoint> pAf = integratorF.buildPreciseInterpolationQuadrature(N, momentums);
				List<PlanePoint> pAx3 = integratorX3.buildPreciseInterpolationQuadrature(N, momentums);
				double preciseFValue = 0;
				double x3Value = 0;

				System.out.println("\nf(x) = sin(x)\n k |     x_k     |     A_k    ");
				for (int i = 0; i < N; i++) {
					System.out.printf(Locale.US, " %d | %11.6f | %11.6f\n", i, pAf.get(i).getX(), pAf.get(i).getY());
					preciseFValue += pAf.get(i).getY() * f.evaluate(pAf.get(i).getX());
					x3Value += pAx3.get(i).getY() * x3.evaluate(pAx3.get(i).getX());
				}
				System.out.printf(Locale.US,
						"\nf(x) = x^%d : %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
						2*N - 1, x3Value, Math.abs(x3Value - momentums.get(2*N - 1)),
						Math.abs(x3Value - momentums.get(2*N - 1)) / x3Value * 100);
				System.out.printf(Locale.US,
						"\nf(x) = sin(x) : %.13f\nАбсолютная погрешность: %.8f\nОтносительная погрешность: %.4f%%\n",
						preciseFValue, Math.abs(preciseFValue - exactValueF),
						Math.abs(preciseFValue - exactValueF) / preciseFValue * 100);

				System.out.println("Хотите ли вы выбрать другие узлы [yes/no]?");
				repeat = input.nextLine();
			}
		}
	}

}
