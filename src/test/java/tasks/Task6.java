package tasks;

import org.math.computational.PlanePoint;
import org.math.computational.Polynomial;
import org.math.computational.functions.Function;
import org.math.computational.functions.TransformR2;
import org.math.computational.ode.FirstOrderIVPSolver;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Scanner;
import java.util.stream.Collectors;

public class Task6 {

	private static Polynomial getTaylorExpansion(double x0, double y0, TransformR2 f) {

		double dy = f.evaluate(x0, y0);
		double ddy = -2*y0*dy;
		double dddy = -2*(Math.pow(dy, 2) + y0*ddy);
		double ddddy = -2*(2*dy*ddy + dy*ddy + y0*dddy);

		return new Polynomial(List.of(y0, dy, ddy / 2, dddy / 6, ddddy / 24));
	}

	public static void main(String[] args) {

		System.out.println("""
                Решение задачи Коши для ОДУ первого порядка, разрешённого относительно производной
                y'(x) = f(x, y)

                1 вариант
                
                f(x, y) = -y^2 + 1
                y(0) = 0
                
                h = 0.1
                N = 10""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)){
			double h = -1;
			while (h < 0) {
				System.out.println("Введите h - расстояние между точками в которых будет вычисляться значение решения");
				h = input.nextDouble();
				if (h < 0) {
					System.out.println("h должно быть больше 0, попробуйте ввести большее h");
				}
			}

			int N = -1;
			while (N < 1) {
				System.out.println("Введите N - количество точек справа от x_0, в которых будет вычисляться значение решения");
				N = input.nextInt();
				if (N < 1) {
					System.out.println("N должно быть не меньше 1, попробуйте ввести большее N");
				}
			}

			double x0 = 0;
			double y0 = 0;
			TransformR2 f = (u, v) -> -1*Math.pow(v, 2) + 1;

			Function phi = t -> (Math.exp(2*t) - 1) / (Math.exp(2*t) + 1); // Solution of given IVP

			List<Double> x = new ArrayList<>();
			for (int i = -2; i <= N; i++) {
				x.add(x0 + i*h);
			}

			List<PlanePoint> exactValues = new ArrayList<>();

			System.out.println("\nТочное решение\n    x_k     |     f(x_k)    ");
			for (int i = 0; i < x.size(); i++) {
				double value = phi.evaluate(x.get(i));
				System.out.printf(Locale.US, "%11.6f | %11.6f\n", x.get(i), value);
				exactValues.add(new PlanePoint(x.get(i), value));
			}

			List<PlanePoint> taylor = new ArrayList<>();
			Polynomial T = getTaylorExpansion(x0, y0, f);
			System.out.println("\nРазложение в ряд Тейлора\n    x_k     |     f(x_k)     |    Абс.погр.");
			for (int i = 0; i < x.size(); i++) {
				taylor.add(new PlanePoint(x.get(i), T.evaluate(x.get(i))));
				System.out.printf(Locale.US, "%11.6f | %14.6f | %12.3e\n",
						taylor.get(i).getX(), taylor.get(i).getY(), Math.abs(taylor.get(i).getY() - exactValues.get(i).getY()));
			}

			FirstOrderIVPSolver solver = new FirstOrderIVPSolver(x0, y0, h, N);

			List<PlanePoint> adams = solver.solveExtrapolationAdams(taylor.stream().limit(5).collect(Collectors.toList()));
			System.out.println("\nЭкстраполяционный метод Адамса 4-го порядка\n    x_k     |     f(x_k)     |    Абс.погр.");
			for (int i = 0; i < adams.size(); i++) {
				System.out.printf(Locale.US, "%11.6f | %14.6f | %12.3e\n",
						adams.get(i).getX(), adams.get(i).getY(), Math.abs(adams.get(i).getY() - exactValues.get(i).getY()));
			}

			List<PlanePoint> rungeKutta = solver.solveRungeKutta(f);
			System.out.println("\nМетод Рунге-Кутты 4-го порядка\n    x_k     |     f(x_k)     |    Абс.погр.");
			for (int i = 0; i < rungeKutta.size(); i++) {
				System.out.printf(Locale.US, "%11.6f | %14.6f | %12.3e\n",
						rungeKutta.get(i).getX(), rungeKutta.get(i).getY(), Math.abs(rungeKutta.get(i).getY() - exactValues.get(i + 2).getY()));
			}

			List<PlanePoint> euler = solver.solveEuler(f);
			System.out.println("\nМетод Эйлера\n    x_k     |     f(x_k)     |    Абс.погр.");
			for (int i = 0; i < euler.size(); i++) {
				System.out.printf(Locale.US, "%11.6f | %14.6f | %12.3e\n",
						euler.get(i).getX(), euler.get(i).getY(), Math.abs(euler.get(i).getY() - exactValues.get(i + 2).getY()));
			}

			List<PlanePoint> euler1 = solver.solveEuler1(f);
			System.out.println("\nМетод Эйлера I\n    x_k     |     f(x_k)     |    Абс.погр.");
			for (int i = 0; i < euler1.size(); i++) {
				System.out.printf(Locale.US, "%11.6f | %14.6f | %12.3e\n",
						euler1.get(i).getX(), euler1.get(i).getY(), Math.abs(euler1.get(i).getY() - exactValues.get(i + 2).getY()));
			}

			List<PlanePoint> euler2 = solver.solveEuler2(f);
			System.out.println("\nМетод Эйлера II\n    x_k     |     f(x_k)     |    Абс.погр.");
			for (int i = 0; i < euler2.size(); i++) {
				System.out.printf(Locale.US, "%11.6f | %14.6f | %12.3e\n",
						euler2.get(i).getX(), euler2.get(i).getY(), Math.abs(euler2.get(i).getY() - exactValues.get(i + 2).getY()));
			}
		}
	}

}
