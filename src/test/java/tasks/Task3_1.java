package tasks;

import org.math.computational.PlanePoint;
import org.math.computational.Polynomial;
import org.math.computational.functions.AlgebraicInterpolator;
import org.math.computational.functions.Function;
import org.math.computational.functions.FunctionRootFinder;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Scanner;

public class Task3_1 {

	public static void main(String[] args) {

		System.out.println("""
                Задача обратного интерполирования
                
                15 вариант
                f(x) = sin(x) - x^2/2
                
                Введите число узлов интерполяции:""");
		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)){

			int M = input.nextInt();

			System.out.println("Введите A - левую границу промежутка для выбора узлов интерполирования: ");
			double A = input.nextDouble();
			System.out.println("Введите B - правую границу промежутка для выбора узлов интерполирования: ");
			double B = input.nextDouble();

			Function f = t -> Math.sin(t) - Math.pow(t, 2) / 2;

			List<PlanePoint> nodes = Utils.divideSegmentEqual(f, A, B, M - 1);
			List<PlanePoint> nodesSwapped = new ArrayList<>();

			System.out.println("Таблица значений в точках для функции f: \n\tx\t | \tf(x)\t");
			for (PlanePoint p : nodes){
				System.out.printf(Locale.US, "%.6f | %.6f\n", p.getX(), p.getY());
				nodesSwapped.add(new PlanePoint(p.getY(), p.getX()));
			}

			String status = "yes";
			AlgebraicInterpolator taskReversedF = new AlgebraicInterpolator(nodesSwapped);
			AlgebraicInterpolator taskF = new AlgebraicInterpolator(nodes);

			while (!(status.equals("no"))) {
				System.out.println("\nВведите значение для которого необходимо найти агрумент: ");
				double F = input.nextDouble();
				int N = Integer.MAX_VALUE;
				while (N > M) {
					System.out.println("\nВведите N - степень интерполяционного многочлена (N<=" + (M - 1) + "): ");
					N = input.nextInt();
					if (N > M) {
						System.out.println("Введено неверное N, попробуйте снова");
					}
				}

				System.out.println("1 способ:");
				double x1 = taskReversedF.interpolate(N, F);
				System.out.println("Значение аргумента: " + x1);
				System.out.println("Значение модуля невязки: " + Math.abs(f.evaluate(x1) - F));

				System.out.println("2 способ: \nВведите точность для нахождения аргумента:");
				double epsilon = input.nextDouble();
				Polynomial p = taskF.getInterpolationPolynomial(N, nodes.get(nodes.size() / 2).getX());
				Function g = t -> p.evaluate(t) - F;
				FunctionRootFinder argumentFinder = new FunctionRootFinder(
						g,
						nodes.get(0).getX(),
						nodes.get(nodes.size() - 1).getX(),
						epsilon
				);

				List<Double> roots = argumentFinder.findRoots(
						(nodes.get(nodes.size() - 1).getX() - nodes.get(0).getX()) / 500
				);

				if (roots.size() > 0) {
					System.out.println("Возможные аргументы при которых функция принимает введённое значение: ");
					for (int i = 0; i < roots.size(); i++) {
						System.out.println("\n" + (i + 1) + " возможный аргумент: " + roots.get(i));
						System.out.println("Значение модуля невязки: " + Math.abs(Math.abs(f.evaluate(roots.get(i)) - F)));
					}
				} else {
					System.out.println("На данном промежутке нет аргументов при которых функция принимала бы "
							+ "заданное значение");
				}

				System.out.println("Хотите ли вы продолжить работу и получить аргумент для другого значения? [yes/no]: ");
				status = input.next();
			}
		}

	}

}
