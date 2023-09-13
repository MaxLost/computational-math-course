package tasks.sem4;

import org.math.computations.PlanePoint;
import org.math.computations.functions.NumericalDifferentiator;
import org.math.computations.functions.Function;

import java.util.List;
import java.util.Locale;
import java.util.Scanner;
import tasks.Utils;

public class Task3_2 {

	public static void main(String[] args) {

		System.out.println("""
                Задача численного дифференцирования
                
                15 вариант
                f(x) = e^(1.5x)
                
                Введите A - левую границу промежутка:""");

		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)) {

			double A = input.nextDouble();

			int M = -1;
			while (M < 2) {
				System.out.println("Введите M - количество значений в таблице: ");
				M = input.nextInt();
				if (M < 2) {
					System.out.println("Слишком малое M, попробуйте ввести большее M");
				}
			}

			double H = -1;
			while (H <= 0) {
				System.out.println("Введите H - расстояние между точками в таблице: ");
				H = input.nextDouble();
				if (H <= 0) {
					System.out.println("Расстояние должно быть положительным, попробуйте ввести его снова");
				}
			}

			Function f = t -> Math.exp(1.5*t);

			List<PlanePoint> nodes = Utils.divideSegmentEqual(f, A, A + (M - 1)*H, M - 1);

			System.out.println("Таблица значений в точках для функции f: \n\tx\t | \tf(x)\t");
			for (PlanePoint p : nodes) {
				System.out.printf(Locale.US, "%.6f | %.6f\n", p.getX(), p.getY());
			}

			NumericalDifferentiator task = new NumericalDifferentiator(nodes);

			List<PlanePoint> dfValues = task.computeDerivatives();
			List<PlanePoint> ddfValues = task.computeSecondDerivatives();

			Function df = t -> 1.5*Math.exp(1.5*t);
			Function ddf = t -> 2.25*Math.exp(1.5*t);

			System.out.println("\nТаблица с результатами вычислений для функции f:\n");
			System.out.println("\t x\t   |\t f(x)\t|\t  f'(x)   | Abs Error f'(x) | Rel Error f'(x) " +
					"|   f''(x)   | Abs Error f''(x) | Rel Error f''(x)");
			System.out.printf(
					Locale.US,
					"%10.6f | %10.3e | %10.3e  | %15.3e | %15.3f |      -     |  \t\t-\t\t  | \t-\t\n",
					nodes.get(0).getX(),
					nodes.get(0).getY(),
					dfValues.get(0).getY(),
					Math.abs(dfValues.get(0).getY() - df.evaluate(nodes.get(0).getX())),
					Math.abs(dfValues.get(0).getY() - df.evaluate(nodes.get(0).getX()))
							/ Math.abs(df.evaluate(nodes.get(0).getX()))
			);

			for (int i = 1; i < nodes.size() - 1; i++) {

				System.out.printf(
						Locale.US,
						"%10.6f | %10.3e | %10.3e  | %15.3e | %15.3f | %10.3e | %16.3e | %10.3f\n",
						nodes.get(i).getX(),
						nodes.get(i).getY(),
						dfValues.get(i).getY(),
						Math.abs(dfValues.get(i).getY() - df.evaluate(nodes.get(i).getX())),
						Math.abs(dfValues.get(i).getY() - df.evaluate(nodes.get(i).getX()))
								/ Math.abs(df.evaluate(nodes.get(i).getX())),
						ddfValues.get(i - 1).getY(),
						Math.abs(ddfValues.get(i - 1).getY() - ddf.evaluate(nodes.get(i).getX())),
						Math.abs(ddfValues.get(i - 1).getY() - ddf.evaluate(nodes.get(i).getX()))
								/ Math.abs(ddf.evaluate(nodes.get(i).getX()))
				);
			}

			System.out.printf(
					Locale.US,
					"%10.6f | %10.3e | %10.3e  | %15.3e | %15.3f |      -     |  \t\t-\t\t  | \t-\t\n",
					nodes.get(M - 1).getX(),
					nodes.get(M - 1).getY(),
					dfValues.get(M - 1).getY(),
					Math.abs(dfValues.get(M - 1).getY() - df.evaluate(nodes.get(M - 1).getX())),
					Math.abs(dfValues.get(M - 1).getY() - df.evaluate(nodes.get(M - 1).getX()))
							/ Math.abs(df.evaluate(nodes.get(M - 1).getX()))
			);

			System.out.println("\nХотите ли вы ввести другие параметры? [yes/no]");
			input.nextLine();
			String status = input.nextLine();
			if (status.equals("yes") || status.equals("YES") || status.equals("Yes")) {
				main(new String[0]);
			}

		}

	}

}
