package tasks.sem4;

import org.math.computations.PlanePoint;
import org.math.computations.functions.AlgebraicInterpolator;
import org.math.computations.functions.Function;

import java.util.*;

import static tasks.Utils.divideSegmentRandom;

public class Task2 {

	public static void main(String[] args) {

		// file for test: src/resources/interpolation_test.txt

		System.out.println("""
                Задача алгебраического интерполирования
                
                15 вариант
                
                Введите число узлов интерполяции:""");
		try (Scanner input = new Scanner(System.in).useLocale(Locale.US)){

			int M = input.nextInt();

			Function f = t -> Math.sin(t) - Math.pow(t, 2) / 2;

			System.out.println("Введите A - левую границу промежутка для выбора узлов интерполирования: ");
			double A = input.nextDouble();
			System.out.println("Введите B - правую границу промежутка для выбора узлов интерполирования: ");
			double B = input.nextDouble();

			List<PlanePoint> nodes = divideSegmentRandom(f, A, B, M);

			/*while (true) {

				System.out.println("Вы хотите задать таблицу из файла? [yes/no]:");
				try {
					String clearance = input.next();
					if (clearance.equals("yes")) {

						while (true) {
							System.out.println("Введите путь от текущей директории (" +
									Paths.get("").toAbsolutePath() + "\\) до файла с таблицей");
							String path = input.next();

							try (Scanner fileReader = new Scanner(Paths.get(path)).useLocale(Locale.US)) {

								for (int i = 0; i < M; i++) {
									double a = fileReader.nextDouble();
									double b = fileReader.nextDouble();
									nodes.add(new PlanePoint(a, b));
								}
								break;

							} catch (IOException e) {
								System.out.println("Ошибка доступа к файлу, попробуйте снова");
							}
						}
						break;
					} else if (clearance.equals("no")) {

						for (int i = 0; i < M; i++) {
							System.out.println("Введите" + (i + 1) + "узел и значение функции в нём:");
							nodes.add(new PlanePoint(input.nextDouble(), input.nextDouble()));
						}
						break;
					} else {
						throw new NoSuchElementException();
					}
				} catch (NoSuchElementException e) {
					System.out.println("Неизвестный ответ, попробуйте снова");
					System.out.println(Arrays.toString(e.getStackTrace()));
				}
			}*/

			System.out.println("Таблица значений в точках для функции f: \n\tx\t | \tf(x)\t");
			for (PlanePoint p : nodes){
				System.out.printf(Locale.US, "%.6f | %.6f\n", p.getX(), p.getY());
			}

			String status = "yes";
			AlgebraicInterpolator task = new AlgebraicInterpolator(nodes);

			while (!(status.equals("no"))) {
				System.out.println("\nВведите точку интерполирования: ");
				double X = input.nextDouble();
				int N = Integer.MAX_VALUE;
				while (N > M) {
					System.out.println("\nВведите N - степень интерполяционного многочлена (N<=" + (M - 1) + "): ");
					N = input.nextInt();
					if (N > M) {
						System.out.println("Введено неверное N, попробуйте снова");
					}
				}

				List<Double> result = task.interpolateLogged(N, X);

				System.out.println("\nЗначение абсолютной фактической погрешности для формы Лагранжа: "
									+ Math.abs(result.get(0) - f.evaluate(X)));
				System.out.println("Значение абсолютной фактической погрешности для формы Ньютона: "
									+ Math.abs(result.get(1) - f.evaluate(X)));

				System.out.println("Хотите ли вы продолжить работу и получить значение для другой точки? [yes/no]: ");
				status = input.next();
			}
		}
	}
}
