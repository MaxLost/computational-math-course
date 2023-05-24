package tasks;

import org.math.computations.functions.Function;
import org.math.computations.functions.Integrator;

import java.util.Locale;
import java.util.Scanner;

import static tasks.Utils.ddfMax;
import static tasks.Utils.dfMax;

public class Task4_3 {

	public static void integrationTest(Function f, Function F, double A, double B, int m, int l) {

		double exactValue = F.evaluate(B) - F.evaluate(A);
		System.out.print("\nТочное значение интеграла: " + exactValue);

		Integrator task = new Integrator(f, A, B);

		double dfMax = dfMax(f, A, B);
		double ddfMax = ddfMax(f, A, B);

		System.out.printf(Locale.US, "\nЗначение h: %5.3e", (B - A) / (l*m));

		double mValue = task.integrate("LR", m);
		double lmValue = task.integrate("LR", l*m);
		double rungeValue = (l * lmValue - mValue) / (l - 1);
		System.out.printf(Locale.US, """
						\n\n\nФормула левых прямоугольников:
						M = %d, L*M = %d
						\tЗначение интеграла при M разбиений: %14.8f
						\tАбсолютная погрешность при M разбиений: %10.3e
						\tОтносительная погрешность при M разбиений: %5.3f%%
						\n\tЗначение интеграла при L*M разбиений: %14.8f
						\tАбсолютная погрешность при L*M разбиений: %10.3e
						\tОтносительная погрешность при L*M разбиений: %5.3f%%
						\n\tУточнённое значение интеграла: %14.8f
						\tАбсолютная погрешность уточнённого значения: %10.3e
						\tОтносительная погрешность уточнённого значения: %5.3f%%""",
				m, l*m, mValue, Math.abs(mValue - exactValue), Math.abs(mValue - exactValue) / Math.abs(mValue) * 100,
				lmValue, Math.abs(lmValue - exactValue), Math.abs(lmValue - exactValue) / Math.abs(lmValue) * 100,
				rungeValue, Math.abs(rungeValue - exactValue), Math.abs(rungeValue - exactValue) / Math.abs(rungeValue) * 100);

		mValue = task.integrate("RR", m);
		lmValue = task.integrate("RR", l * m);
		rungeValue = (l * lmValue - mValue) / (l - 1);
		System.out.printf(Locale.US, """
						\n\n\nФормула правых прямоугольников:
						\tЗначение интеграла при M разбиений: %14.8f
						\tАбсолютная погрешность при M разбиений: %10.3e
						\tОтносительная погрешность при M разбиений: %5.3f%%
						\n\tЗначение интеграла при L*M разбиений: %14.8f
						\tАбсолютная погрешность при L*M разбиений: %10.3e
						\tОтносительная погрешность при L*M разбиений: %5.3f%%
						\n\tУточнённое значение интеграла: %14.8f
						\tАбсолютная погрешность уточнённого значения: %10.3e
						\tОтносительная погрешность уточнённого значения: %5.3f%%""",
				mValue, Math.abs(mValue - exactValue), Math.abs(mValue - exactValue) / Math.abs(mValue) * 100,
				lmValue, Math.abs(lmValue - exactValue), Math.abs(lmValue - exactValue) / Math.abs(lmValue) * 100,
				rungeValue, Math.abs(rungeValue - exactValue), Math.abs(rungeValue - exactValue) / Math.abs(rungeValue) * 100);

		mValue = task.integrate("CR", m);
		lmValue = task.integrate("CR", l*m);
		rungeValue = (Math.pow(l, 2) * lmValue - mValue) / (Math.pow(l, 2) - 1);
		System.out.printf(Locale.US, """
						\n\n\nФормула серединных прямоугольников:
						\tЗначение интеграла при M разбиений: %14.8f
						\tАбсолютная погрешность при M разбиений: %10.3e
						\tОтносительная погрешность при M разбиений: %5.3f%%
						\n\tЗначение интеграла при L*M разбиений: %14.8f
						\tАбсолютная погрешность при L*M разбиений: %10.3e
						\tОтносительная погрешность при L*M разбиений: %5.3f%%
						\n\tУточнённое значение интеграла: %14.8f
						\tАбсолютная погрешность уточнённого значения: %10.3e
						\tОтносительная погрешность уточнённого значения: %5.3f%%""",
				mValue, Math.abs(mValue - exactValue), Math.abs(mValue - exactValue) / Math.abs(mValue) * 100,
				lmValue, Math.abs(lmValue - exactValue), Math.abs(lmValue - exactValue) / Math.abs(lmValue) * 100,
				rungeValue, Math.abs(rungeValue - exactValue), Math.abs(rungeValue - exactValue) / Math.abs(rungeValue) * 100);

		mValue = task.integrate("TR", m);
		lmValue = task.integrate("TR", l * m);
		rungeValue = (Math.pow(l, 2) * lmValue - mValue) / (Math.pow(l, 2) - 1);
		System.out.printf(Locale.US, """
						\n\n\nФормула трапеций:
						\tЗначение интеграла при M разбиений: %14.8f
						\tАбсолютная погрешность при M разбиений: %10.3e
						\tОтносительная погрешность при M разбиений: %5.3f%%
						\n\tЗначение интеграла при L*M разбиений: %14.8f
						\tАбсолютная погрешность при L*M разбиений: %10.3e
						\tОтносительная погрешность при L*M разбиений: %5.3f%%
						\n\tУточнённое значение интеграла: %14.8f
						\tАбсолютная погрешность уточнённого значения: %10.3e
						\tОтносительная погрешность уточнённого значения: %5.3f%%""",
				mValue, Math.abs(mValue - exactValue), Math.abs(mValue - exactValue) / Math.abs(mValue) * 100,
				lmValue, Math.abs(lmValue - exactValue), Math.abs(lmValue - exactValue) / Math.abs(lmValue) * 100,
				rungeValue, Math.abs(rungeValue - exactValue), Math.abs(rungeValue - exactValue) / Math.abs(rungeValue) * 100);

		mValue = task.integrate("SI", m);
		lmValue = task.integrate("SI", l * m);
		rungeValue = (Math.pow(l, 4) * lmValue - mValue) / (Math.pow(l, 4) - 1);
		System.out.printf(Locale.US, """
						\n\n\nФормула Симпсона:
						\tЗначение интеграла при M разбиений: %14.8f
						\tАбсолютная погрешность при M разбиений: %10.3e
						\tОтносительная погрешность при M разбиений: %5.3f%%
						\n\tЗначение интеграла при L*M разбиений: %14.8f
						\tАбсолютная погрешность при L*M разбиений: %10.3e
						\tОтносительная погрешность при L*M разбиений: %5.3f%%
						\n\tУточнённое значение интеграла: %14.8f
						\tАбсолютная погрешность уточнённого значения: %10.3e
						\tОтносительная погрешность уточнённого значения: %5.3f%%""",
				mValue, Math.abs(mValue - exactValue), Math.abs(mValue - exactValue) / Math.abs(mValue) * 100,
				lmValue, Math.abs(lmValue - exactValue), Math.abs(lmValue - exactValue) / Math.abs(lmValue) * 100,
				rungeValue, Math.abs(rungeValue - exactValue), Math.abs(rungeValue - exactValue) / Math.abs(rungeValue) * 100);
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
				Уточнение значения интеграла с помощью метода Рунге
				                
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

			int l = 0;
			while (l < 2) {
				System.out.println("Введите L >= 2 - множитель количества разбиений для вычисления уточнённого значения" +
						"интеграла по методу Рунге-Ромберга");
				l = input.nextInt();
				if (l < 2) {
					System.out.println("M должно быть целым числом больше 1, попробуйте ещё раз.");
				}
			}

			integrationTest(f, F, A, B, m, l);
			//integrationTest(c, C, A, B, m);
			//integrationTest(x, X, A, B, m);
			//integrationTest(x3, X3, A, B, m);
		}

	}
}
