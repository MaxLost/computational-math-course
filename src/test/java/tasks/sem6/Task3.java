package tasks.sem6;

import com.github.sh0nk.matplotlib4j.NumpyUtils;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;
import java.io.IOException;
import java.util.List;
import java.util.Locale;
import org.math.computations.functions.Function;
import org.math.computations.functions.FunctionalR2;
import org.math.computations.functions.FredholmIntegralEquationSolver;

public class Task3 {

    public static void main(String[] args) throws PythonExecutionException, IOException {

        /*
            u(x) - 0.6 * Int[cosh(xy)*u(y), y, 0, 1] = x - 0.6
        */

        FunctionalR2 k = (x, y) -> Math.cosh(x * y);
        Function f = x -> x - 0.6;
        double lambda = 0.6;

        System.out.println("Метод вырожденного ядра:");
        for (int i = 1; i < 11; i++) {

            FredholmIntegralEquationSolver solver = new FredholmIntegralEquationSolver(k, f, lambda,
                0, 1);
            Function u = solver.solveSingularKernel(i);

            List<Double> mesh = NumpyUtils.linspace(0, 1, 1000);
            Function leftSide = x -> u.evaluate(x) - lambda * dotProductL2(y -> Math.cosh(x * y), u, 0, 1);

            double loss = 0.0;
            for (Double x : mesh) {
                loss = Math.max(loss, Math.abs(leftSide.evaluate(x) - f.evaluate(x)));
            }

            System.out.printf(Locale.US, "N = %d | Невязка = %.5e\n", i, loss);
        }

        System.out.println("Метод механических квадратур:");
        int baseAmount = 10;
        for (int i = 0; i < 6; i++) {

            FredholmIntegralEquationSolver solver = new FredholmIntegralEquationSolver(k, f, lambda,
                0, 1);
            Function u = solver.solveQuadrature(baseAmount * (1 << i));

            List<Double> mesh = NumpyUtils.linspace(0, 1, 1000);
            Function leftSide = x -> u.evaluate(x) - lambda * dotProductL2(y -> Math.cosh(x * y), u, 0, 1);

            double loss = 0.0;
            for (Double x : mesh) {
                loss = Math.max(loss, Math.abs(leftSide.evaluate(x) - f.evaluate(x)));
            }

            System.out.printf(Locale.US, "N = %d | Невязка = %.5e\n", baseAmount * (1 << i), loss);
        }

    }

    private static double dotProductL2(Function f, Function g, double lowerBound, double upperBound) {

        double step = (upperBound - lowerBound) / 1000;
        double x = lowerBound + step / 2;

        double value = 0;
        while (x < upperBound) {
            value += f.evaluate(x) * g.evaluate(x);
            x += step;
        }

        return value * step;
    }

}
