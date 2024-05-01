package tasks.sem6;

import com.github.sh0nk.matplotlib4j.NumpyUtils;
import com.github.sh0nk.matplotlib4j.Plot;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Map.Entry;
import org.math.computations.functions.Function;
import org.math.computations.ode.SturmLiouvilleEigenvaluesFinder;

public class Task2 {

    public static void main(String[] args) throws PythonExecutionException, IOException {

        double k = 1.59194;
        double l = 5.42365;

        /*
        u(-1) = u(1) = 0

        p(x) = (k * x + l)
        q(x) = k^2 * (1 / (k * x + l) - k * x)
        */

        Function p = x -> k*x + l;
        Function dp = x -> k;
        Function q = x -> Math.pow(k, 2) * ((1 / (k*x + l)) - k*x);
        double[] boundaryConstraints = new double[] {1, 0, 0, 1, 0, 0};

        SturmLiouvilleEigenvaluesFinder solver = new SturmLiouvilleEigenvaluesFinder(p, q, -1, 1, boundaryConstraints);
        Map<Double, Function> result = solver.find(15, 5);

        Plot plt = Plot.create();

        System.out.printf(
            Locale.US,
            "Max(p) = %f\nMin(p) = %f\n",
            maxFunction(p, -1, 1),
            minFunction(p, -1 , 1)
        );

        System.out.printf(
            Locale.US,
            "Max(q) = %f\nMin(q) = %f\n",
            maxFunction(q, -1, 1),
            minFunction(q, -1 , 1)
        );

        List<Double> mesh = NumpyUtils.linspace(-1, 1, 1000);
        for (Entry<Double, Function> item : result.entrySet()) {
            List<Double> y = new ArrayList<>();
            for (int i = 0; i < mesh.size(); i++) {
                y.add(item.getValue().evaluate(mesh.get(i)));
            }

            plt.plot().add(mesh, y).label(item.getKey().toString());
        }

        plt.legend();
        plt.show();

        for (Entry<Double, Function> item : result.entrySet()) {
            System.out.println(item.getValue());
            Plot plt1 = Plot.create();
            List<Double> lambdau = new ArrayList<>();
            List<Double> Lu = new ArrayList<>();

            for (int i = 0; i < mesh.size(); i++) {
                Lu.add(Lu(item.getValue(), p, q, dp, mesh.get(i)));
                lambdau.add(item.getKey() * item.getValue().evaluate(mesh.get(i)));
            }

            plt1.plot().add(mesh, lambdau).label("λu");
            plt1.plot().add(mesh, Lu).label("Lu");
            plt1.legend();
            plt1.title("λ = " + item.getKey());
            plt1.show();
        }
    }

    private static double Lu(Function u, Function p, Function q, Function dp, double x) {

        double step = 0.001;
        double du = (u.evaluate(x + step) - u.evaluate(x - step)) / (2*step);
        double duu = (u.evaluate(x - step) - 2*u.evaluate(x) + u.evaluate(x + step)) / (step * step);
        //double dp = (p.evaluate(x + step) - p.evaluate(x - step)) / (2*step);
        double Lu = -1*(p.evaluate(x)*duu + dp.evaluate(x) * du) + q.evaluate(x)*u.evaluate(x);

        return Lu;
    }

    private static double minFunction(Function f, double lowerBound, double upperBound) {

        double step = (upperBound - lowerBound) / 1000;
        double x = lowerBound;
        double min = f.evaluate(x);

        while (x < upperBound) {
            x += step;
            min = Math.min(min, f.evaluate(x));
        }
        return min;
    }

    private static double maxFunction(Function f, double lowerBound, double upperBound) {

        double step = (upperBound - lowerBound) / 1000;
        double x = lowerBound;
        double max = f.evaluate(x);

        while (x < upperBound) {
            x += step;
            max = Math.max(max, f.evaluate(x));
        }
        return max;
    }

}
