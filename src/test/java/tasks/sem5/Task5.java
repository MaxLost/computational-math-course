package tasks.sem5;

import java.util.List;
import java.util.Locale;
import org.math.computations.PlanePoint;
import org.math.computations.functions.Function;
import org.math.computations.ode.SecondOrderBVPSolver;

public class Task5 {

    public static void main(String[] args) {

        /*
            ODE:
            - ((4 + x) / (5 + 2x))y'' + (x / 2 - 1)y' + (1 + e^(x / 2))y = 2 + x

            Boundaries:
            y'(-1) = 0
            y'(1) + 2*y(1) = 0
         */

        Function p = x -> (4 + x) / (5 + 2*x);
        Function q = x -> x / 2 - 1;
        Function r = x -> 1 + Math.exp(x / 2);
        Function f = x -> 2 + x;

        SecondOrderBVPSolver solver = new SecondOrderBVPSolver(p, q, r, f,
            new double[]{0, -1, 0, 2, 1, 0});

        List<PlanePoint> firstApproximation
            = solver.solveBvpDiff(-1, 1, 10, 1, false);
        List<PlanePoint> firstApproximationRichardson
            = solver.solveBvpDiff(-1, 1, 20, 1, true);

        List<PlanePoint> secondApproximation
            = solver.solveBvpDiff(-1, 1, 10, 2, false);
        List<PlanePoint> secondApproximationRichardson
            = solver.solveBvpDiff(-1, 1, 20, 2, true);

        System.out.println(
            "\nO(h), n = 10\n    x_k     |     f(x_k)");
        for (PlanePoint planePoint : firstApproximation) {
            System.out.printf(Locale.US, "%11.6f | %14.6f\n",
                planePoint.getX(), planePoint.getY());
        }

        System.out.println(
            "\nO(h), n = 20, with Richardson\n    x_k     |     f(x_k)");
        for (PlanePoint planePoint : firstApproximationRichardson) {
            System.out.printf(Locale.US, "%11.6f | %14.6f\n",
                planePoint.getX(), planePoint.getY());
        }

        System.out.println(
            "\nO(h^2), n = 10\n    x_k     |     f(x_k)");
        for (PlanePoint planePoint : secondApproximation) {
            System.out.printf(Locale.US, "%11.6f | %14.6f\n",
                planePoint.getX(), planePoint.getY());
        }

        System.out.println(
            "\nO(h^2), n = 20, with Richardson\n    x_k     |     f(x_k)");
        for (PlanePoint planePoint : secondApproximationRichardson) {
            System.out.printf(Locale.US, "%11.6f | %14.6f\n",
                planePoint.getX(), planePoint.getY());
        }
    }

}
