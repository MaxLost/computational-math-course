package tasks.sem7;

import static org.math.computations.pde.Utils.generate2dMesh;
import static org.math.computations.pde.Utils.getValuesOnMesh;

import java.util.Locale;
import org.math.computations.PlanePoint;
import org.math.computations.functions.FunctionalR2;
import org.math.computations.pde.EllipticalDirichletSolver;

public class Task1 {

    public static void main(String[] args) {

        // Задача вариант 9 (тест на варианте 6)

        // Lu
        FunctionalR2 p = (x, y) -> 1;
        FunctionalR2 q = (x, y) -> 1;

        // Exact solution
        // FunctionalR2 uExact = (x, y) -> 2 * Math.pow(x, 3) * Math.pow(y, 3);
        FunctionalR2 uExact = (x, y) -> (Math.pow(x, 2) * Math.pow(y, 2)) / 2;

        // Variable
        // FunctionalR2 f = (x, y) -> -1 * (12 * x * y * (Math.pow(x, 2) + Math.pow(y, 2)));
        // FunctionalR2 mu = (x, y) -> uExact.evaluate(x, y);
        FunctionalR2 f = (x, y) -> Math.pow(x, 2) * Math.pow(y, 2);

        class boundaryFunction implements FunctionalR2 {
            @Override
            public double evaluate(double x, double y) {
                if (Math.abs(y - 1) < 10e-6) {
                    return Math.pow(x, 2) / 2;
                } else if (Math.abs(x - 1) < 10e-6) {
                    return Math.pow(y, 2) / 2;
                } else if (Math.abs(x) < 10e-6 || Math.abs(y) < 10e-6) {
                    return 0;
                }
                return 0;
            }
        }
        FunctionalR2 mu = new boundaryFunction();

        double[] boundaries = {0.0, 1.0, 0.0, 1.0};

        EllipticalDirichletSolver solver = new EllipticalDirichletSolver(p, q, f, mu, boundaries, uExact);

        PlanePoint[][] largeMesh = generate2dMesh(boundaries, new int[] {5, 5});
        for (int n = 5; n < 21; n *= 2) {

            System.out.printf("%n%n----- N = %d -----%n%n", n);

            PlanePoint[][] mesh = generate2dMesh(boundaries, new int[] {n, n});
            double[][] uValues = getValuesOnMesh(uExact, mesh);

            if (n == 5 && uExact != null) {
                System.out.printf("%n----- Точное решение -----%n");
                printSolutionOnMesh(uValues, largeMesh);
            }
            System.out.printf("%n----- Метод простой итерации -----%n");
            double[][] x1 = solver.solve("SI", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x1, largeMesh);
            }

            System.out.printf("%n----- Метод итерации с оптимальным параметром -----%n");
            double[][] x2 = solver.solve("IOP", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x2, largeMesh);
            }

            System.out.printf("%n----- Метод Зейделя -----%n");
            double[][] x3 = solver.solve("SE", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x3, largeMesh);
            }

            System.out.printf("%n----- Метод верхней релаксации -----%n");
            double[][] x4 = solver.solve("UR", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x4, largeMesh);
            }

            System.out.printf("%n----- Метод итераций с чебышёвским набором параметров -----%n");
            double[][] x5 = solver.solve("ICP", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x5, largeMesh);
            }

            System.out.printf("%n----- Попеременно-треугольный итерационный метод -----%n");
            double[][] x6 = solver.solve("AT", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x6, largeMesh);
            }

            System.out.printf("%n----- Попеременно треугольный итерационный метод с чебышёвским набором параметров -----%n");
            double[][] x7 = solver.solve("ATCP", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x7, largeMesh);
            }

            System.out.printf("%n----- Метод переменных направлений -----%n");
            double[][] x8 = solver.solve("AD", new int[]{n, n}, 0.01);
            if (n == 5) {
                printSolutionOnMesh(x8, largeMesh);
            }
        }
    }

    private static void printSolutionOnMesh(double[][] values, PlanePoint[][] mesh){

        System.out.printf("%n");
        StringBuilder str = new StringBuilder(" y / x |");
        for (int i = 0; i < mesh.length; i++) {
            str.append(String.format(Locale.US, " %2.2f |", mesh[0][i].getX()));
        }
        System.out.println(str);

        for (int i = 0; i < mesh.length; i++) {
            str = new StringBuilder(String.format(Locale.US, "%2.2f |", mesh[i][0].getX()));
            for (int j = 0; j < mesh[0].length; j++) {
                str.append(String.format(Locale.US, " %2.2f |", values[i][j]));
            }
            System.out.println(str);
        }
        return;
    }

}
