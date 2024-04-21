package tasks.sem6;

import com.github.sh0nk.matplotlib4j.NumpyUtils;
import com.github.sh0nk.matplotlib4j.Plot;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import org.math.computations.functions.Function;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.LinearSystemSolver;
import org.math.computations.matrices.Matrix;
import org.math.computations.ode.SecondOrderBVPSolver;

public class Task1 {

    private static List<Double> readSolution(String path) throws IOException {

        Path file = Paths.get(path);
        List<Double> values = new ArrayList<>();

        try (Scanner input = new Scanner(file)) {
            String line;
            while (!(line = input.nextLine()).equals("end")) {
                values.add(Double.parseDouble(line));
            }
        }

        return values;
    }

    public static void main(String[] args) throws PythonExecutionException, IOException {

        Function p = x -> (5 + x) / (9 - 3*x);
        Function q = x -> (1 - x / 2.0);
        Function r = x -> Math.log(x + 3) / 3;
        Function f = x -> (1 + x) / 3;

        double[] boundaryConstraints = new double[] {-3, -1, 0, 0, 1, 0};

        SecondOrderBVPSolver bvpSolver = new SecondOrderBVPSolver(p, q, r, f, boundaryConstraints);

        List<Double> solution = readSolution("src/resources/sem6/task1_solution.txt");
        List<Double> mesh = NumpyUtils.linspace(-1, 1, solution.size());

        Function n3 = bvpSolver.solveMinSquareMethod(3);
        Function n5 = bvpSolver.solveMinSquareMethod(5);
        Function n7 = bvpSolver.solveMinSquareMethod(7);

        List<Double> n3values = new ArrayList<>();
        List<Double> n5values = new ArrayList<>();
        List<Double> n7values = new ArrayList<>();

        for (Double x : mesh) {
            n3values.add(n3.evaluate(x));
            n5values.add(n5.evaluate(x));
            n7values.add(n7.evaluate(x));
        }

        System.out.println("Метод наименьших квадратов:");
        System.out.println("n = 3 | " + n3);
        System.out.println("n = 5 | " + n5);
        System.out.println("n = 7 | " + n7);
        System.out.println(" ");

        Plot plt = Plot.create();
        plt.plot().add(mesh, solution).label("exact");
        plt.plot().add(mesh, n3values).label("n = 3");
        plt.plot().add(mesh, n5values).label("n = 5");
        plt.plot().add(mesh, n7values).label("n = 7");
        plt.legend();
        plt.title("Метод наименьших квадратов");
        plt.show();
        System.out.println(" ");

        Function n3cc = bvpSolver.solveCollocation(3);
        Function n5cc = bvpSolver.solveCollocation(5);
        Function n7cc = bvpSolver.solveCollocation(7);

        List<Double> n3ccValues = new ArrayList<>();
        List<Double> n5ccValues = new ArrayList<>();
        List<Double> n7ccValues = new ArrayList<>();

        for (Double x : mesh) {
            n3ccValues.add(n3cc.evaluate(x));
            n5ccValues.add(n5cc.evaluate(x));
            n7ccValues.add(n7cc.evaluate(x));
        }

        System.out.println("Метод коллокации с узлами Чебышёва:");
        System.out.println("n = 3 | " + n3cc);
        System.out.println("n = 5 | " + n5cc);
        System.out.println("n = 7 | " + n7cc);
        System.out.println(" ");

        Plot plt1 = Plot.create();
        plt1.plot().add(mesh, solution).label("exact");
        plt1.plot().add(mesh, n3ccValues).label("n = 3");
        plt1.plot().add(mesh, n5ccValues).label("n = 5");
        plt1.plot().add(mesh, n7ccValues).label("n = 7");
        plt1.legend();
        plt1.title("Метод коллокации с узлами Чебышёва");
        plt1.show();

    }

}
