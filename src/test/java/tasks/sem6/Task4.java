package tasks.sem6;

import static org.math.computations.functions.Utils.generateMesh;

import java.util.ArrayList;
import java.util.List;
import org.jzy3d.chart.Chart;
import org.jzy3d.chart.factories.EmulGLChartFactory;
import org.jzy3d.chart.factories.IChartFactory;
import org.jzy3d.colors.Color;
import org.jzy3d.colors.ColorMapper;
import org.jzy3d.colors.colormaps.ColorMapRainbow;
import org.jzy3d.maths.Coord3d;
import org.jzy3d.maths.Range;
import org.jzy3d.plot3d.builder.SurfaceBuilder;
import org.jzy3d.plot3d.primitives.Shape;
import org.jzy3d.plot3d.rendering.canvas.Quality;
import org.math.computations.functions.Function;
import org.math.computations.functions.FunctionalR2;
import org.math.computations.pde.TransferEquationSolver;

public class Task4 {

    public static void main(String[] args) {

        /*
         *  u(x, t)
         *  C - const, C > 0
         *
         *  du/dt + C*du/dx = f(x,t)
         *
         *  f(x,t) = cos(x) + t
         *
         */

        FunctionalR2 f = (x, t) -> Math.cos(x) + t;
        // FunctionalR2 f = (x, t) -> t;
        double C = 1.0;
        List<Double> mesh = generateMesh(0, 2, 200);

        double[] u0 = getInitialApproximation(mesh);

        Range rangeX = new Range(0, 2);
        Range rangeT = new Range(0, 2);

        TransferEquationSolver solver = new TransferEquationSolver(f, C, 0, 2, 0, 2);
        double[][] solution = solver.solveImplicit(u0, 0.01, 0.01);

        System.out.println("DONE");

        List<Coord3d> points = transformTo3D(solution, rangeX, 200, rangeT, 200);

        Shape surface = new SurfaceBuilder().delaunay(points);
        surface
            .setColorMapper(new ColorMapper(new ColorMapRainbow(), surface, new Color(1, 1, 1, .5f)));
        surface.setFaceDisplayed(true);
        surface.setWireframeDisplayed(true);
        surface.setWireframeColor(Color.BLUE);

        IChartFactory factory = new EmulGLChartFactory();

        Chart chart = factory.newChart(Quality.Advanced().setHiDPIEnabled(true));
        chart.getScene().getGraph().add(surface);
        chart.getAxisLayout().setXAxisLabel("Time");
        chart.getAxisLayout().setYAxisLabel("X");
        chart.getAxisLayout().setZAxisLabel("U");
        chart.open();
        chart.addMouse();
        chart.addKeyboard();
    }

    private static double[] getInitialApproximation(List<Double> mesh) {

        double length = mesh.size() * Math.abs(mesh.get(1) - mesh.get(0));
        Function f = x -> Math.sin(x * Math.PI / length);
        double[] result = new double[mesh.size()];

        for (int i = 0; i < mesh.size(); i++) {
            result[i] = f.evaluate(mesh.get(i));
        }

        return result;
    }

    private static List<Coord3d> transformTo3D(double[][] solution, Range x, double stepsX, Range t, double stepsT) {

        List<Coord3d> result = new ArrayList<>();
        double stepX = (x.getMax() - x.getMin()) / stepsX;
        double stepT = (t.getMax() - t.getMin()) / stepsT;
        double pointerT = t.getMin();

        for (int i = 0; i < stepsT; i++) {

            double pointerX = x.getMin();
            for (int j = 0; j < stepsX; j++) {

                result.add(new Coord3d(pointerT, pointerX, solution[i][j]));
                pointerX += stepX;
            }
            pointerT += stepT;
        }

        return result;
    }

}
