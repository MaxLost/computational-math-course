package tasks.sem6;

import static org.math.computations.functions.Utils.generateMesh;

import java.util.List;
import org.jzy3d.chart.Chart;
import org.jzy3d.chart.factories.EmulGLChartFactory;
import org.jzy3d.chart.factories.IChartFactory;
import org.jzy3d.colors.Color;
import org.jzy3d.colors.ColorMapper;
import org.jzy3d.colors.colormaps.ColorMapRainbow;
import org.jzy3d.maths.Coord3d;
import org.jzy3d.maths.Range;
import org.jzy3d.plot3d.builder.Mapper;
import org.jzy3d.plot3d.builder.SurfaceBuilder;
import org.jzy3d.plot3d.builder.concrete.OrthonormalGrid;
import org.jzy3d.plot3d.primitives.Shape;
import org.jzy3d.plot3d.rendering.canvas.Quality;
import org.math.computations.functions.FunctionalR2;
import org.math.computations.pde.TransferEquationSolver;

public class Task4ExactSolutionPlot {

    public static void main(String[] args) {

        double C = 1.0;
        Range rangeX = new Range(0, 2);
        Range rangeT = new Range(0, 2);

        double c1 = 1.0; // Coefficient for initial value

        Mapper mapper = new Mapper(){
            public double f(double t, double x) {
                return (2*Math.pow(C, 2) * ((C*t - x) / C) * (-x / 2) + x*(2*C*t - x)) / (2*Math.pow(C, 2));
            }
        };

        Shape surface = new SurfaceBuilder().orthonormal(new OrthonormalGrid(rangeX, 200, rangeT, 200), mapper);
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

}
