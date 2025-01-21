package org.math.computations.pde;

import org.math.computations.PlanePoint;
import org.math.computations.functions.FunctionalR2;

public class Utils {

    public static PlanePoint[][] generate2dMesh(double[] boundaries, int[] stepAmounts) {

        if (boundaries.length != 4) {
            throw new RuntimeException("Incorrect amount of boundaries. Should be exactly 4.");
        }
        if (stepAmounts.length != 2) {
            throw new RuntimeException("Incorrect amount of step amounts. Should be exactly 2.");
        }

        double stepX = (boundaries[1] - boundaries[0]) / stepAmounts[0];
        double stepY = (boundaries[3] - boundaries[2]) / stepAmounts[1];

        PlanePoint[][] mesh = new PlanePoint[stepAmounts[0] + 1][stepAmounts[1] + 1];
        double pointerX = boundaries[0];

        for (int i = 0; i <= stepAmounts[0]; i++) {
            double pointerY = boundaries[2];
            for (int j = 0; j <= stepAmounts[1]; j++) {
                mesh[i][j] = new PlanePoint(pointerX, pointerY);
                pointerY += stepY;
            }
            pointerX += stepX;
        }

        return mesh;
    }

    public static double[][] getValuesOnMesh(FunctionalR2 f, PlanePoint[][] mesh) {

        double[][] values = new double[mesh.length][mesh[0].length];

        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values[0].length; j++) {
                values[i][j] = f.evaluate(mesh[i][j].getX(), mesh[i][j].getY());
            }
        }

        return values;
    }

}
