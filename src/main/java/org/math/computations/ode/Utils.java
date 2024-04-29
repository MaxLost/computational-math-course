package org.math.computations.ode;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.math.computations.PlanePoint;
import org.math.computations.matrices.DenseMatrix;
import org.math.computations.matrices.Matrix;

public class Utils {

    static Matrix normalizeEuclidean(Matrix matrixA) {

        double mean = org.math.computations.matrices.Utils.euclideanMean(matrixA);
        return matrixA.scalarMultiply(1 / mean);
    }

    static Matrix normalizeVector(Matrix vector) {

        double maxElement = getMaxAbsElementOfVector(vector);

        return vector.scalarMultiply(1 / maxElement);
    }

    static double getMaxAbsElementOfVector(Matrix vector) {
        int[] matrixSize = vector.getSize();

        int maxElementId = 0;
        for (int i = 0; i < matrixSize[0]; i++) {
            if (Math.abs(vector.getElement(0, i)) > vector.getElement(0, maxElementId)) {
                maxElementId = i;
            }
        }

        return vector.getElement(0, maxElementId);
    }

    static List<PlanePoint> generateMesh(double start, double end, int cellAmount) {

        List<PlanePoint> mesh = new ArrayList<>();

        double pointer = start;
        double step = (end - start) / cellAmount;

        mesh.add(new PlanePoint(pointer, 0));
        while (Math.abs(end - pointer) > 1e-6) {
            pointer += step;
            mesh.add(new PlanePoint(pointer, 0));
        }

        return mesh;
    }

    static List<PlanePoint> generateMesh(double start, double end, double step) {

        List<PlanePoint> mesh = new ArrayList<>();

        double pointer = start;

        mesh.add(new PlanePoint(pointer, 0));
        while (Math.abs(end - pointer) > 1e-6) {
            pointer += step;
            mesh.add(new PlanePoint(pointer, 0));
        }

        return mesh;
    }

    static double vectorDotProduct(Matrix vectorA, Matrix vectorB) {

        int vectorSize = vectorA.getSize()[0];
        if (vectorSize != vectorB.getSize()[0]) {
            throw new RuntimeException("Vectors should have same size.");
        }

        double result = 0;
        for (int i = 0; i < vectorSize; i++) {
            result += vectorA.getElement(0, i) * vectorB.getElement(0, i);
        }

        return result;
    }

    static Matrix generateRandomVector(int size) {

        double[][] dataY = new double[size][1];
        Random random = new Random();

        for (int i = 0; i < size; i++) {
            dataY[i][0] = random.nextDouble();
        }

        return new DenseMatrix(size, 1, dataY);
    }
}
