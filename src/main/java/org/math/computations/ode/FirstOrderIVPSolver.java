package org.math.computations.ode;

import org.math.computations.PlanePoint;
import org.math.computations.functions.TransformR2;

import java.util.ArrayList;
import java.util.List;

public class FirstOrderIVPSolver {

	private final double x0;
	private final double y0;
	private final double step;
	private final int N;

	public FirstOrderIVPSolver(double x0, double y0, double step, int pointsAmount) {

		this.x0 = x0;
		this.y0 = y0;
		this.N = pointsAmount;
		this.step = step;
	}

	public List<PlanePoint> solveExtrapolationAdams(List<PlanePoint> fValues, TransformR2 f) {

		double[][] finiteDiff = new double[N + 3][7];
		for (int i = -2; i <= N; i++) {
			finiteDiff[i + 2][0] = x0 + step*i;
		}
		for (int i = 0; i < fValues.size(); i++) {
			finiteDiff[i][1] = fValues.get(i).getY();
			finiteDiff[i][2] = step*f.evaluate(finiteDiff[i][0], finiteDiff[i][1]);
		}
		for (int i = 3; i < 7; i++) {
			for (int j = 0; j < 7 - i; j++) {
				finiteDiff[j][i] = finiteDiff[j + 1][i - 1] - finiteDiff[j][i - 1];
			}
		}

		List<PlanePoint> result = new ArrayList<>(fValues);
		for (int i = 5; i <= N + 2; i++) {
			finiteDiff[i][1] = finiteDiff[i - 1][1] + finiteDiff[i - 1][2] + finiteDiff[i - 2][3] / 2.0 +
					5*finiteDiff[i - 3][4] / 12.0 + 3*finiteDiff[i - 4][5] / 8.0 + 251*finiteDiff[i - 5][6] / 720.0;
			finiteDiff[i][2] = step*f.evaluate(finiteDiff[i][0], finiteDiff[i][1]);
			finiteDiff[i - 1][3] = finiteDiff[i][2] - finiteDiff[i - 1][2];
			finiteDiff[i - 2][4] = finiteDiff[i - 1][3] - finiteDiff[i - 2][3];
			finiteDiff[i - 3][5] = finiteDiff[i - 2][4] - finiteDiff[i - 3][4];
			finiteDiff[i - 4][6] = finiteDiff[i - 3][5] - finiteDiff[i - 4][5];

			result.add(new PlanePoint(finiteDiff[i][0], finiteDiff[i][1]));
		}

		return result;
	}

	public List<PlanePoint> solveRungeKutta(TransformR2 f) {

		double x = x0;
		double y = y0;
		List<PlanePoint> result = new ArrayList<>(List.of(new PlanePoint(x0, y0)));

		for (int i = 1; i <= N; i++) {

			double k1 = step * f.evaluate(x, y);
			double k2 = step * f.evaluate(x + step / 2, y + k1 / 2);
			double k3 = step * f.evaluate(x + step / 2, y + k2 / 2);
			double k4 = step * f.evaluate(x + step, y + k3);

			y += (k1 + 2*k2 + 2*k3 + k4) / 6;
			x += step;
			result.add(new PlanePoint(x, y));
		}

		return result;
	}

	public List<PlanePoint> solveEuler(TransformR2 f) {

		List<PlanePoint> result = new ArrayList<>(List.of(new PlanePoint(x0, y0)));
		double x = x0;
		double y = y0;

		for (int i = 0; i < N; i++) {

			y += step * f.evaluate(x, y);
			x += step;
			result.add(new PlanePoint(x, y));
		}

		return result;
	}

	public List<PlanePoint> solveEuler1(TransformR2 f) {

		List<PlanePoint> result = new ArrayList<>(List.of(new PlanePoint(x0, y0)));
		double x = x0;
		double y = y0;

		for (int i = 0; i < N; i++) {

			double yCenter = y + step * f.evaluate(x, y) / 2;
			y += step * f.evaluate(x + step / 2, yCenter);
			x += step;
			result.add(new PlanePoint(x, y));
		}

		return result;
	}

	public List<PlanePoint> solveEuler2(TransformR2 f) {

		List<PlanePoint> result = new ArrayList<>(List.of(new PlanePoint(x0, y0)));
		double x = x0;
		double y = y0;

		for (int i = 0; i < N; i++) {

			double Y = y + step * f.evaluate(x, y);
			y += step * (f.evaluate(x, y) + f.evaluate(x + step, Y)) / 2;
			x += step;
			result.add(new PlanePoint(x, y));
		}

		return result;
	}


}
