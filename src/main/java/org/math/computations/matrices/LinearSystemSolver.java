package org.math.computations.matrices;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class LinearSystemSolver {

	private final DenseMatrix A;
	private final DenseMatrix B;

	/**
	 * Class for finding unique solution X for linear system of equations: A @ X = B
	 * <br>
	 * Where A - N x N matrix, X - N x 1 vector, B - N x 1 vector
	 * @param A N x N matrix
	 * @param B N x 1 vector
	 */
	public LinearSystemSolver (DenseMatrix A, DenseMatrix B) {
		this.A = A;
		this.B = (DenseMatrix) B.transpose();
	}

	/**
	 * Finds unique solution linear system: A @ X = B
	 * <br>
	 * Where A - N x N matrix and det(A) != 0, X - N x 1 vector, B - N x 1 vector
	 * @return X vector as a list
	 */
	public List<Double> solve() {

		double[][] data = A.toArray();
		double[] b = B.toArray()[0];

		Double[] x = new Double[b.length];

		double det = A.det();

		if (Math.abs(det) < 10e-32){
			throw new RuntimeException("Finding unique solution for this system is impossible because det(A) = 0");
		}

		double[][] M = A.toArray();
		double sum = 0;

		for (int i = 0; i < A.colCount - 1; i++) {
			for (int j = i + 1; j < A.rowCount; j++) {
				for (int k = i + 1; k < A.colCount; k++){
					M[j][k] -= (M[i][k] * M[j][i]) / M[i][i];
				}
			}
		}

		x[0] = b[0];
		for (int i = 1; i < A.colCount; i++) {
			sum = 0;
			for (int j = 0; j < i; j++) {
				sum += M[i][j] * x[j];
			}
			x[i] = b[i] - sum;
		}

		for (int i = A.colCount - 1; i >= 0; i--) {
			sum = 0;
			for (int j = i + 1; j < A.colCount; j++) {
				sum += M[i][j] * x[j];
			}
			x[i] = (x[i] - sum) / M[i][i];
		}

		return Arrays.asList(x);
	}

}
