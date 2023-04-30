package org.math.computational.matrices;

import java.util.ArrayList;
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
		this.A = (DenseMatrix) A.transpose();
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

		List<Double> x = new ArrayList<>();

		double det = A.det();

		if (Math.abs(det) < 10e-16) {
			throw new RuntimeException("Finding unique solution for this system is impossible because det(A) = 0");
		}

		for (int i = 0; i < A.rowCount ; i++) {
			double[] tmp = data[i];
			data[i] = b;
			DenseMatrix T = new DenseMatrix(A.rowCount, A.colCount, data);
			double cramerDet = T.det();
			x.add(cramerDet / det);
			data[i] = tmp;
		}

		return x;
	}

}
