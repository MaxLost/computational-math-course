package org.math.computations.matrices;

/**
 * Class for representation of matrices.
 */
public interface Matrix {
  double EPSILON = 10e-6;

  int[] getSize();

  double getElement(int x, int y);

  /**
   * Single-thread matrix multiplication.
   * <br>
   * (1) A@B = C
   *
   * @param o B matrix in (1)
   * @return result of matrix multiplication, C matrix in (1)
   */
  Matrix mul(Matrix o);

  /**
   * Multi-thread matrix multiplication.
   * <br>
   * (1) A@B = C
   *
   * @param o B matrix in (1)
   * @return result of matrix multiplication, C matrix in (1)
   */
  Matrix dmul(Matrix o);

  Matrix add(Matrix o);

  Matrix scalarMultiply(double n);

  double[][] toArray();

  Matrix copy();
}
