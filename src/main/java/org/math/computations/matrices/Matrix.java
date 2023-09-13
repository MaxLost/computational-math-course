package org.math.computations.matrices;


import java.util.List;

public interface Matrix
{
  double EPSILON = 10e-6;

  int[] getSize();

  double getElement(int x, int y);

  /**
   * Single-thread matrix multiplication
   * <p>
   * (1) A@B = C
   *
   * @param o - B matrix in (1)
   * @return - result of matrix multiplication, C matrix in (1)
   */
  Matrix mul(Matrix o);

  /**
   * Multi-thread matrix multiplication
   * <p>
   * (1) A@B = C
   *
   * @param o - B matrix in (1)
   * @return - result of matrix multiplication, C matrix in (1)
   */
  Matrix dmul(Matrix o);
}
