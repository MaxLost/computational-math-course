package org.math.computations.matrices;


import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;


public class SparseMatrix implements Matrix {

  private final HashMap<Integer, HashMap<Integer, Double>> data;
  private final int rowCount;
  private final int colCount;
  private final int hashCode;

  /**
   * Loads sparse matrix from file
   *
   * @param fileName - path to file with matrix data
   */
  public SparseMatrix(String fileName) {
    Path file = Paths.get("src/resources/" + fileName);
    try (Scanner scanner = new Scanner(file)) {
      ArrayList<String> rows = new ArrayList<>();

      while (scanner.hasNextLine()) {
        rows.add(scanner.nextLine());
      }

      if (rows.size() == 0 || rows.get(0).split(" ").length == 0) {
        this.rowCount = 0;
        this.colCount = 0;
        this.data = null;
      } else {
        this.rowCount = rows.size();
        this.colCount = rows.get(0).split(" ").length;

        HashMap<Integer, HashMap<Integer, Double>> data = new HashMap<>();

        for (int i = 0; i < this.rowCount; i++) {
          String[] numbers = rows.get(i).split(" ");
          for (int j = 0; j < this.colCount; j++) {
            double value = Double.parseDouble(numbers[j]);
            if (Math.abs(value) > EPSILON) {
              data.computeIfAbsent(i, t -> new HashMap<Integer, Double>());
              data.get(i).put(j, value);
            }
          }
        }

        this.data = data;
      }
      this.hashCode = this.hashCode();

    } catch (IOException e) {
      throw new RuntimeException("Cannot open file", e);
    }
  }

  public SparseMatrix(int rowCount, int colCount, HashMap<Integer, HashMap<Integer, Double>> data) {
    this.data = data;
    this.rowCount = rowCount;
    this.colCount = colCount;
    this.hashCode = hashCode();
  }

  @Override
  public double getElement(int x, int y) {
    if (x >= this.colCount || y >= this.rowCount) {
      throw new RuntimeException("Invalid coordinates");
    } else {
      HashMap<Integer, Double> row = this.data.get(y);
      if (row == null) {
        return 0;
      }
      Double value = row.get(x);
      return value == null ? 0 : value;
    }
  }

  @Override
  public int[] getSize() {
    return new int[]{rowCount, colCount};
  }

  /**
   * Single-thread matrix multiplication
   * <p>
   * (1) A@B = C
   *
   * @param o - B matrix in (1)
   * @return - result of matrix multiplication, C matrix in (1)
   */
  @Override
  public Matrix mul(Matrix o) {
    if (o instanceof SparseMatrix) {
      return mulSparse((SparseMatrix) o);
    } else if (o instanceof DenseMatrix) {
      return mulDense((DenseMatrix) o);
    }

    return null;
  }

  private Matrix mulDense(DenseMatrix m) {

    if (this.colCount == m.getSize()[0]) {
      if (this.rowCount == 0 || m.getSize()[1] == 0) {
        return new SparseMatrix(0, 0, null);
      }

      HashMap<Integer, HashMap<Integer, Double>> data = new HashMap<>();

      for (Map.Entry<Integer, HashMap<Integer, Double>> row : this.data.entrySet()) {
        for (int column = 0; column < m.getSize()[1]; column++) {
          double sum = 0;
          for (Map.Entry<Integer, Double> element : row.getValue().entrySet()) {
            sum += element.getValue() * m.getElement(column, element.getKey());
          }
          if (Math.abs(sum) > EPSILON) {
            data.computeIfAbsent(row.getKey(), t -> new HashMap<>());
            data.get(row.getKey()).put(column, sum);
          }
        }
      }
      return new SparseMatrix(this.rowCount, m.getSize()[1], data);

    } else {
      throw new RuntimeException("Unable to multiply matrices due to wrong sizes");
    }
  }

  private Matrix mulSparse(SparseMatrix m) {

    if (this.colCount == m.rowCount) {
      if (this.rowCount == 0 || m.colCount == 0) {
        return new SparseMatrix(0, 0, null);
      }
      HashMap<Integer, HashMap<Integer, Double>> data = new HashMap<>();

      SparseMatrix m1 = (SparseMatrix) m.transpose();

      for (Map.Entry<Integer, HashMap<Integer, Double>> row : this.data.entrySet()) {
        for (Map.Entry<Integer, HashMap<Integer, Double>> column : m1.data.entrySet()) {
          double sum = 0;
          for (Map.Entry<Integer, Double> element : row.getValue().entrySet()) {
            sum += element.getValue() * m1.getElement(element.getKey(), column.getKey());
          }
          if (Math.abs(sum) > EPSILON) {
            data.computeIfAbsent(row.getKey(), t -> new HashMap<>());
            data.get(row.getKey()).put(column.getKey(), sum);
          }
        }
      }
      return new SparseMatrix(this.rowCount, m.colCount, data);

    } else {
      throw new RuntimeException("Unable to multiply matrices due to wrong sizes");
    }
  }

  /**
   * Multi-thread matrix multiplication
   * <p>
   * (1) A@B = C
   *
   * @param o - B matrix in (1)
   * @return - result of matrix multiplication, C matrix in (1)
   */
  @Override
  public Matrix dmul(Matrix o) {

    if (o instanceof SparseMatrix) {
      SparseMatrix m = (SparseMatrix) ((SparseMatrix) o).transpose();
      SparseMatrix n = this;

      if (this.colCount == m.colCount) {
        if (this.rowCount == 0 || m.rowCount == 0) {
          return new SparseMatrix(0, 0, null);
        }

        MulTaskManager task_manager = new MulTaskManager(n.rowCount);
        ConcurrentHashMap<Integer, HashMap<Integer, Double>> data = new ConcurrentHashMap<>();

        class Multiplicator implements Runnable {

          @Override
          public void run() {
            Integer row = task_manager.next();
            while (row != null) {
              if (n.data == null || m.data == null || n.data.get(row) == null) {
                row = task_manager.next();
                continue;
              }
              HashMap<Integer, Double> result = new HashMap<>();

              for (Map.Entry<Integer, HashMap<Integer, Double>> column : m.data.entrySet()) {
                double sum = 0;
                for (Map.Entry<Integer, Double> element : n.data.get(row).entrySet()) {
                  sum += element.getValue() * m.getElement(element.getKey(), column.getKey());
                }
                if (Math.abs(sum) > EPSILON) {
                  result.put(column.getKey(), sum);
                }
              }
              if (!result.isEmpty()) {
                data.put(row, result);
              }
              row = task_manager.next();
            }
          }
        }

        Thread[] threads = new Thread[4];
        for (int i = 0; i < threads.length; i++) {
          threads[i] = new Thread(new Multiplicator());
          threads[i].start();
        }
        for (Thread thread : threads) {
          try {
            thread.join();
          } catch (InterruptedException e) {
            throw new RuntimeException("Multiplication failed! Try again!", e);
          }
        }

        return new SparseMatrix(n.rowCount, m.rowCount, new HashMap<>(data));
      } else {
        throw new RuntimeException("Unable to multiply matrices due to wrong sizes");
      }
    } else {
      return this.mul(o);
    }
  }

  @Override
  public Matrix add(Matrix o) {

    int[] matrixSize = o.getSize();

    if (this.rowCount != matrixSize[0] || this.colCount != matrixSize[1]) {
      throw new RuntimeException("Addition supported only for matrices with equal sizes.");
    }

    double[][] matrixA = this.toArray();

    for (int i = 0; i < matrixSize[0]; i++) {
      for (int j = 0; j < matrixSize[1]; j++) {
        matrixA[i][j] += o.getElement(j, i);
      }
    }

    return new DenseMatrix(matrixSize[0], matrixSize[1], matrixA);
  }

  @Override
  public Matrix scalarMultiply(double n) {

    double[][] matrixA = this.toArray();
    for (int i = 0; i < rowCount; i++) {
      for (int j = 0; j < colCount; j++) {
        matrixA[i][j] *= n;
      }
    }

    return new DenseMatrix(rowCount, colCount, matrixA);
  }

  public Matrix transpose() {
    if (this.rowCount == 0 || this.colCount == 0) {
      return this;
    } else {
      HashMap<Integer, HashMap<Integer, Double>> data = new HashMap<>();
      for (Integer i : this.data.keySet()) {
        for (Integer j : this.data.get(i).keySet()) {
          data.computeIfAbsent(j, t -> new HashMap<Integer, Double>());
          data.get(j).put(i, this.getElement(j, i));
        }
      }
      return new SparseMatrix(this.colCount, this.rowCount, data);
    }
  }

  @Override
  public double[][] toArray() {
    double[][] result = new double[rowCount][colCount];
    for (int i = 0; i < rowCount; i++) {
      for (int j = 0; j < colCount; j++) {
        result[i][j] = this.getElement(j, i);
      }
    }
    return result;
  }

  @Override
  public int hashCode() {

    String caller = String.valueOf((new Throwable().getStackTrace())[1]);
    if (caller.equals("SparseMatrix")) {

      if (this.rowCount == 0 | this.colCount == 0) {
        return 0;
      }

      int a = 0, b = 0;
      for (int i = 0; i < this.colCount; i++) {
        HashMap<Integer, Double> row = this.data.get(i);
        if (row != null) {

          if (row.values().iterator().hasNext()) {
            a += row.values().iterator().next().intValue();
          }
          if (row.values().iterator().hasNext()) {
            b += row.values().iterator().next().intValue();
          }
        }
      }
      return a ^ b;
    }
    return this.hashCode;
  }

  @Override
  public boolean equals(Object o) {

    if (this == o) {
      return true;
    }

    if (o instanceof SparseMatrix) {
      return equalsSparse((SparseMatrix) o);
    }

    if (o instanceof DenseMatrix) {
      return equalsDense((DenseMatrix) o);
    }

    return false;
  }

  private boolean equalsDense(DenseMatrix o) {

    if (this.rowCount != o.getSize()[0] || this.colCount != o.getSize()[1]) {
      return false;
    }

    if (this.rowCount == 0 || this.colCount == 0) {
      return true;
    }

    for (int i = 0; i < this.rowCount; i++) {
      for (int j = 0; j < this.colCount; j++) {
        if (Math.abs(Math.abs(o.getElement(j, i)) - Math.abs(this.getElement(j, i))) > EPSILON) {
          return false;
        }
      }
    }

    return true;
  }

  private boolean equalsSparse(SparseMatrix o) {

    if (this.rowCount != o.rowCount | this.colCount != o.colCount) {
      return false;
    }

    if (this.hashCode() == o.hashCode()) {

      if (this.rowCount == 0 || this.colCount == 0) {
        return true;
      }

      for (int i = 0; i < this.rowCount; i++) {
        for (int j = 0; j < this.colCount; j++) {
          if (Math.abs(Math.abs(o.getElement(j, i)) - Math.abs(this.getElement(j, i))) > EPSILON) {
            return false;
          }
        }
      }

      return true;
    }
    return false;
  }

  @Override
  public String toString() {

    StringBuilder str = new StringBuilder();
    for (int i = 0; i < this.rowCount; i++) {
      for (int j = 0; j < this.colCount; j++) {
        str.append(this.getElement(j, i)).append(" ");
      }
      str.append("\n");
    }
    return str.toString();
  }
}
