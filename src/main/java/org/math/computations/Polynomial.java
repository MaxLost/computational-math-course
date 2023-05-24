package org.math.computations;

import org.math.computations.functions.Function;

import java.util.ArrayList;
import java.util.List;

public class Polynomial implements Function {

	private final ArrayList<Double> coefficients;

	public Polynomial(List<Double> coefficients) {
		this.coefficients = new ArrayList<>(coefficients);
	}

	public ArrayList<Double> getCoefficients() {
		return coefficients;
	}

	public Polynomial multiply(Polynomial g) {

		List<Double> newCoefficients = new ArrayList<>();
		List<Double> gCoefficients = g.getCoefficients();

		for (int i = 0; i < this.coefficients.size() + gCoefficients.size(); i++) {
			newCoefficients.add(0.0);
		}

		for (int i = 0; i < this.coefficients.size(); i++) {
			for (int j = 0; j < gCoefficients.size(); j++) {
				newCoefficients.set(i + j, newCoefficients.get(i + j) + this.coefficients.get(i) * gCoefficients.get(j));
			}
		}

		int deg = 0;
		for (int i = 0; i < newCoefficients.size(); i++) {
			if (newCoefficients.get(i) != null && Math.abs(newCoefficients.get(i)) > 10e-8) {
				deg = i;
			}
		}
		for (int i = newCoefficients.size() - 1; i > deg; i--) {
			newCoefficients.remove(i);
		}

		return new Polynomial(newCoefficients);
	}

	public Polynomial scalarMultiply(double a) {

		List<Double> newCoefficients = new ArrayList<>();

		for (int i = 0; i < this.coefficients.size(); i++) {
			newCoefficients.add(this.coefficients.get(i) * a);
		}

		return new Polynomial(newCoefficients);
	}

	public Polynomial sum(Polynomial g) {

		List<Double> newCoefficients = new ArrayList<>();
		List<Double> gCoefficients = g.getCoefficients();

		int i = 0;

		while (i < Math.min(this.coefficients.size(), gCoefficients.size())) {

			newCoefficients.add(this.coefficients.get(i) + gCoefficients.get(i));
			i++;
		}
		while (i < this.coefficients.size()) {

			newCoefficients.add(this.coefficients.get(i));
			i++;
		}
		while (i < gCoefficients.size()) {

			newCoefficients.add(gCoefficients.get(i));
			i++;
		}

		return new Polynomial(newCoefficients);
	}

	@Override
	public double evaluate(double point) {

		double value = this.coefficients.get(this.coefficients.size() - 1);

		for (int i = this.coefficients.size() - 2; i >= 0; i--) {
			value = value * point + this.coefficients.get(i);
		}

		return value;
	}

	@Override
	public boolean equals(Object o) {

		if (o instanceof Polynomial p) {
			return this.coefficients.equals(p.getCoefficients());
		}

		return false;
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder("");
		for (int i = this.coefficients.size() - 1; i > 0; i--) {
			str.append(this.coefficients.get(i)).append("*x^").append(i).append(" + ");
		}
		str.append(this.coefficients.get(0));
		return str.toString();
	}
}
