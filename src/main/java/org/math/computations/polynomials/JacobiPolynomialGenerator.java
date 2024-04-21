package org.math.computations.polynomials;

import java.util.ArrayList;
import java.util.List;

public class JacobiPolynomialGenerator {

    public static List<Polynomial> getJacobiPolynomials(double alpha, double beta, int n) {
        List<Polynomial> result = new ArrayList<>();
        result.add(new Polynomial(List.of(1.0)));
        result.add(new Polynomial(List.of(alpha + 1 - (alpha + beta + 2) / 2, (alpha + beta + 2) / 2)));

        Polynomial x = new Polynomial(List.of(0.0, 1.0));

        for (int i = 2; i <= n; i++) {
            double a = (2 * (i - 1) + 1 + alpha + beta) * (2 * (i - 1) + 2 + alpha + beta) /
                (2 * i * (i + alpha + beta));
            double b = (Math.pow(alpha, 2) - Math.pow(beta, 2)) * (2 * (i - 1) + 1 + alpha + beta) /
                (2 * i * (2 * (i - 1) + alpha + beta) * (i + alpha + beta));
            double c = (i - 1 + alpha) * (i - 1 + beta) * (2 * i + alpha + beta) /
                (i * (i + alpha + beta) * (2 * (i - 1) + alpha + beta));

            Polynomial p = result.get(i - 1).multiply(new Polynomial(List.of(b, a)))
                .sum(result.get(i - 2).scalarMultiply(-1 * c));
            result.add(p);
        }

        return result;
    }

}
