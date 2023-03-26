package org.math.computational.functions;

import org.math.computational.PlanePoint;

import java.util.ArrayList;
import java.util.List;

public class NumericalDifferentiator {

	private final List<PlanePoint> nodes;
	private final double h;

	public NumericalDifferentiator(List<PlanePoint> nodes) {

		this.nodes = nodes;
		this.h = nodes.get(1).getX() - nodes.get(0).getX();
	}

	public List<PlanePoint> computeDerivatives () {

		List<PlanePoint> result = new ArrayList<>();

		result.add(
				new PlanePoint(
						nodes.get(0).getX(),
						(-3*nodes.get(0).getY() + 4*nodes.get(1).getY() - nodes.get(2).getY()) / (2*h)
				)
		);

		for (int i = 1; i < nodes.size() - 1; i++) {
			result.add(
					new PlanePoint(
							nodes.get(i).getX(),
							(nodes.get(i + 1).getY() - nodes.get(i - 1).getY()) / (2*h)
					)
			);
		}

		int k = nodes.size() - 1;

		result.add(
				new PlanePoint(
						nodes.get(k).getX(),
						(3*nodes.get(k).getY() - 4*nodes.get(k - 1).getY() + nodes.get(k - 2).getY()) / (2*h)
				)
		);

		return result;
	}

	public List<PlanePoint> computeSecondDerivatives () {

		List<PlanePoint> result = new ArrayList<>();

		for (int i = 1; i < nodes.size() - 1; i++) {
			result.add(
					new PlanePoint(
							nodes.get(i).getX(),
							(nodes.get(i + 1).getY() - 2*nodes.get(i).getY() + nodes.get(i - 1).getY()) / (h * h)
					)
			);
		}

		return result;
	}

}
