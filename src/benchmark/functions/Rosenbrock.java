package benchmark.functions;

import benchmark.Function;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author thang.tb153544
 */
public class Rosenbrock extends Function {

	@Override
	public double getValue(double[] x) {
		double v[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i] - bias[i];
		}

		double sum = 0;
		int dim = x.length;
		for (int i = 0; i < dim - 1; i++) {
			double xi = v[i];
			double xnext = v[i + 1];
			double _new = 100 * (xnext - xi * xi) * (xnext - xi * xi) + (xi - 1) * (xi - 1);
			sum = sum + _new;
		}
		return sum;
	}

}
