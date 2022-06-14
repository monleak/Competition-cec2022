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
public class Schwefel extends Function {

	@Override
	public double getValue(double[] x) {
		double v[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i] - bias[i];
		}

		int dim = x.length;
		double sum = 0;
		for (int i = 0; i < x.length; i++) {
			sum = sum + v[i] * Math.sin(Math.sqrt(Math.abs(v[i])));
		}
		return 418.9829 * dim - sum;

	}

}
