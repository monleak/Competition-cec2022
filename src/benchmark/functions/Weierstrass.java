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
public class Weierstrass extends Function {

	@Override
	public double getValue(double[] x) {
		double v[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i] - bias[i];
		}

		double a = 0.5;
		double b = 3;
		int kmax = 20;
		double obj = 0;
		int D = x.length;
		for (int i = 1; i <= D; i++) {
			for (int k = 0; k <= kmax; k++) {
				obj = obj + Math.pow(a, k) * Math.cos(2 * Math.PI * Math.pow(b, k) * (v[i - 1] + 0.5));
			}
		}
		for (int k = 0; k <= kmax; k++) {
			obj = obj - D * Math.pow(a, k) * Math.cos(2 * Math.PI * Math.pow(b, k) * 0.5);
		}
		return obj;
	}

}
