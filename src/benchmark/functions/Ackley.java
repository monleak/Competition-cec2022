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
public class Ackley extends Function {

	@Override
	public double getValue(double[] x) {
		double v[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			v[i] = x[i] - bias[i];
		}

		double sum1 = 0;
		double sum2 = 0;
		for (int i = 0; i < x.length; i++) {
			sum1 += v[i] * v[i];
			sum2 += Math.cos(2 * Math.PI * v[i]);
		}
		double avgsum1 = sum1 / x.length;
		double avgsum2 = sum2 / x.length;
		return -20 * Math.exp(-0.2 * Math.sqrt(avgsum1)) - Math.exp(avgsum2) + 20 + Math.exp(1);
	}

}
