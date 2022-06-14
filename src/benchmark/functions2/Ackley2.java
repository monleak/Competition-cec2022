package benchmark.functions2;

import benchmark.Function;

public class Ackley2 extends Function {

	@Override
	public double getValue(double[] x) {
		double vars[] = new double[x.length];
		for (int i = 0; i < x.length; i++) {
			vars[i] = x[i] - bias[i];
		}
		double v[] = new double[x.length];
		for (int row = 0; row < x.length; row++) {
			v[row] = 0.0;
			for (int col = 0; col < x.length; col++) {
				v[row] += matrix[row][col] * vars[col];
			}
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
