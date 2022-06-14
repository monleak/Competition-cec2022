package benchmark.functions2;

import benchmark.Function;

public class Griewank2 extends Function {

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
		double sum2 = 1;
		for (int i = 0; i < x.length; i++) {
			sum1 += v[i] * v[i];
			sum2 = sum2 * Math.cos(v[i] / (Math.sqrt(i + 1)));
		}
		return 1 + 1.0 / 4000 * sum1 - sum2;
	}

}
