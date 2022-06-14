package benchmark.functions2;

import benchmark.Function;

public class Schwefel2 extends Function {

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
		int dim = x.length;
		double sum = 0;
		for (int i = 0; i < x.length; i++) {
			sum = sum + v[i] * Math.sin(Math.sqrt(Math.abs(v[i])));
		}
		return 418.9829 * dim - sum;

	}

}
