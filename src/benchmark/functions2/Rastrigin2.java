package benchmark.functions2;

import benchmark.Function;

public class Rastrigin2 extends Function {

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
		double obj = dim * 10;
		for (int i = 0; i < dim; i++) {
			obj = obj + (v[i] * v[i] - 10 * (Math.cos(2 * Math.PI * v[i])));
		}
		return obj;
	}

}
