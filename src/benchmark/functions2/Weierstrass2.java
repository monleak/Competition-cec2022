package benchmark.functions2;

import benchmark.Function;

public class Weierstrass2 extends Function {

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
