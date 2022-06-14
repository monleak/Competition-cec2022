package benchmark.functions2;

import benchmark.Function;

public class Sphere2 extends Function {
	
	public Sphere2() {
		
	}
	
	public Sphere2(int dim, double[] bias, double[][] matrix, double[] lb, double[] ub) {
		super(dim, bias, matrix, lb, ub);
	}

	@Override
	public double getValue(double[] x) {
		double vars[] = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            vars[i] = x[i] - bias[i];
        }
        double v[] = new double[x.length];
        for (int row = 0; row < x.length; row++) {
            v[row] =0.0;
            for (int col = 0; col < x.length; col++) {
                    v[row]+= matrix[row][col]*vars[col];
            }
        }
        double sum = 0;
        for(double d : v){
            sum += d*d;
        }
        return sum;
	}

}
