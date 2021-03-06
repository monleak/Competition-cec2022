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
public class Griewank extends Function {

	@Override
    public double getValue(double[] x) {
        double v[] = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            v[i] = x[i] - bias[i];
        }

        double sum1 =0;
        double sum2 = 1;
        for (int i = 0; i < x.length; i++) {
            sum1+= v[i]*v[i];
            sum2 = sum2 * Math.cos(v[i]/(Math.sqrt(i+1)));
        }
        return 1+1.0/4000*sum1-sum2;
    }

    
}
