package util;

import basic.Params;

public class Utils {

    public static double cauchy_g(double mu, double gamma) {
        return mu + gamma * Math.tan(Math.PI * (Params.rand.nextDouble() - 0.5));
    }

    public static double gauss(double mu, double sigma) {
        return mu + sigma * Math.sqrt(-2.0 * Math.log(Params.rand.nextDouble()))
                * Math.sin(2.0 * Math.PI * Params.rand.nextDouble());
    }

    public static double[] ones(int n, double scale) {
        double[] matrix = new double[n];
        for (int i = 0; i < n; i++) {
            matrix[i] = 1 * scale;
        }
        return matrix;
    }

    public static double[][] i_matrix(int n) {
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = (int) 1;
                } else {
                    matrix[i][j] = matrix[j][i] = (int) 0;
                }
            }

        }
        return matrix;
    }

}
