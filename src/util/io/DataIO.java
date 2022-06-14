package util.io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * @author cuonglv.hust@gmail.com
 * @date 25/02/2021
 *
 */
public class DataIO {

    /**
     * Read the shift file
     *
     * @param file
     * @param dim
     * @return the bias if reading file successfully, return null in the otherwise
     */
    public static double[] readBias(String file, int dim) {
        double[] bias = null;
        try {
            Scanner in = new Scanner(new File(file));
            bias = new double[dim];
            for (int i = 0; i < dim; i++) {
                bias[i] = in.nextDouble();
            }
            in.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return bias;
    }

    /**
     * Read the rotation matrix
     *
     * @param file
     * @param dim
     * @return the matrix if reading file successfully, return null in the otherwise
     */
    public static double[][] readMatrix(String file, int dim) {
        double[][] matrix = null;
        try {
            Scanner in = new Scanner(new File(file));
            matrix = new double[dim][dim];
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    matrix[i][j] = in.nextDouble();
                }
            }
            in.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return matrix;
    }

    public static void main(String[] args) {
        double bias[] = DataIO.readBias("data\\WCCI-Competition\\SO-Manytask-Benchmarks\\Tasks\\benchmark_1\\bias_1",
                50);
        double[][] matrix = DataIO
                .readMatrix("data\\WCCI-Competition\\SO-Manytask-Benchmarks\\Tasks\\benchmark_1\\matrix_1", 50);

        for (double d : bias) {
            System.out.println(d);
        }

        for (double[] row : matrix) {
            for (double item : row) {
                System.out.print(item + "\t");
            }
            System.out.println();
        }
    }

}
