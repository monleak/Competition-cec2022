package benchmark;

import java.util.ArrayList;

import basic.Params;
import basic.Problem;
import basic.Task;
import benchmark.functions.Ackley;
import benchmark.functions.Griewank;
import benchmark.functions.Rastrigin;
import benchmark.functions.Rosenbrock;
import benchmark.functions.Schwefel;
import benchmark.functions.Sphere;
import benchmark.functions.Weierstrass;
import benchmark.functions2.Ackley2;
import benchmark.functions2.Griewank2;
import benchmark.functions2.Rastrigin2;
import benchmark.functions2.Rosenbrock2;
import benchmark.functions2.Schwefel2;
import benchmark.functions2.Sphere2;
import benchmark.functions2.Weierstrass2;
import util.io.DataIO;

/**
 * @author cuonglv.hust@gmail.com
 * @date 25/02/2021
 *
 */
public class ProblemConstructor {

    public static ArrayList<Problem> getComplexProblems() {
        ArrayList<Problem> problems = new ArrayList<Problem>();
        final int problems_num = 10;

        for (int index = 1; index <= problems_num; index++) {
            Problem prob = new Problem();
            int dim = 50;
            Task task_1;
            Task task_2;
            switch (index) {
                case 1:
                    int func_id = 6;
                    Function func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    Function func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 2:
                    func_id = 7;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 3:
                    func_id = 17;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 4:
                    func_id = 13;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 5:
                    func_id = 15;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 6:
                    func_id = 21;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 7:
                    func_id = 22;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 8:
                    func_id = 5;
                    func1 = new Cec14_func(func_id, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(func_id, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 9:
                    func1 = new Cec14_func(11, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(16, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                case 10:
                    func1 = new Cec14_func(20, index, 1, dim);
                    func1.dim = dim;
                    func1.LB = ones(dim, -100);
                    func1.UB = ones(dim, 100);
                    task_1 = new Task(1, func1);

                    func2 = new Cec14_func(21, index, 2, dim);
                    func2.dim = dim;
                    func2.LB = ones(dim, -100);
                    func2.UB = ones(dim, 100);
                    task_2 = new Task(2, func2);

                    prob.addTask(task_1);
                    prob.addTask(task_2);
                    break;
                default:

            }
            ;

            problems.add(prob);
        }

        System.out.println("Complex benchmark constructed successfully!");
        return problems;

    }

    public static ArrayList<Problem> get10TasksBenchmark() {
        ArrayList<Problem> list_prob = new ArrayList<Problem>();
        Problem prob = new Problem();
        int task_size = 10;

        for (int task_id = 1; task_id <= task_size; task_id++) {
            Task task = new Task(task_id);
            Function func;

            switch (task_id) {
                case 1: // sphere 1
                    func = new Sphere();
                    func.dim = 50;
                    func.LB = ones(func.dim, -100);
                    func.UB = ones(func.dim, 100);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = 0;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 2:
                    func = new Sphere();
                    func.dim = 50;
                    func.LB = ones(func.dim, -100);
                    func.UB = ones(func.dim, 100);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = 80;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 3:
                    func = new Sphere();
                    func.dim = 50;
                    func.LB = ones(func.dim, -100);
                    func.UB = ones(func.dim, 100);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = -80;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 4:
                    func = new Weierstrass();
                    func.dim = 25;
                    func.LB = ones(func.dim, -0.5);
                    func.UB = ones(func.dim, 0.5);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = -0.4;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 5:
                    func = new Rosenbrock();
                    func.dim = 50;
                    func.LB = ones(func.dim, -50);
                    func.UB = ones(func.dim, 50);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = -1;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 6:
                    func = new Ackley();
                    func.dim = 50;
                    func.LB = ones(func.dim, -50);
                    func.UB = ones(func.dim, 50);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = 40;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 7:
                    func = new Weierstrass();
                    func.dim = 50;
                    func.LB = ones(func.dim, -0.5);
                    func.UB = ones(func.dim, 0.5);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = -0.4;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 8:
                    func = new Schwefel();
                    func.dim = 50;
                    func.LB = ones(func.dim, -500);
                    func.UB = ones(func.dim, 500);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);
                    for (int i = 0; i < func.dim; i++) {
                        func.bias[i] = 0;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 9:
                    func = new Griewank();
                    func.dim = 50;
                    func.LB = ones(func.dim, -100);
                    func.UB = ones(func.dim, 100);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);

                    int k = func.dim / 2;
                    for (int i = 0; i < k; i++) {
                        func.bias[i] = -80;
                    }
                    for (int i = k; i < func.dim; i++) {
                        func.bias[i] = 80;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                case 10:
                    func = new Rastrigin();
                    func.dim = 50;
                    func.LB = ones(func.dim, -50);
                    func.UB = ones(func.dim, 50);
                    func.bias = new double[func.dim];
                    func.matrix = i_matrix(func.dim);

                    int m = func.dim / 2;
                    for (int i = 0; i < m; i++) {
                        func.bias[i] = 40;
                    }
                    for (int i = m; i < func.dim; i++) {
                        func.bias[i] = -40;
                    }

                    task.setFunction(func);
                    prob.addTask(task);
                    break;
                default:
                    System.err.println("Invalid function");
            }
            ;
        }

        list_prob.add(prob);

        System.out.println("10 tasks benchmark constructed successfully!");

        return list_prob;
    }

    public static ArrayList<Problem> get50TasksBenchmark() {
        ArrayList<Problem> problems = new ArrayList<Problem>();
        final int problems_num = 10;

        for (int index = 1; index <= problems_num; index++) {
            Problem problem = new Problem();
            int tasks_num = 50;
            int dim = 50;
            int choice_functions[] = null;
            switch (index) {
                case 1:
                    choice_functions = new int[] { 1 };
                    break;
                case 2:
                    choice_functions = new int[] { 2 };
                    break;
                case 3:
                    choice_functions = new int[] { 4 };
                    break;
                case 4:
                    choice_functions = new int[] { 1, 2, 3 };
                    break;
                case 5:
                    choice_functions = new int[] { 4, 5, 6 };
                    break;
                case 6:
                    choice_functions = new int[] { 2, 5, 7 };
                    break;
                case 7:
                    choice_functions = new int[] { 3, 4, 6 };
                    break;
                case 8:
                    choice_functions = new int[] { 2, 3, 4, 5, 6 };
                    break;
                case 9:
                    choice_functions = new int[] { 2, 3, 4, 5, 6, 7 };
                    break;
                case 10:
                    choice_functions = new int[] { 3, 4, 5, 6, 7 };
                    break;
                default:
                    System.out.println("Invalid input: ID should be in [1,10]");
            }
            ;

            for (int task_id = 1; task_id <= tasks_num; task_id++) {
                int func_id = choice_functions[(task_id - 1) % choice_functions.length];

                String file_dir = Params.MANY_TASKS_BENCHMARKS_PATH + index;
                String file_matrix = file_dir + "/matrix_" + task_id;
                String file_bias = file_dir + "/bias_" + task_id;

                // read rotation matrix and shift file
                double matrix[][] = DataIO.readMatrix(file_matrix, dim);
                double bias[] = DataIO.readBias(file_bias, dim);

                Task task = new Task(task_id);
                Function func;
                switch (func_id) {
                    case 1:
                        func = new Sphere2();
                        func.dim = dim;
                        func.LB = ones(dim, -100);
                        func.UB = ones(dim, 100);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    case 2:
                        func = new Rosenbrock2();
                        func.dim = dim;
                        func.LB = ones(dim, -50);
                        func.UB = ones(dim, 50);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    case 3:
                        func = new Ackley2();
                        func.dim = dim;
                        func.LB = ones(dim, -50);
                        func.UB = ones(dim, 50);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    case 4:
                        func = new Rastrigin2();
                        func.dim = dim;
                        func.LB = ones(dim, -50);
                        func.UB = ones(dim, 50);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    case 5:
                        func = new Griewank2();
                        func.dim = dim;
                        func.LB = ones(dim, -100);
                        func.UB = ones(dim, 100);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    case 6:
                        func = new Weierstrass2();
                        func.dim = dim;
                        func.LB = ones(dim, -0.5);
                        func.UB = ones(dim, 0.5);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    case 7:
                        func = new Schwefel2();
                        func.dim = dim;
                        func.LB = ones(dim, -500);
                        func.UB = ones(dim, 500);
                        func.matrix = matrix;
                        func.bias = bias;

                        task.setFunction(func);
                        problem.addTask(task);
                        break;
                    default:
                        System.out.println("Invalid function");
                }
                ;
            }

            problems.add(problem);
        }
        System.out.println("50 tasks benchmark constructed successfully!");
        return problems;
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
