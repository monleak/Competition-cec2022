package ls;

import basic.Individual;
import basic.Task;

/**
 * The strategy of Davies, Swann, and Campey with Gram-Schmidt orthogonalization.
 *
 * https://www.researchgate.net/publication/220690578_Evolution_and_Optimum_Seeking
 *
 * @author cuonglv.hust@gmail.com
 *
 */
public class DSCG implements LocalSearch {

    private final double INIT_STEP_SIZE = 0.02;
    private final double EPSILON = 1E-8;
    private final int EVALS_PER_LINE_SEARCH = 50;

    @Override
    public Individual search(Individual start_point, int fes, Task task) {
        return runDSCG(start_point, fes, task);
    }

    private Individual runDSCG(Individual start_point, int fes, Task task) {
        double s = INIT_STEP_SIZE;
        int evals_per_linesearch = EVALS_PER_LINE_SEARCH;
        int dim = start_point.getDim();

        Individual result = new Individual(start_point.getChromosome());
        result.setFitness(start_point.getFitness());

        Individual[] x = new Individual[dim + 2];
        for (int i = 0; i < x.length; i++) {
            x[i] = new Individual(dim);
        }

        for (int i = 0; i < dim; i++) {
            x[0].setGene(i, start_point.getGene(i));
        }
        x[0].setFitness(start_point.getFitness());

        int direct = 1, evals = 0;
        double[][] v = new double[dim + 1][dim];
        double[][] a = new double[dim][dim];
        for (int i = 0; i <= dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (i == j) {
                    v[i][j] = 1;
                } else {
                    v[i][j] = 0;
                }
            }
        }

        while (true) {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < dim; j++) {
                    a[i][j] = 0;
                }
            }

            // line search
            while (evals < fes - evals_per_linesearch) {
                evals += lineSearch(x[direct - 1], x[direct], this.EVALS_PER_LINE_SEARCH, task, s, v[direct - 1]);

                // x_fit[direct - 1] >= x_fit[direct]
                for (int i = 1; i <= direct; i++) {
                    for (int j = 0; j < dim; j++) {
                        a[i - 1][j] += x[direct].getGene(j) - x[direct - 1].getGene(j);
                    }
                }

                // update best
                if (result.getFitness() > x[direct].getFitness()) {
                    result.setFitness(x[direct].getFitness());
                    result.setChromosome(x[direct].getChromosome().clone());
                }

                if (direct < dim) {
                    direct++;
                } else {
                    break;
                }
            }

            if (evals >= fes || direct < dim) {
                break;
            }

            // Eventually one more line search
            double z[] = new double[dim];
            double norm_z = 0;
            for (int i = 0; i < dim; i++) {
                z[i] = x[dim].getGene(i) - x[0].getGene(i);
                norm_z += z[i] * z[i];
            }
            norm_z = Math.sqrt(norm_z);

            if (norm_z == 0) {
                for (int i = 0; i < dim; i++) {
                    x[dim + 1].setGene(i, x[dim].getGene(i));
                }
                x[dim + 1].setFitness(x[dim].getFitness());

                // Termination criterion
                s *= 0.1;
                if (s <= EPSILON) {
                    // end the search
                    break;
                } else {
                    // next loop
                    direct = 1;

                    for (int i = 0; i < dim; i++) {
                        x[0].setGene(i, x[dim + 1].getGene(i));
                    }
                    x[0].setFitness(x[dim + 1].getFitness());
                    continue;
                }
            } else {
                for (int i = 0; i < dim; i++) {
                    v[dim][i] = z[i] / norm_z;
                }
                direct = dim + 1;

                int rest_eval = fes - evals;
                int overall_ls_eval;
                if (rest_eval < evals_per_linesearch) {
                    overall_ls_eval = rest_eval;
                } else {
                    overall_ls_eval = evals_per_linesearch;
                }

                evals += lineSearch(x[direct - 1], x[direct], overall_ls_eval, task, s, v[direct - 1]);

                // update best
                if (result.getFitness() > x[direct].getFitness()) {
                    result.setFitness(x[direct].getFitness());
                    result.setChromosome(x[direct].getChromosome().clone());
                }

                // Check appropriateness of step length
                norm_z = 0;
                for (int i = 0; i < dim; i++) {
                    double tmp = x[direct].getGene(i) - x[0].getGene(i);
                    norm_z += tmp * tmp;
                }
                norm_z = Math.sqrt(norm_z);

                if (norm_z < s) {
                    // Termination criterion
                    s *= 0.1;
                    if (s <= EPSILON) {
                        // end the search
                        break;
                    } else {
                        // next loop
                        direct = 1;

                        for (int i = 0; i < dim; i++) {
                            x[0].setGene(i, x[dim + 1].getGene(i));
                        }
                        x[0].setFitness(x[dim + 1].getFitness());
                        continue;
                    }
                } else {
                    // Orthogonalization
                    // v = GramSchmidtOrthogonalization(v, a);

                    direct = 2;
                    for (int i = 0; i < dim; i++) {
                        x[0].setGene(i, x[dim].getGene(i));
                        x[1].setGene(i, x[dim + 1].getGene(i));
                    }
                    x[0].setFitness(x[dim].getFitness());
                    x[1].setFitness(x[dim + 1].getFitness());

                    continue;
                }
            }
        }

        if (result.getFitness() > x[dim + 1].getFitness()) {
            result = x[dim + 1];
        }

        return result;
    }

    private int lineSearch(Individual start_point, Individual result, int fes, Task task, double step_size,
                           double[] v) {
        int evals = 0;
        int dim = start_point.getDim();
        double s = step_size;
        boolean change;
        boolean interpolation_flag = false;
        Individual x0 = new Individual(dim);
        Individual x = new Individual(dim);

        for (int i = 0; i < dim; i++) {
            x0.setGene(i, start_point.getGene(i));
            x.setGene(i, x0.getGene(i) + s * v[i]);
        }

        x0.setFitness(start_point.getFitness());
        x.setFitness(task.calculateFitnessValue(x.getChromosome()));
        evals++;

        // For Lagrangian quadratic interpolation
        double F[] = new double[3];
        Individual interpolation_points[] = new Individual[3];
        for (int i = 0; i < 3; i++) {
            interpolation_points[i] = new Individual(dim);
        }

        // x1, x2 of for the Lagrangian quadratic interpolation
        for (int i = 0; i < dim; i++) {
            interpolation_points[0].setGene(i, x0.getGene(i));
            interpolation_points[1].setGene(i, x.getGene(i));
        }
        F[0] = x0.getFitness();
        F[1] = x.getFitness();

        // Step backward
        if (x.getFitness() > x0.getFitness()) {
            for (int i = 0; i < dim; i++) {
                x.setGene(i, x.getGene(i) - 2 * s * v[i]);
            }
            s = -s;

            x.setFitness(task.calculateFitnessValue(x.getChromosome()));
            evals++;

            if (x.getFitness() <= x0.getFitness()) {
                change = true;

                // update x1 and x2 for the Lagrangian quadratic interpolation
                for (int i = 0; i < dim; i++) {
                    interpolation_points[0].setGene(i, x0.getGene(i));
                    interpolation_points[1].setGene(i, x.getGene(i));
                }
                F[0] = x0.getFitness();
                F[1] = x.getFitness();
            } else {
                change = false;
                interpolation_flag = true;

                // x1, x2, x3 for the Lagrangian quadratic interpolation
                for (int i = 0; i < dim; i++) {
                    interpolation_points[2].setGene(i, interpolation_points[1].getGene(i));
                    interpolation_points[1].setGene(i, interpolation_points[0].getGene(i));
                    interpolation_points[0].setGene(i, x.getGene(i));
                }
                F[2] = F[1];
                F[1] = F[0];
                F[0] = x.getFitness();
            }
        } else {
            // activate further steps
            change = true;
        }

        // Further steps
        while (change) {
            s *= 2;

            for (int i = 0; i < dim; i++) {
                x0.setGene(i, x.getGene(i));
                x.setGene(i, x0.getGene(i) + s * v[i]);
            }

            x0.setFitness(x.getFitness());
            x.setFitness(task.calculateFitnessValue(x.getChromosome()));
            evals++;

            if (x.getFitness() < x0.getFitness()) {
                // update x1 and x2 for the Lagrangian quadratic interpolation
                for (int i = 0; i < dim; i++) {
                    interpolation_points[0].setGene(i, x0.getGene(i));
                    interpolation_points[1].setGene(i, x.getGene(i));
                }
                F[0] = x0.getFitness();
                F[1] = x.getFitness();
            } else {
                change = false;
                interpolation_flag = true;

                // x3 = x
                for (int i = 0; i < dim; i++) {
                    interpolation_points[2].setGene(i, x.getGene(i));
                }
                F[2] = x.getFitness();

                // generate x = x0 + 0.5s
                s *= 0.5;
                for (int i = 0; i < dim; i++) {
                    x.setGene(i, x0.getGene(i) + s * v[i]);
                }
                x.setFitness(task.calculateFitnessValue(x.getChromosome()));
                evals++;

                // reject the one which is furthest from the point that has the smallest value
                // of the objective function
                if (x.getFitness() > F[1]) {
                    // x2 is smallest
                    // reject x3, new x3 = x
                    for (int i = 0; i < dim; i++) {
                        interpolation_points[2].setGene(i, x.getGene(i));
                    }
                    F[2] = x.getFitness();
                } else {
                    // x is smallest
                    // reject current x1, new x1 = x2, new x2 = x
                    for (int i = 0; i < dim; i++) {
                        interpolation_points[0].setGene(i, interpolation_points[1].getGene(i));
                        interpolation_points[1].setGene(i, x.getGene(i));
                    }
                    F[0] = F[1];
                    F[1] = x.getFitness();
                }
            }

            if (evals >= fes - 2) {
                change = false;
            }
        }

        // Lagrangian quadratic interpolation
        if (interpolation_flag && (F[0] - 2 * F[1] + F[2] != 0)) {
            // x = x2 + Lagrangian quadratic interpolation
            for (int i = 0; i < dim; i++) {
                x.setGene(i, interpolation_points[1].getGene(i) + s * (F[0] - F[2]) / (2.0 * (F[0] - 2 * F[1] + F[2])));
            }
            x.setFitness(task.calculateFitnessValue(x.getChromosome()));
            evals++;

            if (x.getFitness() < F[1]) {
                // best found = Lagrangian quadratic interpolation point
                result.setChromosome(x.getChromosome().clone());
                result.setFitness(x.getFitness());
            } else {
                // best found = x2
                result.setChromosome(interpolation_points[1].getChromosome().clone());
                result.setFitness(F[1]);
            }
        } else {
            // best found = x2
            result.setChromosome(interpolation_points[1].getChromosome().clone());
            result.setFitness(F[1]);
        }

        return evals;
    }

    private double[][] GramSchmidtOrthogonalization(double[][] v, double[][] a) {
        int n = v[0].length;
        double[][] new_v = new double[v.length][n];

        double norm0 = 0;
        for (int i = 0; i < n; i++) {
            norm0 += a[0][i] * a[0][i];
        }
        norm0 = Math.sqrt(norm0);
        if (norm0 != 0) {
            for (int i = 0; i < n; i++) {
                new_v[0][i] = a[0][i] / norm0;
            }
        } else {
            for (int i = 0; i < n; i++) {
                new_v[0][i] = v[0][i];
            }
        }

        for (int i = 1; i < n; i++) {
            double[] w = new double[n];
            double[] vectemp = new double[n];

            for (int j = 0; j < i; j++) {
                double temp = 0;
                for (int d = 0; d < n; d++) {
                    temp += a[i][d] * new_v[j][d];
                }

                for (int d = 0; d < n; d++) {
                    vectemp[d] += temp * new_v[j][d];
                }
            }

            double norm = 0;
            for (int d = 0; d < n; d++) {
                w[d] = a[i][d] - vectemp[d];
                norm += w[d] * w[d];
            }
            norm = Math.sqrt(norm);

            if (norm != 0) {
                for (int d = 0; d < n; d++) {
                    new_v[i][d] = w[d] / norm;
                }
            } else {
                for (int d = 0; d < n; d++) {
                    new_v[i][d] = v[i][d];
                }
            }
        }

        return new_v;
    }

//	public static void main(String[] args) {
//		Params.rand = new Random(0);
//		Params.countEvals = 0;
//		Params.maxEvals = 100000;
//
//		ProblemManager pm = new ProblemManager(ProblemConstructor.get50TasksBenchmark());
//		Task task = pm.getProblem(0).getTask(1);
//		int dim = pm.getProblem(0).DIM;
//
//		Individual x0 = new Individual(dim);
//		x0.randomInit();
//		x0.setFitness(task.calculateFitnessValue(x0.getChromosome()));
//		System.out.println("Init point: f = " + x0.getFitness());
//
//		DSCG solver = new DSCG();
//		Individual x = solver.search(x0, 100000, task);
//		System.out.println("Result: f = " + x.getFitness());
//		System.out.println(Params.countEvals);
//	}
}
