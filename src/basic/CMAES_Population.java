package basic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;

import util.Utils;
import RWS.RWS;

public class CMAES_Population {
    private Problem problem; //bộ benchmark (ID 1==>10)
    private int task_index; //tác vụ thứ i
    private ArrayList<Individual> individuals; //Dân số hiện tại.


    private int lambda; //Số lượng con cái mỗi lần lặp
    private int mu; //Số lượng cá thể được chọn để tái tổ hợp
    private double sigma; //Độ lệch chuẩn tổng thể (step size)
    private double mueff; //phương sai hiệu quả
    private double ccov; //learning rate
    private double ccovsep; //learning rate khi chế độ chẩn đoán đang hoạt động
    private double chiN; // Kỳ vọng của || N (0, I) ||.
    private double cs; //Tham số tích lũy kích thước bước.
    private double cc; //Tham số tích lũy
    private double damps; //Giảm chấn cho kích thước bước
    private double[] weights; //Trọng lượng để tái tổ hợp
    private double[] diagD; //Hệ số tỷ lệ
    private double[] xmean; //Trung tâm hiện tại của bản phân phối.
    private double[] pc; //Con đường tiến hóa
    private double[] ps; //Con đường tiến hóa liên hợp
    private double[][] B; //Hệ tọa độ
    private double[][] C; //Ma trận hiệp phương sai hiện tại

    private int lastEigenupdate; //Lần lặp cuối cùng là quá trình phân hủy giá trị eigenvalue đã được tính toán.
    private Individual best; //Cá thể tốt nhất được tìm thấy cho đến nay
    private int N; //Số chiều không gian tìm kiếm
    private int iteration; //Đếm số lần lặp

    /**
     * if true, perform consistency checks to ensure the algorithm remains
     * numerically stable.
     */
    private boolean checkConsistency;

    /**
     * The number of iterations in which only the covariance diagonal is used. This
     * enhancement helps speed up the algorithm when there are many decision
     * variables. Set to {@code 0} to always use the full covariance matrix.
     */
    private int diagonalIterations;

    public CMAES_Population(Problem problem, int task_idx, boolean checkConsistency) {
        this.problem = problem;
        this.task_index = task_idx;
        this.individuals = new ArrayList<Individual>();

        this.N = problem.DIM;
        this.lambda = (int) (4 + 3 * Math.floor(Math.log(problem.DIM)));
        this.checkConsistency = checkConsistency;
    }

    public void initialize() {
        sigma = 0.5;
        diagonalIterations = 150 * N / lambda;
        diagD = new double[N];
        pc = new double[N];
        ps = new double[N];
        B = new double[N][N];
        C = new double[N][N];

        for (int i = 0; i < N; i++) {
            pc[i] = 0;
            ps[i] = 0;
            diagD[i] = 1;

            for (int j = 0; j < N; j++) {
                B[i][j] = 0;
            }

            for (int j = 0; j < i; j++) {
                C[i][j] = 0;
            }

            B[i][i] = 1;
            C[i][i] = diagD[i] * diagD[i];
        }

        xmean = new double[N];
        for (int i = 0; i < N; i++) {
            xmean[i] = Params.rand.nextDouble();
        }

        chiN = Math.sqrt(N) * (1.0 - 1.0 / (4.0 * N) + 1.0 / (21.0 * N * N));
        mu = (int) Math.floor(lambda / 2.0);
        weights = new double[mu];

        for (int i = 0; i < mu; i++) {
            weights[i] = Math.log(mu + 1) - Math.log(i + 1);
        }

        double sum = Utils.sum(weights);

        for (int i = 0; i < mu; i++) {
            weights[i] /= sum;
        }

        double sumSq = Utils.sumSq(weights);

        mueff = 1.0 / sumSq;

        cs = (mueff + 2) / (N + mueff + 3);
        damps = (1 + 2 * Math.max(0, Math.sqrt((mueff - 1.0) / (N + 1)) - 1)) + cs;
        cc = 4.0 / (N + 4.0);
        ccov = 2.0 / (N + 1.41) / (N + 1.41) / mueff
                + (1 - (1.0 / mueff)) * Math.min(1, (2 * mueff - 1) / (mueff + (N + 2) * (N + 2)));
        ccovsep = Math.min(1, ccov * (N + 1.5) / 3.0);

    }

    /**
     * Performs eigenvalue decomposition to update B and diagD.
     */
    public void eigenDecomposition() {
        lastEigenupdate = iteration;

        if (diagonalIterations >= iteration) {
            for (int i = 0; i < N; i++) {
                diagD[i] = Math.sqrt(C[i][i]);
            }
        } else {
            // set B <- C
            for (int i = 0; i < N; i++) {
                for (int j = 0; j <= i; j++) {
                    B[i][j] = B[j][i] = C[i][j];
                }
            }

            // eigenvalue decomposition
            double[] offdiag = new double[N];
            Utils.tred2(N, B, diagD, offdiag);
            Utils.tql2(N, diagD, offdiag, B);

//			if (checkConsistency) {
//				Utils.checkEigenSystem(N, C, diagD, B);
//			}

            // assign diagD to eigenvalue square roots
            for (int i = 0; i < N; i++) {
                if (diagD[i] < 0) { // numerical problem?
                    System.err.println("an eigenvalue has become negative");
                    diagD[i] = 0;
                }

                diagD[i] = Math.sqrt(diagD[i]);
            }
        }
    }

    /**
     * Test and correct any numerical issues.
     */
    public void testAndCorrectNumerics() {
        // flat fitness, test is function values are identical
        if (individuals.size() > 0) {
            sortPop();

            if (individuals.get(0).getFitness() == individuals.get(Math.min(lambda - 1, lambda / 2 + 1) - 1).getFitness()) {
                System.err.println("flat fitness landscape, consider reformulation of fitness, step size increased");
                sigma *= Math.exp(0.2 + cs / damps);
            }
        }

        // align (renormalize) scale C (and consequently sigma)
        double fac = 1.0;

        if (Utils.max(diagD) < 1e-6) {
            fac = 1.0 / Utils.max(diagD);
        } else if (Utils.min(diagD) > 1e4) {
            fac = 1.0 / Utils.min(diagD);
        }

        if (fac != 1.0) {
            sigma /= fac;

            for (int i = 0; i < N; i++) {
                pc[i] *= fac;
                diagD[i] *= fac;

                for (int j = 0; j <= i; j++) {
                    C[i][j] *= fac * fac;
                }
            }
        }
    }
    public void sortPop() {
        //Sắp xếp lại các cá thể trong quần thể
        this.individuals.sort(new Comparator<Individual>() {
            @Override
            public int compare(Individual o1, Individual o2) {
                return Double.valueOf(o1.getFitness()).compareTo(o2.getFitness());
            }
        });
    }
    public void samplePopulation() {
        iteration++;
//        Individual best = null;

        if ((iteration - lastEigenupdate) > 1.0 / ccov / N / 5.0) {
            eigenDecomposition();
        }

        if (checkConsistency) {
            testAndCorrectNumerics();
        }

        individuals.clear();
        if (this.diagonalIterations >= iteration) {
            for (int i = 0; i < this.lambda; i++) {
                double[] x = new double[N];
                for (int j = 0; j < N; j++) {
                    do {
                        x[j] = xmean[j] + sigma * diagD[j] * Params.rand.nextGaussian();
                    } while (x[j] > 1 || x[j] < 0);
                }

                Individual indiv = new Individual(x);
                indiv.setFitness(problem.getTask(this.task_index).calculateFitnessValue(indiv.getChromosome()));
                individuals.add(indiv);

                if (best == null || best.getFitness() > indiv.getFitness()) {
                    best = indiv;
                }
            }
        } else {
            for (int i = 0; i < this.lambda; i++) {
                boolean feasible;
                double[] x = new double[N];

                do {
                    feasible = true;
                    double[] artmp = new double[N];
                    for (int j = 0; j < N; j++) {
                        artmp[j] = diagD[j] * Params.rand.nextGaussian();
                    }

                    for (int j = 0; j < N; j++) {
                        double sum = 0.0;
                        for (int k = 0; k < N; k++) {
                            sum += B[j][k] * artmp[k];
                        }

                        x[j] = xmean[j] + sigma * sum;

                        if (x[j] > 1 || x[j] < 0) {
                            feasible = false;
                            break;
                        }
                    }
                } while (!feasible);

                Individual indiv = new Individual(x);
                indiv.setFitness(problem.getTask(this.task_index).calculateFitnessValue(indiv.getChromosome()));
                individuals.add(indiv);

                if (best == null || best.getFitness() > indiv.getFitness()) {
                    best = indiv;
                }
            }
        }
    }

    public void updateDistribution() {
        double[] xold = xmean.clone();
        double[] BDz = new double[N];
        double[] artmp = new double[N];

        sortPop();
        // calculate xmean and BDz
        for (int i = 0; i < N; i++) {
            xmean[i] = 0;

            for (int j = 0; j < mu; j++) {
                if(j<individuals.size())
                    xmean[i] += weights[j] * individuals.get(j).getGene(i);
            }

            BDz[i] = Math.sqrt(mueff) * (xmean[i] - xold[i]) / sigma;
        }

        // cumulation for sigma (ps) using B*z
        if (diagonalIterations >= iteration) {
            // given B=I we have B*z = z = D^-1 BDz
            for (int i = 0; i < N; i++) {
                ps[i] = (1.0 - cs) * ps[i] + Math.sqrt(cs * (2.0 - cs)) * BDz[i] / diagD[i];
            }
        } else {
            for (int i = 0; i < N; i++) {
                double sum = 0.0;

                for (int j = 0; j < N; j++) {
                    sum += B[j][i] * BDz[j];
                }

                artmp[i] = sum / diagD[i];
            }

            for (int i = 0; i < N; i++) {
                double sum = 0.0;

                for (int j = 0; j < N; j++) {
                    sum += B[i][j] * artmp[j];
                }

                ps[i] = (1.0 - cs) * ps[i] + Math.sqrt(cs * (2.0 - cs)) * sum;
            }
        }

        // calculate norm(ps)^2
        double psxps = 0;
        for (int i = 0; i < N; i++) {
            psxps += ps[i] * ps[i];
        }

        // cumulation for covariance matrix (pc) using B*D*z
        int hsig = 0;
        if (Math.sqrt(psxps) / Math.sqrt(1.0 - Math.pow(1.0 - cs, 2.0 * iteration)) / chiN < 1.4 + 2.0 / (N + 1)) {
            hsig = 1;
        }

        for (int i = 0; i < N; i++) {
            pc[i] = (1.0 - cc) * pc[i] + hsig * Math.sqrt(cc * (2.0 - cc)) * BDz[i];
        }

        // update of C
        for (int i = 0; i < N; i++) {
            for (int j = (diagonalIterations >= iteration ? i : 0); j <= i; j++) {
                C[i][j] = (1.0 - (diagonalIterations >= iteration ? ccovsep : ccov)) * C[i][j]
                        + ccov * (1.0 / mueff) * (pc[i] * pc[j] + (1 - hsig) * cc * (2.0 - cc) * C[i][j]);

                for (int k = 0; k < mu; k++) {
                    if(k<individuals.size())
                        C[i][j] += ccov * (1 - 1.0 / mueff) * weights[k] * (individuals.get(k).getGene(i) - xold[i]) * (individuals.get(k).getGene(j) - xold[j]) / sigma / sigma;
                }
            }
        }

        // update of sigma
        sigma *= Math.exp(((Math.sqrt(psxps) / chiN) - 1) * cs / damps);
    }
    public Individual getBest() {
        return best;
    }
    public void setBest(Individual c){
        best = c;
    }
    public ArrayList<Individual> getIndividuals() {
        return this.individuals;
    }
//    public Individual run() {
//        Params.countEvals = 0;
//
//        best = samplePopulation();
//        System.out.println(iteration + "\t" + best.getFitness());
//
//        while (Params.countEvals < Params.maxEvals) {
//            updateDistribution();
//            Individual indiv = samplePopulation();
//
//            if (indiv.getFitness() < best.getFitness()) {
//                best = indiv;
//            }
//
//            System.out.println(iteration + "\t" + best.getFitness());
//        }
//
//        return best;
//    }
}
