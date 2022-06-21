package core;

import java.util.ArrayList;
import java.util.Arrays;

import basic.*;
import ls.DSCG;
import ls.LocalSearch;
import util.Utils;
import util.crossover.SBX;

/**
 * Multitasking L-SHADE with local search
 * @author cuonglv.hust@gmail.com
 * @date 25/5/2021
 */
public class MLSHADE_LS {

    private Problem prob;

    /**
     *  One population for each task
     */
    //Danh sách các quần thể
    private ArrayList<LSHADE_Population> pop;
    private ArrayList<CMAES_Population> pop1;

    private double[] FCP;

    public MLSHADE_LS(Problem prob) {
        this.prob = prob;
        this.pop = new ArrayList<LSHADE_Population>();
        this.pop1 = new ArrayList<CMAES_Population>();
        for (int i = 0; i < prob.TASKS_NUM; i++) {
            //Khởi tạo quần thể cho từng tác vụ
            LSHADE_Population sub_pop = new LSHADE_Population(prob, i,false);
            this.pop.add(sub_pop);
            CMAES_Population sub_pop1 = new CMAES_Population(prob,i,false);
            this.pop1.add(sub_pop1);
        }

        this.FCP = new double[prob.TASKS_NUM];
        Arrays.fill(FCP,0.5);
    }

     //LSHADE without inter-task interaction
    public void run1(double[][] mem_f) {
        //Không có tương tác giữa các tác vụ khác nhau
        Params.countEvals = 0;
        Params.maxEvals = prob.TASKS_NUM * Params.MAX_EVALS_PER_TASK;

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            pop.get(i).randomInit(Params.MAX_POP_SIZE);
            pop.get(i).update();
            System.out.println(pop.get(i).getBest().getFitness());
        }

        int gen = 0;
        while (Params.countEvals < Params.maxEvals) {
            System.out.print(gen + ":\t" + pop.get(0).getIndividuals().size() + "\t");
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                ArrayList<Individual> offspring = new ArrayList<Individual>();
                for (Individual indiv : pop.get(i).getIndividuals()) {
                    offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));
                }
                pop.get(i).getIndividuals().clear();
                pop.get(i).getIndividuals().addAll(offspring);
                pop.get(i).update();
                System.out.print(pop.get(i).getBest().getFitness() + "\t");
            }
            System.out.println();

            gen++;
        }
    }
    public double[] run2(double[][] mem_f) {
        Params.countEvals = 0;
        Params.maxEvals = prob.TASKS_NUM * Params.MAX_EVALS_PER_TASK;
        double[] result = new double[prob.TASKS_NUM];

        int recordCounter = 0; //Biến đếm số lần ghi lại kết quả
        int evalsPerRecord = Params.maxEvals / Params.recordsNum;

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            // Khởi tạo quần thể của từng tác vụ
            pop.get(i).randomInit(Params.MAX_POP_SIZE);
            pop.get(i).update();
        }

        if (Params.countEvals >= (recordCounter + 1) * evalsPerRecord) {
            //Lấy giá trị tốt nhất của từng thế hệ
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                mem_f[i][recordCounter] = pop.get(i).getBest().getFitness();
            }
            recordCounter++;
        }

        int gen = 0; //biến đếm thế hệ
        boolean stop = false;
        while (Params.countEvals < Params.maxEvals && !stop) {
            stop = true;
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                if (pop.get(i).getBest().getFitness() < Params.EPSILON) {
                    //nếu đạt giá trị tối ưu ==> bỏ qua k chạy nữa
                    continue;
                }
                stop = false;

                // reproduction
                ArrayList<Individual> offspring = new ArrayList<Individual>();
                for (Individual indiv : pop.get(i).getIndividuals()) {
                    int t = Params.rand.nextInt(prob.TASKS_NUM);
                    if (t == i) {
                        // lai ghép cùng tác vụ sử dụng LSHADE
                        offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));

                    } else {
                        double rmp;
                        if (pop.get(i).getBest_partner() == t) {
                            rmp = 1;
                        } else {
                            double mu_rmp = pop.get(i).getRmp(t);
                            do {
                                rmp = Utils.gauss(mu_rmp, 0.1);
                            } while (rmp <= 0 || rmp > 1);
                        }

                        if (Params.rand.nextDouble() <= rmp) {
                            // inter-task crossover
                            Individual indiv2 = pop.get(t).getRandomIndiv();
                            ArrayList<Individual> child = this.inter_crossover(indiv, indiv2);

                            // select best among 2 offspring
                            for (Individual c : child) {
                                c.setFitness(prob.getTask(i).calculateFitnessValue(c.getChromosome()));
                            }
                            Individual survival = child.get(0);
                            if (survival.getFitness() > child.get(1).getFitness()) {
                                survival = child.get(1);
                            }

                            double delta_f = indiv.getFitness() - survival.getFitness();
                            if (delta_f == 0) {
                                offspring.add(survival);
                            } else if (delta_f > 0) {
                                pop.get(i).addSuccessRmp(t, rmp);
                                pop.get(i).addDiffFInterX(t, delta_f);
                                offspring.add(survival);
                            } else {
                                offspring.add(indiv);
                            }
                        } else {
                            // intra
                            offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));
                        }
                    }
                }

                pop.get(i).getIndividuals().clear();
                pop.get(i).getIndividuals().addAll(offspring);

                // update RMP, F, CR, population size
                pop.get(i).update();

                // local search on the best solution using DSCG method
                if (gen % 50 == 0) {
                    LocalSearch ls = new DSCG();
                    Individual neib = ls.search(pop.get(i).getIndividual(0), 2000, prob.getTask(i));
                    if (neib.getFitness() < pop.get(i).getIndividual(0).getFitness()) {
                        pop.get(i).getIndividual(0).setFitness(neib.getFitness());
                        pop.get(i).getIndividual(0).setChromosome(neib.getChromosome().clone());
                    }
                }

                while (Params.countEvals >= (recordCounter + 1) * evalsPerRecord) {
                    for (int j = 0; j < prob.TASKS_NUM; j++) {
                        mem_f[j][recordCounter] = pop.get(j).getBest().getFitness();
                    }
                    recordCounter++;
                }
            }



            gen++;
//			for (int i=0; i<prob.TASKS_NUM; i++) {
//				System.out.print(pop.get(i).getBest().getFitness() + "\t");
//			}
//			System.out.println();
        }

        while (recordCounter < Params.recordsNum) {
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                mem_f[i][recordCounter] = pop.get(i).getBest().getFitness();
            }
            recordCounter++;
        }

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            System.out.println(i + 1 + ":\t" + String.format("%.6f", pop.get(i).getBest().getFitness()));
            result[i] = pop.get(i).getBest().getFitness();
        }
        System.out.println("-------------------------");

        return result;
    }

    public double[] run3(double[][] mem_f){
        //2 quẩn thể - CMAES và LSHADE
        Params.countEvals = 0;
        Params.maxEvals = prob.TASKS_NUM * Params.MAX_EVALS_PER_TASK;
        double[] result = new double[prob.TASKS_NUM];

        int recordCounter = 0; //Biến đếm số lần ghi lại kết quả
        int evalsPerRecord = Params.maxEvals / Params.recordsNum;

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            // Khởi tạo quần thể của từng tác vụ
//            pop1.get(i).initialize();
//            pop1.get(i).updateDistribution();

            pop.get(i).randomInit(Params.MAX_POP_SIZE);
            pop.get(i).update();

            pop1.get(i).initialize();
            pop1.get(i).getIndividuals().clear();
            pop1.get(i).getIndividuals().addAll(pop.get(i).getIndividuals());
            pop1.get(i).setBest(pop.get(i).getBest());
            pop1.get(i).updateDistribution();
        }

        if (Params.countEvals >= (recordCounter + 1) * evalsPerRecord) {
            //Lấy giá trị tốt nhất của từng thế hệ
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                mem_f[i][recordCounter] = pop.get(i).getBest().getFitness();
            }
            recordCounter++;
        }

        int gen = 0; //biến đếm thế hệ
        boolean stop = false;
        int[] initPop = new int[prob.TASKS_NUM];
        while (Params.countEvals < Params.maxEvals && !stop) {
            stop = true;
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                if (pop.get(i).getBest().getFitness() < Params.EPSILON || pop1.get(i).getBest().getFitness() < Params.EPSILON) {
                    pop.get(i).getIndividuals().addAll(pop1.get(i).getIndividuals());
                    pop.get(i).update();
                    if(pop1.get(i).getBest().getFitness() < pop.get(i).getBest().getFitness()){
                        pop.get(i).getIndividuals().get(0).setChromosome(pop1.get(i).getBest().getChromosome());
                        pop.get(i).getIndividuals().get(0).setFitness(pop1.get(i).getBest().getFitness());
                    }
                    //nếu đạt giá trị tối ưu ==> bỏ qua k chạy nữa
                    continue;
                }
                stop = false;

                // reproduction
                ArrayList<Individual> offspring = new ArrayList<Individual>();
                if (Params.countEvals < 0.6 * Params.maxEvals) {
                    pop1.get(i).updateDistribution();
                    pop1.get(i).samplePopulation();
//                    System.out.println(Params.countEvals);
                }else{
                    if(initPop[i] == 0){
                        pop.get(i).getIndividuals().addAll(pop1.get(i).getIndividuals());
                        pop.get(i).update();
                        if(pop1.get(i).getBest().getFitness() < pop.get(i).getBest().getFitness()){
                            pop.get(i).getIndividuals().get(0).setChromosome(pop1.get(i).getBest().getChromosome());
                            pop.get(i).getIndividuals().get(0).setFitness(pop1.get(i).getBest().getFitness());
                        }
                        initPop[i]++;
                    }
                    for (Individual indiv : pop.get(i).getIndividuals()) {
                        int t = Params.rand.nextInt(prob.TASKS_NUM);
                        if (t == i) {
                            // lai ghép cùng tác vụ sử dụng LSHADE
                            offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));
                        } else {
                            double rmp;
                            if (pop.get(i).getBest_partner() == t) {
                                rmp = 1;
                            } else {
                                double mu_rmp = pop.get(i).getRmp(t);
                                do {
                                    rmp = Utils.gauss(mu_rmp, 0.1);
                                } while (rmp <= 0 || rmp > 1);
                            }

                            if (Params.rand.nextDouble() <= rmp) {
                                // inter-task crossover
                                Individual indiv2 = pop.get(t).getRandomIndiv();
                                ArrayList<Individual> child = this.inter_crossover(indiv, indiv2);

                                // select best among 2 offspring
                                for (Individual c : child) {
                                    c.setFitness(prob.getTask(i).calculateFitnessValue(c.getChromosome()));
                                }
                                Individual survival = child.get(0);
                                if (survival.getFitness() > child.get(1).getFitness()) {
                                    survival = child.get(1);
                                }

                                double delta_f = indiv.getFitness() - survival.getFitness();
                                if (delta_f == 0) {
                                    offspring.add(survival);
                                } else if (delta_f > 0) {
                                    pop.get(i).addSuccessRmp(t, rmp);
                                    pop.get(i).addDiffFInterX(t, delta_f);
                                    offspring.add(survival);
                                } else {
                                    offspring.add(indiv);
                                }
                            } else {
                                // intra
                                offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));
                            }
                        }
                    }

                    pop.get(i).getIndividuals().clear();
                    pop.get(i).getIndividuals().addAll(offspring);
                    // update RMP, F, CR, population size
                    pop.get(i).update();

                    // local search on the best solution using DSCG method
                    if (gen % 50 == 0) {
                        LocalSearch ls = new DSCG();
                        Individual neib = ls.search(pop.get(i).getIndividual(0), 2000, prob.getTask(i));
                        if (neib.getFitness() < pop.get(i).getIndividual(0).getFitness()) {
                            pop.get(i).getIndividual(0).setFitness(neib.getFitness());
                            pop.get(i).getIndividual(0).setChromosome(neib.getChromosome().clone());
                        }
                    }

                    if (gen % 20 == 0) {
                        pop1.get(i).getIndividuals().addAll(pop.get(i).getTop());
                        pop1.get(i).updateDistribution();
                        pop1.get(i).samplePopulation();
                        if(pop1.get(i).getBest().getFitness() < pop.get(i).getBest().getFitness()){
                            pop.get(i).getIndividuals().get(0).setChromosome(pop1.get(i).getBest().getChromosome());
                            pop.get(i).getIndividuals().get(0).setFitness(pop1.get(i).getBest().getFitness());
                        }
//                        pop.get(i).getArchive().addAll(pop1.get(i).getIndividuals());
//                        pop.get(i).update();
                        for (int k = 0; k < pop.get(i).getIndividuals().size()*Params.ARC_RATE*0.3; k++) {
                            if(k<pop1.get(i).getIndividuals().size()){
                                pop.get(i).getArchive().remove(Params.rand.nextInt(pop.get(i).getArchive().size()));
                                pop.get(i).getArchive().add(pop1.get(i).getIndividuals().get(k));
                            }
                        }
                    }
                }

                while (Params.countEvals >= (recordCounter + 1) * evalsPerRecord) {
                    for (int j = 0; j < prob.TASKS_NUM; j++) {
                        mem_f[j][recordCounter] = pop.get(j).getBest().getFitness();
                    }
                    recordCounter++;
                }
            }
            gen++;
        }
        while (recordCounter < Params.recordsNum) {
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                mem_f[i][recordCounter] = pop.get(i).getBest().getFitness();
            }
            recordCounter++;
        }

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            System.out.println(i + 1 + ":\t" + String.format("%.6f", pop.get(i).getBest().getFitness()));
            result[i] = pop.get(i).getBest().getFitness();
        }
        System.out.println("-------------------------");

        return result;
    }
    public double[] run4(double[][] mem_f){
        //2 quẩn thể - CMAES và LSHADE
        Params.countEvals = 0;
        Params.maxEvals = prob.TASKS_NUM * Params.MAX_EVALS_PER_TASK;
        double[] result = new double[prob.TASKS_NUM];

        int recordCounter = 0; //Biến đếm số lần ghi lại kết quả
        int evalsPerRecord = Params.maxEvals / Params.recordsNum;

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            // Khởi tạo quần thể của từng tác vụ
            pop.get(i).randomInit(Params.MAX_POP_SIZE);
            pop.get(i).update();

            pop1.get(i).initialize();
            pop1.get(i).getIndividuals().clear();
            pop1.get(i).getIndividuals().addAll(pop.get(i).getIndividuals());
            pop1.get(i).setBest(pop.get(i).getBest());
            pop1.get(i).updateDistribution();
        }

        if (Params.countEvals >= (recordCounter + 1) * evalsPerRecord) {
            //Lấy giá trị tốt nhất của từng thế hệ
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                mem_f[i][recordCounter] = pop.get(i).getBest().getFitness();
            }
            recordCounter++;
        }

        int gen = 0; //biến đếm thế hệ
        boolean stop = false;
        while (Params.countEvals < Params.maxEvals && !stop) {
            stop = true;
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                if (pop.get(i).getBest().getFitness() < Params.EPSILON) {
                    //nếu đạt giá trị tối ưu ==> bỏ qua k chạy nữa
                    continue;
                }
                stop = false;
                System.out.println(pop.get(1).getBest().getFitness() +" "+pop1.get(1).getBest().getFitness());
                if(Params.rand.nextDouble() < FCP[i]){
                    //run CMAES
                    pop1.get(i).updateDistribution();
                    pop1.get(i).samplePopulation();

                    if(pop1.get(i).getBest().getFitness() < pop.get(i).getBest().getFitness()){
                        pop.get(i).sortPop();
                        pop.get(i).getIndividuals().get(0).setFitness(pop1.get(i).getBest().getFitness());
                        pop.get(i).getIndividuals().get(0).setChromosome(pop1.get(i).getBest().getChromosome());
                    }
                }else {
                    // reproduction
                    ArrayList<Individual> offspring = new ArrayList<Individual>();
                    for (Individual indiv : pop.get(i).getIndividuals()) {
                        int t = Params.rand.nextInt(prob.TASKS_NUM);
                        if (t == i) {
                            // lai ghép cùng tác vụ sử dụng LSHADE
                            offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));
                        } else {
                            double rmp;
                            if (pop.get(i).getBest_partner() == t) {
                                rmp = 1;
                            } else {
                                double mu_rmp = pop.get(i).getRmp(t);
                                do {
                                    rmp = Utils.gauss(mu_rmp, 0.1);
                                } while (rmp <= 0 || rmp > 1);
                            }

                            if (Params.rand.nextDouble() <= rmp) {
                                // inter-task crossover
                                Individual indiv2 = pop.get(t).getRandomIndiv();
                                ArrayList<Individual> child = this.inter_crossover(indiv, indiv2);

                                // select best among 2 offspring
                                for (Individual c : child) {
                                    c.setFitness(prob.getTask(i).calculateFitnessValue(c.getChromosome()));
                                }
                                Individual survival = child.get(0);
                                if (survival.getFitness() > child.get(1).getFitness()) {
                                    survival = child.get(1);
                                }

                                double delta_f = indiv.getFitness() - survival.getFitness();
                                if (delta_f == 0) {
                                    offspring.add(survival);
                                } else if (delta_f > 0) {
                                    pop.get(i).addSuccessRmp(t, rmp);
                                    pop.get(i).addDiffFInterX(t, delta_f);
                                    offspring.add(survival);
                                } else {
                                    offspring.add(indiv);
                                }
                            } else {
                                // intra
                                offspring.add(pop.get(i).operateCurrentToPBest1BinWithArchive(indiv));
                            }
                        }
                    }

                    pop.get(i).getIndividuals().clear();
                    pop.get(i).getIndividuals().addAll(offspring);
                    // update RMP, F, CR, population size
                    pop.get(i).update();

                    if(pop.get(i).getBest().getFitness() < pop1.get(i).getBest().getFitness()){
                        pop1.get(i).getIndividuals().addAll(pop.get(i).getIndividuals());
                        pop1.get(i).updateDistribution();
                        pop1.get(i).getBest().setFitness(pop.get(i).getBest().getFitness());
                        pop1.get(i).getBest().setChromosome(pop.get(i).getBest().getChromosome());
                    }
                }

                // local search on the best solution using DSCG method
                if (gen % 50 == 0) {
                    LocalSearch ls = new DSCG();
                    Individual neib = ls.search(pop.get(i).getIndividual(0), 2000, prob.getTask(i));
                    if (neib.getFitness() < pop.get(i).getIndividual(0).getFitness()) {
                        pop.get(i).getIndividual(0).setFitness(neib.getFitness());
                        pop.get(i).getIndividual(0).setChromosome(neib.getChromosome().clone());
                    }
                }

                while (Params.countEvals >= (recordCounter + 1) * evalsPerRecord) {
                    for (int j = 0; j < prob.TASKS_NUM; j++) {
                        mem_f[j][recordCounter] = pop.get(j).getBest().getFitness();
                    }
                    recordCounter++;
                }
            }
            gen++;
        }
        while (recordCounter < Params.recordsNum) {
            for (int i = 0; i < prob.TASKS_NUM; i++) {
                mem_f[i][recordCounter] = pop.get(i).getBest().getFitness();
            }
            recordCounter++;
        }

        for (int i = 0; i < prob.TASKS_NUM; i++) {
            System.out.println(i + 1 + ":\t" + String.format("%.6f", pop.get(i).getBest().getFitness()));
            result[i] = pop.get(i).getBest().getFitness();
        }
        System.out.println("-------------------------");

        return result;
    }

    private ArrayList<Individual> inter_crossover(Individual p1, Individual p2) {
        ArrayList<double[]> offspring_gene = SBX.generateOffspring(p1.getChromosome(), p2.getChromosome());
        for (int i=0;i< prob.DIM;i++){
            if(Params.rand.nextDouble()<0.5){
                double temp = offspring_gene.get(0)[i];
                offspring_gene.get(0)[i] = offspring_gene.get(1)[i];
                offspring_gene.get(1)[i] = temp;
            }
        }
        ArrayList<Individual> offspring = new ArrayList<Individual>();
        offspring.add(new Individual(offspring_gene.get(0)));
        offspring.add(new Individual(offspring_gene.get(1)));
        return offspring;
    }
}
