package basic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;

import util.Utils;
import RWS.RWS;

public class LSHADE_Population {
    //Mỗi tác vụ có 1 quần thể
    private Problem problem; //bộ benchmark (ID 1==>10)
    private int task_index; //tác vụ thứ i
    private ArrayList<Individual> individuals; //Danh sách các cá thể
    private Individual best; //Cá thể tốt nhất trong quần thể

    private ArrayList<Individual> archive; //
    private ArrayList<Individual> top; //Chứa n cá thể tốt nhất của quần thể

    private double[] mem_cr; //Bộ nhớ lịch sử thành công CR, F
    private double[] mem_f;

    private ArrayList<Double> s_cr; //  Lưu trữ các giá trị thành công từng thế hệ
    private ArrayList<Double> s_f; //   Sử dụng để cập nhật bộ nhớ lịch sử thành công sau mỗi thế hệ


    private ArrayList<Double> diff_f; //chênh lệch fitness của cá thể con và cha
                                    //Sử dụng để tính độ ảnh hưởng của index i đến việc cập nhật mem_cr,mem_f

    private double[] PWI; //vecto quán tính toàn dân số
    private ArrayList<Double[]> mv;
    private int NI; // Đếm số con thành công

    private int mem_pos;
    private int gen;

    private double rmp[]; //Bộ nhớ lịch sử thành công rmp
    private ArrayList<Double>[] s_rmp; //   Lưu trữ các giá trị thành công từng thế hệ
                                    //      Sử dụng để cập nhật bộ nhớ lịch sử thành công sau mỗi thế hệ
    private ArrayList<Double>[] diff_f_inter_x; //Sử dụng để tính độ ảnh hưởng của index i đến việc cập nhật rmp[]
    private int best_partner;
    private int count_evals; //Số lượng đánh giá

    public LSHADE_Population(Problem problem, int task_idx) {
        this.problem = problem;
        this.task_index = task_idx;
        this.individuals = new ArrayList<Individual>();

        mem_cr = new double[Params.H];
        mem_f = new double[Params.H];
        s_cr = new ArrayList<Double>();
        s_f = new ArrayList<Double>();
        diff_f = new ArrayList<Double>();
        mem_pos = 0;
        PWI = new double[problem.DIM];
        mv = new ArrayList<Double[]>();
        NI = 0;

        for (int i = 0; i < Params.H; i++) {
            mem_cr[i] = 0.5;
            mem_f[i] = 0.5;
        }

        archive = new ArrayList<Individual>();
        top = new ArrayList<Individual>();
        gen = 0;

        rmp = new double[problem.TASKS_NUM];
        s_rmp = new ArrayList[problem.TASKS_NUM];
        diff_f_inter_x = new ArrayList[problem.TASKS_NUM];
        this.best_partner = -1;

        for (int i = 0; i < problem.TASKS_NUM; i++) {
            rmp[i] = Params.INIT_RMP;
            s_rmp[i] = new ArrayList<Double>();
            diff_f_inter_x[i] = new ArrayList<Double>();
        }
        count_evals = 0;
    }

    public void randomInit(int size) {
        this.individuals.clear();
        //size - số lượng cá thể khởi tạo
        for (int i = 0; i < size; i++) {
            Individual indiv = new Individual(problem.DIM);
            indiv.randomInit();
            indiv.setFitness(problem.getTask(this.task_index).calculateFitnessValue(indiv.getChromosome()));
            count_evals++;
            this.individuals.add(indiv);
        }
    }

    public void update() {
        gen++;

        // update F, CR memory
        if (s_cr.size() > 0) {
            mem_cr[mem_pos] = 0;
            mem_f[mem_pos] = 0;
            double temp_sum_cr = 0;
            double temp_sum_f = 0;
            double sum_diff = 0;

            for (double d : diff_f) {
                sum_diff += d;
            }

            for (int i = 0; i < s_cr.size(); i++) {
                double weight = diff_f.get(i) / sum_diff;

                mem_f[mem_pos] += weight * s_f.get(i) * s_f.get(i);
                temp_sum_f += weight * s_f.get(i);

                mem_cr[mem_pos] += weight * s_cr.get(i) * s_cr.get(i);
                temp_sum_cr += weight * s_cr.get(i);
            }

            mem_f[mem_pos] /= temp_sum_f;

            if (temp_sum_cr == 0 || mem_cr[mem_pos] == -1)
                mem_cr[mem_pos] = -1;
            else
                mem_cr[mem_pos] /= temp_sum_cr;

            mem_pos++;
            if (mem_pos >= Params.H) {
                mem_pos = 0;
            }

            s_cr.clear();
            s_f.clear();
            diff_f.clear();
        }

        // update pop size
        int plan_pop_size = (int) Math.round(Params.MAX_POP_SIZE
                + (Params.MIN_POP_SIZE - Params.MAX_POP_SIZE) * (1.0 * Params.countEvals/Params.maxEvals));
        if (plan_pop_size < Params.MIN_POP_SIZE) {
            plan_pop_size = Params.MIN_POP_SIZE;
        }

        sortPop();
        if (best == null || best.getFitness() > this.individuals.get(0).getFitness()) {
            best = this.individuals.get(0);
        }

        if (this.individuals.size() > plan_pop_size) {
            while (this.individuals.size() > plan_pop_size) {
                //Xóa các cá thể tệ
                this.individuals.remove(this.individuals.size() - 1);
            }

            // resize the archive size
            int arc_size = (int) (this.individuals.size() * Params.ARC_RATE);
            while (archive.size() > arc_size) {
                //Xóa ngẫu nhiên cá thể
                archive.remove(Params.rand.nextInt(archive.size()));
            }
        }

        // update pbest list
        top.clear();
        int pbest_size = (int) (Params.BEST_RATE * this.individuals.size());
        pbest_size = Math.max(pbest_size,2); //tối thiểu phải bằng 2
        for (int i = 0; i < pbest_size; i++) {
            //Chứa n cá thể tốt nhất của quần thể
            top.add(this.individuals.get(i));
        }

        // update RMP
        best_partner = -1;
        double maxRmp = 0;
        for (int i = 0; i < problem.TASKS_NUM; i++) {
            if (i != this.task_index) {
                double good_mean = 0;
                if (s_rmp[i].size() > 0) {
                    double sum = 0;
                    for (double d : diff_f_inter_x[i]) {
                        sum += d;
                    }

                    double val1 = 0, val2 = 0, w;
                    for (int k = 0; k < s_rmp[i].size(); k++) {
                        w = diff_f_inter_x[i].get(k) / sum;
                        val1 += w * s_rmp[i].get(k) * s_rmp[i].get(k);
                        val2 += w * s_rmp[i].get(k);
                    }
                    good_mean = val1 / val2;

                    if (good_mean > rmp[i] && good_mean > maxRmp) {
                        maxRmp = good_mean;
                        best_partner = i;
                    }
                }

                double c1 = s_rmp[i].size() > 0 ? 1.0 : 1.0 - Params.C;
                rmp[i] = c1 * rmp[i] + Params.C * good_mean;
                rmp[i] = Math.max(0.01, Math.min(1, rmp[i]));

                s_rmp[i].clear();
                diff_f_inter_x[i].clear();
            }
        }
        //update PWI
        if(NI > 0.1*individuals.size()){
            for (int i=0;i<PWI.length;i++){
                for (int j=0;j<mv.size();j++){
                    PWI[i]+=mv.get(j)[i]/NI;
                }
            }
        }else {
            for (int i=0;i<PWI.length;i++){
                PWI[i]=0;
            }
        }
        mv.clear();
        NI=0;
    }
    public void updatePop1() {
        //Không giảm quy mô quần thể
        gen++;

        // update F, CR memory
        if (s_cr.size() > 0) {
            mem_cr[mem_pos] = 0;
            mem_f[mem_pos] = 0;
            double temp_sum_cr = 0;
            double temp_sum_f = 0;
            double sum_diff = 0;

            for (double d : diff_f) {
                sum_diff += d;
            }

            for (int i = 0; i < s_cr.size(); i++) {
                double weight = diff_f.get(i) / sum_diff;

                mem_f[mem_pos] += weight * s_f.get(i) * s_f.get(i);
                temp_sum_f += weight * s_f.get(i);

                mem_cr[mem_pos] += weight * s_cr.get(i) * s_cr.get(i);
                temp_sum_cr += weight * s_cr.get(i);
            }

            mem_f[mem_pos] /= temp_sum_f;

            if (temp_sum_cr == 0 || mem_cr[mem_pos] == -1)
                mem_cr[mem_pos] = -1;
            else
                mem_cr[mem_pos] /= temp_sum_cr;

            mem_pos++;
            if (mem_pos >= Params.H) {
                mem_pos = 0;
            }

            s_cr.clear();
            s_f.clear();
            diff_f.clear();
        }


        sortPop();
        if (best == null || best.getFitness() > this.individuals.get(0).getFitness()) {
            best = this.individuals.get(0);
        }

        // update pbest list
        top.clear();
        int pbest_size = (int) (Params.BEST_RATE * this.individuals.size());
        pbest_size = Math.max(pbest_size,2); //tối thiểu phải bằng 2
        for (int i = 0; i < pbest_size; i++) {
            //Chứa n cá thể tốt nhất của quần thể
            top.add(this.individuals.get(i));
        }

        // update RMP
        best_partner = -1;
        double maxRmp = 0;
        for (int i = 0; i < problem.TASKS_NUM; i++) {
            if (i != this.task_index) {
                double good_mean = 0;
                if (s_rmp[i].size() > 0) {
                    double sum = 0;
                    for (double d : diff_f_inter_x[i]) {
                        sum += d;
                    }

                    double val1 = 0, val2 = 0, w;
                    for (int k = 0; k < s_rmp[i].size(); k++) {
                        w = diff_f_inter_x[i].get(k) / sum;
                        val1 += w * s_rmp[i].get(k) * s_rmp[i].get(k);
                        val2 += w * s_rmp[i].get(k);
                    }
                    good_mean = val1 / val2;

                    if (good_mean > rmp[i] && good_mean > maxRmp) {
                        maxRmp = good_mean;
                        best_partner = i;
                    }
                }

                double c1 = s_rmp[i].size() > 0 ? 1.0 : 1.0 - Params.C;
                rmp[i] = c1 * rmp[i] + Params.C * good_mean;
                rmp[i] = Math.max(0.01, Math.min(1, rmp[i]));

                s_rmp[i].clear();
                diff_f_inter_x[i].clear();
            }
        }
    }

    public Individual operateCurrentToPBest1BinWithArchive(Individual current) {
        Individual r1, r2, pbest;

        int rand_pos = Params.rand.nextInt(Params.H);
        double mu_cr = mem_cr[rand_pos];
        double mu_f = mem_f[rand_pos];
        double cr, f;

        if (mu_cr == -1) {
            cr = 0;
        } else {
            cr = Utils.gauss(mu_cr, 0.1);
            if (cr > 1)
                cr = 1;
            else if (cr < 0)
                cr = 0;
        }

        do {
            f = Utils.cauchy_g(mu_f, 0.1);
        } while (f <= 0);
        if (f > 1)
            f = 1;

        int dem=0; // Biến đếm số lần chọn (Để tránh trường hợp lặp vô hạn)
        do {
            if(dem<10){
                pbest = top.get(RWS.runRWS(top));
                dem++;
            }else pbest = top.get(Params.rand.nextInt(top.size()));
        } while (pbest.getID() == current.getID());
        dem=0;
        do {
            if(dem<10){
                r1 = individuals.get(RWS.runRWS(individuals));
                dem++;
            }else r1 = individuals.get(Params.rand.nextInt(individuals.size()));
        } while (r1.getID() == current.getID() || r1.getID() == pbest.getID());

        dem = 0;
        if (archive.size() > 0 && Params.rand.nextDouble() <= archive.size() / (1.0 * archive.size() + individuals.size())) {
            r2 = archive.get(Params.rand.nextInt(archive.size()));
        } else {
            do {
                if(dem<10){
                    r2 = individuals.get(RWS.runRWS(individuals));
                    dem++;
                }else r2 = individuals.get(Params.rand.nextInt(individuals.size()));
            } while (r2.getID() == current.getID() || r2.getID() == r1.getID() || r2.getID() == pbest.getID());
        }

        int j_rand = Params.rand.nextInt(problem.DIM);
        Individual child = new Individual(problem.DIM);
        for (int j = 0; j < problem.DIM; j++) {
            if (Params.rand.nextDouble() <= cr || j == j_rand) {
                child.setGene(j, current.getGene(j)
                        + f * (pbest.getGene(j) - current.getGene(j) + r1.getGene(j) - r2.getGene(j) + PWI[j]));

                // bound handling
                if (child.getGene(j) > 1) {
                    child.setGene(j, (current.getGene(j) + 1) / 2.0);
                } else if (child.getGene(j) < 0) {
                    child.setGene(j, (current.getGene(j) + 0) / 2.0);
                }
            } else {
                child.setGene(j, current.getGene(j));
            }
        }

        child.setFitness(problem.getTask(task_index).calculateFitnessValue(child.getChromosome()));
        count_evals++;
        if (child.getFitness() == current.getFitness()) {

            Double[] Temp_mv = new Double[child.getChromosome().length];
            for(int i=0;i<child.getChromosome().length;i++){
                Temp_mv[i] = child.getGene(i) - current.getGene(i);
            }
            mv.add(Temp_mv);
            NI++;

            return child;
        } else if (child.getFitness() < current.getFitness()) {
            s_cr.add(cr);
            s_f.add(f);
            diff_f.add(current.getFitness() - child.getFitness());

            if (archive.size() < Params.ARC_RATE * this.individuals.size()) {
                archive.add(current);
            } else {
                archive.remove(Params.rand.nextInt(archive.size()));
                archive.add(current);
            }

            Double[] Temp_mv = new Double[child.getChromosome().length];
            for(int i=0;i<child.getChromosome().length;i++){
                Temp_mv[i] = child.getGene(i) - current.getGene(i);
            }
            mv.add(Temp_mv);
            NI++;

            return child;
        } else {
            return current;
        }
    }

//    public Individual runCMA_ES(Individual current){
//
//    }

    public Individual operateRand1(Individual current) {
        Individual r1, r2, r3;

        int rand_pos = Params.rand.nextInt(Params.H);
        double mu_cr = mem_cr[rand_pos];
        double mu_f = mem_f[rand_pos];
        double cr, f;

        if (mu_cr == -1) {
            cr = 0;
        } else {
            cr = Utils.gauss(mu_cr, 0.1);
            if (cr > 1)
                cr = 1;
            else if (cr < 0)
                cr = 0;
        }

        do {
            f = Utils.cauchy_g(mu_f, 0.1);
        } while (f <= 0);
        if (f > 1)
            f = 1;

        do {
            r1 = individuals.get(Params.rand.nextInt(individuals.size()));
        } while (r1.getID() == current.getID());

        do {
            r2 = individuals.get(Params.rand.nextInt(individuals.size()));
        } while (r2.getID() == current.getID() || r1.getID() == r2.getID());

        do {
            r3 = individuals.get(Params.rand.nextInt(individuals.size()));
        } while (r3.getID() == current.getID() || r3.getID() == r1.getID() || r3.getID() == r2.getID());

        int j_rand = Params.rand.nextInt(problem.DIM);
        Individual child = new Individual(problem.DIM);
        for (int j = 0; j < problem.DIM; j++) {
            if (Params.rand.nextDouble() <= cr || j == j_rand) {
                child.setGene(j, r1.getGene(j)
                        + f * (r2.getGene(j) - r3.getGene(j)));

                // bound handling
                if (child.getGene(j) > 1) {
                    child.setGene(j, (current.getGene(j) + 1) / 2.0);
                } else if (child.getGene(j) < 0) {
                    child.setGene(j, (current.getGene(j) + 0) / 2.0);
                }
            } else {
                child.setGene(j, current.getGene(j));
            }
        }

        child.setFitness(problem.getTask(task_index).calculateFitnessValue(child.getChromosome()));
        count_evals++;
        if (child.getFitness() == current.getFitness()) {
            return child;
        } else if (child.getFitness() < current.getFitness()) {
            s_cr.add(cr);
            s_f.add(f);
            diff_f.add(current.getFitness() - child.getFitness());

            if (archive.size() < Params.ARC_RATE * this.individuals.size()) {
                archive.add(current);
            } else {
                archive.remove(Params.rand.nextInt(archive.size()));
                archive.add(current);
            }
            return child;
        } else {
            return current;
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

    public double getRmp(int task) {
        return rmp[task];
    }

    public Individual getRandomIndiv() {
        return this.individuals.get(Params.rand.nextInt(this.individuals.size()));
    }

    public int getBest_partner() {
        return best_partner;
    }

    public ArrayList<Individual> getIndividuals() {
        return this.individuals;
    }

    public Individual getIndividual(int index) {
        return this.individuals.get(index);
    }

    public void addIndividuals(Collection<Individual> indivs) {
        this.individuals.addAll(indivs);
    }

    public void addSuccessRmp(int partner, double rmp) {
        this.s_rmp[partner].add(rmp);
    }

    public void addDiffFInterX(int partner, double diff) {
        this.diff_f_inter_x[partner].add(diff);
    }

    public ArrayList<Individual> getTop() {
        return top;
    }
    public ArrayList<Individual> getArchive(){
        return archive;
    }
    public Individual getBest() {
        return best;
    }
    public int getCount_evals(){
        return count_evals;
    }

}
