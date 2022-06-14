package main;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Random;

import basic.Params;
import basic.ProblemManager;
import benchmark.ProblemConstructor;
import core.MLSHADE_LS;

public class MainMany {

    public static void main(String[] args) throws IOException {
        Params.recordsNum = 1000; //Số lần ghi lại kết quả

        ProblemManager pm = new ProblemManager(ProblemConstructor.get50TasksBenchmark());
        for (int i=1;i<=1;i++/*int i = 0; i < pm.problems.size(); i++*/) {
            //Chạy 10 bộ benchmark
            int taskNum = pm.getProblem(i).TASKS_NUM; //số lượng tác vụ của bộ benchmark ID = i
            System.out.println("=============== Problem " + (i+1) + " =================");
            double[][][] mem_f = new double[Params.REPT][taskNum][Params.recordsNum];
            double mean[] = new double[taskNum];

            for (int seed = 0; seed < Params.REPT; seed++) {
                //Chạy 30 lần
                Params.rand = new Random(seed);
                System.out.println("Program running with seed = " + seed);
                MLSHADE_LS solver = new MLSHADE_LS(pm.getProblem(i));
                double[] result = solver.run2(mem_f[seed]);

                for (int t = 0; t < taskNum; t++) {
                    mean[t] += result[t] / Params.REPT;
                }
            }
            System.out.println("Mean result:");
            for (int t = 0; t < taskNum; t++) {
                System.out.println("Task " + t + ":\t" + String.format("%.6f", mean[t]));
            }

            // print results
            String rootFolder = "results/";
            File dir = new File(rootFolder);
            if (!dir.exists()) {
                dir.mkdir();
            }

            String dirType = rootFolder + "MTO-Many";
            dir = new File(dirType);
            if (!dir.exists()) {
                dir.mkdir();
            }

            int evalsPerRecord = Params.maxEvals / Params.recordsNum;
            String fitnessFile = dirType + "/MTOMSO_P" + (i + 1) + ".txt";
            DataOutputStream outFit = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fitnessFile)));
            for (int r = 0; r < Params.recordsNum; r++) {
                outFit.writeBytes((r + 1) * evalsPerRecord + ", ");
                for (int seed = 0; seed < Params.REPT; seed++) {
                    for (int task = 0; task < taskNum; task++) {
                        if (seed < Params.REPT - 1 || task < taskNum - 1) {
                            outFit.writeBytes(String.format("%.6f", mem_f[seed][task][r]) + ", ");
                        } else {
                            outFit.writeBytes("" + String.format("%.6f", mem_f[seed][task][r]) + "\n");
                        }
                    }
                }
            }
            outFit.close();

            String dirMean = dirType + "/Mean";
            dir = new File(dirMean);
            if (!dir.exists()) {
                dir.mkdir();
            }

            String meanFile = dirMean + "/MTOMSO_P" + (i + 1) + ".txt";
            DataOutputStream outMean = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(meanFile)));
            for (int task = 0; task < taskNum; task++) {
                outMean.writeBytes("" + (task + 1) + ", " + String.format("%.6f", mean[task]) + "\n");
            }
            outMean.close();
        }
    }
}
