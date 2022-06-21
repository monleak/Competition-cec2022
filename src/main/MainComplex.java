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

public class MainComplex {

    public static void main(String[] args) throws IOException {
        Params.recordsNum = 100;

        ProblemManager pm = new ProblemManager(ProblemConstructor.getComplexProblems());
        for (int i = 0; i < pm.problems.size(); i++) {
            int taskNum = pm.getProblem(i).TASKS_NUM;
            System.out.println("=============== Problem " + (i+1) + " =================");
            double[][][] mem_f = new double[Params.REPT][taskNum][Params.recordsNum];
            double mean[] = new double[taskNum]; //Lưu kết quả trung bình sau 30 lần chạy

            for (int seed = 0; seed < Params.REPT; seed++) {
                Params.rand = new Random(seed);
                System.out.println("Program running with seed = " + seed);
                MLSHADE_LS solver = new MLSHADE_LS(pm.getProblem(i));
                double[] result = solver.run3(mem_f[seed]);

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

            String dirType = rootFolder + "MTO-Complex";
            dir = new File(dirType);
            if (!dir.exists()) {
                dir.mkdir();
            }

            int evalsPerRecord = Params.maxEvals / Params.recordsNum;
            String fitnessFile = dirType + "/MTOSOO_P" + (i + 1) + ".txt";
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

            String meanFile = dirMean + "/MTOSOO_P" + (i + 1) + ".txt";
            DataOutputStream outMean = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(meanFile)));
            for (int task = 0; task < taskNum; task++) {
                outMean.writeBytes("" + (task + 1) + ", " + String.format("%.6f", mean[task]) + "\n");
            }
            outMean.close();
        }
    }
}
