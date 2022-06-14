package RWS;

import basic.Individual;
import basic.Params;

import java.util.ArrayList;
//RWS: https://www.researchgate.net/publication/359964152_RWS-L-SHADE_An_Effective_L-SHADE_Algorithm_Incorporation_Roulette_Wheel_Selection_Strategy_for_Numerical_Optimisation
public class RWS {
    public static int runRWS(ArrayList<Individual> Individuals){
        double F[] = new double[Individuals.size()];
        double P[] = new double[Individuals.size()];
        double maxF=0,sumF=0;
        for(int i=0;i<Individuals.size();i++){
            if(maxF<Individuals.get(i).getFitness()){
                maxF = Individuals.get(i).getFitness();
            }
        }
        for (int i=0;i<Individuals.size();i++){
            F[i]=maxF-Individuals.get(i).getFitness();
            sumF+=F[i];
        }
        for (int i=0;i<Individuals.size();i++){
            F[i] = F[i]/sumF;
        }
        P[0]=F[0];
        for (int i=1;i<Individuals.size();i++){
            P[i]=P[i-1]+F[i];
        }
        double r = Params.rand.nextDouble();
        for (int i=0;i<Individuals.size();i++){
            if(i==0){
                if(r<=P[i]) return i;
            }else {
                if(r>P[i-1] && r<=P[i]) return i;
            }
        }
        return Params.rand.nextInt(Individuals.size());
    }
}
