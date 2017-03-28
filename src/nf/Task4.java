package nf;

import MyMath.Function;
import MyMath.Integral;
import nf.functions.BetaFunction;
import nf.functions.BetaTimesFFunction;
import nf.functions.FFunction;

public class Task4 {
    public static void go(){
        int n = 3;
        double a = 0;
        double b = 1;
        int accuracy = 10;
        double eps = Math.pow(0.1, accuracy);

        FFunction f = new FFunction();
        double[] fi = new double[n];
        for(int i = 0; i < n; i++){
            BetaTimesFFunction betaTimesF = new BetaTimesFFunction(i);
            fi[i] = Integral.value(betaTimesF, a, b, eps);
        }
    }
}
