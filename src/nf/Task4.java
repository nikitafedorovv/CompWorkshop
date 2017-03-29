package nf;

import Jama.Matrix;
import MyMath.Integral;
import nf.functions.AlphaTimesBetaFunction;
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
        Matrix fi = new Matrix(n, 1);
        for(int i = 0; i < n; i++){
            fi.set(i, 0, Integral.value((new BetaTimesFFunction(i)), a, b, eps));
        }
        Matrix A = new Matrix(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                A.set(i, j, Integral.value(new AlphaTimesBetaFunction(i, j), a, b, eps));
            }
        }
        Matrix A2 = new Matrix(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == j){
                    A2.set(i, j, 1 - A.get(i, j));
                } else {
                    A2.set(i, j, - A.get(i, j));
                }
            }
        }
        Matrix c = A2.solve(fi);
        c.print(10,10);
    }
}
