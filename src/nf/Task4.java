package nf;

import Jama.Matrix;
import MyMath.Integral;
import nf.functions.*;

public class Task4 {

    private static int n1 = 3; //число слагаемых при разложении в ряд тейлора
    private static double a = 0;
    private static double b = 1;
    private static int accuracy = 12;
    private static double eps = Math.pow(0.1, accuracy + 2);
    private static double eps1 = Math.pow(0.1, 5); //для второго
    private static int n = 10; //по скольки точкам будем сравнивать методы 1 и 2


    private static double u1(double x, Matrix c){
        FFunction f = new FFunction();
        double answer = f.value(x);
        for(int i = 0; i < n1; i++){
            answer += c.get(i, 0) * (new AlphaFunction(i)).value(x);
        }
        return answer;
    }

    private static double u2(double x, Matrix A, Matrix xj, Matrix u){
        int n = A.getRowDimension();
        FFunction f = new FFunction();
        HFunction H = new HFunction();
        double answer = f.value(x);
        for(int j = 0; j < n; j++){
            answer += A.get(j,0) * H.value(x, xj.get(j,0)) * u.get(j, 0);
        }
        return answer;
    }

    private static Matrix method1(){
        Matrix fi = new Matrix(n1, 1);
        for(int i = 0; i < n1; i++){
            fi.set(i, 0, Integral.value((new BetaTimesFFunction(i)), a, b, eps));
        }
        Matrix A = new Matrix(n1, n1);
        for(int i = 0; i < n1; i++){
            for(int j = 0; j < n1; j++){
                A.set(i, j, Integral.value(new AlphaTimesBetaFunction(i, j), a, b, eps));
            }
        }
        Matrix system = new Matrix(n1, n1);
        for(int i = 0; i < n1; i++){
            for(int j = 0; j < n1; j++){
                if(i == j){
                    system.set(i, j, 1 - A.get(i, j));
                } else {
                    system.set(i, j, - A.get(i, j));
                }
            }
        }
        Matrix c = system.solve(fi);

        return c;
    }

    private static Matrix[] method2cycle(int n){
        Matrix x;
        Matrix A;
        Matrix u;

        double step = (b - a) / n;
        x = new Matrix(n,1);
        for (int i = 0; i < n; i++) {
            x.set(i,0,a + i * step);
        }

        A = new Matrix(n,1);
        A.set(0,0, step / 2);
        A.set(n - 1,0,step / 2);
        for (int i = 1; i < n - 1; i++) {
            A.set(i,0, step);
        }

        Matrix system = new Matrix(n, n);
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                if (k == j) {
                    system.set(k, j, 1 - A.get(j,0) * (new HFunction().value(x.get(k,0), x.get(j,0))));
                } else {
                    system.set(k, j, -A.get(j,0) * (new HFunction().value(x.get(k,0), x.get(j,0))));
                }
            }
        }

        FFunction fFunction = new FFunction();

        Matrix f = new Matrix(n, 1);
        for (int k = 0; k < n; k++) {
            f.set(k, 0, fFunction.value(x.get(k,0)));
        }

        u = system.solve(f);

        Matrix uInFirstXValue = new Matrix(1,1);
        Matrix uInSecondXValue = new Matrix(1,1);
        Matrix uInThirdValue = new Matrix(1,1);
        uInFirstXValue.set(0,0, u2(x.get(0,0), A, x, u));
        uInSecondXValue.set(0,0, u2(x.get(1,0), A, x, u));
        uInThirdValue.set(0,0, u2(x.get(2,0), A, x, u));

        Matrix[] answer = new Matrix[6];

        answer[0] = uInFirstXValue;
        answer[1] = uInSecondXValue;
        answer[2] = uInThirdValue;
        answer[3] = x;
        answer[4] = A;
        answer[5] = u;

        return answer;
    }

    private static Matrix[] method2(){

        int n = 4;

        Matrix[] before;
        Matrix[] after = method2cycle(n);
        do {
            n *= 2;
            before = after;
            after = method2cycle(n);
        } while (
                        (Math.abs(before[0].get(0,0) - after[0].get(0,0)) < eps1)
                                &
                        (Math.abs(before[1].get(0,0) - after[1].get(0,0)) < eps1)
                                &
                        (Math.abs(before[2].get(0,0) - after[2].get(0,0)) < eps1)
                );

        return after;
    }

    public static void go(){
        Matrix c = method1();

        Matrix[] method2 = method2();

        Matrix x = method2[1];
        Matrix A = method2[2];
        Matrix u = method2[3];

        double step = (b - a) / n;
        double maxDifference = 0;
        for(int i = 0; i < n; i++){
            double currentX = a + step * i;
            double difference = Math.abs(u1(currentX, c) - u2(currentX, A, x, u));
            if (difference > maxDifference) {
                maxDifference = difference;
            }
        }
        System.out.println("maxDifference: " + maxDifference);
    }
}
