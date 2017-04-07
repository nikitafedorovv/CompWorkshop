package nf;

import Jama.Matrix;
import MyMath.Integral;
import nf.functions.*;

public class Task4 {

    private static int n1; //число слагаемых при разложении в ряд тейлора
    private static double a = 0;
    private static double b = 1;
    private static int accuracy = 12;
    private static double eps = Math.pow(0.1, accuracy + 2); //по сути, это для подсчета интеграла в первом
    private static double eps1 = 0.00001; //для второго //5!
    private static int n = 20; //по скольки точкам будем сравнивать методы 1 и 2

    private static double u1(double x, Matrix c){
        double answer = new FFunction().value(x);
        for(int i = 0; i < n1; i++){
            answer += c.get(i, 0) * new AlphaFunction(i).value(x);
        }
        return answer;
    }

    private static double u2(double x, Matrix A, Matrix xj, Matrix u){
        int n = A.getRowDimension();
        double answer = new FFunction().value(x);
        for(int j = 0; j < n; j++){
            answer += A.get(j,0) * new HFunction().value(x, xj.get(j,0)) * u.get(j, 0);
            //System.out.println("a = "+A.get(j,0)+", h = "+H.value(x, xj.get(j,0))+", u = "+u.get(j, 0));
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
                A.set(i, j, Integral.value(new AlphaTimesBetaFunction(j, i), a, b, eps));
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
        A.set(0,0, step * 0.5);
        A.set(n - 1,0,step * 0.5);
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

        Matrix f = new Matrix(n, 1);
        for (int k = 0; k < n; k++) {
            f.set(k, 0, new FFunction().value(x.get(k,0)));
        }

        u = system.solve(f);

        Matrix uAtStartXValue = new Matrix(1,1);
        Matrix uAtCenterXValue = new Matrix(1,1);
        Matrix uAtlastXValue = new Matrix(1,1);
        uAtStartXValue.set(0,0, u2(a, A, x, u));
        uAtCenterXValue.set(0,0, u2((a + b) * 0.5, A, x, u));
        uAtlastXValue.set(0,0, u2(b, A, x, u));

        Matrix[] answer = new Matrix[6];

        answer[0] = uAtStartXValue;
        answer[1] = uAtCenterXValue;
        answer[2] = uAtlastXValue;
        answer[3] = x;
        answer[4] = A;
        answer[5] = u;

        return answer;
    }

    private static Matrix[] method2(){

        int n = 3;

        Matrix[] before;
        Matrix[] after = method2cycle(n);
        do {
            n *= 2;
            before = after;
            after = method2cycle(n);

//            System.out.println(n +" "+after[0].get(0,0));
//            System.out.println(n +" "+ after[1].get(0,0));
//            System.out.println(n +" "+ after[2].get(0,0));

        } while (
                        (Math.abs(before[0].get(0,0) - after[0].get(0,0)) >= eps1)
                                ||
                        (Math.abs(before[1].get(0,0) - after[1].get(0,0)) >= eps1)
                                ||
                        (Math.abs(before[2].get(0,0) - after[2].get(0,0)) >= eps1)
                );
        System.out.println("В методе 2 (мех. кв.) остановка при n = " + n);
        return after;
    }

    public static void go(){
        n1 = 3;
        Matrix c3 = method1();

        n1 = 4;
        Matrix c4 = method1();

        Matrix[] method2 = method2();

        Matrix x = method2[3];
        Matrix A = method2[4];
        Matrix u = method2[5];

        //double maxDifference = 0;
        for(int i = 0; i <= n; i++){
            double currentX = (a * (n - i) + b * i) / n;
            n1 = 3;
            double first3 = u1(currentX, c3);
            n1 = 4;
            double first4 = u1(currentX, c4);
            double second = u2(currentX, A, x, u);
            //double difference = Math.abs(first3 - second);
            System.out.println("Точка: " + currentX);
            System.out.println("Значение решения первым методом (зам. ядра) в этой точке, в тейлоре три слагаемых: " + first3);
            System.out.println("Значение решения первым методом (зам. ядра) в этой точке, в тейлоре четыре слагаемых: " + first4);
            System.out.println("Значение решения вторым методом (мех. кв.) в этой точке: " + second);
            //System.out.println("Разница: "+difference);

//            if (difference > maxDifference) {
//                maxDifference = difference;
//            }
        }
        //System.out.println("maxDifference: " + maxDifference);
    }
}