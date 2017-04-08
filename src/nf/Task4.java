package nf;

import Jama.Matrix;
import MyMath.Factorial;

public class Task4 {

    private static int n1; //число слагаемых при разложении в ряд тейлора в первом методе
    private static double a = 0;
    private static double b = 1;
    private static double eps = 0.0001; //для второго метода
    private static int n = 20; //по скольки точкам будем сравнивать методы 1 и 2

    private static double epsInt = 0.000000000001; //для интегралов (метод трапеции)


    private static double hFunction(double x, double y){
        return -1.0 / 3.0 * Math.exp( -x * y);
    }

    private static double fFunction(double x){
        return 0.5 * (2 - x + x * x);
    }

    public static double alphaFunction(double x, int index){
        byte sign = (byte)(2 * (index % 2) - 1);
        return sign * (1.0 / 3.0) * Math.pow(x, index) / Factorial.get(index);
    }

    public static double betaFunction(double y, int index){
        return Math.pow(y, index);
    }

    public static double betaTimesFFunction(double y, int index){
        return betaFunction(y, index) * fFunction(y);
    }

    public static double alphaTimesBetaFunction(double y, int i, int j){
        return alphaFunction(y, j) * betaFunction(y, i);
    }

    public static double integralBetaTimesF(int i, double a, double b, double eps){
        double newValue = 0;
        double oldValue = eps;
        long n = 1;
        while(Math.abs(newValue - oldValue) >= eps){
            oldValue = newValue;
            double step = (b - a) / n;
            double point = a;
            double sum = 0;
            while(point < b){
                sum += (betaTimesFFunction(point, i) + betaTimesFFunction(point + step, i)) * 0.5 * step;
                point += step;
            }
            newValue = sum;
            n *= 2;
        }
        return newValue;
    }

    public static double integralAlphaTimesBeta(int i, int j, double a, double b, double eps){
        double newValue = 0;
        double oldValue = eps;
        long n = 1;
        while(Math.abs(newValue - oldValue) >= eps){
            oldValue = newValue;
            double step = (b - a) / n;
            double point = a;
            double sum = 0;
            while(point < b){
                sum += (alphaTimesBetaFunction(point, i, j) + alphaTimesBetaFunction(point + step, i, j)) * 0.5 * step;
                point += step;
            }
            newValue = sum;
            n *= 2;
        }
        return newValue;
    }

    private static double u1(double x, Matrix c){
        double answer = fFunction(x);
        for(int i = 0; i < n1; i++){
            answer += c.get(i, 0) * alphaFunction(x, i);
        }
        return answer;
    }

    private static double u2(double x, Matrix A, Matrix xj, Matrix u){
        int n = A.getRowDimension();
        double answer = fFunction(x);
        for(int j = 0; j < n; j++){
            answer += A.get(j,0) * hFunction(x, xj.get(j,0)) * u.get(j, 0);
        }
        return answer;
    }

    private static Matrix method1(){

        Matrix fi = new Matrix(n1, 1);
        for(int i = 0; i < n1; i++){
            fi.set(i, 0, integralBetaTimesF(i, a, b, epsInt));
        }
        Matrix A = new Matrix(n1, n1);
        for(int i = 0; i < n1; i++){
            for(int j = 0; j < n1; j++){
                A.set(i, j, integralAlphaTimesBeta(i, j, a, b, epsInt));
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
                    system.set(k, j, 1 - A.get(j,0) * hFunction(x.get(k,0), x.get(j,0)));
                } else {
                    system.set(k, j, -A.get(j,0) * hFunction(x.get(k,0), x.get(j,0)));
                }
            }
        }

        Matrix f = new Matrix(n, 1);
        for (int k = 0; k < n; k++) {
            f.set(k, 0, fFunction(x.get(k,0)));
        }

        u = system.solve(f);

        Matrix uAtStartXValue = new Matrix(1,1);
        Matrix uAtCenterXValue = new Matrix(1,1);
        Matrix uAtLastXValue = new Matrix(1,1);
        uAtStartXValue.set(0,0, u2(a, A, x, u));
        uAtCenterXValue.set(0,0, u2((a + b) * 0.5, A, x, u));
        uAtLastXValue.set(0,0, u2(b, A, x, u));

        Matrix[] answer = new Matrix[6];

        answer[0] = uAtStartXValue;
        answer[1] = uAtCenterXValue;
        answer[2] = uAtLastXValue;
        answer[3] = x;
        answer[4] = A;
        answer[5] = u;

        return answer;
    }

    private static Matrix[] method2(){

        System.out.println();

        int n = 4;

        Matrix[] before;
        Matrix[] after = method2cycle(n);
        do {
            n *= 2;
            before = after;
            after = method2cycle(n);

        } while (
                        (Math.abs(before[0].get(0,0) - after[0].get(0,0)) >= eps)
                                ||
                        (Math.abs(before[1].get(0,0) - after[1].get(0,0)) >= eps)
                                ||
                        (Math.abs(before[2].get(0,0) - after[2].get(0,0)) >= eps)
                );
        System.out.println("В методе 2 (мех. кв.) остановка при n = " + n);
        return after;
    }

    public static void go(){
        System.out.println("\nРешаем методом 1, 3 слагаемых в тейлоре...");
        n1 = 3;
        Matrix c3 = method1();

        System.out.println("Решаем методом 1, 4 слагаемых в тейлоре...");
        n1 = 4;
        Matrix c4 = method1();

        System.out.println("Решаем методом 2 (eps = " + eps + ")...");
        Matrix[] method2 = method2();

        Matrix x = method2[3];
        Matrix A = method2[4];
        Matrix u = method2[5];

        double maxDifferenceTaylor = 0;
        double maxDifference = 0;
        double differenceTaylor;
        double difference;
        System.out.println("\nТочка      Зам. ядра 3            Зам. ядра 4            Мех. кв");
        System.out.println("__________________________________________________________________________");


        for(int i = 0; i <= n; i++){
            double currentX = (a * (n - i) + b * i) / n;
            n1 = 3;
            double first3 = u1(currentX, c3);
            n1 = 4;
            double first4 = u1(currentX, c4);
            double second = u2(currentX, A, x, u);

            differenceTaylor = Math.abs(first3 - first4);
            difference = Math.abs(first4 - second);

            System.out.printf("%.3f      ", currentX);
            System.out.printf("%.15f      ", first3);
            System.out.printf("%.15f      ", first4);
            System.out.printf("%.15f      \n", second);

            if (differenceTaylor > maxDifferenceTaylor) {
                maxDifferenceTaylor = differenceTaylor;
            }

            if (difference > maxDifference) {
                maxDifference = difference;
            }
        }
        System.out.println("__________________________________________________________________________");
        System.out.println("\nМаксимальная разница в этих точках между решениями первым методом " +
                "(зам. ядра, разное количество слагаемых " +
                "при разложении в ряд тейлора): " + maxDifferenceTaylor);
        System.out.println("Максимальная разница в этих точках между решением первым методом " +
                "(зам. ядра, 4 слагаемых) и вторым (мех. кв.): " + maxDifference);
    }
}