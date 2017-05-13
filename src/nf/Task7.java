package nf;

import Jama.Matrix;

public class Task7 {
    private static double alpha = 3.6;
    private static int n = 10;
    private static float eps = 0.0000000001f; //для интегралов
    private static int nIntegr = 8192;

    private static double a = -1, b = 1;

    private static double p(double x){
        return 1;
    }
    private static double q(double x){
        return -(x * x + alpha);
    }
    private static double r(double x){
        return -2 * x;
    }
    private static double f(double x){
        return 2 * (3 * x * x - alpha) / Math.pow((x * x + alpha), 3) + 2 * x / (alpha + 1);
    }
    private static double solution(double x){
        return 1 / (x * x + alpha) - 1 / (alpha + 1);
    }

    private static double chebyshev(int i){
        return Math.cos(Math.PI * (2 * i - 1) / (double)(2 * n));
    }

    private static double P(double x, int i, int k){
        if (k == 0){
            if (i == 0) return 1;
            if (i == 1) return 2 * x;
            if (i > 1) return ((i+1)*(2*i+1)*P(x,i-1,k)*x-(i+1)*i*P(x,i-2,k))/((i+2)*i);
        }
        if (k == 1){
            if (i == 0) return 0;
            if (i == 1) return 2;
            if (i > 1) return ((i+1)*(2*i+1)*(P(x,i-1,0)+x*P(x,i-1,1))-(i+1)*i*P(x,i-2,1))/((i+2)*i);
        }
        if (k == 2){
            if ((i == 0) || (i == 1)) return 0;
            if (i > 1) return ((i+1)*(2*i+1)*(2*P(x,i-1,1)+x*P(x,i-1,2))-(i+1)*i*P(x,i-2,2))/((i+2)*i);
        }

        return 0;
    }

    private static double w(double x, int i, int k){
        if (k == 0) return (1 - x * x) * P(x, i - 1, 0);
        if (k == 1) return -2 * x * P(x, i - 1, 0) + (1 - x * x) * P(x, i - 1, 1);
        if (k == 2) return -2 * P(x, i - 1, 0) - 4 * x * P(x, i - 1, 1) + (1 - x * x) * P(x, i - 1, 2);

        return 0;
    }

    private static double LWi(double x, int i){
        return p(x) * w(x, i, 2) + q(x) * w(x, i, 1) + r(x) * w(x, i, 0);
    }

    private static double approximateSolution(double x, Matrix c){
        double sum = 0;
        for(int i = 0; i < n; i++){
            sum += c.get(i, 0) * w(x, i + 1, 0);
        }
        return sum;
    }

    private static void collocation(){
        System.out.println("\nМетод коллокации:");
        Matrix system = new Matrix(n, n);
        Matrix b = new Matrix(n, 1);
        for(int j = 0; j < n; j++){
            b.set(j, 0, f(chebyshev(j + 1)));
            for(int i = 0; i < n; i ++){
                system.set(j, i, LWi(chebyshev(j + 1),i + 1));
            }
        }
        Matrix c = system.solve(b);

        System.out.println(" тчк       Точное решение           Численное      Модуль разности");
        System.out.println("-0.5: " + solution(-0.5) + " " + approximateSolution(-0.5, c) + " " + Math.abs(solution(-0.5) - approximateSolution(-0.5, c)));
        System.out.println("   0: " + solution(0) + " " + approximateSolution(0, c) + " " + Math.abs(solution(0) - approximateSolution(0, c)));
        System.out.println(" 0.5: " + solution(0.5) + " " + approximateSolution(0.5, c) + " " + Math.abs(solution(0.5) - approximateSolution(0.5, c)));
    }


    //теперь готовимся к голеркину:

    private static double wjfFunction(double x, int j){ //функция под правым интегралом
        return w(x, j, 0) * f(x);
    }
    private static double integralWjf(int j){ //сам правый интеграл. считаем по формуле Симпсона
        double newValue = 0;
        double oldValue = eps;
        long n = nIntegr;
        while(Math.abs(newValue - oldValue) >= eps) {
            oldValue = newValue;
            double h = (b - a) / (double) n / 2;
            newValue = wjfFunction(a, j) + wjfFunction(b, j);
            for (int k = 2; k <= 2 * n - 2; k += 2)
                newValue += 2 * wjfFunction(a + k * h, j);
            for (int k = 1; k <= 2 * n - 1; k += 2)
                newValue += 4 * wjfFunction(a + k * h, j);
            newValue *= h / 3;
            n *= 2;
        }
        return newValue;
    }
    private static double wjLFunction(double x, int i, int j){ //функция под левым интегралом (тот что под суммой)
        return w(x, j,0) * LWi(x, i);
    }
    private static double integralWjL(int i, int j){ //сам левый интеграл (что под суммой). считаем по формуле Симпсона
        double newValue = 0;
        double oldValue = eps;
        long n = nIntegr;
        while(Math.abs(newValue - oldValue) >= eps) {
            oldValue = newValue;
            double h = (b - a) / (double) n / 2;
            newValue = wjLFunction(a, i, j) + wjLFunction(b, i, j);
            for (int k = 2; k <= 2 * n - 2; k += 2)
                newValue += 2 * wjLFunction(a + k * h, i, j);
            for (int k = 1; k <= 2 * n - 1; k += 2)
                newValue += 4 * wjLFunction(a + k * h, i, j);
            newValue *= h / 3;
            n *= 2;
        }
        return newValue;
    }

    private static void galerkin(){
        System.out.println("\nМетод Галеркина:");
        Matrix system = new Matrix(n, n);
        Matrix b = new Matrix(n, 1);
        for(int j = 0; j < n; j++){
            b.set(j, 0, integralWjf(j + 1));
            for(int i = 0; i < n; i ++){
                system.set(j, i, integralWjL(i + 1, j + 1));
            }
        }
        Matrix c = system.solve(b);

        System.out.println(" тчк       Точное решение           Численное      Модуль разности");
        System.out.println("-0.5: " + solution(-0.5) + " " + approximateSolution(-0.5, c) + " " + Math.abs(solution(-0.5) - approximateSolution(-0.5, c)));
        System.out.println("   0: " + solution(0) + " " + approximateSolution(0, c) + " " + Math.abs(solution(0) - approximateSolution(0, c)));
        System.out.println(" 0.5: " + solution(0.5) + " " + approximateSolution(0.5, c) + " " + Math.abs(solution(0.5) - approximateSolution(0.5, c)));
    }

    public static void go(){
        collocation();
        galerkin();
    }
}
