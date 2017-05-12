package nf;

import Jama.Matrix;

public class Task7 {
    private static double alpha = 3.6;
    private static int n = 10; //для метода коллокации

    //коэффициенты граничных условий
    private static double a1 = 1;
    private static double a2 = 0;
    private static double b1 = 1;
    private static double b2 = 0;

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

    private static double chebyshev(int i){ //i=1..n
        return Math.cos(Math.PI * (2 * i - 1) / (double)(2 * n));
    }

    private static double P(double x, int i, int k){ //Л i=0..
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

        System.out.println("Чёт с индексами i/k намутил в P(double x, int i, int k) :(" + i + " " + k);
        return 0;
    }

    private static double w(double x, int i, int k){ //НП i=1..n
        //i++;
        if (k == 0) return (1 - x * x) * P(x, i - 1, 0);
        if (k == 1) return -2 * x * P(x, i - 1, 0) + (1 - x * x) * P(x, i - 1, 1);
        if (k == 2) return -2 * P(x, i - 1, 0) - 4 * x * P(x, i - 1, 1) + (1 - x * x) * P(x, i - 1, 2);

        System.out.println("Чёт с производной намутил в w(double x, int i, int k) :(");
        return 0;
    }

    private static double LWi(double x, int i){//i=1..n
        return p(x) * w(x, i, 2) + q(x) * w(x, i, 1) + r(x) * w(x, i, 0);
    }

    private static double collocSolution(double x, Matrix c){
        double sum = 0;
        for(int i = 0; i < n; i++){
            sum += c.get(i, 0) * w(x, i + 1, 0);
        }
        return sum;
    }

    private static void collocation(){
        Matrix system = new Matrix(n, n);
        Matrix b = new Matrix(n, 1);
        for(int j = 0; j < n; j++){
            b.set(j, 0, f(chebyshev(j + 1)));
            for(int i = 0; i < n; i ++){
                system.set(j, i, LWi(chebyshev(j + 1),i + 1));
            }
        }
        Matrix c = system.solve(b);

        System.out.println("-0.5: " + solution(-0.5) + " " + collocSolution(-0.5, c) + " " + Math.abs(solution(-0.5) - collocSolution(-0.5, c)));
        System.out.println("   0: " + solution(0) + " " + collocSolution(0, c) + " " + Math.abs(solution(0) - collocSolution(0, c)));
        System.out.println(" 0.5: " + solution(0.5) + " " + collocSolution(0.5, c) + " " + Math.abs(solution(0.5) - collocSolution(0.5, c)));

    }

    private static void galerkin(){

    }

    public static void go(){
        collocation();
        galerkin();
    }
}
