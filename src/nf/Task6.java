package nf;

import Jama.Matrix;

public class Task6 {
    private static float h = 0.1f;
    private static float tau1 = 0.004f;
    private static float tau2 = 0.02f;
    private static float T = 0.02f; //верхняя граница для t
    private static int N = (int)(1 / h);
    private static int M1 = (int)(T / tau1);
    private static int M2 = (int)(T / tau2);


    private static double exp(double x){
        return Math.exp(x);
    }
    private static double sin(double x){
        return Math.sin(x);
    }
    private static double abs(double x){
        return Math.abs(x);
    }

    private static double solution(double x, double t){
        return exp(0.25 * t) * sin(0.5 * x) - 2 * exp(-t) * (1 - x * x);
        //return x * exp(-t) - x * x * exp(-t) + exp(-0.25 * t) * Math.cos(0.5 * x);
    }
    private static double f(double x, double t){
        return 0.5 * exp(0.25 * t) * sin(0.5 * x) - 2 * exp(-t) * (1 + x * x);
        //return exp(-t) * (x * x - x + 2);
    }
    private static double u0(double x){
        return sin(0.5 * x) - 2 * (1 - x * x);
        //return Math.cos(0.5 * x) + (1 - x) * x;
    }
    private static double u1(double t){
        return -2 * exp(-t);
        //return exp(-0.25 * t);
    }
    private static double u2(double t){
        return exp(0.25 * t) * sin(0.5);
        //return exp(-0.25 * t) * Math.cos(0.5);
    }


    private static void explicitScheme(){
        System.out.println("Явная разностная схема");

        double[] x = new double[N + 1];
        double[] t = new double[M1 + 1];
        for(int i = 0; i <= N; i++){
            x[i] = i * h;
        }
        for(int j = 0; j <= M1; j ++){
            t[j] = j * tau1;
        }

        double[][] u = new double[M1 + 1][N + 1];

        //заполняем левый и правый столбики благодаря граничным условиям
        for(int j = 0; j <= M1; j++){
            u[j][0] = u1(t[j]); //левый
            u[j][N] = u2(t[j]); //правый
        }

        //заполняем нижнюю строчку благодаря граничным условиям
        for(int i = 1; i <= N - 1; i++){
            u[0][i] = u0(x[i]);
        }

        double max = 0;
        //заполняем остальное по формуле явной схемы и сразу считаем погрешность
        for(int j = 0; j <= M1 - 1; j++){
            for(int i = 1; i <= N - 1; i++){
                u[j + 1][i] = (1 - 2 * tau1 / (h * h)) * u[j][i] + (tau1 / (h * h)) * (u[j][i - 1] + u[j][i + 1]) + tau1 * f(x[i], t[j]);
                if(max < abs(u[j + 1][i] - solution(x[i], t[j + 1]))) {
                    max = abs(u[j + 1][i] - solution(x[i], t[j + 1]));
                }
            }
        }

        System.out.println("Погрешность: " + max + "\n");

    }

    private static double[] threeDiagSolve(double[][] system, double[] b){
        double[] uSol = new double[N - 1];  //решение системы
        double[] mi = new double[N - 1]; //прогоночные коэффициенты
        double[] ki = new double[N - 1]; //прогоночные коэффициенты
        //первые прогоночные коэффициенты
        mi[0] = -system[0][1] / system[0][0];
        ki[0] = b[0] / system[0][0];
        //остальные прогоночные коэффициенты
        for(int i = 1; i < N - 2; i++){
            mi[i] = -system[i][i + 1] / (system[i][i - 1] * mi[i - 1] + system[i][i]);
            ki[i] = (b[i] - system[i][i - 1] * ki[i - 1]) / (system[i][i - 1] * mi[i - 1] + system[i][i]);
        }
        //зная коэффициенты, находим решение
        uSol[N - 2] = (b[N - 2] - system[N - 2][N - 3] * ki[N - 3]) / (system[N - 2][N - 3] * mi[N - 3] + system[N - 2][N - 2]);
        for(int i = N - 3; i >= 0; i--){
            uSol[i] = mi[i] * uSol[i + 1] + ki[i];
        }
        return uSol;
    }

    private static void implicitScheme(){
        System.out.println("Неявная разностная схема");

        double[] x = new double[N + 1];
        double[] t = new double[M2 + 1];
        for(int i = 0; i <= N; i++){
            x[i] = i * h;
        }
        for(int j = 0; j <= M2; j ++){
            t[j] = j * tau2;
        }

        double[][] u = new double[M2 + 1][N + 1];
        double[][] F = new double[M2 + 1][N]; //from 1 to N-1, 0 не используем. т е индексы как в теории

        //заполняем левый и правый столбики благодаря граничным условиям
        for(int j = 0; j <= M2; j++){
            u[j][0] = u1(t[j]); //левый
            u[j][N] = u2(t[j]); //правый
        }

        //заполняем нижнюю строчку благодаря граничным условиям
        for(int i = 1; i <= N - 1; i++){
            u[0][i] = u0(x[i]);
        }

        double max = 0;
        for(int j = 0; j <= M2 - 1; j++){
            double[][] A = new double[N - 1][N - 1]; //трёхдиаг. система, которую будем решать
            double[] b = new double[N - 1]; //вектор (Au = b), b[i] = -F[j][i]
            for(int i = 1; i <= N - 1; i++){
                F[j][i] = (1 - tau2 / (h * h)) * u[j][i] + (tau2 / (2 * h * h)) * (u[j][i - 1] + u[j][i + 1]) + tau2 * f(x[i], t[j]);
                b[i - 1] = -F[j][i];
            }

            double aa = tau2 / (2 * h * h);
            double bb = -(1 + tau2 / (h * h));
            double cc = tau2 / (2 * h * h);

//            double aa = 1;
//            double bb = -3;
//            double cc = 1;

            b[0] = -F[j][1] - aa * u[j + 1][0];
            A[0][0] = bb;
            A[0][1] = cc;
            b[N - 2] = -F[j][N - 1] - cc * u[j + 1][N];
            A[N - 2][N - 3] = aa;
            A[N - 2][N - 2] = bb;
            for(int k = 1; k < N - 2; k++){
                A[k][k - 1] = aa;
                A[k][k] = bb;
                A[k][k + 1] = cc;
            }

            double[] uSol = threeDiagSolve(A, b);

            //заполняем пустующие u[j + 1][] и сразу считаем погрешность
            for(int i = 1; i < N - 1; i ++){
                u[j + 1][i] = uSol[i - 1];
                if(max < abs(u[j + 1][i] - solution(x[i], t[j + 1]))) {
                    max = abs(u[j + 1][i] - solution(x[i], t[j + 1]));
                }
                //System.out.println(u[j + 1][i] + " " + solution(x[i], t[j + 1]));
            }
        }
        System.out.println("Погрешность: " + max);
    }

    public static void go(){
        explicitScheme();
        implicitScheme();
    }
}
