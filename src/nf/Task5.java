package nf;

import java.util.Scanner;

public class Task5 {
    private static float alpha = 0.7f;
    private static float a = 0;
    private static float b = 1;
    private static float a0 = 1;
    private static float a1 = -2;
    private static float A = 1 / alpha;
    private static float b0 = 1;
    private static float b1 = 0;
    private static float B = 1 / (1 + alpha);

    private static int n = 100;

    private static double h = (b - a) / n; //шаг

    private static double p(double x) {
        return - (x * x + alpha);
    }
    private static double q(double x) {
        return - 2 * x;
    }
    private static double f(double x) {
        return 2 * (3 * x * x - alpha) / Math.pow((x * x + alpha), 3);
    }
    private static double solution(double x) {
        return 1 / (x * x + alpha);
    }

    // Метод прогонки с основной сеткой
    private static void basicGrid(){
        System.out.println("_____________________________________________________________________________\n");
        System.out.println("Основная сетка:");
        double[] ai = new double[n]; //коэффициенты для системы уравнений - трехдиагональной матрицы
        double[] bi = new double[n]; //коэффициенты для трехдиагональной матрицы
        double[] ci = new double[n]; //коэффициенты для трехдиагональной матрицы
        double[] mi = new double[n]; //прогоночные коэффициенты
        double[] ki = new double[n]; //прогоночные коэффициенты
        double[] y = new double[n + 1];  //приближенные значения решения
        // коэффициенты для 2х уравнений из краевых условий
        double
                p1 = a1 / (a1 - a0 * h),
                p2 = b1 / (b0 * h + b1),
                d1 = -A * h / (a1 - a0 * h),
                d2 = B * h / (b0 * h + b1);
        // находим коэффициенты для остальных уравнений системы
        for (int i = 1; i < n; i++) {
            ai[i] = 1 - h * p(a + i * h) / 2;
            bi[i] = 1 + h * p(a + i * h) / 2;
            ci[i] = 2 - h * h * q(a + i * h);
        }
        //первые прогоночные коэффициенты
        mi[0] = p1;
        ki[0] = d1;

        //остальные прогоночные коэффициенты
        for (int i = 1; i < n; i++) {
            mi[i] = bi[i] / (ci[i] - ai[i] * mi[i - 1]);
            ki[i] = (ai[i] * ki[i - 1] - h * h * f(a + i * h)) / (ci[i] - ai[i] * mi[i - 1]);
        }
        //Зная коэффициенты, находим решение
        y[n] = (p2 * ki[n - 1] + d2) / (1 - p2 * mi[n - 1]);

        for (int i = n - 1; i >= 0; i--)
            y[i] = mi[i] * y[i + 1] + ki[i];

        //выводим значения решения и точного решения в узлах
        System.out.println("Точное значение | Приближенное значение");
        for (int i = 0; i < n + 1; i++)
            System.out.printf("  %.10f  |  %.10f \n", solution(a + i * h), y[i]);

        //ищем погрешность — максимум модулей разности
        double maxabs = Math.abs(solution(a) - y[0]), max = solution(a) - y[0];
        for (int i = 1; i < n + 1; i++)
            if (Math.abs(solution(a + i * h) - y[i]) > maxabs) {
                maxabs = Math.abs(solution(a + i * h) - y[i]);
                max = solution(a + i * h) - y[i];
            }

        System.out.println("\nПогрешность: " + max);

//        //ищем среднюю погрешность
//        double sum = 0;
//        for(int i = 0; i < n + 1; i++){
//            sum += Math.abs(solution(a - h / 2) - y[0]);
//        }
//        System.out.println("Средняя погрешность: " + sum / (n + 1));
    }

    // метод прогонки с использованием расширенной сетки (сдвинутой)
    private static void extendedGrid() {
        System.out.println("\n_____________________________________________________________________________\n");
        System.out.println("Сдвинутая сетка:");

        double[] ai = new double[n + 1]; //коэффициенты для системы уравнений - трехдиагональной матрицы
        double[] bi = new double[n + 1]; //коэффициенты для трехдиагональной матрицы
        double[] ci = new double[n + 1]; //коэффициенты для трехдиагональной матрицы
        double[] mi = new double[n + 1]; //прогоночные коэффициенты
        double[] ki = new double[n + 1]; //прогоночные коэффициенты
        double[] y = new double[n + 2];  //приближенные значения решения
        //коэффициенты для уравнений, которые получаются из граничных условий
        double
                p1 = (2 * a1 + a0 * h) / (2 * a1 - a0 * h),
                p2 = (2 * b1 - b0 * h) / (b0 * h + 2 * b1),
                d1 = -2 * A * h / (2 * a1 - a0 * h),
                d2 = 2 * B * h / (b0 * h + 2 * b1);
        //коэффициенты для остальных уравнений системы
        for (int i = 1; i < n + 1; i++) {
            ai[i] = 1 - h * p(a - h / 2 + i * h) / 2;
            bi[i] = 1 + h * p(a - h / 2 + i * h) / 2;
            ci[i] = 2 - h * h * q(a - h / 2 + i * h);
        }
        //коэффициенты прогонки
        mi[0] = p1;
        ki[0] = d1;
        for (int i = 1; i < n + 1; i++) {
            mi[i] = bi[i] / (ci[i] - ai[i] * mi[i - 1]);
            ki[i] = (ai[i] * ki[i - 1] - h * h * f(a - h / 2 + i * h)) / (ci[i] - ai[i] * mi[i - 1]);
        }
        //обратный прогон - зная коэффициенты, находим решение
        y[n + 1] = (p2 * ki[n] + d2) / (1 - p2 * mi[n]);
        for (int i = n; i >= 0; i--)
            y[i] = mi[i] * y[i + 1] + ki[i];

        //выводим значения решения и точного решения в узлах
        System.out.println("Точное значение | Приближенное значение");
        for (int i = 0; i < n + 2; i++)
            System.out.printf("  %.10f  |  %.10f \n", solution(a - h / 2 + i * h), y[i]);

        //погрешность - максимум модулей разности
        double maxa = Math.abs(solution(a - h / 2) - y[0]), max = solution(a - h / 2) - y[0];
        for (int i = 1; i < n + 1; i++) {
            if (Math.abs(solution(a - h / 2 + i * h) - y[i]) > maxa) {
                maxa = Math.abs(solution(a - h / 2 + i * h) - y[i]);
                max = solution(a - h / 2 + i * h) - y[i];
            }
        }
        System.out.println("\nПогрешность: " + max);
//        double sum = 0;
//        for(int i = 0; i < n + 1; i++){
//            sum += Math.abs(solution(a - h / 2) - y[0]);
//        }
//        System.out.println("Средняя погрешность: " + sum / (n + 1));
    }

    public static void go() {
        System.out.println("Разностный метод решения краевой задачи для обыкновенного дифференциального уравнения второго порядка. Метод прогонки.");
        basicGrid();
        extendedGrid();
    }
}