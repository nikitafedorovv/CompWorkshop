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

    private static double p(double x){
        return - (x * x + alpha);
    }
    private static double q(double x){
        return - 2 * x;
    }
    private static double f(double x){
        return 2 * (3 * x * x - alpha) / Math.pow((x * x + alpha), 3);
    }
    private static double solution(double x){
        return 1 / (x * x + alpha);
    }

    // Метод прогонки с основной сеткой
    private static void basicGrid(){
        double h = (b - a) / n; //шаг
        double[] ai = new double[n + 2];
        double[] bi = new double[n + 2];
        double[] ci = new double[n + 2];
        double[] mi = new double[n + 2];
        double[] ki = new double[n + 2];
        double[] y = new double[n + 2];
        // коэффициенты для 2х уравнений из краевых условий
        double p1 = a1 / (a1 - a0 * h), p2 = b1 / (b0 * h + b1), d1 = -A * h / (a1 - a0 * h), d2 = B * h / (b0 * h + b1);
        // находим коэффициенты для системы уравнений
        for (int i = 1; i < n; i++)
        {
            ai[i] = 1 - h * p(a + i * h) / 2;
            bi[i] = 1 + h * p(a + i * h) / 2;
            ci[i] = 2 - h * h * q(a + i * h);
        }
        //прогоночные коэффициенты
        mi[0] = p1;
        ki[0] = d1;

        for (int i = 1; i < n; i++)
        {
            mi[i] = bi[i] / (ci[i] - ai[i] * mi[i - 1]);
            ki[i] = (ai[i] * ki[i - 1] - h * h * f(a + i * h)) / (ci[i] - ai[i] * mi[i - 1]);
        }
        //Обратный прогон, то есть, зная коэффициенты, находим решение
        y[n] = (p2 * ki[n - 1] + d2) / (1 - p2 * mi[n - 1]);

        for (int i = n - 1; i >= 0; i--)
            y[i] = mi[i] * y[i + 1] + ki[i];

        //находим значения решения и точного решения в узлах
        for (int i = 0; i < n + 1; i++)
            System.out.println(solution(a + i * h) + " " + y[i]);

        double maxa = Math.abs(solution(a) - y[0]), max = solution(a) - y[0];

        //ищем максимальную погрешность
        for (int i = 1; i < n + 1; i++)
            if (Math.abs(solution(a + i * h) - y[i]) > maxa)
            {
                maxa = Math.abs(solution(a + i * h) - y[i]);
                max = solution(a + i * h) - y[i];
            }

        System.out.println("Максимальная погрешность: " + max);
    }

    // метод прогонки с использованием расширенной сетки (сдвинутой)
    private static void extendedGrid()
    {
        System.out.println("Расширенная сетка:");
        //шаг
        double h = (b - a) / n;
        double[] ai = new double[n + 2];
        double[] bi = new double[n + 2];
        double[] ci = new double[n + 2];
        double[] mi = new double[n + 2];
        double[] ki = new double[n + 2];
        double[] y = new double[n + 2];
        //коэффициенты для уравнений, которые получаются из граничных условий
        double p1 = (2 * a1 + a0 * h) / (2 * a1 - a0 * h), p2 = (2 * b1 - b0 * h) / (b0 * h + 2 * b1), d1 = -2 * A * h / (2 * a1 - a0 * h), d2 = 2 * B * h / (b0 * h + 2 * b1);
        //остальные уравнения системы
        for (int i = 1; i < n + 1; i++)
        {
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
        //обратный прогон
        y[n + 1] = (p2 * ki[n] + d2) / (1 - p2 * mi[n]);
        for (int i = n; i >= 0; i--)
            y[i] = mi[i] * y[i + 1] + ki[i];
        //значения в узлах
        for (int i = 0; i < n + 2; i++)
            System.out.println(solution(a - h / 2 + i * h) + " " + y[i] );
        //максимальная погрешность
        double maxa = Math.abs(solution(a - h / 2) - y[0]), max = solution(a - h / 2) - y[0];
        for (int i = 1; i < n + 1; i++) {
            if (Math.abs(solution(a - h / 2 + i * h) - y[i]) > maxa) {
                maxa = Math.abs(solution(a - h / 2 + i * h) - y[i]);
                max = solution(a - h / 2 + i * h) - y[i];
            }
        }
        System.out.println("Максимальная погрешность: " + max);
    }

    public static void go() {
        System.out.println("Разностный метод решения краевой задачи для обыкновенного дифференциального уравнения второго порядка. Метод прогонки.");
        int n;
        System.out.println("Введите n: ");
        Scanner sc = new Scanner(System.in);
        n = sc.nextInt();
        basicGrid();
        extendedGrid();
    }
}