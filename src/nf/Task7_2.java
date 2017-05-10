package nf;

import java.util.Scanner;

public class Task7_2 {
    private static float alpha = 3.6f;
    private static double[][] c = new double[100][1];
    private static int[][] per = new int[100][1];
    private final static double pi = Math.PI;

    // полином Якоби:
    // заполняем согласно формулам:
    private static double P(int n, double x) {
        if (n == 0)
            return 1;
        else if (n == 1)
            return 2 * x;
        else
            return ((n + 1) * (2 * n + 1) * x * P(n - 1, x) - (n + 1) * n * P(n - 2, x)) / ((n + 2) * n);
    }

    // первая производная полинома Якоби
    private static double derP1(int n, double x) {
        if (n == 0)
            return 0;
        else if (n == 1)
            return 2;
        else
            return ((n + 1) * (2 * n + 1)) / ((n + 2) * n) * (P(n - 1, x) + x * derP1(n - 1, x)) - (n + 1) / (n + 2) * derP1(n - 2, x);
    }

    // вторая производная полинома Якоби
    private static double derP2(int n, double x) {
        if (n == 0)
            return 0;
        else if (n == 1)
            return 0;
        else
            return ((n + 1) * (2 * n + 1)) / ((n + 2) * n) * (2 * derP1(n - 1, x) + x * derP2(n - 1, x)) - (n + 1) / (n + 2) * derP2(n - 2, x);
    }

    //
    private static double q(double x) {
        return -(x * x + 3.6);
    }

    private static double r(double x) {
        return -2 * x;
    }

    // f
    private static double f(double x) {
        return 2 * (3 * x * x - 3.6) / Math.pow((x * x + 3.6), 3) + 2 * x / 4.6;
    }

    // точное решение
    private static double sol(double x) {
        return 1.0 / (x * x + 3.6) - 1.0 / 4.6;
    }

    // L(Wi)
    private static double L(int n, double x) {
        return -2 * P(n - 1, x) - 4 * x * derP1(n - 1, x) + (1 - x * x) * derP2(n - 1, x) + q(x) * (-2 * x * P(n - 1, x) + (1 - x * x) * derP1(n - 1, x)) + r(x) * (1 - x * x) * P(n - 1, x);
    }

    // координатная сетка
    private static double w(int n, double x) {
        return (1 - x * x) * P(n - 1, x);
    }

    /*
    Следующие 3 метода из предыдущей задачи, для решения системы линейных уравнений
    */
    private static int MaxTabX(double a[][], int n, int x) {
        int maxx = x, maxy = x;
        for (int i = x; i < n; i++)
            for (int j = x; j < n; j++)
            {
                if (a[i][j] > a[maxx][maxy])
                {
                    maxx = i;
                    maxy = j;
                }
            }
        return maxx;
    }

    private static int MaxTabY(double a[][], int n, int x) {
        int maxx = x, maxy = x;
        for (int i = x; i < n; i++)
            for (int j = x; j < n; j++) {
                if (a[i][j] > a[maxx][maxy]) {
                    maxx = i;
                    maxy = j;
                }
            }
        return maxy;
    }

    private static void mainElTab(double a[][], double b[][], int n) {
        for (int i = 0; i < 100; i++)
            per[i][0] = i;
        boolean t = false, t1 = true;
        for (int i = 0; i < n; i++) {
            int b1 = MaxTabX(a, n, i);
            int b2 = MaxTabY(a, n, i);
            for (int j = i; j < n; j++) {
                double y = a[b1][j];
                a[b1][j] = a[i][j];
                a[i][j] = y;
            }
            double y = b[b1][0];
            b[b1][0] = b[i][0];
            b[i][0] = y;
            for (int j = 0; j < n; j++) {
                double y0 = a[j][b2];
                a[j][b2] = a[j][i];
                a[j][i] = y0;
            }
            int y1 = per[i][0];
            per[i][0] = per[b2][0];
            per[b2][0] = y1;
            if (a[i][i] != 0) {
                double q = 1.0 / a[i][i];
                for (int j = i; j < n; j++)
                    a[i][j] = a[i][j] * q;
                b[i][0] = b[i][0] * q;
                for (int i1 = i + 1; i1 < n; i1++) {
                    double p = a[i1][i];
                    for (int j = i; j < n; j++) {
                        a[i1][j] = a[i1][j] - a[i][j] * p;
                    }
                    b[i1][0] = b[i1][0] - b[i][0] * p;
                }
            }
        }
        if (a[n - 1][n - 1] == 0 && b[n - 1][0] != 0) {
            t1 = false;
        } else if (a[n - 1][n - 1] == 0 && b[n - 1][0] == 0) {
            c[per[n - 1][0]][0] = 1.0;
            t = true;
            t1 = true;
        } else {
            c[per[n - 1][0]][0] = b[n - 1][0];
            t1 = true;
            t = false;
        } for (int i = n - 2; i >= 0; i--) {
            if (a[i][i] != 0) {
                c[per[i][0]][0] = b[i][0];
                for (int j = i + 1; j < n; j++)
                    c[per[i][0]][0] -= a[i][j] * c[per[j][0]][0];
            } else {
                double a1 = b[i][0];
                for (int j = i + 1; j < n; j++)
                    a1 -= a[i][j] * c[per[j][0]][0];
                if (a1 == 0) {
                    c[per[i][0]][0] = 1.0;
                    t = true;
                    t1 = true;
                } else {
                    t1 = false;
                }
            }
        }
	/*cout << "per: ";
	for (int i = 0; i < n; i++)
	cout << per[i][0] + 1 << ' ';
	cout << endl;*/
    }


    // 1 способ: Метод коллокации
    private static void collocation(int n)
    {
        System.out.println("Решение методом коллокации:");
        // инициализация
        double[][] A = new double[100][100];
        double[][] b = new double[100][1];

        for (int i = 0; i < 100; i++)
            for (int j = 0; j < 100; j++)
                A[i][j] = 0;
        for (int i = 0; i < 100; i++)
            b[i][0] = 0;
        // считаем в точке == КОРЕНЬ МН-НА ЧЕБЫШЕВА
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = L(j + 1, Math.cos(pi * (2.0 * (i + 1) - 1) / (2 * n)));

        for (int i = 0; i < n; i++)
            b[i][0] = f(Math.cos(pi * (2.0 * (i + 1) - 1) / (2 * n)));

	/*for (int i = 0; i < n; i++)
	{
	for (int j = 0; j < n; j++)
	cout << A[i][j] << ' ';
	cout << endl;
	}
	for (int i = 0; i < n; i++)
	cout << b[i][0] << endl;*/

        mainElTab(A, b, n);

	/*for (int i = 0; i < n; i++)
	cout << "c[" << i + 1 << "] = " << setprecision(15) << c[i][0] << endl;*/

        double a = 0.0;
        for (int i = 0; i < n; i++)
            a += c[i][0] * w(i + 1, -0.5);
        System.out.println("Приближенное значение в (.) (-0.5): " + a);
        System.out.println("Точное значение в (.) (-0.5): " + sol(-0.5));
        System.out.println("Погрешность в (.) (-0.5): " + Math.abs(a - sol(-0.5)));

        a = 0.0;
        for (int i = 0; i < n; i++)
            a += c[i][0] * w(i + 1, 0.0);
        System.out.println("Приближенное значение в (.) (0): " + a);
        System.out.println("Точное значение в (.) (0): " + sol(0.0));
        System.out.println("Погрешность в (.) (0): " + Math.abs(a - sol(0)));
        a = 0.0;
        for (int i = 0; i < n; i++)
            a += c[i][0] * w(i + 1, 0.5);
        System.out.println("Приближенное значение в (.) (0.5): " + a);
        System.out.println("Точное значение в (.) (0.5): " + sol(0.5));
        System.out.println("Погрешность в (.) (0.5): " + Math.abs(a - sol(0.5)));
    }

    //считаем интеграл (метод Симпсона)
    private static double integration(int i1, int j1, boolean f1) {
        double h = 2.0 / 200;
        double a = 0.0;
        if (f1) {
            for (int i = 0; i <= 200; i++) {
                if (i == 0 || i == 200)
                    a += w(i1, i * h - 1) * f(i * h - 1);
                else if (i % 2 == 1)
                    a += 4 * w(i1, i * h - 1) * f(i * h - 1);
                else if (i % 2 == 0)
                    a += 2 * w(i1, i * h - 1) * f(i * h - 1);
            }
            a *= h / 3;
        } else {
            for (int i = 0; i <= 200; i++) {
                if (i == 0 || i == 200)
                    a += w(i1, i * h - 1) * L(j1, i * h - 1);
                else if (i % 2 == 1)
                    a += 4 * w(i1, i * h - 1) * L(j1, i * h - 1);
                else if (i % 2 == 0)
                    a += 2 * w(i1, i * h - 1) * L(j1, i * h - 1);
            }
            a *= h / 3;
        }
        return a;
    }

    //определитель матрицы
    private static double determinant(double mas[][], int n) {
        if (n == 2) {
            return mas[0][0] * mas[1][1] - mas[0][1] * mas[1][0];
        } else {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                double mas1[][] = new double[100][100];
                int i2 = 0, j2 = 0;
                for (int i1 = 1; i1 < n; i1++) {
                    for (int j1 = 0; j1 < n; j1++)
                        if (j1 != j) {
                            mas1[i2][j2] = mas[i1][j1];
                            j2++;
                        }
                    i2++;
                    j2 = 0;
                }
                if (j % 2 == 0)
                    sum += mas[0][j] * determinant(mas1, n - 1);
                else
                    sum -= mas[0][j] * determinant(mas1, n - 1);
            }
            return sum;
        }
    }

    private static void Gal(int n) {
        System.out.println("Решение методов Галеркина:");
        //инициализация
        double[][] A = new double[100][100];
        double[][] b = new double[100][1];
        for (int i = 0; i < 100; i++)
            c[i][0] = 0;
        for (int i = 0; i < 100; i++)
            for (int j = 0; j < 100; j++)
                A[i][j] = 0;
//        for (int i = 0; i < 100; i++)
//            b[i][0];  //!!!!!!!!!!??????!!!!!!!!!!!!???!!!!!!!!!!!!!!!!!!!!!!????!!!!!!!!!!!!!!!?????!!!!!!!!!!!!!!!!!

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                A[i][j] = integration(i + 1, j + 1, false);
        for (int i = 0; i < n; i++)
            b[i][0] = integration(i + 1, 0, true);
        double z = determinant(A, n);
        if (z == 0.0) {
            System.out.println("Нет решений");
            return;
        }
        //решаем
        mainElTab(A, b, n);
	/*double a = 0.0;
	for (int i = 0; i < n; i++)
	a += c[i][0] * (1 - 0.25) * P(i, -0.5);
	double y = abs(a - sol(-0.5));
	a = 0.0;
	for (int i = 0; i < n; i++)
	a += c[i][0] * P(i, 0.0);
	if (abs(a - sol(0.0)) > y)
	y = abs(a - sol(0.0));
	a = 0.0;
	for (int i = 0; i < n; i++)
	a += c[i][0] * (1 - 0.25) * P(i, 0.5);
	if (abs(a - sol(0.5)) > y)
	y = abs(a - sol(0.5));
	while (y > e)
	{
	n++;
	for (int i = 0; i < 100; i++)
	c[i][0] = 0;
	for (int i = 0; i < 100; i++)
	for (int j = 0; j < 100; j++)
	A[i][j] = 0;
	for (int i = 0; i < 100; i++)
	b[i][0];
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	A[i][j] = integration(i + 1, j + 1, false);
	for (int i = 0; i < n; i++)
	b[i][0] = integration(i + 1, 0, true);
	double z = determinant(A, n);
	if (z == 0.0)
	{
	cout << "can't solve by this method" << endl;
	return;
	}
	mainElTab(A, b, n);
	double a = 0.0;
	for (int i = 0; i < n; i++)
	a += c[i][0] * (1 - 0.25) * P(i, -0.5);
	double y = abs(a - sol(-0.5));
	a = 0.0;
	for (int i = 0; i < n; i++)
	a += c[i][0] * P(i, 0.0);
	if (abs(a - sol(0.0)) > y)
	y = abs(a - sol(0.0));
	a = 0.0;
	for (int i = 0; i < n; i++)
	a += c[i][0] * (1 - 0.25) * P(i, 0.5);
	if (abs(a - sol(0.5)) > y)
	y = abs(a - sol(0.5));
	}*/
        double a = 0.0;
        for (int i = 0; i < n; i++)
            a += c[i][0] * w(i + 1, -0.5);
        System.out.println("Приближенное значение в (.) (-0.5): " + a);
        System.out.println("Точное значение в (.) (-0.5): " + sol(-0.5));
        System.out.println("Погрешность в (.) (-0.5): " + Math.abs(a - sol(-0.5)));
        a = 0.0;
        for (int i = 0; i < n; i++)
            a += c[i][0] * w(i + 1, 0.0);
        System.out.println("Приближенное значение в (.) (0): " + a);
        System.out.println("Точное значение в (.) (0): " + sol(0.0));
        System.out.println("Погрешность в (.) (0): " + Math.abs(a - sol(0)));
        a = 0.0;
        for (int i = 0; i < n; i++)
            a += c[i][0] * w(i + 1, 0.5);
        System.out.println("Приближенное значение в (.) (0.5): " + a);
        System.out.println("Точное значение в (.) (0.5): " + sol(0.5));
        System.out.println("Погрешность в (.) (0.5): " + Math.abs(a - sol(0.5)));
    }

    public static void go() {
        System.out.println("Методы решения ОДУ в аналитическом виде:");
        int n;
        double e;
        System.out.println("Кол-во функций: ");
        Scanner sc = new Scanner(System.in);
        n = sc.nextInt();
        collocation(n);
        //cout << "enter the accuracy: "; cin >> e;
        Gal(n);
        //cout << setprecision(15) << integration(1, 1, true) << endl;
        //cout << w(1, -0.98) * f(-0.98) * 4 << endl;
    }
}
