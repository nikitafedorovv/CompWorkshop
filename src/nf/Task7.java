package nf;

import java.util.Scanner;

public class Task7 {

    private static void swap(double[][] A, int i, int j, int k, int p){
        double c = A[i][j];
        A[i][j] = A[k][p];
        A[k][p] = c;
    }

    //решение
    private static double F(double x, int alpha) {
        return 1.0/(x*x+alpha)-1.0/(alpha+1);
    }
    //полином Якоби считаем
    private static double P(double x, int i, int k) {
        if (k==0) {
            if (i==0)
                return 1;
            if (i == 1)
                return 2*x;
            else
                return ((i+1)*(2*i+1)*P(x,i-1,k)*x-(i+1)*i*P(x,i-2,k))/(double)(i+2)/(double)(i);
        }

        if (k==1) {
            if (i==0)
                return 0;
            if (i==1)
                return 2;
            else
                return ((i+1)*(2*i+1)*(P(x,i-1,0)+x*P(x,i-1,1))-(i+1)*i*P(x,i-2,1))/(i+2)/i;
        }

        else { //k=0/1/2
            if (i==0)
                return 0;
            if (i==1)
                return 0;
            else
                return ((i+1)*(2*i+1)*(2*P(x,i-1,1)+x*P(x,i-1,2))-(i+1)*i*P(x,i-2,2))/(i+2)/i;
        }
    }

    //сетка
    private static double w( double x, int i, int k) { //функциональная сетка, k-порядок производной
        if (k==0)
            return (1-x*x)*P(x,i,k);
        if (k==1)
            return -2*x*P(x,i,0)+(1-x*x)*P(x,i,1);
        else
            return -2*P(x,i,0)-2*x*P(x,i,1)-2*x*P(x,i,1)+(1-x*x)*P(x,i,2);
    }

    //Для решения системы линейных уравнений
    private static int find_max_stol(double[][] A,int n,int k)
    {
        int max_st=k;
        double max=Math.abs(A[k][k]);
        for (int i=k+1;i<n;i++)
            if (max<Math.abs(A[i][k]))
            {
                max_st=i;
                max=Math.abs(A[i][k]);
            }
        //cout <<max_st<<endl;
        return max_st;
    }

    private static void main_stol(double A[][], int n, double x[])
    {
        for (int k=0;k<n;k++)
        {
            int max_st=find_max_stol(A,n,k);
            for (int j=0;j<=n;j++)
                swap(A,k,j,max_st,j);

            if (A[k][k]!=0)
            {
                double mult = A[k][k];
                for (int j=k;j<=n;j++)
                    A[k][j]=A[k][j]/mult;
            }

            for (int i=k+1;i<n;i++)
            {
                double mult = A[i][k];
                for (int j=k;j<=n;j++)
                    A[i][j]=A[i][j]-A[k][j]*mult;
            }

		/*for (int i=0;i<n;i++)
		{
		for (int j=0;j<=n;j++)
		cout<< A[i][j]<<' ';
		cout <<endl<<endl;
		}*/
        }

        for (int i=n-1;i>=0;i--)
        {
            double tmp=0;
            for (int j=i+1;j<n;j++)
                tmp+=A[i][j]*x[j];
            x[i]=A[i][n]-tmp;
        }
    }


    private static double L(int i, int j, double x[],double p[], double q[],double r[])
    {
        return p[i]*w(x[i],j,2)+q[i]*w(x[i],j,1)+r[i]*w(x[i],j,0);
    }


    //наше уравнение, которое нам дано
    private static double LL(int i, double x,int alpha)
    {
        return 1*w(x,i,2)-(x*x+alpha)*w(x,i,1)-2*x*w(x,i,0);
    }

    //вычисление интегралов для метода Галеркина
    private static double MyFunction(double x, int i, int j,int alpha, int who)
    {
        if (who==1)
            return w(x,j,0)*LL(i,x,alpha);
        else
            return (2*(3*x*x-alpha)/(x*x+alpha)/(x*x+alpha)/(x*x+alpha)+2*x/(alpha+1))*w(x,j,0);
    }

    //метод Симпсона для решения интеграла
    private static double SimpsonIntegration(int n, int i,int j, int alpha,int who) //для поиска коэффициентов матриц
    {
        int a=-1;
        int b=1;
        n=100;
        double h = ((double)b - (double)a) / (double)n/2;
        double IntegSum = MyFunction(a,i,j,alpha,who)+MyFunction(b,i,j,alpha,who);
        for (int k = 2; k<=2*n-2; k+=2 )
            IntegSum+=2*MyFunction(a+k*h,i,j,alpha,who);
        for (int k = 1; k<=2*n-1; k+=2 )
            IntegSum+=4*MyFunction(a+k*h,i,j,alpha,who);
        IntegSum*=h/3;
        return IntegSum;
    }


    private static void Galyorkin_method(int n, int alpha,double y[]) {
        double A[][]=new double[n][];
        for (int i=0;i<n;i++)
            A[i]=new double[n+1];

        for (int j=0;j<n;j++)
            for (int i=0;i<n;i++)
                A[j][i]=SimpsonIntegration(n,i,j,alpha,1);

        for (int i=0;i<n;i++)
            A[i][n]=SimpsonIntegration(n,i,i,alpha,2);

        double c[] = new double[n];
        main_stol(A,n,c);

        double result[]=new double[3];
        System.out.println("Метод Галеркина: ");
        for (int i=0;i<3;i++)
        {
            double sum=0;
            for (int j=0;j<n;j++)
                sum+=c[j]*w(-0.5+0.5*i,j,0);
            System.out.println("y[" + (-0.5+0.5*i) + "] = " + sum);
            result[i]=y[i]-sum;
        }

        System.out.println("Погрешность:");
        for (int i=0;i<3;i++)
            System.out.println(result[i] + " ");
    }

    public static void go() {
        System.out.println("Методы решения ОДУ в аналитическом виде");
        System.out.println("Введите n");
        int n;
        Scanner sc = new Scanner(System.in);
        n = sc.nextInt();
        System.out.println("[a, b] = ["+(-1)+", "+1+"]");

        int a=-1;
        int b=1;
        int alpha=4;

        double  x[] = new double[n];
        double pi=3.141592653589793;
        // густая сетка узлов (корни многочлена Чебышева 1 порядка)
        for (int i=0;i<n;i++)
            x[i]=Math.cos((2*i+1)*pi/2.0/(double)(n));

        double p[]=new double[n];
        for (int i=0;i<n;i++) p[i]=1;

        double f[]=new double[n];
        for (int i=0;i<n;i++)
            f[i]=2*(3*x[i]*x[i]-alpha)/(x[i]*x[i]+alpha)/(x[i]*x[i]+alpha)/(x[i]*x[i]+alpha)+2*x[i]/(alpha+1);

        double r[]=new double[n];
        for (int i=0;i<n;i++)
            r[i]=-2*x[i];

        double q[]=new double[n];
        for (int i=0;i<n;i++)
            q[i]=-(x[i]*x[i]+alpha);

        double A[][]=new double[n][];
        for (int i=0;i<n;i++)
            A[i]=new double[n+1];

        //заполняем матрицу А и f, и решаем методом главного элемента по столбцу (таким образом решение - с(i)
	    double c[]=new double[n];
        for (int i=0;i<n;i++)
            for (int j=0;j<n;j++)
                A[i][j]=L(i,j,x,p,q,r);

        for (int i=0;i<n;i++)
            A[i][n]=f[i];

        main_stol(A,n,c);


	    double y[]=new double[3];

        System.out.println("Точное решение:");
        for (int i=0;i<3;i++)
        {
            y[i]=F(-0.5+0.5*i,alpha);
            System.out.println("y[" + (-0.5+0.5*i) + "] = "+y[i]);
        }


        System.out.println("Метод коллокаций: ");
	    double result[]=new double[3];
        for (int i=0;i<3;i++)
        {
            double sum=0;
            for (int j=0;j<n;j++)
                sum+=c[j]*w(-0.5+0.5*i,j,0);
            System.out.println("y[" + (-0.5+0.5*i) + "] = "+sum);
            result[i]=y[i]-sum;
        }

        System.out.println("Погрешность:");
        for (int i=0;i<3;i++)
            System.out.println(result[i]+" ");

        //метод Галеркина
        Galyorkin_method(n,alpha,y);
    }

}