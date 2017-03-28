package MyMath;

public class Integral {

    public static double value(Function f, double a, double b, double eps){ //методом трапеций
        if(a > b){
            double c = a;
            a = b;
            b = c;
        }
        double newValue = 0;
        double oldValue = eps;
        long n = 1;
        while(Math.abs(newValue - oldValue) >= eps){
            oldValue = newValue;
            double step = (b - a)/n;
            double point = a;
            double sum = 0;
            while(point < b){
                sum += (f.value(point) + f.value(point + step)) * 0.5 * step;
                point += step;
            }
            newValue = sum;
            n *= 2;
        }
        return newValue;
    }
}
