package nf.functions;

public class HFunction implements MyMath.Function2 {
    public double value(double x, double y){
        return -1.0 / 3.0 * Math.exp( -x * y);
    }
}
