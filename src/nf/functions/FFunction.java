package nf.functions;

public class FFunction implements MyMath.Function {
    public double value(double x){
        return 0.5 * (2 - x + x * x);
    }
}
