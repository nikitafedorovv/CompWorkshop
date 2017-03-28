package nf.functions;

/**
 * Created by nikitafedorov on 29/03/2017.
 */
public class ExpFunction implements MyMath.Function {
    public double value(double x){
        return Math.exp(x);
    }
}
