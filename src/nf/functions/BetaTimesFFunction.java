package nf.functions;

/**
 * Created by nikitafedorov on 29/03/2017.
 */
public class BetaTimesFFunction implements MyMath.Function {
    private int index;
    public double value(double y){
        BetaFunction beta = new BetaFunction(index);
        FFunction f = new FFunction();
        return beta.value(y) * f.value(y);
    }

    public int getIndex() {
        return index;
    }

    public BetaTimesFFunction(int index){
        this.index = index;
    }
}
