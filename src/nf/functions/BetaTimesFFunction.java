package nf.functions;

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
