package nf.functions;

public class BetaTimesFFunction implements MyMath.Function {
    private int index;
    public double value(double y){
        return new BetaFunction(index).value(y) * new FFunction().value(y);
    }

    public int getIndex() {
        return index;
    }

    public BetaTimesFFunction(int index){
        this.index = index;
    }
}
