package nf.functions;

public class BetaFunction implements MyMath.Function{
    private int index;
    public double value(double y){
        return Math.pow(y, index);
    }
    public int getIndex() {
        return index;
    }
    public BetaFunction(int index){
        this.index = index;
    }

}
