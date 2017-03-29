package nf.functions;

public class AlphaTimesBetaFunction implements MyMath.Function {
    private int i;
    private int j;
    public double value(double y){
        AlphaFunction alpha = new AlphaFunction(i);
        BetaFunction beta = new BetaFunction(j);
        return alpha.value(y) * beta.value(y);
    }

    public int getI() {
        return i;
    }

    public int getJ() {
        return j;
    }

    public AlphaTimesBetaFunction(int i, int j){
        this.i = i;
        this.j = j;
    }
}
