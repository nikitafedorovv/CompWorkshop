package nf.functions;

public class AlphaTimesBetaFunction implements MyMath.Function {
    private int i;
    private int j;
    public double value(double y){
        return new AlphaFunction(j).value(y) * new BetaFunction(i).value(y);
    }

    public int getI() {
        return i;
    }

    public int getJ() {
        return j;
    }

    public AlphaTimesBetaFunction(int j, int i){
        this.i = i;
        this.j = j;
    }
}
