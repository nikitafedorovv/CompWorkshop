package nf.functions;

import MyMath.Factorial;

public class AlphaFunction implements MyMath.Function{
    private int index;
    public double value(double x){
        byte sign = (byte)(2 * (index % 2) - 1);
        return sign * (1.0 / 3.0) * Math.pow(x, index) / Factorial.get(index);
    }
    public int getIndex() {
        return index;
    }

    public AlphaFunction(int index){
        this.index = index;
    }

}
