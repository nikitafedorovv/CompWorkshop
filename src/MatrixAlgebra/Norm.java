package MatrixAlgebra;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * Created by nikitafedorov on 19/03/2017.
 */
public class Norm {
    public static double norm1(Matrix matrix){
        double max = 0;
        for(int i = 0; i < matrix.getColumnDimension(); i++){
            double sum = 0;
            for(int j = 0; j < matrix.getRowDimension(); j++){
                sum += Math.abs(matrix.get(i,j));
            }
            if (sum > max) {
                max = sum;
            }
        }
        return max;
    }

    public static double norm2(Matrix matrix){
        EigenvalueDecomposition eig = matrix.times(matrix.transpose()).eig();

        double[] eigIm = eig.getImagEigenvalues();
        double[] eigRe = eig.getRealEigenvalues();

        double max = 0;
        for(int i = 0; i < eigRe.length; i++){
            double abs = Math.abs(eigRe[i] * eigRe[i] + eigIm[i] * eigIm[i]);
            if (abs > max) {
                max = abs;
            }
        }
        return Math.sqrt(max);
    }

    public static double normInf(Matrix matrix){
        double max = 0;
        for(int i = 0; i < matrix.getRowDimension(); i++){
            double sum = 0;
            for(int j = 0; j < matrix.getColumnDimension(); j++){
                sum += Math.abs(matrix.get(i,j));
            }
            if (sum > max){
                max = sum;
            }
        }
        return max;
    }
}
