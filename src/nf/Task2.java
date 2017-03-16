package nf;

import Jama.Matrix;

public class Task2 {

    public static Matrix gaussSolve(Matrix system){
        system = new Matrix(system.getArrayCopy());
        for(int k = 0; k < system.getRowDimension(); k++){
            for(int i = k; i < system.getRowDimension(); i++){
                if (system.get(i, k) == 0) { //если ноль, то складываем с чем-нибудь из низа, не нулевым в этом месте
                    int v = i + 1;
                    while(system.get(v,k) == 0){
                        v++;
                    }
                    for(int m = 0; m < system.getColumnDimension(); m++){
                        system.getArray()[i][m] += system.get(v,m);
                    }
                }
                double valueToDivide = system.get(i, k);
                for(int l = 0; l < system.getColumnDimension(); l++){
                    system.getArray()[i][l] /= valueToDivide;
                }
            }
            for(int j = k + 1; j < system.getRowDimension(); j++){

                if (system.get(j, k) != 0) { //если ноль, то отнимать уже ничего и не надо
                    for(int l = 0; l < system.getColumnDimension(); l++){
                        system.getArray()[j][l] -= system.get(k,l);
                    }
                }
            }
        }

        double[] answerAsArray = new double[system.getRowDimension()];

        for(int i = system.getRowDimension() - 1; i >= 0 ; i--){
            double sum = 0;
            for(int j = i + 1; j < system.getRowDimension(); j++) {
                sum += system.get(i,j) * answerAsArray[j];
            }
            answerAsArray[i] = system.get(i, system.getColumnDimension() - 1) - sum;
        }

        Matrix answer = new Matrix(1, system.getRowDimension());
        answer.getArray()[0] = answerAsArray;
        answer = answer.transpose();
        return answer;
    }

    public static Matrix[] lu(Matrix matrix){
        int n = matrix.getColumnDimension();
        double[][] u = new double[n][n];
        double[][] l = new double[n][n];

        for(int j = 0; j < n; j++){
            u[0][j] = matrix.get(0, j);
        }

        for(int j = 0; j < n; j++){
            l[j][0] = matrix.get(j,0) / u[0][0];
        }

        for(int i = 1; i < n; i++){
            for(int j = i; j < n; j++){
                double sum = 0;
                for(int k = 0; k <= i - 1; k++){
                    sum += l[i][k] * u[k][j];
                }
                u[i][j] = matrix.get(i, j) - sum;
            }
            for(int j = i; j < n; j++){
                double sum = 0;
                for(int k = 0; k <= i - 1; k++){
                    sum += l[j][k] * u[k][i];
                }
                l[j][i] = (matrix.get(j, i) - sum)/u[i][i];
            }
        }
        Matrix[] result = new Matrix[2];
        result[0] = new Matrix(l);
        result[1] = new Matrix(u);
        return result;
    }

    public static Matrix luSolve(Matrix[] lu, Matrix b){

        double[] y = new double[b.getRowDimension()];

        for(int i = 0; i < b.getRowDimension(); i++){

            double sum = 0;

            for(int j = 0; j < i; j++) {
                sum += lu[0].get(i, j) * y[j];
            }

            y[i] = b.get(i,0) - sum;
        }

        double[] x = new double[y.length];

        for(int i = y.length - 1; i >= 0 ; i--){
            double sum = 0;
            for(int j = i + 1; j < y.length; j++) {
                sum += lu[1].get(i, j) * x[j];
            }
            x[i] = (y[i] - sum)/lu[1].get(i, i);
        }

        Matrix answer = new Matrix(1, x.length);
        answer.getArray()[0] = x;
        answer = answer.transpose();
        return answer;
    }

    public static void main(String[] args) {
        double[][] systemAsArray = {{0.00008164, -0.4772, 4.7292, 5.07181},
                                    {0.7564    , -0.4772, 1.9592, 2.51058},
                                    {0.00010364,  1.1228, 1.4092, 3.02007}};
        Matrix system = new Matrix(systemAsArray);

        Matrix gaussS = gaussSolve(system);
        Matrix[] lu = lu(system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2));
        Matrix luS = luSolve(lu, system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1));

//        gaussS.print(10,10);
//        luS.print(10,10);

        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(gaussS)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(30,30);

        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(luS)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(30,30);
    }
}
