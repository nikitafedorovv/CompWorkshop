package nf;

import Jama.Matrix;
import java.util.Scanner;

public class Task3 {

    private static Matrix normirY(Matrix y1){
        int indexOfGreaterCoord = 0;
        double max = 0;
        for(int i = 0; i < y1.getRowDimension(); i++){
            if(Math.abs(y1.get(i, 0)) > max){
                max = Math.abs(y1.get(i, 0));
                indexOfGreaterCoord = i;
            }
        }
        y1 = y1.times(1 / y1.get(indexOfGreaterCoord,0));
        return y1;
    }

    private static double lambdaMaxMethod1(Matrix A, int accuracy){
        Matrix y0, y1 = new Matrix(A.getColumnDimension(),1,1);
        double eps = 1/Math.pow(10, accuracy);
        double lambdaMax0, lambdaMax1 = 0;
        long counter = 0;
        do {
            lambdaMax0 = lambdaMax1;
            y0 = y1;
            y1 = A.times(y1);
            lambdaMax1 = y1.get(0,0) / y0.get(0,0);
            y1 = normirY(y1);
            counter++;
        } while (Math.abs(lambdaMax1 - lambdaMax0) >= eps);

        System.out.println("\nMax lambda by method 1 (power) counted in " + counter + " steps.");

        return lambdaMax1;
    }

    private static double lambdaMaxMethod2(Matrix A, int accuracy){
        Matrix y0, y1 = Matrix.random(A.getColumnDimension(), 1);
        long counter = 0;
        double eps = 1/Math.pow(10, accuracy);
        double lambdaMax0, lambdaMax1 = 0;
        do {
            lambdaMax0 = lambdaMax1;
            y0 = y1;
            y1 = A.times(y1);
            lambdaMax1 = y1.transpose().times(y0).get(0,0) / y0.transpose().times(y0).get(0,0);
            y1 = normirY(y1);
            counter++;
        } while (Math.abs(lambdaMax0 - lambdaMax1) >= eps);

        System.out.println("\nMax lambda by method 2 (scalar) counted in " + counter + " steps.");

        return lambdaMax1;
    }

    private static double eigenvaluesSumMethod1(Matrix A, int accuracy){
        double lambdaMax = lambdaMaxMethod1(A, accuracy);

        Matrix ETimesLambdaMax = new Matrix(A.getRowDimension(), A.getColumnDimension());
        for(int i = 0; i < A.getRowDimension(); i++){
            ETimesLambdaMax.set(i, i, lambdaMax);
        }
        Matrix B = A.minus(ETimesLambdaMax);

        return lambdaMaxMethod1(B, accuracy) + 2 * lambdaMax;
    }

    private static double eigenvaluesSumMethod2(Matrix A, int accuracy){
        double lambdaMax = lambdaMaxMethod2(A, accuracy);

        Matrix ETimesLambdaMax = new Matrix(A.getRowDimension(), A.getColumnDimension());
        for(int i = 0; i < A.getRowDimension(); i++){
            ETimesLambdaMax.set(i, i, lambdaMax);
        }
        Matrix B = A.minus(ETimesLambdaMax);

        return lambdaMaxMethod2(B, accuracy) + 2 * lambdaMax;
    }

    private static Matrix[] diagPredomHB(Matrix A, Matrix b) {
        Matrix D = new Matrix(A.getRowDimension(), A.getColumnDimension());
        Matrix E = new Matrix(A.getRowDimension(), A.getColumnDimension());
        for(int i = 0; i < A.getRowDimension(); i++){
            D.set(i, i, A.get(i,i));
            E.set(i, i, 1);
        }

        Matrix H = E.minus(D.inverse().times(A));
        Matrix g = D.inverse().times(b);

        Matrix[] answer = new Matrix[2];
        answer[0] = H;
        answer[1] = g;
        return answer;
    }

    private static Matrix[] hermitAndPosHB(Matrix A, Matrix b, double eigenvaluesSum) {

        Matrix E = new Matrix(A.getRowDimension(), A.getColumnDimension());
        for(int i = 0; i < A.getRowDimension(); i++){
            E.set(i, i, 1);
        }

        Matrix H = E.minus(A.times(2 / eigenvaluesSum));
        Matrix g = b.times(2 / eigenvaluesSum);

        Matrix[] answer = new Matrix[2];
        answer[0] = H;
        answer[1] = g;
        return answer;
    }

    private static void iterate(Matrix[] HAndG, Matrix A, Matrix b, Matrix x0, int accuracy){
        Matrix H = HAndG[0];
        Matrix g = HAndG[1];

        double normH = MatrixAlgebra.Norm.normInf(H);
        double coeff = 1;
        if(normH >= 0.5) {
            coeff = normH / (1 - normH);
        }

        long counter = 0;
        Matrix x1, x2 = new Matrix(x0.getArrayCopy());
        double eps = 1/Math.pow(10, accuracy);
        do {
            x1 = x2;
            x2 = H.times(x1).plus(g);
            counter++;
        } while(MatrixAlgebra.Norm.normInf(x2.minus(x1)) >= eps);

        System.out.print("is in error by:");
        (b.minus(A.times(x2))).print(accuracy + 5, accuracy + 5);
        System.out.println("Done in " + counter + " steps");
    }

    private static void iterateSeidel(Matrix A, Matrix b, Matrix x0, int accuracy) {
        Matrix[] HAndG = diagPredomHB(A, b);
        Matrix H = HAndG[0], g = HAndG[1];
        double eps = 1/Math.pow(10, accuracy);

        long counter = 0;
        Matrix x1, x2 = new Matrix(x0.getArrayCopy());
        do{
            x1 = new Matrix(x2.getArrayCopy());
            for(int i = 0; i < x2.getRowDimension(); i++){
                double sum = 0;
                for(int j = 0; j < x2.getRowDimension(); j++){
                    sum += H.get(i,j) * x2.get(j,0);
                }
                //С предыдущими/следующими координатами тут как раз очень удобно. Не нужно на две суммы разбивать
                x2.set(i,0, sum + g.get(i,0));
            }
            counter++;
        } while(MatrixAlgebra.Norm.normInf(x1.minus(x2)) >= eps);

        System.out.print("is in error by:");
        (b.minus(A.times(x2))).print(accuracy + 5, accuracy + 5);
        System.out.println("Done in " + counter + " steps");
    }

    public static void go(){

        double[][] systemAsArray = {{2.20219, 0.33266, 0.16768, 2.17523},
                                    {0.33266, 3.17137, 0.54055, 6.58335},
                                    {0.16768, 0.54055, 4.92343, 6.36904}};
        Matrix x0 = Matrix.random(systemAsArray.length, 1);

        System.out.print("Enter accuracy (integer): ");
        int accuracy = (new Scanner(System.in)).nextInt();

        Matrix system = new Matrix(systemAsArray);
        Matrix A = system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2);
        Matrix b = system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1);

        System.out.print("\n\n\nDiagPredom matrix type solution ");
        iterate(diagPredomHB(A, b), A, b, x0, accuracy);

        System.out.println("\n");
        eigenvaluesSumMethod1(A, accuracy);

//        System.out.print("\n\n\nHermit matrix type solution ");
//        iterate(hermitAndPosHB(A, b, eigenvaluesSumMethod1(A, accuracy)), A, b, x0, accuracy);
        System.out.print("\n\n\nHermit matrix type solution ");
        iterate(hermitAndPosHB(A, b, eigenvaluesSumMethod2(A, accuracy)), A, b, x0, accuracy);
        System.out.print("\n\n\nSeidel solution ");
        iterateSeidel(A, b, x0, accuracy);
    }
}
