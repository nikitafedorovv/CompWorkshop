package nf;

import Jama.Matrix;

public class Task3 {

    private static Matrix[] type1(Matrix A, Matrix b) {
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

    private static Matrix[] type2(Matrix A, Matrix b) {
//        Matrix D = new Matrix(A.getRowDimension(), A.getColumnDimension());
//        Matrix E = new Matrix(A.getRowDimension(), A.getColumnDimension());
//        for(int i = 0; i < A.getRowDimension(); i++){
//            D.set(i, i, A.get(i,i));
//            E.set(i, i, 1);
//        }

//        Matrix H;
//        Matrix g;
//
        Matrix[] answer = new Matrix[2];
//        answer[0] = H;
//        answer[1] = g;
        return answer;
    }

    private static void iterate(Matrix[] HAndG, Matrix x0, int accuracy){
        Matrix H = HAndG[0];
        Matrix g = HAndG[1];

        if(MatrixAlgebra.Norm.normInf(H) >= 0.5) {
            System.out.println("WRONG WAY, see 'count' METHOD'");
        }

        long counter = 0;
        Matrix x2 = new Matrix(x0.getArrayCopy());
        Matrix x1;
        double eps = 1/Math.pow(10, accuracy);
        do {
            x1 = x2;
            x2 = H.times(x1).plus(g);
            counter++;
        } while(MatrixAlgebra.Norm.normInf(x2.minus(x1)) >= eps);

        System.out.println("\nType1 solution:");
        x2.print(accuracy + 5, accuracy + 5);
        System.out.println("Done in " + counter + " steps");
    }

    public static void go(){

        double[][] systemAsArray = {{2.20219, 0.33266, 0.16768, 2.17523},
                                    {0.33266, 3.17137, 0.54055, 6.58335},
                                    {0.16768, 0.54055, 4.92343, 6.36904}};
        Matrix x0 = Matrix.random(systemAsArray.length, 1);
        int accuracy = 10;

        Matrix system = new Matrix(systemAsArray);
        Matrix A = system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2);
        Matrix b = system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1);

        System.out.println("\nSolution by LU-decomposition:");
        A.solve(b).print(accuracy + 5, accuracy + 5);

        iterate(type1(A, b), x0, accuracy);
//        iterate(type2(A, b), x0, accuracy);


    }
}
