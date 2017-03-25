package nf;

import Jama.Matrix;

public class Task2 {

    public static Matrix gaussSolve(Matrix system){
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

    private static int maxElementInRow(Matrix m){
        double max = 0;
        int indexMax = 0;
        for(int i = 0; i < m.getColumnDimension() - 1; i++){ //ищем максимальный элемент в первой строке
            double current = Math.abs(m.get(0,i));
            if (current > max){
                indexMax = i;
                max = current;
            }
        }
        return indexMax;
    }

    private static int maxElementInColumn(Matrix m){
        double max = 0;
        int indexMax = 0;
        for(int i = 0; i < m.getRowDimension(); i++){ //ищем максимальный элемент в первом столбце
            double current = Math.abs(m.get(i,0));
            if (current > max){
                indexMax = i;
                max = current;
            }
        }
        return indexMax;
    }

    private static int[] maxElement(Matrix m){
        int[] indexMax = new int[2];
        double max = 0;
        for(int i = 0; i < m.getRowDimension(); i++){
            for(int j = 0; j < m.getColumnDimension() - 1; j++){
                if (Math.abs(m.get(i,j)) > max){
                    max = Math.abs(m.get(i,j));
                    indexMax[0] = i;
                    indexMax[1] = j;
                }
            }
        }
        return indexMax;
    }

    private static Matrix swapCols(Matrix m, int indexMax){
        if (indexMax > 0){ //меняем местами первый и с максимальным элементом столбцы
            double temp;
            for(int i = 0; i < m.getColumnDimension() - 1; i++){
                temp = m.get(i, 0);
                m.set(i,0, m.get(i, indexMax));
                m.set(i, indexMax, temp);
            }
        }
        return m;
    }

    private static Matrix swapRows(Matrix m, int indexMax){
        if (indexMax > 0){ //меняем местами первую и с максимальным элементом строки
            double temp;
            for(int i = 0; i < m.getColumnDimension(); i++){
                temp = m.get(0, i);
                m.set(0,i, m.get(indexMax, i));
                m.set(indexMax, i, temp);
            }
        }
        return m;
    }

    private static Matrix divideFirstRowByDouble(Matrix m, double max){
        if (max != 1) {
            for (int i = 0; i < m.getColumnDimension(); i++) {  //делим первую строку на максимальный
                m.set(0, i, m.get(0, i) / max);
            }
        }
        return m;
    }

    private static Matrix zerosInFirstColumn(Matrix m){
        for(int i = 1; i < m.getRowDimension(); i++){ //делаем нули в первом столбце под единицей
            double firstElement = m.get(i, 0);
            for(int j = 0; j < m.getColumnDimension(); j++){
                m.set(i, j,m.get(i, j) - m.get(0, j) * firstElement);
            }
        }
        return m;
    }

    private static Matrix countAnswer(Matrix m, int indexMax){
        double[] answerAsArray = new double[m.getRowDimension()];

        if(m.getRowDimension() > 1){
            Matrix x = gaussModifMaxInRow(m.getMatrix(1, m.getRowDimension() - 1, 1, m.getColumnDimension() - 1));

            double sum = 0;
            for(int i = 1; i < m.getColumnDimension() - 1; i++){
                sum += m.get(0, i) * x.get(i - 1, 0);
            }
            double x0 = m.get(0, m.getColumnDimension() - 1) - sum; //посчитали новую координату решения

            if (indexMax != 0){
                answerAsArray[0] = x.get(indexMax - 1,0);  //правильно заполняем ответ, с учетом того, что двигали столбцы
                for(int i = 1; i < indexMax; i++){
                    answerAsArray[i] = x.get(i - 1, 0);
                }
                answerAsArray[indexMax] = x0;
                for(int i = indexMax + 1; i < m.getColumnDimension() - 1; i++){
                    answerAsArray[i] = x.get(i - 1, 0);
                }
            } else { //не двигали столбцы
                for(int i = 1; i < m.getColumnDimension() - 1; i++){
                    answerAsArray[i] = x.get(i - 1, 0);
                }
                answerAsArray[0] = x0;
            }

        } else {
            answerAsArray[0] = m.get(0,1);
        }

        Matrix answer = new Matrix(1, m.getColumnDimension() - 1);
        answer.getArray()[0] = answerAsArray;
        answer = answer.transpose();
        return answer;
    }

    private static Matrix countAnswer(Matrix m){
        return countAnswer(m, 0);
    }

    public static Matrix gaussModifMaxInRow(Matrix system){
        int indexMax = maxElementInRow(system);
        double max = system.get(0,indexMax);
        swapCols(system,indexMax);
        divideFirstRowByDouble(system, max);
        zerosInFirstColumn(system);
        return countAnswer(system, indexMax);
    }

    public static Matrix gaussModifMaxInColumn(Matrix system){
        int indexMax = maxElementInColumn(system);
        double max = system.get(indexMax,0);
        swapRows(system,indexMax);
        divideFirstRowByDouble(system, max);
        zerosInFirstColumn(system);
        return countAnswer(system);
    }

    public static Matrix gaussModifMax(Matrix system){

        int[] indexMax = maxElement(system);
        double max = system.get(indexMax[0],indexMax[1]);
        swapRows(system, indexMax[0]);
        swapCols(system, indexMax[1]);
        divideFirstRowByDouble(system, max);
        zerosInFirstColumn(system);
        return countAnswer(system, indexMax[1]);
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

    public static void go() {
        double[][] systemAsArray = {{0.00008164, -0.4772, 4.7292, 5.07181},
                                    {0.7564    , -0.4772, 1.9592, 2.51058},
                                    {0.00010364,  1.1228, 1.4092, 3.02007}};

        double[][] system2AsArray = {{3, -4, 2, 5},
                                     {2,  3, 5, 5},
                                     {6,  3, 1, 4}};

        double[][] eeAsArray = new double[systemAsArray.length][systemAsArray.length];

        for(int i = 0; i < eeAsArray.length; i++){
            eeAsArray[i][eeAsArray.length - i - 1] = 1;
        }

        Matrix ee = new Matrix(eeAsArray);

        Matrix system = new Matrix(systemAsArray);

        int symb = 31;

        Matrix gaussS = ee.times(gaussSolve(new Matrix(system.getArrayCopy())));
        Matrix gaussMaxInRow = ee.times(gaussModifMaxInRow(new Matrix(system.getArrayCopy())));
        Matrix gaussMaxInColumn = ee.times(gaussModifMaxInColumn(new Matrix(system.getArrayCopy())));
        Matrix gaussMax = ee.times(gaussModifMax(new Matrix(system.getArray())));
        Matrix[] lu = lu(system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2));
        Matrix luS = luSolve(lu, system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1));

//        gaussS.print(10,10);
//        luS.print(10,10);

        System.out.print("Gauss:");
        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(gaussS)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(symb,symb);
//        gaussS.print(symb,symb);

        System.out.print("Gauss modification - maximum in row:");
        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(gaussMaxInRow)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(symb,symb);
//        gaussMaxInRow.print(symb,symb);

        System.out.print("Gauss modification - maximum in column:");
        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(gaussMaxInColumn)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(symb,symb);
//        gaussMaxInColumn.print(symb,symb);

        System.out.print("Gauss modification - maximum in matrix:");
        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(gaussMax)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(symb,symb);
//        gaussMax.print(symb,symb);

        System.out.print("LU:");
        system.getMatrix(0, system.getRowDimension() - 1, 0, system.getColumnDimension() - 2)
                .times(luS)
                .minus(system.getMatrix(0, system.getRowDimension() - 1, system.getColumnDimension() - 1, system.getColumnDimension() - 1))
                .print(symb,symb);
//        luS.print(symb,symb);
    }
}