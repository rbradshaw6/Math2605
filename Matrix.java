import java.util.Scanner;

public class Matrix {
    protected double[][] arrayVersion;
    protected double[][] Matrix1;
    protected double[][] Matrix2;
    protected double[][] resultantMatrix;
    protected Matrix transposeMatrix;
    protected double[] Lvalues;
    protected Matrix L;
    protected Matrix U;
    protected double[][] identity;

    private int matrixRows;
    private int matrixColumns;

    public Matrix(double[][] matrix) { //assumes you pass in a consistent matrix
        this.arrayVersion = matrix; //for passing into different methods
        this.matrixRows = matrix.length;
        this.matrixColumns = matrix[0].length;
        this.Matrix1 = new double[matrixRows][matrixColumns];
        for (int i = 0; i < matrixRows; i++)
            for (int j = 0; j < matrixColumns; j++)
                    this.Matrix1[i][j] = matrix[i][j];

        this.resultantMatrix = new double[matrixRows][matrixColumns];
        for (int z = 0; z < matrixRows; z++) {
            for(int zz = 0; zz < matrixColumns; zz++) {
                resultantMatrix[z][zz] = 0.0;

            }
        }
    }


    public Matrix(int rows, int columns) { //assumes you pass in a consistent matrix
        Scanner scan = new Scanner(System.in);

        double[][] Matrix1 = new double[rows][columns];
        this.matrixRows = rows;
        this.matrixColumns = columns;

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                System.out.println((i + 1) + " row" + " | " + (j + 1) + " column value?");
                Matrix1[i][j] = scan.nextDouble();
            }
        }

        this.arrayVersion = Matrix1;
        this.Matrix1 = Matrix1;
        this.display();

        this.resultantMatrix = new double[matrixRows][matrixColumns];
        for (int z = 0; z < matrixRows; z++) {
            for(int zz = 0; zz < matrixColumns; zz++) {
                resultantMatrix[z][zz] = 0.0;
            }
        }

    }

    public Matrix createIdentity() {
        double[][] identityz = new double[this.matrixColumns][this.matrixColumns];
        for (int i = 0; i < matrixColumns; i++) {
            for (int j = 0; j < matrixColumns; j++) {
                if (i == j) {
                    identityz[i][j] = 1.0;
                } else {
                    identityz[i][j] = 0.0;
                }
            }
        }
        double[][] a = identityz;

        Matrix g = new Matrix(a);
        this.identity = a;
        return g;

    }

     public Matrix transpose() {
            double[][] transposeArray = new double[this.matrixColumns][this.matrixRows];

            for (int i = 0; i < matrixRows; i++) {
                for (int j = 0; j < matrixColumns; j++) {
                    transposeArray[j][i] = this.Matrix1[i][j];
                }
            }


            this.transposeMatrix = new Matrix(transposeArray);

            System.out.println(" Transpose Matrix:\n");

            transposeMatrix.display();

            return transposeMatrix;
       }

       public Vector matrixVectorMult(Vector v) {
           double[][] mPlaceholder = new double[v.length][1];
           for (int i = i; i < v.length; i++) {
               mPlaceholder[i][1] = v.array[i];
           }
           Matrix result = matrixMultiplier(new Matrix(mPlaceholder));
           double[] newVector = new double[v.length];
           for (int i = 0; i < v.length; i++) {
               newVector[i] = result.arrayVersion[i][1];
           }
           return new Vector(newVector)
       }

     //http://www.programiz.com/python-programming/examples/multiply-matrix used to troubleshoot triple for loop
    public Matrix matrixMultiplier(Matrix secondMatrix) {
        if (this.matrixColumns != secondMatrix.matrixRows) {

            System.out.println("Number of columns of 1st matrix do not match the rows of 2nd.");
            return null;
        } else {
            this.resultantMatrix = new double[this.matrixRows][secondMatrix.matrixColumns];
            for (int z = 0; z < matrixRows; z++) {
                for(int zz = 0; zz < secondMatrix.matrixColumns; zz++) {
                    resultantMatrix[z][zz] = 0.0;
                }
            }

            for (int i = 0; i < this.matrixRows; i++) {
                for (int j = 0; j < secondMatrix.matrixColumns; j++) {
                    for(int z = 0; z < secondMatrix.matrixRows; z++) {
                        resultantMatrix[i][j] += (this.Matrix1[i][z] * secondMatrix.Matrix1[z][j]);
                    }
                }
            }


            System.out.println("\n Multiplied Matrix");
            Matrix multipliedMatrix = new Matrix(resultantMatrix);

            multipliedMatrix.display();
            return multipliedMatrix;

        }
    }

    public void QR() {
        Matrix identity = this.createIdentity();

    }

    public Matrix multiplyByConstant(double constant) {
        Matrix copyMatrix = new Matrix(this.arrayVersion);
        double[][] temp = copyMatrix.arrayVersion;

        for (int i = 0; i < matrixRows; i++) {
            for (int j = 0; j < matrixColumns; j++) {
                temp[i][j] *= constant;
            }
        }

        Matrix newMatrix = new Matrix(temp);
        newMatrix.display();

        return newMatrix;
    }

    private void display() {
        for (int i1 = 0; i1 < this.matrixRows; i1++) { //display
            for (int j1 = 0; j1 < this.matrixColumns; j1++) {
                System.out.print(this.Matrix1[i1][j1] + "\t");
            }
            System.out.print("\n");
        }
    }

/*  Experimental code section */
public Matrix REF() {
        double[][] matrix = this.arrayVersion;


        int test = 0;
        int recursL = test;
        int minsize = matrixRows < matrixColumns ? matrixRows : matrixColumns;

        for (int i = 0; i < minsize-1; i++) {
            int max = i;
            for (int j = i + 1; j < matrixRows-1; j++) {
                if (Math.abs(matrix[j][i]) > Math.abs(matrix[max][i])) {
                    max = j;
                }
            }
            if (i != max) {
                double[] pl = matrix[i];
                matrix[i] = matrix[max];
                matrix[max] = pl;
            }

            if (Math.abs(matrix[i][i]) <= 1/(Math.pow(10, 10))) {
                return null;
            }

            for (int k = i + 1; k < matrixRows; k++) {
                double temp = (matrix[k][i])/(matrix[i][i]);
                //Lvalues[recursL] = temp;
                //System.out.println(Lvalues[recursL]);
                //recursL++;
                //System.out.println(temp);
                for (int z = i; z < matrixColumns; z++) {
                    matrix[k][z] -= matrix[i][z] * temp;
                }
            }
        }

        for (int i1 = matrixRows - 1; i1 > -1; i1--) {
            double temp2 = matrix[i1][i1];
//          System.out.println(temp2);
            for(int j = 0; j < i1; j++) {

            }
            //matrix[i1][i1] /= temp2;


            for (int k = matrixRows; k < matrixColumns; k++) {
                matrix[i1][k] /= temp2;
            }
        }



        Matrix newMatrix = new Matrix(matrix);
        return newMatrix;
    }

    public Matrix inverse() {
        if (matrixRows != matrixColumns) { throw new IllegalArgumentException("Matrix must be square!"); }
        double[][] newArr = new double[matrixRows][matrixColumns * 2];
        for (int i = 0; i < matrixRows; i++) {
            for (int j = 0; j < matrixColumns; j++) {
                newArr[i][j] = Matrix1[i][j];
            }
        }
        for (int i = 0; i< newArr.length; i++) {
            for (int j = matrixColumns; j < newArr[0].length; j++) {
                if (i == j - matrixColumns) {
                    newArr[i][j] = 1;
                } else {
                    newArr[i][j] = 0;
                }
            }
        }

        //modified REF
        double[][] matrix = newArr;


        int minsize = matrix.length < matrix[0].length ? matrix.length : matrix[0].length;

        for (int i = 0; i < minsize-1; i++) {
            int max = i;
            for (int j = i + 1; j < matrix.length - 1; j++) {
                if (Math.abs(matrix[j][i]) > Math.abs(matrix[max][i])) {
                    max = j;
                }
            }
            if (i != max) {
                double[] pl = matrix[i];
                matrix[i] = matrix[max];
                matrix[max] = pl;
            }

            if (Math.abs(matrix[i][i]) <= 1/(Math.pow(10, 10))) {
                return null;
            }

            for (int k = 0; k < matrix.length; k++) {
                if (k != i) {
                    double temp = (matrix[k][i])/(matrix[i][i]);
                    for (int z = i; z < matrix[0].length; z++) {
                        matrix[k][z] -= matrix[i][z] * temp;
                    }
                } else {
                    double temp = ((matrix[i][i] - 1) / matrix[i][i]);
                    for (int z = i; z < matrix[0].length; z++) {
                        matrix[k][z] -= matrix[i][z] * temp;
                    }
                }
            }
        }

        for (int i1 = matrix.length - 1; i1 > -1; i1--) {
            double temp2 = matrix[i1][i1];
//          System.out.println(temp2);
            for(int j = 0; j < i1; j++) {

            }
            //matrix[i1][i1] /= temp2;


            for (int k = matrix.length; k < matrix[0].length; k++) {
                matrix[i1][k] /= temp2;
            }
        }
        double[][] n = new double[matrixRows][matrixColumns];
        for (int i = 0; i < matrixRows; i++) {
            for (int j = 0; j < matrixColumns; j++) {
                n[i][j] = matrix[i][j + matrixColumns];
            }
        }
        Matrix neo = new Matrix(n);
        return neo;
    }

    public void lu_decomposition() {
        double[][] matrix = this.arrayVersion;
        double[][] Larr = new double[matrixRows][matrixRows];

        for (int i = 0; i < Larr[0].length; i++) {
            for (int j = 0; j < Larr.length; j++) {
                if (i == j) {
                    Larr[i][j] = 1.0;
                } else {
                    Larr[i][j] = 0.0;
                }
            }
        }
        Matrix l = new Matrix(Larr);

        this.Lvalues = new double[100];

        int test = 0;
        int recursL = test;

        for (int i = 0; i < matrixColumns-1; i++) {

            double[][] Ltemp = new double[matrixRows][matrixRows];

            for (int i2 = 0; i2 < Ltemp[0].length; i2++) {
                for (int j2 = 0; j2 < Ltemp.length; j2++) {
                    if (i2 == j2) {
                        Ltemp[i2][j2] = 1.0;
                    } else {
                        Ltemp[i2][j2] = 0.0;
                    }
                }
            }

            int max = i;
            for (int j = i + 1; j < matrixRows-1; j++) {
                if (Math.abs(matrix[j][i]) > Math.abs(matrix[max][i])) {
                    max = j;
                }
            }
            if (false /*i != max*/) {
                double[] pl = matrix[i];
                matrix[i] = matrix[max];
                matrix[max] = pl;

                // switch columns in L matrix
                for (int j = 0; j < Ltemp.length; j++) {
                    double pl2 = Ltemp[j][i];
                    Ltemp[j][i] = Ltemp[j][max];
                    Ltemp[j][max] = pl2;
                }
                i--;
            } else {

                if (Math.abs(matrix[i][i]) <= 1/(Math.pow(10, 10))) {
                    throw new IllegalArgumentException("LU not possible");
                }

                for (int k = i + 1; k < matrixRows; k++) {
                    double temp = (matrix[k][i])/(matrix[i][i]);
                    //Lvalues[recursL] = temp;
                    //System.out.println(Lvalues[recursL]);
                    //recursL++;
                    Ltemp[k][i] = -1 * temp;
                    //System.out.println(temp);
                    for (int z = i; z < matrixColumns; z++) {
                        matrix[k][z] -= matrix[i][z] * temp;
                    }
                }
            }
            Matrix Lt = new Matrix(Ltemp);
            l = Lt.matrixMultiplier(l);
        }

//        for (int i1 = matrixRows - 1; i1 > -1; i1--) {
//            double temp2 = matrix[i1][i1];
////          System.out.println(temp2);
//            for(int j = 0; j < i1; j++) {
//
//            }
//            //matrix[i1][i1] /= temp2;
//
//
//            for (int k = matrixRows; k < matrixColumns; k++) {
//                matrix[i1][k] /= temp2;
//            }
//        }
        Matrix u = new Matrix(matrix);

        System.out.println("U: ");
        u.display();
        System.out.println("L: ");
        l.inverse().display();

    }

    public double determinant() {
        if (matrixRows != matrixColumns) {
            throw new IllegalArgumentException("Matrix must be square!");
        }
        if (matrixRows == 1) {
            return arrayVersion[0][0];
        } else if (matrixRows == 2) {
            return ((arrayVersion[0][0] * arrayVersion[1][1]) - (arrayVersion[1][0] * arrayVersion[0][1]));
        } else {
            double sum = 0.0;
            int i = 0;
            for (int j=0; j < matrixRows; j++) {
                double[][] temp = new double[matrixRows-1][matrixColumns-1];
                for (int i2=0, i3=0; i2 < matrixRows; i2++) {
                    for (int j2=0, j3=0; j2 < matrixColumns; j2++) {
                        if (i2 != i && j2 != j) {
                            temp[i3][j3] = arrayVersion[i2][j2];
                            j3++;
                        }
                    }
                    if (i2 != i) {
                        i3++;
                    }
                }
                Matrix tempm = new Matrix(temp);
                sum += Math.pow(-1, i + j) * tempm.determinant() * arrayVersion[i][j];
            }
            return sum;
        }
    }

    public static void main(String[] args) {
        double[][] m1 = {{1, 2, 4}, {3, 8, 14}, {2, 6, 14}};
        //double[][] m2 = {{1, -2, -2, -3}, {3, -9, 0, -9}, {-1, 2, 4, 7}, {-3, -6, 26, 2}};
        Matrix m = new Matrix(m1);
        //m.lu_decomposition();
//        double[][] m3 = {{1, 2}, {3, 4}, {3, 2}};
//        double[][] m4 = {{4, -1}, {5, 0}};
//        Matrix matrix3 = new Matrix(m3);
//        Matrix matrix4 = new Matrix(m4);
//        matrix3.matrixMultiplier(matrix4);
        // double[][] m5 = {{1, 0, 0}, {2, 1, 0}, {3, 2, 1}};
        // Matrix matrix5 = new Matrix(m5);
        // Matrix matrix6 = matrix5.inverse();
        // System.out.println("***");
        // matrix6.display();
        m.display();
        System.out.println(m.determinant());
    }

}
