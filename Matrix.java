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

     //http://www.programiz.com/python-programming/examples/multiply-matrix used to troubleshoot triple for loop
    public Matrix matrixMultiplier(Matrix secondMatrix) {
        if (this.matrixColumns != secondMatrix.matrixRows) {

            System.out.println("Number of columns of 1st matrix do not match the rows of 2nd.");
            return null;
        } else {
            this.resultantMatrix = new double[this.matrixRows][secondMatrix.matrixColumns];
            for (int z = 0; z < matrixRows; z++) {
                for(int zz = 0; zz < matrixRows; zz++) {
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