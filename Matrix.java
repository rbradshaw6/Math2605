import java.util.Scanner;
import java.util.Random;

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

            return transposeMatrix;
       }

       public Vector matrixVectorMult(Vector v) {
           double[][] mPlaceholder = new double[v.length][1];
           for (int i = 0; i < v.length; i++) {
               mPlaceholder[i][0] = v.vector[i];
           }
           Matrix result = matrixMultiplier(new Matrix(mPlaceholder));
           double[] newVector = new double[v.length];
           for (int i = 0; i < v.length; i++) {
               newVector[i] = result.arrayVersion[i][0];
           }
           return new Vector(newVector);
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


            Matrix multipliedMatrix = new Matrix(resultantMatrix);
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
            System.out.println();
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
    
    
    //http://www.iaa.ncku.edu.tw/~dychiang/lab/program/mohr3d/source/Jama%5CQRDecomposition.html 
    public void qr_fact_househ() {
        double[] rd = new double[matrixColumns];
        double[][] QR = Matrix1;
        for (int i = 0; i < matrixColumns; i++) {
           double k = 0;
           for(int j = i; j < matrixRows; j++) {
               k = Math.hypot(k, QR[j][i]);
           }
           
           if (k != 0.0) {
               if (QR[i][i] < 0) {
                   k = -k;
                }
                for (int j = i; j < matrixRows; j++) {
                   QR[j][i] /= k;
                }
                QR[i][i] += 1.0;
              
               
               for (int l = i+1; l < matrixColumns; l++) {
                   double s = 0.0; 
                   for (int j = i; j < matrixRows; j++) {
                      s += QR[j][i]*QR[j][l];
                   }
                   s = -s/QR[i][i];
                   for (int j = i; j < matrixRows; j++) {
                      QR[j][l] += s*QR[j][i];
                   }
                }
            }
           rd[i] = -k;
        }

        double[][] Q = new double[matrixRows][matrixColumns];
        for (int k = matrixColumns-1; k >= 0; k--) {
           for (int i = 0; i < matrixRows; i++) {
              Q[i][k] = 0.0;
           }
           Q[k][k] = 1.0;
           for (int j = k; j < matrixColumns; j++) {
              if (QR[k][k] != 0) {
                 double s = 0.0;
                 for (int i = k; i < matrixRows; i++) {
                    s += QR[i][k]*Q[i][j];
                 }
                 s = -s/QR[k][k];
                 for (int i = k; i < matrixRows; i++) {
                    Q[i][j] += s*QR[i][k];
                 }
              }
           }
        }
        Q = arrayErrorFix(Q);
        System.out.println("Q is:");
        Matrix X = new Matrix(Q);
        X.display();
        
       
        double[][] R = new double[matrixColumns][matrixColumns];
        for (int i = 0; i < matrixColumns; i++) {
           for (int j = 0; j < matrixColumns; j++) {
              if (i < j) {
                 R[i][j] = QR[i][j];
              } else if (i == j) {
                 R[i][j] = rd[i];
              } else {
                 R[i][j] = 0.0;
              }
           }
        }
        System.out.println("R is:");
        R = arrayErrorFix(R);
        X = new Matrix(R);
        X.display();
           
    }
    //http://stackoverflow.com/questions/13438073/qr-decomposition-algorithm-using-givens-rotations
    public void qr_fact_givens() {
        int m = matrixRows;
        int n = matrixColumns;
        Matrix Q = createIdentity();
        Matrix R = this;
        
        for (int j = 0; j < n; j++) {
            for(int i = (m - 1); i > j; i--) {
                Matrix G = createIdentity();
                double c = givensC(R.Matrix1[i - 1][j], R.Matrix1[i][j]);
                double s = givensS(R.Matrix1[i - 1][j], R.Matrix1[i][j]);
                G.Matrix1[i-1][i-1] = c;
                G.Matrix1[i-1][i] = -s;
                G.Matrix1[i][i-1] = s;
                G.Matrix1[i][i] = c;
                Q = Q.matrixMultiplier(G);
                R = G.matrixMultiplier(R);
                
                
            }

        }
        
        System.out.println("Q is:");
        Q.display();
        System.out.println("R is:");
        R.display();
    }
    
    
    private double givensC(double a, double b) {
        if (b == 0) {
            return 1;
        } else {
            if (Math.abs(b) > Math.abs(a)) {
                double r = a / b;
                double s = 1 / Math.sqrt(1 + Math.pow(r, 2));
                return s*r;
                
            } else {
                double r = b / a;
                return 1/ Math.sqrt(1 + Math.pow(r, 2));
            }
           
        }
    }

    private double givensS(double a, double b) {
        if (b == 0) {
            return 0;
        } else {
            if (Math.abs(b) > Math.abs(a)) {
                double r = a / b;
                return 1 / Math.sqrt(1 + Math.pow(r, 2));
                
            } else {
                double r = b / a;
                double c = 1/ Math.sqrt(1 + Math.pow(r, 2));
                return c * r;
            }
           
        }
    }

    public void qr_given() {
 
        int m = matrixRows;
        int n = matrixColumns;
        Matrix Q = createIdentity();
        Matrix Gn = createIdentity();
        Matrix A = new Matrix(Matrix1);
         
         
        // The for loops that begin the Givens rotation matrices.
 
        for (int i=0; i<n; i++) {
                for (int j=(n-1); j>i; j--) {                                       
                     
                    double a = Matrix1[j-1][i];
                    double b = Matrix1[j][i];   
                    double c = a/(Math.sqrt(a*a+b*b));
                    double s = -b/(Math.sqrt(a*a+b*b));
                     
                    Gn.Matrix1[j][j] = c;
                    Gn.Matrix1[j][j-1] = s;
                    Gn.Matrix1[j-1][j] = -s;
                    Gn.Matrix1[j-1][j-1] = c;           
 
                    A = Gn.matrixMultiplier(A);
               
                     
                    Q = Gn.matrixMultiplier(Q);
 
                    // Turning the Gn matrix back into the identity.
                     
                    Gn = Gn.createIdentity();
                                
                }//end j    
        }//end i
         
         
        System.out.println("Q:");
        Q.display();
         
        System.out.println("R:");
        A.display();
         
        Q = Q.transpose();
        Matrix answer = Q.matrixMultiplier(A);
        System.out.println("Did it work? Q x R:");
        answer.display();
        /*
        // Calculating the maximum norm of QR-H.
         
        answer = subtractMatrix(answer, Hilbert);
        double maxNorm = normOfInfinity(answer);
        System.out.println();
        System.out.println("The maximum norm of QR - H:");
        System.out.println(maxNorm);
        */
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
    
    public double getNorm() {
    	Matrix a = this;
    	Matrix b = this.transpose();
    	
    	Matrix c = a.matrixMultiplier(b);
    	double yolo = c.trace();
    	double yolosqrt = Math.sqrt(yolo);
    	System.out.println("Norm of matrix:" + yolosqrt);
    	return yolosqrt;
    	
    }
    
    public double trace() {
    	double[][] matrix = this.arrayVersion;
    	double sum = 0;
    	for (int i = 0; i < matrixRows; i++) {
    		for (int j = 0; j < matrixColumns; j++) {
    			if (i == j) {
    				sum += matrix[i][j];
    			}
    		}
    	}
    	return sum;
    }
    
    public void power_method() {
    	Matrix matrix = this;
    	double[] randomVector = new double[matrixRows];
    	
    	Random rand = new Random(); //creating the first random vector
    	for (int i = 0; i < matrixRows; i++) {
    		double randomValue = rand.nextDouble();
    		randomVector[i] = randomValue;
    	}
    	
    	Vector[] steps = new Vector[10001];
    	
    	double[] u = randomVector;
    	Vector uu = new Vector(u);
    	Vector mult = matrix.matrixVectorMult(uu);
    	double[] g = mult.vector;
    	steps[0] = mult.multiplyByConstant(1/g[0]);
    	
    	for (int i = 0; i < 10000; i++) {
    			Vector lol = matrix.matrixVectorMult(steps[i]);
    			double[] lolz = lol.vector;
    			Vector what = lol.multiplyByConstant(1/steps[i].vector[0]);
    			steps[i + 1] = what;
    	}
    	
    	steps[10000].display();
    	
    	/*

    	Vector v = new Vector(randomVector);
    	Vector iterate = matrix.matrixVectorMult(v);
    	double[] b = iterate.vector;
    	*/
    	
    	/*
    	for (int i = 0; i < 10000; i++) {    		
    		double mag = iterate.magnitude();
    		
    		for (int j = 0; j < b.length; j++) {
    			b[j] /= b[0];
    		}
    	}
    	
    	Vector newV = new Vector(b);
    	newV.display();
    	*///////////////
    }
    
    public double[] getEigenvalues() {
    	return null;
    }
    
    private double[][] arrayErrorFix(double[][] array) {
        for(int i = 0; i < array.length; i++) {
            for(int j = 0; j < array[0].length; j++) {
                if(array[i][j] != 0) {
                    array[i][j] *= -1;
                }
            }
        }
        
        return array;
    }

    public static void main(String[] args) {
        double[][] m1 = {{1, 1}, {1, 2}};
       // double[] b = {1,1};
       // Vector g = new Vector(b);
        //System.out.println(g.magnitude());
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
        m.qr_given();
        //m.display();
        //m.power_method();
    }
 

}
