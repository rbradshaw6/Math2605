

public class ConvoCodes {

    public static Vector encodeY0(Vector plaintext) {
        int length = plaintext.length;
        double[][] encryptionArray = new double[length + 3][length];
        for (int i = 0; i < length + 3; i++) {
            if (i < length) {
                encryptionArray[i][i] = 1;
            }
            if (i - 2 < length && i - 2 >= 0) {
                encryptionArray[i][i-2] = 1;
            }
            if (i - 3 < length && i - 3 >= 0) {
                encryptionArray[i][i-3] = 1;
            }
        }
        Matrix encryptionMatrix = new Matrix(encryptionArray);
        Vector v = encryptionMatrix.matrixVectorMult(plaintext);
        for (int i = 0; i < v.length; i++) {
            v.vector[i] = v.vector[i] % 2;
        }
        return v;
    }

    public static Vector encodeY1(Vector plaintext) {
        int length = plaintext.length;
        double[][] encryptionArray = new double[length + 3][length];
        for (int i = 0; i < length + 3; i++) {
            if (i < length) {
                encryptionArray[i][i] = 1;
            }
            if (i - 1 < length && i - 1 >= 0) {
                encryptionArray[i][i-2] = 1;
            }
            if (i - 3 < length && i - 3 >= 0) {
                encryptionArray[i][i-3] = 1;
            }
        }
        Matrix encryptionMatrix = new Matrix(encryptionArray);
        Vector v = encryptionMatrix.matrixVectorMult(plaintext);
        for (int i = 0; i < v.length; i++) {
            v.vector[i] = v.vector[i] % 2;
        }
        return v;
    }

    public static Vector jacobi(Matrix a, Vector y, Vector x0, double tol) {
        //TODO Use jacobi iterations to decrypt the vector. after every
        // iteration, mod the vector by 2
        // will use y stream to decrypt.
        if (a.matrixRows != a.matrixColumns) { throw new IllegalArgumentException("Matrix must be square!"); }
        int length = a.matrixRows;

        double[][] dMatrix = new double[length][length];
        double[][] rMatrix = new double[length][length];
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i==j) {
                    dMatrix[i][j] = a.arrayVersion[i][j];
                } else {
                    rMatrix[i][j] = a.arrayVersion[i][j];
                }
            }
        }
        Matrix d = new Matrix(dMatrix);
        Matrix r = new Matrix(rMatrix);
        int i = 0;

        while (i < 1000 && ConvoCodes.calcError(a, x0, y) > tol) {
            Vector rx = r.matrixVectorMult(x0);
            x0 = d.inverse().matrixVectorMult(y.subtract(rx));
            for (int k = 0; k < x0.length; k++) {
                x0.vector[k] = Math.abs(x0.vector[k] % 2);
            }

            // double[] newx = new double[x.length];
            // for (int j = 0; j < newx.length; j++) {
            //     double sum = 0;
            //     for (int k = 0; k < encryptionArray[j].length; k++) {
            //         if (j != k) {
            //             sum += encryptionArray[j][k] * x.vector[j]
            //         }
            //     }
            //     newx[j] = (y.vector[j] - sum) / encryptionArray[j][j];
            // }
            // x = new Vector(newx);


            i++;
        }
        double[] finalArr = new double[length-3];
        for (int j = 0; j < finalArr.length; j++) {
            finalArr[j] = x0.vector[j];
        }
        x0 = new Vector(finalArr);
        return x0;
    }

    public static Vector gauss_seidel(Matrix a, Vector y, Vector x0, double tol) {
        if (a.matrixRows != a.matrixColumns) { throw new IllegalArgumentException("Matrix must be square!"); }
        int length = a.matrixRows;

        double[][] lMatrix = new double[length][length];
        double[][] uMatrix = new double[length][length];
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i>=j) {
                    lMatrix[i][j] = a.arrayVersion[i][j];
                } else {
                    uMatrix[i][j] = a.arrayVersion[i][j];
                }
            }
        }
        Matrix l = new Matrix(lMatrix);
        Matrix u = new Matrix(uMatrix);
        int i = 0;

        //---

        while (i < 1000 && ConvoCodes.calcError(a, x0, y) > tol) {
            Vector ux = u.matrixVectorMult(x0);
            x0 = l.inverse().matrixVectorMult(y.subtract(ux));
            for (int k = 0; k < x0.length; k++) {
                x0.vector[k] = Math.abs(x0.vector[k] % 2);
            }
            i++;
        }
        double[] finalArr = new double[length-3];
        for (int j = 0; j < finalArr.length; j++) {
            finalArr[j] = x0.vector[j];
        }
        x0 = new Vector(finalArr);
        return x0;
    }

    public static double calcError(Matrix a, Vector x, Vector b) {
        Vector b1 = a.matrixVectorMult(x);
        for (int i=0; i < b1.vector.length; i++) {
            b1.vector[i] = b1.vector[i] % 2;
        }
        Vector c = b.subtract(b1);
        return c.norm();
    }

    public static void main(String[] args) {
        double[] data = {1, 0, 1, 1, 0};
        Vector plaintext = new Vector(data);
        plaintext.display();
        System.out.println("\n\n");
        Vector y0 = ConvoCodes.encodeY0(plaintext);
        y0.display();
        System.out.println("\n\n");

        int length = y0.length;
        double[][] encryptionArray = new double[length][length];
        for (int i = 0; i < length; i++) {
            if (i < length) {
                encryptionArray[i][i] = 1;
            }
            if (i - 2 < length && i - 2 >= 0) {
                encryptionArray[i][i-2] = 1;
            }
            if (i - 3 < length && i - 3 >= 0) {
                encryptionArray[i][i-3] = 1;
            }
        }
        double[] guess = {0,0,0,0,0,0,0,0};
        Vector guess0 = new Vector(guess);
        Vector x = ConvoCodes.gauss_seidel(new Matrix(encryptionArray), y0, guess0, Math.pow(10, -8));
        x.display();
    }
}
