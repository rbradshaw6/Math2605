

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

    public static Vector jacobi(Vector y0) {
        //TODO Use jacobi iterations to decrypt the vector. after every
        // iteration, mod the vector by 2
        // will use y0 stream to decrypt.

        int length = y0.length;
        double[] xArr = new double[length];
        Vector x = new Vector(xArr);

        // Make D and R Vectors from the x -> y0 encryption matrix.
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
        double[][] dMatrix = new double[length][length];
        double[][] rMatrix = new double[length][length];
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if (i==j) {
                    dMatrix[i][j] = encryptionArray[i][j];
                } else {
                    rMatrix[i][j] = encryptionArray[i][j];
                }
            }
        }
        Matrix d = new Matrix(dMatrix);
        Matrix r = new Matrix(rMatrix);
        // d and r matrices have been made. now iterate
        int i = 0;

        while (i < 1000 && ConvoCodes.calcError(new Matrix(encryptionArray), x, y0) > Math.pow(10, -8)) {
            Vector rx = r.matrixVectorMult(x);
            x = d.inverse().matrixVectorMult(y0.subtract(rx));
            for (int k = 0; k < x.length; k++) {
                x.vector[k] = Math.abs(x.vector[k] % 2);
            }

            // double[] newx = new double[x.length];
            // for (int j = 0; j < newx.length; j++) {
            //     double sum = 0;
            //     for (int k = 0; k < encryptionArray[j].length; k++) {
            //         if (j != k) {
            //             sum += encryptionArray[j][k] * x.vector[j]
            //         }
            //     }
            //     newx[j] = (y0.vector[j] - sum) / encryptionArray[j][j];
            // }
            // x = new Vector(newx);


            i++;
        }
        double[] finalArr = new double[length-3];
        for (int j = 0; j < finalArr.length; j++) {
            finalArr[j] = x.vector[j];
        }
        x = new Vector(finalArr);
        return x;
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
        Vector x = ConvoCodes.jacobi(y0);
        x.display();
    }
}
