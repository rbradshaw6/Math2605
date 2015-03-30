

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
}
