

public class Tests {

    public static void main(String[] args) {
        // Tests go here.
        
        //2x2 test for QR
        double[][] m1 = {{1, 1}, {1, 2}};
        Matrix m = new Matrix(m1);
        m.qr_fact_givens();
        m.qr_fact_househ();
        
        //2x2 test for LU
        double[][] m3 = {{1, 1}, {1, 2}};
        Matrix m2 = new Matrix(m1);
        m2.lu_fact();
        
        
    }
}
