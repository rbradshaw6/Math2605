

public class Leslie {


    public static double power_method(Matrix a, double tol, Vector u) {
        double diff = 0;
        int i = 0;
        while ((i < 10000 && diff > tol) || i == 0) {
            Vector newu = a.matrixVectorMult(u).multiplyByConstant(1 / u.vector[0]);
            //newu.display();
            diff = Math.abs(u.vector[0] - newu.vector[0]);
            u = newu;
            i++;
        }
        System.out.println("Iterations (max 10000):");
        System.out.println(i);
        System.out.println("Error:");
        System.out.println(diff);
        if (i == 10000) {
            System.out.println("Tolerance not reached. After 10000 iterations, error is:" + diff);
        }
        System.out.println("Final alpha value:");
        return u.vector[0];
    }

    public static Vector calcPop(Matrix a, Vector x, int decades) {
        for (int i = 0; i < decades; i++) {
            x = a.matrixVectorMult(x);
        }
        System.out.println("Population after " + decades + " deades:");
        x.display();
        return x;
    }

    public static void main(String[] args) {
        double[][] a = {{0, 0.6, 1.1, .9, .1, 0, 0, 0, 0},
                {.7, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, .85, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, .9, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, .9, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, .88, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, .8, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, .77, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, .40, 0}};
        Matrix leslie = new Matrix(a);
        double tol = Math.pow(10, -8);
        double[] x = {2.1, 1.9, 1.8, 2.1, 2.0, 1.7, 1.2, 0.9, 0.5};
        double[] bControl = {5.188, 4.445, 1.25, 1.453, 1.458, 1.663, 1.408, 1.047, 0.37};
        Vector u0 = new Vector(x);
        System.out.println("Power method:");
        double alpha = Leslie.power_method(leslie, tol, u0);

        calcPop(leslie, u0, 1);
        calcPop(leslie, u0, 2);
        calcPop(leslie, u0, 3);
        calcPop(leslie, u0, 4);
        calcPop(leslie, u0, 5);


    }
}
