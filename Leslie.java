

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
        System.out.println("Iterations:");
        System.out.println(i);
        System.out.println("Error:");
        System.out.println(diff);
        if (i == 10000) {
            throw new IllegalArgumentException("Tolerance not reached. After 1000 iterations, error is:" + diff);
        }
        return u.vector[0];
    }

    public static Vector calcPop(Matrix a, Vector x, int years) {
        for (int i = 0; i < years; i++) {
            x = a.matrixVectorMult(x);
        }
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
        Vector u0 = new Vector(bControl);
        System.out.println("Power method:");
        double alpha = Leslie.power_method(leslie, tol, u0);
        System.out.println("Alpha:");
        System.out.println(alpha);

        System.out.println("\nPopulation in 2010:");
        calcPop(leslie, u0, 1).display();
        System.out.println("\nPopulation in 2020:");
        calcPop(leslie, u0, 2).display();
        System.out.println("\nPopulation in 2030:");
        calcPop(leslie, u0, 3).display();
        System.out.println("\nPopulation in 2040:");
        calcPop(leslie, u0, 4).display();
        System.out.println("\nPopulation in 2050:");
        calcPop(leslie, u0, 5).display();


    }
}
