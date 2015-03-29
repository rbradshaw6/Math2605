
/**
* Okay look. Vectors and matrices are IMMUTABLE (because Bradshaw started it
* that way, and I'm too lazy to change it). That means that all ops - mult,
* add, etc., create NEW objects and return them as the answer. They do not
* change the underlying matrix.

* Vectors are n x 1, meaning n rows and 1 column. This is important for
* vector matrix multiplication, so be careful.
*/
public class Vector {
    protected double[] vector;
    protected int length;

    public Vector(double[] vector) {
        this.length = vector.length;
        this.vector = vector;
    }

    public double dotProduct(Vector v) {
        if (v.length != length) {
            throw new IllegalArgumentException("Vectors must have same length.");
        }
        double sum = 0;
        for (int i = 0; i < length; i++) {
            sum += vector[i] * v.vector[i];
        }
    }

    public Vector add(Vector v) {
        if (v.length != length) {
            throw new IllegalArgumentException("Vectors must have same length.");
        }
        double[] sum = new double[length];
        for (int i = 0; i < length; i++) {
            sum[i] = vector[i] + v.vector[i];
        }
    }
}
