
README PART 2:

Example tests are in the main method of ConvoCodes.java

-------------------------

To encode a stream:

First create a Vector object from the stream.

Then, to get Y0 and Y1, run:

ConvoCodes.encodeY0(Vector plaintext)

and

ConvoCodes.encodeY1(Vector plaintext)

these will return the encoded vector, as well as printing them out to the console.

-------------------------


To run gauss_seidel:

run gauss_seidel(Matrix a, Vector y, Vector x0, double tol)
where y is the input stream and x0 is the initial guess.
