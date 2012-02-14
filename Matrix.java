class Matrix {

    //Class variables
    private int rows;
    private int cols;
    
    private double[][] data;

    //Constructor to create the empty matrix
    public Matrix(int rowsIn, int colsIn){
	//Iterator variables
	int cCount, rCount;

	rows = rowsIn;
	cols = colsIn; 

	data = new double[rows][cols];	

    }

    //Constructor to fill matrix with a double array
    public Matrix(int rowsIn, int colsIn, double[][] fill){
	//Iterator variables
	int cCount, rCount;

	rows = rowsIn;
	cols = colsIn; 

	data = fill.clone();

    }

    //Constructor to add a single array
    public Matrix(int rowsIn, int colsIn, double[] fill){
	//Iterator variables
	int cCount, rCount;

	rows = rowsIn;
	cols = colsIn; 

	data = new double[rows][cols];

	//If there is only one column, fill the column
	if (cols == 1){
	    for (rCount = 0; rCount < rows; ++rCount){
		data[rCount][0] = fill[rCount];
	    }  
	}
	else{ //Else, fill the first row
	    data[0] = fill;	    
	}
    }

    //Function to return an item
    public double getItem (int row, int col){   
	return data[row][col];
    }

    //Function to set an item
    public double setItem (int row, int col, double value){
	data[row][col] = value;
	return value;
    }

    //Function to return a row
    public Matrix row (int rowIn){

	double[] rowVals = new double[cols];
	
	System.arraycopy(data[rowIn], 0, rowVals, 0, data[rowIn].length);

	return new Matrix(1, cols, rowVals); 
    }

    //Function to return a column
    public Matrix col (int colIn){

	int rCount;
	
	double[] colVals = new double[rows];
	
	for (rCount = 0; rCount < rows; ++rCount){
	    colVals[rCount] = data[rCount][colIn];
	}
	
	return new Matrix(rows, 1, colVals); 
    }	

    //Function to return the matrix as a list
    //One row at a time
    public double[] as_list (){
	int rCount, cCount, counter = 0;
	double[] retVal = new double[cols*rows];
	
	for (rCount = 0; rCount < rows; rCount++){
	    for (cCount = 0; cCount < cols; cCount++){
		retVal[counter] = data[rCount][cCount];
		counter++;
	    }
	} 

	return retVal;
    }

    //Function to return the array as a list
    public String toString(){
	double[] retVal = this.as_list();

	return java.util.Arrays.toString(retVal);
    }

    //Function to create the identity matrix
    public Matrix identity(int rows){
	int rCount;
	Matrix retMatrix = new Matrix(rows, rows);

	for (rCount = 0; rCount < rows; ++rCount){
	    retMatrix.setItem(rCount, rCount, 1);
	}

	return retMatrix;
    }

    //Lays the values along the diagonal
    public Matrix diagonal(double[] vals){
	int rCount;
	Matrix retMatrix = new Matrix(vals.length, vals.length);

	for (rCount = 0; rCount < vals.length; ++rCount){
	    retMatrix.setItem(rCount, rCount, vals[rCount]);
	}

	return retMatrix;
    }

    //Adds to matrices together provided they have the same dimensions
    public Matrix add (Matrix matA, Matrix matB){
	
	int rCount, cCount;
	Matrix matC = null;

	if (checkDimensions(matA, matB)){
	    
	    matC = new Matrix(matA.rows, matB.rows);

	    for (rCount = 0; rCount < matC.rows; rCount++){
		for (cCount = 0; cCount < matC.cols; cCount++){
		    matC.setItem(rCount, cCount, matA.getItem(rCount, cCount) + matB.getItem(rCount, cCount));
		}
	    }	   

	}

	return matC;
    }

    public double scalar_mult(Matrix other) {
  if (cols == 1 && other.cols == 1 && rows == other.rows) {
      double sum = 0.;
      for (int r = 0; r < rows; r++)
          sum += getItem(r, 0)*other.getItem(r, 0);

      return sum;
  }
  else throw new ArithmeticException("Incorrect matrix dimensions");
    }

    public static Matrix mult(Matrix A, double[] B_values) {
  double[][] B_values_fill = new double[B_values.length][1];
  for (int i = 0; i < B_values.length; i++)
      B_values_fill[i][0] = B_values[i];

  Matrix B = new Matrix(B_values.length, 1, B_values_fill);
  System.out.println("B:");
  B.printMatrix();
  System.out.println();
  return mult(A, B);
    }

    /*
     * M = Matrix(A.rows,B.cols)
     * for r in xrange(A.rows):
     *     for c in xrange(B.cols):
     *         for k in xrange(A.cols):
     *             M[r,c] += A[r,k]*B[k,c]
     */
    public static Matrix mult(Matrix A, Matrix B) {
  if (A.cols != B.rows)
      throw new ArithmeticException("Incorrect matrix dimensions");

        Matrix result = new Matrix(A.rows, B.cols);
  for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < B.cols; c++) {
    for (int k = 0; k < A.cols; k++) {
        double left  = A.getItem(r, k);
        double right = B.getItem(k, c);
        double old = result.getItem(r, c);
        result.setItem(r, c, old+left*right);
    }
      }
  }
  return result;
    }

    //Helper function to ensure the dimensions are the same
    public boolean checkDimensions(Matrix matA, Matrix matB){
	if ((matA.cols == matB.cols) && (matA.rows == matB.rows)){
	   return true;
	 }
	else{
	    return false;
	} 
    }

    //Transpose
    public static Matrix transpose (Matrix matA){

	int rCount, cCount;
	Matrix matT = new Matrix(matA.cols, matA.rows);

	for (rCount = 0; rCount < matT.rows; rCount++){
	  for (cCount = 0; cCount < matT.cols; cCount++){
	      matT.setItem(rCount, cCount, matA.getItem(cCount, rCount));
	  }
	}	   

	return matT;
    }

    //A function to print the matrix
    public void printMatrix(){
	int cCount, rCount;
	String currLine;

	for (rCount = 0; rCount < rows; ++rCount){
	    
	    currLine = "["; 
	    
	    for (cCount = 0; cCount < cols; ++cCount){
		currLine += data[rCount][cCount] + " ";
	    }

	    currLine += "]";

	    System.out.println(currLine);
	}

    }    

    public static void main(String[] args){

	//Using this main class for testing
	Matrix testMatrix = new Matrix(6, 6);
	double testVal = 45.4543243;

	testMatrix.setItem(1, 1, testVal);
	testMatrix.setItem(1, 2, testVal);

	Matrix newMatrix = testMatrix.col(1);

	double[] vals = new double[4];
	vals[0] = 0;
	vals[1] = 1;
	vals[2] = 2;
	vals[3] = 3;

	testMatrix = newMatrix.diagonal(vals);
	newMatrix = newMatrix.identity(4);
	
	Matrix threeMatrix = testMatrix.add(testMatrix, newMatrix);

	threeMatrix.printMatrix();

  double[][] matA_values = {
      { 14, 9,  3  },
      { 2,  11, 15 },
      { 0,  12, 17 },
      { 5,  2,  3  }
  };

  double[][] matB_values = {
      { 12, 25 },
      { 9,  10 },
      { 8,  5  }
  };

  Matrix matA = new Matrix(4, 3, matA_values);
  Matrix matB = new Matrix(3, 2, matB_values);
  mult(matA, matB).printMatrix();

  double[] matC_values = { 2, 3, 4 };
  mult(matA, matC_values).printMatrix();	
    }
}
