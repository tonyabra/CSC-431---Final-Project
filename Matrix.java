class Matrix {

    //Class variables
    private int rows;
    private int cols;
    
    private double[][] data;

    //Constructor to fill matrix with value
    public Matrix(int rowsIn, int colsIn, double fill){
	//Iterator variables
	int cCount, rCount;

	rows = rowsIn;
	cols = colsIn; 
	
	data = new float[rows][cols];

	for (rCount = 0; rCount < rows; ++rCount){
	   for (cCount = 0; cCount < cols; ++cCount){
	       data[rCount][cCount] = fill;
	    }
	}		
    }

    //Constructor to fill matrix with value
    public Matrix(int rowsIn, int colsIn, double fill){
	//Iterator variables
	int cCount, rCount;

	rows = rowsIn;
	cols = colsIn; 
	
	data = new float[rows][cols];
	
    }

    //Function to return an item
    public double getItem (int row, int col){   
	return data[row][col];
    }

    //Function to set an item
    public double setItem (int row, int col, float value){
	data[row][col] = value;
	return value;
    }

    //Function to return a row
    public double
    

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

	Matrix testMatrix = new Matrix(2, 2, 1);

	testMatrix.printMatrix();
    }
}