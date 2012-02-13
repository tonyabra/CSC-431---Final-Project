class Matrix {

    //Class variables
    private int rows;
    private int cols;
    
    private float[][] data;

    public Matrix(int rowsIn, int colsIn, float fill){
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

    //Function to return an item
    public float getItem (int row, int col){   
	return data[row][col];
    }

    //Function to set an item
    public float setItem (int row, int col, float value){
	data[row][col] = value;
	return value;
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

	Matrix testMatrix = new Matrix(2, 2, 1);

	testMatrix.printMatrix();
    }
}