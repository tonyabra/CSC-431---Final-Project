// CSC431 Final Project 
// Members: Tony Abraham, Heidi Cochran, Larcen Razi, Mike Storm

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


public class Numeric {
	
	//This function returns the square root for the supplied value
	public static double sqrt(double x){
        if (x>=0)return Math.sqrt(Math.abs(x));
        else throw new UnsupportedOperationException();
	}

	//This function is the overloaded condition number function
	//This first implementation accepts a function and a value
	public static double condition_number(Function f, double x){
		return new Derivative(f).execute(x) *x/f.execute(x);
	}
	
	//This is the second overloaded condition number function
	//This function accepts a matrix and determines the condition number
	public static double condition_number(Matrix m){
		double val;
		
		if(m instanceof Matrix){
			Matrix tempM = new Matrix(m.rowCount(), m.colCount());
			for(int i = 0; i < m.rowCount(); i++){
				for(int k = 0; k < m.colCount(); k++){
					val = m.getItem(i, k);
					tempM.setItem(i, k, 1/val);
				}
			}
			
			return norm(m)*norm(tempM);
		}
		  else throw new UnsupportedOperationException();
	}

	//This is the exponential function that raises a matrix to a power
	public static Matrix exp(Matrix x, double...params){
		double ap = 1* 10^-6;
		double rp = 1*10^-4;
		double ns=40;
		if(params.length>0) {
			ap=params[0];
		}
		if(params.length>1){
			rp=params[1];
		}
		if(params.length>2){
			ns=params[2];
		}
	       Matrix t =  Matrix.identity(x.colCount());
	       Matrix s =t;
	       for(int k=1; k<=ns;k++)
	       {	    	   
	           t = Matrix.div(Matrix.mult(t, x), k); // next term
	           s = Matrix.add(s , t);   //add next term
	           if((norm(t)<Math.max(ap,norm(s)*rp)))
	        		   { 
	        	  
	       throw new ArithmeticException("no convergence"); 	   	                    
	       }   
	    } 
	    return s;
	}
	
	//This function checks whether or not the supplied Matrix is positive definite
	public boolean is_positive_definite(Matrix A){
	    if(! is_symmetric(A))
	        {
	        return false; 
	        }
	    try{
	        Cholesky(A);
	        return true;
	      
	       }catch(Exception e)
	           {
	    	   return false;
	           }
	}

	
	//This function determines if the Matrix is symmetric
	public static boolean is_symmetric(Matrix A) {
        int N = A.colCount();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (A.getItem(i,j) != A.getItem(j,i)) return false;

            }
         
        }
        return true;
	}	
	
	//This function determines if the Matrix is close to symmetric
    public static boolean is_almost_symmetric(Matrix matrix){	
	double ap = 0.000001;
	double rp = 0.0001;
		
		return (is_almost_symmetric(matrix, ap, rp));
    }
	
    //This function determines if the matrix is close to symmetric while adding 
    //the accuracy values
    public static boolean is_almost_symmetric(Matrix matrix, double ap, double rp) {
    	
	if (matrix.rowCount() != matrix.colCount())
	    return false;

	for (int r = 0; r < matrix.rowCount()-1; r++) {
	    for (int c = 0; c < r; c++) {
	        double delta = Math.abs(matrix.getItem(r, c) - matrix.getItem(c, r));
		if (delta > ap && delta > Math.max(Math.abs(matrix.getItem(r, c)), Math.abs(matrix.getItem(c, r)))*rp)
		    return false;
	    }
	}
	return true;
    }
    
    public static Matrix FitLeastSquares(double[][] listOfPoints, MyFunctionInterface[] listOfFunctions){
    	int sizeOfPoints = listOfPoints.length;
    	int sizeOfFunctions = listOfFunctions.length;
    	MyFunctionInterface functionObj = null;
    	double funcReturn = 0;
    	double weight = 0;
    	Matrix A = new Matrix(sizeOfPoints, sizeOfFunctions);
    	Matrix b = new Matrix(sizeOfPoints, 0);//Check if this OK
    	Matrix c = null;
    	int nbrOfRows = A.rowCount();
    	for(int i=0; i< nbrOfRows; i++){
    	if((listOfPoints[i]).length > 2){
    	weight = 1/listOfPoints[i][2];
    	}else{
    	weight = 1;
    	}
    	b.setItem(i,0,weight*listOfPoints[i][0]);
    	for(int j=0; j<A.colCount(); j++){
    	/*get the passed function, which is an implementation of MyFunctionInterface*/
    	functionObj = listOfFunctions[j];
    	funcReturn = functionObj.function(listOfPoints[i][0]);
    	A.setItem(i,j,weight*funcReturn );
    	}
    	}
    	c = Matrix.rDiv( (Matrix.mult(Matrix.transpose(A),A)),1);

    	return c;
    	}


    	interface MyFunctionInterface {
    	    double function(double x);
    	}
    	class MyFunctionImpl implements MyFunctionInterface {
    	public double function(double x){
    	return x*x;
    	}
    	}

        static class MyFunctionImplOne implements MyFunctionInterface {
        	public double function(double x){
        	return x*x;
        	}
        	}

        	static class MyFunctionImplTwo implements MyFunctionInterface {
        	public double function(double x){
        	return x*x*x;
        	}
        	}

    	


    //Self-explanatory function to determine how close to zero Matrix values are
    public static boolean is_almost_zero(Matrix matrix) {
      return is_almost_zero(matrix, 0.0000001);
    }

    public static boolean is_almost_zero(Matrix matrix, double ap) {
      return is_almost_zero(matrix, ap, 0.0000001);
    }

    public static boolean is_almost_zero(Matrix matrix, double ap, double rp) {
	double delta;
    for (int r = 0; r < matrix.rowCount(); r++) {
	    for (int c = 0; c < matrix.colCount(); c++) {
		delta = Math.abs(matrix.getItem(r, c) - matrix.getItem(c, r));
		if (delta > ap && delta > Math.max(Math.abs(matrix.getItem(r, c)), Math.abs(matrix.getItem(c, r)))*rp)
		    return false;
	    }
	}
    
	return true;
    }
    

    //The norm function, overloaded to accept just a Matrix
    public static double norm(Matrix matrix) {
        return norm(matrix, 1);
      }

      /*
       * def norm(A,p=1):
       *   if isinstance(A,(list,tuple)):
       *     return sum(x**p for x in A)**(1.0/p)
       *   elif isinstance(A,Matrix):
       *     if A.rowCount()==1 or A.colCount()==1:
       *        return sum(norm(A[r,c])**p \
       *          for r in xrange(A.rowCount()) \
       *          for c in xrange(A.colCount()))**(1.0/p)
       *     elif p==1:
       *       return max([sum(norm(A[r,c]) \
       *         for r in xrange(A.rowCount())) \
       *         for c in xrange(A.colCount())])
       *     else:
       *       raise NotImplementedError
       *   else:
       *     return abs(A)
       
    	public static double norm(Matrix matrix, int p) {
        if (matrix.rowCount() == 1 || matrix.colCount() == 1) {
          double n = 0.;
          for (int r = 0; r < matrix.rowCount(); r++)
            for (int c = 0; c < matrix.colCount(); c++)
              n += Math.pow(matrix.getItem(r, c), p);
          return Math.pow(n, (1./((double)p)));
        }
        else if (p == 1) {
          List<Double> norms = new ArrayList<Double>();
          for (int r = 0; r < matrix.rowCount(); r++)
            for (int c = 0; c < matrix.colCount(); c++)
              norms.add(norm(matrix.getItem(r, c)));
          Collections.sort(norms);
          return norms.get(norms.size()-1);
        }
        else
          throw new RuntimeException("Not implemented");
      }
      
      */

    //The function will determine the matrix norm
    public static double norm(Matrix m, int... params){
    	int p = 1;
		if (params.length > 0) p=params[0];
    	double retSum =0;
		if(m.rowCount() == 1 || m.colCount() ==1){
			for(int r=0;r<m.rowCount();r++){
				for(int c=0; c<m.colCount(); c++){
					retSum += Math.pow(Math.abs(m.getItem(r,c)),p);
				}
			}
			retSum = Math.pow(retSum, 1.0/p);
		}else if(p==1){
			for(int r=0;r<m.rowCount();r++){
				for(int c=0; c< m.colCount(); c++){
					retSum += Math.abs(m.getItem(r,c));
				}
			}			
		}
		return retSum;
	}
    
     //The norm function for lists
      public static double norm(double[] list) {
        double[] copy = (double[])list.clone();
        Arrays.sort(copy);
        return copy[copy.length-1];
      }
      
      //The norm function for values
      public static double norm(double val) {
        return Math.abs(val);
      }
      
      //The Cholesky function for Matrices
      public static Matrix Cholesky(Matrix m) {
    	  double temp, p;
    	  	if(!is_almost_symmetric(m)){
    	  		throw new ArithmeticException("not symmetric");
    	  	}
    	  
    	  	Matrix l = (Matrix) m.clone();
    	  	
    	  	for (int k=0; k< m.colCount()-1; k++){
    	  		if(l.getItem(k,k) <=0){
    	  			throw new ArithmeticException("not positive definitive");}
    	  		p = Math.sqrt(l.getItem(k,k));
    	  		l.setItem(k,k, p);
    	  		for(int i=k+1;i< l.rowCount();i++){
    	  			p = l.getItem(i,k) / p;
    	  		}
    	  		for(int j=k+1; j< l.rowCount();j++){
    	  			p= l.getItem(j,k);
    	  			for(int i=k+1; i< l.colCount();i++){
    	  				temp = l.getItem(i,j);
    	  				temp-= p*l.getItem(i,k);
    	  				l.setItem(i,j,temp);
    	  			}
    	  		}
    	  	}
    	  	
    	  	return l;
      }

      //The following overloaded functions handle the fixed point algorithm
      public static double solve_fixed_point(Function f, double x) {
          return solve_fixed_point(f, x, 0.0000001);
        }

        public static double solve_fixed_point(Function f, double x, double ap) {
          return solve_fixed_point(f, x, ap, 0.00001);
        }

        public static double solve_fixed_point(Function f, double x, double ap, double rp) {
          return solve_fixed_point(f, x, ap, rp, 100);
        }

        public static double solve_fixed_point(Function f_arg, double x, double ap, double rp, int ns) {
          final Function f = f_arg;
          Function g = new Function() {
            public double execute(double x) {
              return f.execute(x)+x;
            }
          };
          Function Dg = new Derivative(g);
          for (int k = 0; k < ns; k++) {
            if (Math.abs(Dg.execute(x)) >= 1)
              throw new RuntimeException("error D(g)(x)>=1");
            double x_old = x;
            x = g.execute(x);
            if (k > 2 && norm(x_old-x) < Math.max(ap, norm(x)*rp))
              return x;
          }
          throw new RuntimeException("no convergence");
        }


        //The solve bisection functions overloaded to allow different levels of specificity
        public static double solve_bisection(Function f, double a, double b) {
          return solve_bisection(f, a, b, 0.0000001);
        }

        public static double solve_bisection(Function f, double a, double b, double ap) {
          return solve_bisection(f, a, b, ap, 0.00001);
        }

        public static double solve_bisection(Function f, double a, double b, double ap, double rp) {
          return solve_bisection(f, a, b, ap, rp, 100);
        }

        public static double solve_bisection(Function f, double a, double b, double ap, double rp, int ns) {
          double fa = f.execute(a);
          double fb = f.execute(b);
          if (fa == 0)
            return a;
          else if (fb == 0)
            return b;
          else if (fa*fb > 0)
            throw new RuntimeException("f(a) and f(b) must have opposite sign");

          for (int k = 0; k < ns; k++) {
            double x = (a+b)/2.;
            double fx = f.execute(x);
            if (fx == 0 || norm(b-a) < Math.max(ap, norm(x)*rp))
              return x;
            else if (fx * fa < 0) {
              b = x;
              fb = fx;
            }
            else {
              a = x;
              fa = fx;
            }
          }
          throw new RuntimeException("no convergence");
        }

      //This function is the basic Newton method for solving a function
	public static double solve_newton(Function f, double x, double ap, double rp, int ns){
		
		double fx, Dfx;
		double x_old;
		
		for (int k = 0; k < ns; k++){
			fx = f.execute(x);
			Dfx = new Derivative(f).execute(x);
			
			if (Math.abs(Dfx) < ap){
				throw new ArithmeticException("Unstable Solution");
			}
		
			x_old = x; 
			x = x - fx/Dfx;
			
			if ((k > 2) && (Math.abs(x-x_old) < Math.max(ap, (Math.abs(x)*rp)))){ 
				return x;
			}
		}
		
		throw new ArithmeticException("No Convergence.");
	}
	
	//This is the overloaded newton function
	public static double solve_newton(Function f, double x){
		double ap= 1e-6;
		double rp= 1e-4;
		int ns = 20;
		
		return solve_newton(f, x, ap, rp, ns);
	}
	
	
	//This function handles the overloaded solve secant values
	public static double solve_secant(Function f, double x, double ap, double rp, int ns){
		double fx, Dfx;
		double x_old, fx_old;

		fx = f.execute(x);
		Dfx = new Derivative(f).execute(x);
		
		for (int k = 0; k < ns; k++){
			
			if (Math.abs(Dfx) < ap){
				throw new ArithmeticException("Unstable Solution");
			}

			x_old = x; 
			fx_old = fx; 
			x = x - fx/Dfx;
			
			if ((k > 2) && (Math.abs(x-x_old) < Math.max(ap, (Math.abs(x)*rp)))){ 
				return x;
			}			
			
			fx = f.execute(x);
			Dfx = (fx -fx_old)/(x-x_old);
		}
		
		throw new ArithmeticException("No Convergence.");
	}
	
	public static double solve_secant(Function f, double x){
		double ap= 1e-6;
		double rp= 1e-4;
		int ns = 20;
		
		return solve_secant(f, x, ap, rp, ns);
	}
	
	
	//This is the stabilized newton solver
	public static double solve_newton_stabilized(Function f, double a, double b, double ap, double rp, int ns){

		double fa = f.execute(a);
		double fb = f.execute(b);
		double x, fx, Dfx, x_old, fx_old; 
		
		if (fa == 0){
			return a;
		}
		else if (fb == 0){
			return b;
		}
		else if (fa * fb > 0){
			throw new ArithmeticException("f(a) and f(b) must have opposite signs");			
		}
		
		x = (a+b)/2;
		fx = f.execute(x);
		Dfx = new Derivative(f).execute(x);
		
		for (int k = 0; k < ns; k++){
			x_old = x;
			fx_old = fx;
			
			if (Math.abs(Dfx) > ap){
				x = x - fx/Dfx;
			}
			
			if ((x == x_old) || (x < a) || (x > b)){
				x = (a+b)/2;
			}
			
			fx = f.execute(x);
			
			if ((fx == 0) || (Math.abs(x-x_old) < Math.max(ap, (Math.abs(x)*rp)))){ 
				return x;
			}
			
			Dfx = (fx - fx_old)/(x-x_old);
			if (fx * fa < 0){
				b = x; 
				fb = fx;
			}
			else{
				a = x; 
				fa = fx;
			}	
		}
		throw new ArithmeticException("No Convergence.");		
	}
	
	public static double solve_newton_stabilized(Function f, double a, double b){
		double ap= 1e-6;
		double rp= 1e-4;
		int ns = 20;

		return solve_newton_stabilized(f, a, b, ap, rp, ns);
	}
	
	
	//The following are the bisectional optimization overloaded functions
    public static double optimize_bisection(Function f, double a, double b) {
        return optimize_bisection(f, a, b, 0.0000001);
      }
      
      public static double optimize_bisection(Function f, double a, double b, double ap) {
        return optimize_bisection(f, a, b, ap, 0.00001);
      }

      public static double optimize_bisection(Function f, double a, double b, double ap, double rp) {
        return optimize_bisection(f, a, b, ap, rp, 100);
      }
      
      public static double optimize_bisection(Function f, double a, double b, double ap, double rp, int ns) {
        double Dfa = new Derivative(f).execute(a);
        double Dfb = new Derivative(f).execute(b);
        if (Dfa == 0)
          return a;
        else if (Dfb == 0)
          return b;
        else if (Dfa*Dfb > 0)
          throw new RuntimeException("First derivatives at f(a) and f(b) must have opposite sign");

        for (int k = 0; k < ns; k++) {
          double x = (a+b)/2;
          double Dfx = new Derivative(f).execute(x);
          if (Dfx == 0 || norm(b-a) < Math.max(ap, norm(x)*rp))
            return x;
          else if (Dfx * Dfa < 0) {
            b = x;
            Dfb = Dfx;
          }
          else {
            a = x;
            Dfa = Dfx;
          }
        }
        throw new RuntimeException("No convergence");
      }
	
      
      //This is the optimized newton functions overloaded for convenience
      public static double optimize_newton(Function f, double x) {
          return optimize_newton(f, x, 0.0000001);
        } 

        public static double optimize_newton(Function f, double x, double ap) {
          return optimize_newton(f, x, ap, 0.00001);
        } 

        public static double optimize_newton(Function f, double x, double ap, double rp) {
          return optimize_newton(f, x, ap, rp, 20);
        } 

        public static double optimize_newton(Function f, double x, double ap, double rp, int ns) {
          for (int k = 0; k < ns; k++) {
            double Dfx = new Derivative(f).execute(x);
            double DDfx = new Derivative(f).executeDD(x);
            if (norm(DDfx) < ap)
              throw new RuntimeException("unstable solution");
            double x_old = x;
            x = x-Dfx/DDfx;
            if (norm(x-x_old) < Math.max(ap, norm(x)*rp))
              return x;
          }
          throw new RuntimeException("no convergence");
        }   

        //This is the optimized secant functions
        public static double optimize_secant(Function f, double x) {
            return optimize_secant(f, x, 0.0000001);
          }
          
          public static double optimize_secant(Function f, double x, double ap) {
            return optimize_secant(f, x, ap, 0.00001);
          }

          public static double optimize_secant(Function f, double x, double ap, double rp) {
            return optimize_secant(f, x, ap, rp, 100);
          }

          public static double optimize_secant(Function f, double x, double ap, double rp, int ns) {
            double fx = f.execute(x);
            double Dfx = new Derivative(f).execute(x);
            double DDfx = new Derivative(f).executeDD(x);
            for (int k = 0; k < ns; k++) {
              if (Dfx == 0)
                return x;
              else if (norm(DDfx) < ap)
                throw new RuntimeException("unstable solution");

              double x_old = x;
              double Dfx_old = Dfx;
              x = x-Dfx/DDfx;

              if (norm(x-x_old) < Math.max(ap, norm(x)*rp))
                return x;
              fx = f.execute(x);
              Dfx = new Derivative(f).execute(x);
              DDfx = (Dfx - Dfx_old)/(x-x_old);
            }
            throw new RuntimeException("no convergence");
          }        
        

        //Functions that handle the optimized stabilized newton 
        public static double optimize_newton_stabilized(Function f, double a, double b, double ap, double rp, int ns){

		double Dfa = new Derivative(f).execute(a);
		double Dfb = new Derivative(f).execute(b);
		double x, fx, Dfx, DDfx, x_old, fx_old, Dfx_old; 
		
		if (Dfa == 0){
			return a;
		}
		else if (Dfb == 0){
			return b;
		}
		else if (Dfa * Dfb > 0){
			throw new ArithmeticException("Df(a) and Df(b) must have opposite signs");			
		}
		
		x = (a+b)/2;
		fx = f.execute(x);
		Dfx = new Derivative(f).execute(x);
		DDfx = new Derivative(f).executeDD(x);
		
		for (int k = 0; k < ns; k++){
			if (Dfx == 0){
				return x;
			}
			
			x_old = x;
			fx_old = fx;
			Dfx_old = Dfx;
		
			if (Math.abs(DDfx) > ap){
				x = x - Dfx/DDfx;
			}
			
			if ((x == x_old) || (x < a) || (x > b)){
				x = (a+b)/2;
			}
			
			if (Math.abs(x-x_old) < Math.max(ap, (Math.abs(x)*rp))){ 
				return x;
			}
			
			fx = f.execute(x);
			Dfx = (fx - fx_old)/(x-x_old);
			DDfx = (Dfx-Dfx_old)/(x-x_old);
			
			if (Dfx * Dfa < 0){
				b = x; 
				Dfb = Dfx;
			}
			else{
				a = x;
				Dfa = Dfx;
			}
		}
		
		throw new ArithmeticException("No Convergence.");
	}
	
	public static double optimize_newton_stabilized(Function f, double a, double b){
		double ap= 1e-6;
		double rp= 1e-4;
		int ns = 20;

		return optimize_newton_stabilized(f, a, b, ap, rp, ns);
	}
 
	//The optimized golden search functions
	 public static double optimize_golden_search(Function f, double a, double b) {
	      return optimize_golden_search(f, a, b, 0.0000001);
	    }
	    
	    public static double optimize_golden_search(Function f, double a, double b, double ap) {
	      return optimize_golden_search(f, a, b, ap, 0.00001);
	    }

	    public static double optimize_golden_search(Function f, double a, double b, double ap, double rp) {
	      return optimize_golden_search(f, a, b, ap, rp, 100);
	    }
	    
	    public static double optimize_golden_search(Function f, double a, double b, double ap, double rp, int ns) {
	      double tau = (Math.sqrt(5.)-1.)/2.;
	      double x1 = a+(1.-tau)*(b-a);
	      double x2 = a+tau*(b-a);
	      double fa = f.execute(a);
	      double f1 = f.execute(x1);
	      double f2 = f.execute(x2);
	      double fb = f.execute(b);
	      for (int k = 0; k < ns; k++) {
	        if (f1 > f2) {
	          a = x1;
	          fa = f1;
	          x1 = x2;
	          f1 = f2;
	          x2 = a+tau*(b-a);
	          f2 = f.execute(x2);
	        }
	        else {
	          b = x2;
	          fb = f2;
	          x2 = x1;
	          f2 = f1;
	          x1 = a+(1.-tau)*(b-a);
	          f1 = f.execute(x1);
	        }
	        if (k > 2 && norm(b-a) > Math.max(ap, norm(b)*rp))
	          return b;
	      }
	      throw new RuntimeException("no convergence");
	    }	
	
	  public static void main(String[] args){

		  //This function is for testing 
		  
		  Function f = new Function(){
			    public double execute(double a) {
				      return a*a-5.0*a; 
				    }		  
		  };
		  
		  Matrix testMatrix = new Matrix(6, 6);
	       double testVal = 45.4543243;
	       //System.out.println(condition_number(f, 1)); //0.74999...
	        Matrix A = new Matrix(2,2);
	    	A.setItem(0, 0, 1);
	    	A.setItem(0, 1, 2);
	    	A.setItem(1, 0, 3);
	    	A.setItem(1, 1, 4);
	        System.out.println(condition_number(A)); //21.0
	    	System.out.println(exp(A)); //[[51.96..., 74.73...], [112.10..., 164.07...]]
	    	Matrix B = new Matrix(3,3);
	    	B.setItem(0, 0, 4);
	    	B.setItem(0, 1, 2);
	    	B.setItem(0, 2, 1);
	    	B.setItem(1, 0, 2);
	    	B.setItem(1, 1, 9);
	    	B.setItem(1, 2, 3);
	    	B.setItem(2, 0, 1);
	    	B.setItem(2, 1, 3);
	    	B.setItem(2, 2, 16);
	        Matrix L = Cholesky(B);
	        Matrix lt = Matrix.transpose(L);
	        Matrix mult = Matrix.mult(L, lt);
	        Matrix result = Matrix.subtract(lt, mult);

	        System.out.println(is_almost_zero(result)); //True
		  
		  
		  f = new Function(){
			    public double execute(double a) {
				      return (a-2.)*(a-5.)/10;
				    }		  
		  };
		  
		  System.out.println(solve_fixed_point(f, 1.0, 0));
		  
		  f = new Function() {
			    public double execute(double a) {
			      return (a-2.)*(a-5.);
			    }
			  };
		  		  
		   System.out.println(solve_bisection(f,1.0,3.0));  
		   System.out.println(solve_newton(f, 1.0)); 
		   System.out.println(solve_secant(f, 1.0)); 
       	   System.out.println(solve_newton_stabilized(f, 1.0, 3.0));
       	   System.out.println(optimize_bisection(f, 2.0, 5.0));
       	   System.out.println(optimize_newton(f,3.0));
       	   System.out.println(optimize_secant(f,3.0));
       	   System.out.println(optimize_newton_stabilized(f, 2.0, 5.0));
       	   System.out.println(optimize_golden_search(f,2.0,5.0));
       	   
           MyFunctionImplOne functionOne = new MyFunctionImplOne();
           MyFunctionImplTwo functionTwo = new MyFunctionImplTwo();

//           MyFunctionInterface[] listOfFunctions = {functionOne, functionTwo };
//           double[][] listOfPoints = {{20,12,1},{13,14,1},{50,82,1}};
//           Matrix coef = FitLeastSquares(listOfPoints, listOfFunctions);
//           System.out.println("coef result is...------------:");
//           coef.printMatrix();
           
	    }
}


