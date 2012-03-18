// CSC431 Final Project 
// Members: Tony Abraham, Heidi Cochran, Larcen Razi, Mike Storm
//
//This class provides a means of finding the derivative for a user-supplied class. 
//It has the derivative and second derivative functions supplied by Numeric.py
public class Derivative implements Function {
    protected Function f;
    protected double h;

    public Derivative(Function f) {
      this.f = f;
      this.h = 0.0000001;
    }

    public Derivative(Function f, double h) {
      this.f = f;
      this.h = h;
    }

    public double execute(double a) {
      return (f.execute(a+h)-f.execute(a-h))/2/h;
    }
    
    public double executeDD(double a){
    	return (f.execute(a+h)-2.0*f.execute(a)+f.execute(a-h))/(h*h);
    }
}
