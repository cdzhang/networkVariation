package quorumSensing;

public class Components {
    public double Ae;
    public double P;
    public double r;
    public double R;
    public double i;
    public double I;
    public double A;
    Components(){
    }
    Components(double P, double r, double R, double i, double I, double A){
        this.P = P;
        this.r = r;
        this.R = R;
        this.i = i;
        this.I = I;
        this.A = A;
    }
    public void setValues(double P, double r, double R, double i, double I, double A){
        this.P = P;
        this.r = r;
        this.R = R;
        this.i = i;
        this.I = I;
        this.A = A;
    }
    public void setAe(double Ae){
        this.Ae = Ae;
    }
   /* public void print(){
        System.out.format("%7s,%7s,%7s,%7s,%7s,%7s,%7s","Ae","P","r","R","i","I","A");
        System.out.format("\n%5.5f,%5.5f,%5.5f,%5.5f,%5.5f,%5.5f,%5.5f\n",Ae,P,r,R,i,I,A);
    }*/
    public double[] getX0(){
        double[] X0 = {P,r,R,i,I,A};
        return X0;
    }
}
