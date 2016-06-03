package quorumSensing;

public class Basal {
    private double[] 
    B={0.000136684747792146,
    0.0102216510333074,
    0.0125089032299558,
    0.00102313950697868,
    0.000149798525769270,
    4.53717629423953,
    0.00391076449604863,
    0.00613883704187208,
    0.321845638538543,
    3.30000063093836e-06,
    0.0155942612460342,
    0.000102402928451709,
    0.000116934896039553,
    0.0136516614435384,
    0.00614417570710256,
    13.3984086347505,
    0.408357478338906
    };//parameters used in this network
    private double k_RA = B[0];
    private double d_P = B[1]; 
    private double k_r = B[2]; 
    private double d_R = B[3]; 
    private double r_0 = B[4]; 
    private double K_r = B[5]; 
    private double V_r = B[6]; 
    private double d_r = B[7]; 
    private double V_A = B[8]; 
    private double d_A = B[9];
    private double k_i = B[10];
    private double d_I = B[11];
    private double i_0 = B[12];
    private double V_i = B[13];
    private double d_i = B[14];
    private double K_i = B[15];
    public double sgm = B[16];
    double Ae = 23.561563490726638;
    public static void main(String[] arg){
    	Basal b = new Basal();
    	//System.out.println( b.f());
    	b.run();
    }
    public void run(){
    	for(double P=0;P<30;P=P+1)
    		System.out.println(P+"\t"+fP(P)+"\t"+0);
    	Ae = 0;
    	for(double P=0;P<30;P=P+1)
    		System.out.println(P+"\t"+fP(P)+"\t");    
    	Ae = 15;
    	for(double P=0;P<30;P=P+1)
    		System.out.println(P+"\t"+fP(P)+"\t");    
    	Ae = 30;
    	for(double P=0;P<30;P=P+1)
    		System.out.println(P+"\t"+fP(P)+"\t");    
    }
    private double fP(double P){
        double r = r_0/d_r + (V_r/d_r)*P/(K_r+P);
        double R = k_r*r/d_R;
        double i = i_0/d_i + (V_i/d_i)*P/(K_i+P);
        double I = k_i*i/d_I;
        double A = d_P*P/(k_RA*R);
        return V_A*I - d_A*A + sgm*(Ae-A);
    }
    private double dfP(double P){
        double dP = 1e-3;
        double df = (fP(P+dP)-fP(P))/dP;
        return df;
    }
    private double f(double Ae){//p satisfies fp'(p)=0 && fp''(p)>0, this function
    	//returns fp(p);
    	this.Ae = Ae;
    	double p1 = 1e-5;
    	while(dfP(p1)>=0)
    		p1 = 0.9*p1;
    	double p2 = 1.1*p1;
    	while(dfP(p2)<0)
    		p2 = 1.1*p2;
    	double ERR = 1e-5;
    	double p = (p1+p2)/2;
    	double dp = dfP(p);
    	while(Math.abs(dp)>ERR){
    		if(dp<0)
    			p1 = p;
    		else
    			p2 = p;
    		p = (p1+p2)/2;
    		dp = dfP(p);
    	}
    	return fP(p);
    }
    public double f(){//return the Ae that makes f(Ae)=0
    	double A1 = 0.1;
    	while(f(A1)>0)
    		A1 = 0.9*A1;
    	double A2 = 10;
    	while(f(A2)<0)
    		A2 = 1.1*A2;
    	double A = (A1+A2)/2;
    	double fA = f(A);
    	double ERR = 1e-5;
    	while(Math.abs(fA)>ERR){
    		if(fA<0)
    			A1=A;
    		else
    			A2=A;
    		A = (A1+A2)/2;
    		fA = f(A);
    	}
    	return A;
    }
}
