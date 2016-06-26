package quorumSensing;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Random;

public class Basal {
	private double[] B = {
		0.00014,
		0.01,
		0.0128,
		0.001,
		1.5e-4,
		4.5,
		4e-3,
		6.0e-3,
		0.33,
		3.3e-6,
		1.6e-2,
		1e-4,
		1.17e-4,
		1.4e-2,
		6.0e-3,
		13.2,
		0.4
	};
	double maxA = 200;
	double maxP = 500;
	private int N = B.length;
    private double[] 
    BS={0.000136684747792146,
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
    };
    
    private double[] Bhigh = {
    		1.60948513585537E-4,
    		0.00824384217096337,
    		0.014715292670677753,
    		8.243842170963541E-4,
    		1.6112667271061965E-4,
    		4.12961264699604,
    		0.004322200963413728,
    		0.00494630530257811,
    		0.3785760192443947,
    		3.2999952208858548E-6,
    		0.01835520093306161,
    		8.277277460529322E-5,
    		1.254833772246951E-4,
    		0.015114858269182104,
    		0.004966366476317545,
    		12.077339775230532,
    		0.3297543856891649
    };
    //parameters used in this network*/
	/*    	double[] b={0.9053749590054542,1.0920412475199646,0.9197116626070527,
    	   1.0857409599284473,0.9966876982013599,1.0582875635478466,0.9140012511190845,
    	   1.0275544131511434,1.005819792187465,1.0375555196397994,0.9384849339627492,
    	   1.0865349544071143,0.9027890233743637,0.9124597704671469,1.0303501229416336,
    	   0.9944017050538064,1.0973286676043663};*///this set of parameters don't have switch behavior
    private double k_RA,d_P,k_r,d_R,r_0,K_r,V_r,d_r,V_A,d_A,k_i,d_I,i_0,V_i,d_i,K_i,sgm;
    Basal(){
    	assignParameters();
    }
    Basal(double[] B){
    	this.B = B;
    }
    private double Ae = 23.561563490726638;

    public static void main(String[] arg) throws FileNotFoundException{
    	Basal b = new Basal();
    	//System.out.println( b.f());
    	b.run();
    }
    public void run() throws FileNotFoundException{
    	//sample();
    	//plot fp.png
    	
    	double[] b={0.9053749590054542,1.0920412475199646,0.9197116626070527,
    	    	   1.0857409599284473,0.9966876982013599,1.0582875635478466,0.9140012511190845,
    	    	   1.0275544131511434,1.005819792187465,1.0375555196397994,0.9384849339627492,
    	    	   1.0865349544071143,0.9027890233743637,0.9124597704671469,1.0303501229416336,
    	    	   0.9944017050538064,1.0973286676043663};
    	double[] B2 = new double[N];
    	for(int i=0;i<N;i++)
    		B2[i] = b[i]*B[i];
    	print(B2);
    	//assignParameters(B2);
    	//print(maxCurAe());
    	//print(lxdso());
    	/*double[] 
    		    B1={0.000136684747792146,
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
    		    };
    	assignParameters(B1);
    	print(AeThreshold());
    	B1[0] = B1[0]*1.01;
    	assignParameters(B1);
    	print(AeThreshold());
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = B1[i]/B[i];
    	}
    	double[] g = gradient(b);
    	print1(g);*/
    	/*sample("samples1.txt",0.01,10000);
    	sample("samples2.txt",0.02,10000);
    	sample("samples3.txt",0.03,10000);
    	sample("samples4.txt",0.04,10000);
    	sample("samples5.txt",0.05,10000);
    	sample("samples6.txt",0.06,10000);
    	sample("samples7.txt",0.07,10000);
    	sample("samples8.txt",0.08,10000);
    	sample("samples9.txt",0.09,10000);
    	sample("samples10.txt",0.1,10000);*/
    	//normalSample("normalSample.txt",0.205,100000);
    	//print(AeThreshold());
    	//double[] A={5,15,30};
    	//plotfP(B,A);
    	//print(AeThreshold());
    	/*for(double p=0;p<30;p=p+0.1){
    		Ae = 10;
    		double f10 = fP(p);
    		Ae = 50;
    		double f15 = fP(p);
    		Ae = 100;
    		double f30 = fP(p);
    		print(p+","+f10+","+f15+","+f30);
    	}*/
    	//plotAe();
    	/*double[] b3 = {0.8918911445507934,0.8862346818298884,1.6852800621207324
    			,0.8404437297327034,1.0395578269897432,1.3531874877338934,
    			1.4039010432289332,0.7303463320627843,1.2476386661855678,1.001115111497602,
    			1.3443218878507148,1.341018053850371,0.9611265643266024,
    			0.874026874871655,1.0459633473926349,0.8347142287354439,0.8568066184994949};
    	double[] B3 = new double[N];
    	for(int i=0;i<N;i++){
    		B3[i] = b3[i] * B[i];
    	}
    	assignParameters(B3);
    	//plotAe();
    	double[] A = {5};
    	//plotfP(B3,A);
    	//print(bottom(0));
    	plotAe();*/
    	//set3();
    }
    public void set3(){
    	double A = 0;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double bt = bottom(A,b);
    	double dr = 1e-3;
    	while(bt < 2){
    		double[] g = bottomGradient(A,b);
    		for(int i=0;i<N;i++)
    			b[i] = b[i] + g[i]*dr;
    		bt = bottom(A,b);
    	}
    	print(bt);
    	print(b);
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B[i];
    	print(BB);
    }
    public void plotfP(double[] B1,double[] A){
    	assignParameters(B1);
    	for(double p=0;p<0.5;p=p+0.001){
    		String s = p+"";
    		for(double Ae:A){
    			this.Ae=Ae;
    			s = s+"\t"+fP(p);
    		}
    		print(s);
    	}
    }
    private void plotAe(){
    	/* double[] B1={0.000136684747792146,
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
    		    };
    	assignParameters(B1);*/
    	for(double Ae=0.0;Ae<=20.1;Ae=Ae+0.1){
    		LinkedList<Double> sl = solveStationary(Ae);
    		String s = Ae + "";
    		for(double so:sl){
    			s += "\t" + so;
    		}
    		print(s);
    	}
    }
    private void sample(String file,double range,int s) throws FileNotFoundException{
    	FileOutputStream fi = new FileOutputStream(file);
		PrintWriter out = new PrintWriter(fi,true);
    	for(int i=1;i<s;i++){
    		double Ai = oneSample(range);
    		out.println(Ai);
        	if(i%100==0)
        		System.out.println(i);
    	}
    	out.close();
    }
    private double oneSample(double range){
    	Random rd = new Random();
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1 - range + 2*range*rd.nextDouble();
    	}
    	/*String bs = "";
    	for(double bi:b)
    		bs +=  bi+"\t";
    	System.out.println(bs);*/
    	return AeThreshold(b);
    }
    private void normalSample(String file,double std,int s) throws FileNotFoundException{
    	FileOutputStream fi = new FileOutputStream(file);
		PrintWriter out = new PrintWriter(fi,true);
    	for(int i=1;i<s;i++){
    		double Ai = oneNormalSample(std);
    		out.println(Ai);
        	//if(i%100==0)
        		//System.out.println(i);
    	}
    	out.close();
    }
    private double oneNormalSample(double std){
    	Random rd = new Random();
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = rd.nextGaussian() * std + 1;
    		while(b[i] <=0 ){
    			b[i] = rd.nextGaussian() * std + 1;
    		}
    	}
    	print(b);
    	return AeThreshold(b);
    }
    private double fP(double P){
        double r = r_0/d_r + (V_r/d_r)*P/(K_r+P);
        double R = k_r*r/d_R;
        double i = i_0/d_i + (V_i/d_i)*P/(K_i+P);
        double I = k_i*i/d_I;
        double A = d_P*P/(k_RA*R);
        return V_A*I - d_A*A + sgm*(Ae-A);
    }
    private double dfP(double P,double dP){//for test purpose, test if dP is enough
        double df = (fP(P+dP)-fP(P))/dP;
        return df;
    }
    private double dfP(double P){
        double dP = 1e-4;
        double df = (fP(P+dP)-fP(P))/dP;
        return df;
    }
    private double ddfP(double P){
        double dP = 1e-4;
        double dfp = (dfP(P+dP)-dfP(P))/dP;
        return dfp;
    }
    
    private Double bottom(double Ae){//p satisfies fp'(p)=0 && fp''(p)>0, this function
    	this.Ae = Ae;
    	double p1 = 1e-5;
    	while(dfP(p1)>=0)
    		p1 = 0.9*p1;
    	double p2 = 1.01*p1;
    	while(dfP(p2)<0){
    		p2 = 1.01*p2;
    		if(p2>maxP)
    			return null;
    	}
    	double ERR = 1e-8;
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
    	//System.out.println("p="+p);
    	return fP(p);
    }
    private Double bottom(double Ae,double[] b){
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = b[i]*B[i];
    	}
    	assignParameters(BB);
    	return bottom(Ae);
    }
    private double[] bottomGradient(double A,double[] b){
    	double[] g = new double[N];
    	double db=1e-3;
    	double bt = bottom(A,b);
    	for(int i=0;i<N;i++){
        	double[] bi = (double[]) b.clone();
        	bi[i] = b[i] + db;
        	g[i] = (bottom(A,bi)-bt)/db;
    	}
    	return g;
    }

    private double peak(double Ae){
    	this.Ae = Ae;
    	double p1 = 300;
    	double p2 = p1;
    	while(dfP(p2)<0)
    		p2 = 0.9*p2;
    	double ERR = 1e-8;
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
    public double AeThreshold(){//return the Ae that makes f(Ae)=0
    	double A1 = 0.1;
    	Double b1 = bottom(A1);
    	if(b1==null){
    		print("123");
    		return maxCurAe();
    	}
    	while(bottom(A1)>0)
    		A1 = 0.9*A1;
    	double A2 = 10;
    	while(bottom(A2)<0)
    		A2 = 1.01*A2;
    	double A = (A1+A2)/2;
    	double fA = bottom(A);
    	double ERR = 1e-5;
    	while(Math.abs(fA)>ERR){
    		if(fA<0)
    			A1=A;
    		else
    			A2=A;
    		A = (A1+A2)/2;
    		fA = bottom(A);
    	}
    	return A;
    }
    public double AeThreshold(double[] b){
    	double[] B1 = new double[N];
    	for(int i=0;i<N;i++){
    		B1[i] = B[i]*b[i];
    	}
    	assignParameters(B1);
    	return AeThreshold();
    }
    private double[] gradient(double[] b){
    	double Ab = AeThreshold(b);
    	double db = 1e-3;
    	double[] g = new double[N];
    	for(int i=0;i<N;i++){
    		double[] bi = (double[]) b.clone();
    		bi[i] = b[i] + db;
    		g[i] = (AeThreshold(bi)-Ab)/Ab/db;
    	}
    	return g;
    }
    
    public LinkedList<Double> solveStationary(double Ae){
    	this.Ae = Ae;
        LinkedList<Double> solutions = new LinkedList<Double>();
        //high value, fP(h2) must be negative
        double h2 = (k_i*(i_0+V_i)/(d_I*d_i)+sgm*Ae)/(sgm*d_P*d_R*d_r/(k_RA*k_r*(r_0+V_r)))+1;
        //find all the solutions and store them in variable "solutions"
        double h1 = 0.9*h2;
        double f2 = fP(h2);
        double err = 1e-5;
        double s = h2;
        boolean hasMoreSolution = true;
        while(Math.abs(h2) > err){
            f2 = fP(h2);
            h1 = 0.9*h2;
            while(fP(h1)*f2 > 0){
                h1 = 0.9*h1;
                if(Math.abs(h1) < err){
                    hasMoreSolution = false;
                    break;
                }
            }
            if(hasMoreSolution == false){
                break;
            }else{
                s = bisection(h1,h2);
                solutions.add(s);
                h2 = 0.99*s;
            }
        }
        return solutions;
    }
    public double so(double Ae){
    	LinkedList<Double> dd = solveStationary(Ae);
    	return dd.get(0);
    }
    public double dso(double Ae){
    	double dA = 0.1;
    	double ds = (so(Ae+dA)-so(Ae))/dA;
    	return ds;
    }
    public double ddso(double Ae){
    	double dA = 0.1;
    	double dds = (dso(Ae+dA)-dso(Ae))/dA;
    	return dds;
    }
    double lxdso(){
    	double A1 = 0.1;
    	double A2 = A1*1.1;
    	while(ddso(A2)>=0)
    		A2 = A2*1.1;
    	double ERR = 1e-3;
    	double A = 0;
    	while(A2-A1>ERR){
    		A = (A1+A2)/2;
    		if(ddso(A)>0)
    			A1 = A;
    		else
    			A2 = A;
    	}
    	return A;
    }
    double curvature(double A){
    	double d = dso(A);
    	double dd = ddso(A);
    	return dd / Math.pow(1+d*d, 1.5);
    }
    double maxCurAe(){
    	double Am = -1;
    	double cAm = -1;
    	for(double A=0;A<=maxA;A++){
    		double cA = curvature(A);
    		if(cAm < cA){
    			Am = A;
    			cAm = cA;
    		}else if(cAm - cA > 0.01)
    			break;
    	}
    	return Am;
    }
    public double bisection(double x1,double x2){
        if(fP(x1)*fP(x2) >0){
            System.out.println("fP(x1,Ae)*(fP(x2,Ae) >0");
            return -1;
        }
        double ERR = 1e-7;
        double x = (x1+x2)/2;
        while (Math.abs(fP(x)) > ERR){
            if(fP(x)*fP(x1) >=0)
                x1 = x;
            else
                x2 = x;
            x = (x1+x2)/2;
        }
        return x;
    }
    private void trace(){
    	double[] b0 = new double[N];
    	for(int i=0;i<N;i++){
    		b0[i] = 1;
    	}
    	double dx = 0.001;
    	double A0 = AeThreshold(b0);
    	double Ab = A0;
    	double[] b = (double[]) b0.clone();
    	while(Ab <= 2*A0){
    		double[] g = gradient(b);
    		for(int i=0;i<N;i++){
    			b[i] = b[i] + g[i]*dx;
    		}
    		Ab = AeThreshold(b);
    	}
    	for(int i=0;i<N;i++){
    		System.out.println(b0[i]+","+b[i]+","+(Math.abs(b[i]-b0[i])/b0[i]));
    	}
    	System.out.println(A0+","+Ab);
    }
    public void sensitivity(double[] B1){
    	assignParameters(B1);
    	
    }
    public void test(){
    	//test bottom(Ae)
    	/*
    	Ae = 10;
    	double ba = bottom(Ae);
    	double P = 0.4831226169142312;
    	System.out.println(fP(P));
    	System.out.println(dfP(P));
    	System.out.println(ba);
    	double dp = 1e-3;
    	System.out.println((dfP(P+dp)-dfP(P))/dp);
    	*/
    	//test AeThreshold()
    	/*double range = 0.05;
    	Random rd = new Random();
    	double[] B1 = new double[N];
    	for(int i=0;i<N;i++){
    		B1[i] = (1 - range + 2*range*rd.nextDouble())*B[i];
    	}
    	assignParameters(B1);
    	double A = AeThreshold();
    	System.out.println(bottom(A*0.95)+","+bottom(A)+","+bottom(A*1.05));*/
    	//test solveStationary(Ae)
    	/*
    	 *     	//solve steady state
    	double[] Aes = {10,16.62,30};
    	for(double Ae:Aes){
    		LinkedList<Double> solutions = solveStationary(Ae);
    		String s = "";
    		for(double so:solutions){
    			s += so + ","+fP(so)+"\t";
    		}
    		print(s);
    		this.Ae = Ae;
    	}
    	 */
    }
    private void assignParameters(){
    	k_RA = B[0];
    	d_P = B[1]; 
    	k_r = B[2]; 
    	d_R = B[3]; 
    	r_0 = B[4]; 
    	K_r = B[5]; 
    	V_r = B[6]; 
    	d_r = B[7]; 
    	V_A = B[8]; 
    	d_A = B[9];
    	k_i = B[10];
    	d_I = B[11];
    	i_0 = B[12];
    	V_i = B[13];
    	d_i = B[14];
    	K_i = B[15];
    	sgm = B[16];
    }
    private void assignParameters(double[] B){
    	k_RA = B[0];
    	d_P = B[1]; 
    	k_r = B[2]; 
    	d_R = B[3]; 
    	r_0 = B[4]; 
    	K_r = B[5]; 
    	V_r = B[6]; 
    	d_r = B[7]; 
    	V_A = B[8]; 
    	d_A = B[9];
    	k_i = B[10];
    	d_I = B[11];
    	i_0 = B[12];
    	V_i = B[13];
    	d_i = B[14];
    	K_i = B[15];
    	sgm = B[16];
    }
    private void print(Object s){
    	System.out.println(s+"");
    }
    private void print(double[] d){
    	String s = "";
    	int L = d.length;
    	for(int i = 0; i < L;i++){
    		if(i!=L-1)
    			s = s+d[i]+",";
    		else
    			s = s+d[i]+":";
    	}
    	print(s);
    }
    private void print1(double[] d){
    	for(double d1:d)
    		print(d1);
    }
}
