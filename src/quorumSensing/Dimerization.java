package quorumSensing;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.Random;

public class Dimerization {
	private double[] B1 = {
			4.326E-5
			,0.00502
			,0.01785
			,4.698E-4
			,1.872E-4
			,1.086
			,0.03043
			,4.302E-4
			,0.001684
			,3.0E-6
			,0.01632
			,3.72E-4
			,6.291E-5
			,0.1546
			,2.548E-4
			,27.87
			,0.3988
			,4.044E-6
			,0.01567
	};
	private double[] B2 ={//base line
			4.668E-5
			,0.004623
			,0.01926
			,4.327E-4
			,1.875E-4
			,1.054
			,0.03283
			,3.962E-4
			,0.001817
			,3.0E-6
			,0.01761
			,3.426E-4
			,6.301E-5
			,0.01668
			,2.347E-4
			,25.8
			,0.3673
			,4.363E-6
			,0.01443
	};
	private double[] B3 = {
			6.273E-5
			,0.002761
			,0.02588
			,2.584E-4
			,2.714E-4
			,0.5973
			,0.04412
			,2.366E-4
			,0.00214
			,3.0E-6
			,0.02074
			,2.311E-4
			,7.951E-5
			,0.0156
			,1.583E-4
			,27.624
			,0.2193
			,5.864E-6
			,0.008619
			};
	private double[] B4={
			6.160247505039312E-5
			,9.22167484175848E-4
			,0.025414985781329582
			,7.837730154909747E-5
			,2.5757099613948507E-4
			,0.9725494874521773
			,0.03462621495095533
			,1.24664309372095E-6
			,4.031796374038114E-4
			,3.180925829880949E-6
			,0.0011843036212531833
			,5.447569045059144E-4
			,6.254014005294303E-5
			,0.0027199148661958394
			,3.6417397665888996E-4
			,37.571795672012065
			,0.0485655883641746
			,2.9124754200183724E-6
			,0.019854682278356926
	};
	double maxA = 200;
	double maxP = 10000;
	double pdfm = 0;
	private int N = B1.length;
	private double k_RA,d_P,k_r,d_R,r_0,K_r,V_r,d_r,V_A,d_A,k_i,d_I,i_0,V_i,d_i,K_i,sgm,k_M,d_M;
	double lowMax = 12;
	double highMin = 40;
	double highMax = 5000;
	double peak0 = -5;
	double[] Bfit = B1.clone();
    Dimerization(){
    	assignParameters();
    }
    Dimerization(double[] B){
    	this.B1 = B;
    }
    private double Ae = 0;
    public static void main(String[] arg) throws FileNotFoundException{
    	Dimerization b = new Dimerization();
    	b.run();
    }
    public void run() throws FileNotFoundException{
/*    	double[] b = ones(N);
    	B1=B2;
    	assignParameters(B1);
    	plotAe();
    	Ae=0;
    	print(allComponents());
    	Ae=30;
    	print(allComponents());
    	double[] B = set2();
    	print1(B);
    	B1=B;
    	assignParameters(B1);
    	plotAe();*/
    	//print1(gradient(b));
    	//double[] B = set1();
    	//assignParameters(B2);
    	//double[] A = {0,15,30};
    	//plotfP(B2,A);
    	//plotdfP(B2,A);
    	Ae=30;
    	print(allComponents());
    }
    void restrictValues(double[] B){
    	if(B[16]>0.45) B[16] = 0.45;
    	if(B[16]<0.05) B[16] = 0.05;
    	if(B[10]<1.6e-2) B[10] = 1.6e-2;
    	if(B[2]<1.5e-2) B[2] = 1.6e-2;
    }
    void restrict(double a,double[] B, double[] B0){
    	for(int i=0;i<N;i++){
    		double ai = (B[i]-B0[i])/B0[i];
    		if(ai>a)
    			B[i] = B0[i]*(1+a);
    		else if(ai<-a)
    			B[i] = B0[i]*(1-a);
    	}
    }
  /*  public boolean fit(int K){
    	double[] BB = Bfit.clone();
    	boolean needFit = false;
    	for(int i=0;i<K;i++){
    		if(!isASimilar(BB)){
    			BB = AediffNextPoint(BB);
    			if(!checkB(BB)){
    				print("1");
    			}
    			restrictValues(BB);
    			needFit = true;
    		}
    	}
    	for(int i=0;i<K;i++){
    		if(!isLow(BB)){
    			Ae = 0;
    			BB = nextPoint(BB,1e-8,lowMax);
    			if(!checkB(BB)){
    				print("2");
    			}
    			restrictValues(BB);
    	    	assignParameters(BB);
    	    	double min = minConcentration();
    	    	print(min);
    	    	needFit = true;
    		}
    	}
    	for(int i=0;i<K;i++){
    		if(!isHigh(BB)){
    			Ae = 40;
    			BB = nextPoint(BB,highMin,highMax);
    			if(!checkB(BB)){
    				print("3");
    			}
    			restrictValues(BB);
    	    	assignParameters(BB);
    	    	double min = minConcentration();
    	    	double max = maxConcentration();
    	    	print(min+","+max);
    	    	needFit = true;
    	    	
    		}
    	}
    	for(int i=0;i<K;i++){
    		if(!isAe0(BB)){
    			BB = peakNextPoint(BB);
    			restrictValues(BB);
    			assignParameters(BB);
    			if(!checkB(BB)){
    				print("4");
    			}
    			double p = peak(0);
    			print(p);
    			needFit = true;
    		}
    	}
    	print1(BB);
    	Bfit = BB;
    	return needFit;
    }*/
    double[] fit(double[] B){
    	double tB = totalScore(B);
    	double ERR = 1e-3;
    	int i=0;
    	while(tB > ERR){
    		print("tB="+tB);
    		B = nextPoint(B);
    		tB = totalScore(B);
    		i++;
    	}
    	return B;
    }
    double[] nextPoint(double[] B){
    	double tB = totalScore(B);
    	double db = 1e-3;
    	double delta=5e-2;
    	double[] g = new double[N];
    	for(int i=0;i<N;i++){
    		double[] Bi = B.clone();
    		Bi[i] = B[i]*(1+db);
    		double tBi = totalScore(Bi);
    		g[i] = (tBi-tB)/db;
    		double gi = Math.abs(g[i]);
    		if(gi>1 && delta > 1e-3/gi)
    			delta = 1e-3/gi;
    	}
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = B[i]*(1-delta*g[i]);
    	}
    	//restrictValues(BB);
    	return BB;
    }
    double[] AediffNextPoint(double[] B){
    	double adb = Aediff(B);
    	print(adb);
    	double db = 1e-3;
    	double delta = 1e-3;
    	double[] g = new double[N];
    	for(int i=N-3;i<N;i++){
    		double[] Bi = (double[]) B.clone();
    		Bi[i] = B[i]*(1+db);
    		double adbi = Aediff(Bi);
    		double gi = (adbi-adb)/db;
    		g[i] = gi;
    		double di = 1e-3/Math.abs(gi);
    		if(delta>di)
    			delta = di;
    	}
    	double[] BB = new double[N];
    	if(adb>0)
    		delta = -1*delta;
    	for(int i=0;i<N;i++)
    		BB[i] = B[i]*(1+delta*g[i]);
    	return BB;
    }
    double Aediff(double[] B){
    	assignParameters(B);
    	double[] ac = allComponents();
    	double A = ac[5];
    	return (A-Ae)/Ae;
    }
    //criteria
    double isASimilar(double[] B){
    	Ae = 40;
    	double ab = Aediff(B);
    	if(Math.abs(ab)>0.2){
    		print("A and Ae differ " + ab);
    		return (Math.abs(ab)-0.2)/0.2;
    	}
    	return 0;
    }
    boolean isThreshold(){
    	return true;
    }
    double isLow(double[] B){
    	double Ae1 = Ae;
    	Ae = 0;
    	assignParameters(B);
    	double max = maxConcentration();
    	Ae = Ae1;
    	if(max > lowMax){
    		print("When Ae = 0, max="+max);
    		return (max-lowMax)/lowMax;
    	}
    	return 0;
    }
    double isHigh(double[] B){
    	double Ae1 = Ae;
    	Ae = 40;
    	assignParameters(B);
    	double min = minConcentration();
    	double max = maxConcentration();
    	Ae = Ae1;
    	double x1=0,x2=0;
    	if(min < highMin){
    		print("min="+min);
    		x1 = (highMin-min)/highMin;
    	}
    	if(max>highMax){
    		print("max="+max);
    		x2 = (max-highMax)/highMax;
    	}
    	return x1 + x2;
    }
    double isAe0(double[] B){
    	assignParameters(B);
    	double p = peak(0);
    	if(p>peak0){
    		return (peak0-p)/peak0;
    	}
    	if(p<peak0)
    		return (p-peak0)/peak0;
    	return 0;
    }
    double isAe(double[] B){
    	assignParameters(B);
    	double AT = AeThreshold();
    	if(AT > 21)
    		return (AT-21)/21;
    	if(AT < 15)
    		return (15-AT)/15;
    	return 0;
    }
    public double totalScore(double[] B){
    	return isASimilar(B)+isLow(B)+isHigh(B)+isAe0(B);
    }
    //fit
 
    boolean checkB(double[] B){
    	for(double bi:B){
    		if(bi<=0)
    			return false;
    	}
    	return true;
    }


    double[] fitAe(double[] B){
    	while(Math.abs(Aediff(B))>0.1){
    		B = AediffNextPoint(B);
    	}
    	return B;
    }
    double[] allComponents(double P){
    	if(Math.abs(fP(P))>1e-5) return null;
    	double M = k_M/d_M*P*P;
        double r = r_0/d_r + (V_r/d_r)*M/(K_r+M);
        double R = k_r*r/d_R;
        double i = i_0/d_i + (V_i/d_i)*M/(K_i+M);
        double I = k_i*i/d_I;
        double A = d_P*P/(k_RA*R);
        double[] com = {M,r,R,i,I,A,P};
        return com;
    }
    double[] allComponents(){
    	LinkedList<Double> so =  solveStationary(Ae);
    	
    	double P = so.getFirst();
    	if(Ae<5){
    		P = so.getLast();
    	}
    	return allComponents(P);
    }
    double minConcentration(){
    	double[] ap = allComponents();
    	double a=-1;
    	int aN = ap.length;
    	for(int i=0;i<aN;i++){
    		if(a<0||ap[i]<a){
    			a=ap[i];
    		}
    	}
    	return a;
    }
    double maxConcentration(){
    	double[] ap = allComponents();
    	double a=-1;
    	int aN = ap.length;
    	for(int i=0;i<aN;i++){
    		if(a<0||ap[i]>a){
    			a=ap[i];
    		}
    	}
    	return a;
    }
    double[] peakNextPoint(double[] B){
    	assignParameters(B);
    	double pB = peak(0);
    	double db = 1e-3;
    	double delta = 1e-2;
    	double[] g = new double[N];
    	for(int i=0;i<N;i++){
    		double[] Bi = (double[]) B.clone();
    		Bi[i] = B[i]*(1+db);
    		assignParameters(Bi);
    		double pBi = peak(0);
    		double gi = (pBi - pB)/db;
    		g[i] = gi;
    		double di = 1e-3/Math.abs(gi);
    		if(delta > di)
    			delta = di;
    	}
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = B[i]*(1-g[i]*delta);
    	}
    	return BB;
    }
    double[] nextPointMin(double[] B){
    	assignParameters(B);
    	double mB = minConcentration();
    	double db = 1e-3;
    	double delta = 1e-2;
    	double[] newB = new double[N];
    	double[] g = new double[N];
    	for(int i=0;i<N;i++){
    		double[] Bi = (double[]) B.clone();
    		Bi[i] = B[i]*(1+db);
    		assignParameters(Bi);
    		double mBi = minConcentration();
    		double gi = (mBi-mB)/db;
    		g[i] = gi;
    		double di = 1e-3/Math.abs(gi);
    		if(delta>di)
    			delta=di;
    	}
    	for(int i=0;i<N;i++){
    		newB[i] = (1+delta*g[i])*B[i];
    	}
    	return newB;
    }
    double[] nextPointMax(double[] B){
    	assignParameters(B);
    	double mB = maxConcentration();
    	double db = 1e-3;
    	double delta = 1e-2;
    	double[] newB = new double[N];
    	double[] g = new double[N];
    	for(int i=0;i<N;i++){
    		double[] Bi = (double[]) B.clone();
    		Bi[i] = B[i]*(1+db);
    		assignParameters(Bi);
    		double mBi = maxConcentration();
    		double gi = (mBi-mB)/db;
    		//print("gi="+gi);
    		g[i] = gi;
    		double di = 1e-3/Math.abs(gi);
    		if(delta>di)
    			delta=di;
    	}
    	for(int i=0;i<N;i++){
    		newB[i] = (1-delta*g[i])*B[i];
    	}
    	return newB;
    }
    double[] nextPoint(double[] B,double minc,double maxc){
    	assignParameters(B);
    	double minB = minConcentration();
    	double maxB = maxConcentration();
    	double[] BB;
		if(minB/minc < maxc/maxB){
			BB = nextPointMin(B);
		}else
			BB = nextPointMax(B);
		return BB;
    }

    public double[] set5(double minc, double maxc){
    	double[] b = new double[N];
    	double[] BB = (double[]) B1.clone();
    	double min = minConcentration();
    	double max = maxConcentration();
    	while( max > maxc || min < minc){
    		BB = nextPoint(BB,minc,maxc);
    		min = minConcentration();
    		max = maxConcentration();
    		print(min+","+max);
    	}
    	for(int i=0;i<N;i++)
    		b[i] = BB[i]/B1[i];
    	print(b);
    	print(BB);
    	return BB;
    }
    public double[] set5(double A1,double minc1,double maxc1,double A2,double minc2,double maxc2){
    	double[] BB = (double[]) B1.clone();
    	Ae = A1;
    	double min1 = minConcentration();
    	double max1 = maxConcentration();
    	Ae = A2;
    	double min2 = minConcentration();
    	double max2 = maxConcentration();
    	int i=0;
    	while(min1<minc1||max1>maxc1||min2<minc2||max2>maxc2){
    		if(i%2==0){
    			Ae = A1;
    			BB = nextPoint(BB,minc1,maxc1);
    		}else{
    			Ae = A2;
    			BB = nextPoint(BB,minc2,maxc2);
    		}
        	Ae = A1;
        	min1 = minConcentration();
        	max1 = maxConcentration();
        	Ae = A2;
        	min2 = minConcentration();
        	max2 = maxConcentration();
        	print(min1+","+max1+","+min2+","+max2);
        	i++;
        	if(i%100==0){
        		print(i);
        		print(BB);
        	}
    	}
    	return BB;
    }


    public void run2(){
/*    	//set4();
    	
    	plotfP(B4,A);
    	double[] b = new double[N];
    	for(int i=0;i<N;i++)
    		b[i] = B4[i]/B1[i];
    	//print(peakDfp(10,b));
    	
*/    
/*    	assignParameters(B4);
    	print(maxCurAe());
    	print(minCurAe());
		print(AeThreshold());
		//plotAe();
		for(double A=0;A<25;A=A+0.1){
			print(A+"\t"+dso(A)+"\t"+ddso(A)+"\t"+curvature(A));
		}*/
    	double[] A = {0,15,30};
    	//plotdfP(B4,A);
    	assignParameters(B1);
    	peakDfp(0);
    	
    	print(pdfm);
    	for(double A1=10;A1<20;A1=A1+0.1){
    		Ae=A1;
    		print(A1+"\t"+fP(pdfm));
    	}
    	}
    private void plotAe(){
    	for(double Ae=0.0;Ae<=40.1;Ae=Ae+0.1){
    		LinkedList<Double> sl = solveStationary(Ae);
    		String s = Ae + "";
    		for(double so:sl){
    			s += "\t" + so;
    		}
    		print(s);
    	}
    }
    public void plotfP(double[] B,double[] A){
    	assignParameters(B);
    	for(double p=0;p<5000;p=p+1){
    		String s = p+"";
    		for(double Ae:A){
    			this.Ae=Ae;
    			s = s+"\t"+fP(p);
    			//print(fP(21652.9));
    		}
    		print(s);
    	}
    }
    public void plotdfP(double[] B,double[] A){
    	assignParameters(B);
    	for(double p=0;p<500;p=p+1){
    		String s = p+"";
    		for(double Ae:A){
    			this.Ae=Ae;
    			s = s+"\t"+dfP(p);
    		}
    		print(s);
    	}
    }
    private double fP(double P){
    	double M = k_M/d_M*P*P;
        double r = r_0/d_r + (V_r/d_r)*M/(K_r+M);
        double R = k_r*r/d_R;
        double i = i_0/d_i + (V_i/d_i)*M/(K_i+M);
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
    		BB[i] = b[i]*B1[i];
    	}
    	assignParameters(BB);
    	return bottom(Ae);
    }
    private double peak(double Ae){
    	this.Ae = Ae;
    	double p1 = 50000;
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
    private double peakDfp(double Ae){
    	this.Ae = Ae;
    	double p1 = 0.1;
    	while(ddfP(p1)<=0)
    		p1=p1*0.9;
    	double p2 = p1;
    	while(ddfP(p2)>=0)
    		p2 = 1.01*p2;
    	double ERR = 1e-6;
    	double p = (p1+p2)/2;
    	double dp = ddfP(p);
    	while(Math.abs(dp)>ERR){
    		if(dp<0)
    			p2 = p;
    		else
    			p1 = p;
    		p = (p1+p2)/2;
    		dp = ddfP(p);
    	}
    	pdfm = p;
    	//print("p="+p);
    	return dfP(p);
    }
    private double peakDfp(double Ae, double[] b){
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = b[i]*B1[i];
    	}
    	assignParameters(BB);
    	return peakDfp(Ae);
    }
    private double peak(double Ae, double[] b){
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = b[i]*B1[i];
    	}
    	assignParameters(BB);
    	return peak(Ae);
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
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = B1[i]*b[i];
    	}
    	assignParameters(BB);
    	return AeThreshold();
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
    private double[] peakGradient(double A,double[] b){
    	double[] g = new double[N];
    	double db=1e-3;
    	double bp = peak(A,b);
    	for(int i=0;i<N;i++){
        	double[] bi = (double[]) b.clone();
        	bi[i] = b[i] + db;
        	g[i] = (peak(A,bi)-bp)/db;
    	}
    	return g;
    }
    private double[] peakDfpGradient(double A,double[] b){
    	double[] g = new double[N];
    	double db=1e-3;
    	double bp = peakDfp(A,b);
    	for(int i=0;i<N;i++){
        	double[] bi = (double[]) b.clone();
        	bi[i] = b[i] + db;
        	g[i] = (peakDfp(A,bi)-bp)/db;
    	}
    	return g;
    }
   
    private double[] gradient(double[] b){
    	double Ab = AeThreshold(b);
    	//print(Ab);
    	double db = 1e-3;
    	double[] g = new double[N];
    	for(int i=0;i<N;i++){
    		double[] bi = (double[]) b.clone();
    		bi[i] = b[i] + db;
    		double gi = AeThreshold(bi);
    		g[i] = (gi-Ab)/Ab/db;
    	}
    	return g;
    }
    public double[] set3(){
    	double A = 0;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double[] b0 = ones(N);
    	double a = 0.5;
    	double bt = bottom(A,b);
    	double dr = 1e-1;
    	while(bt < 0.001){
    		double[] g = bottomGradient(A,b);
    		for(int i=0;i<N;i++){
    			if(Math.abs(b[i]-1)>a*0.99)
    				continue;
    			double gi = Math.abs(g[i]);
    			double di = 1e-2/gi;
    			if(dr>di) dr=di;
    		}
    		for(int i=0;i<N;i++){
    			b[i] = b[i] + g[i]*dr;
    		}
    		restrict(a,b,b0);
    		bt = bottom(A,b);
    		print("dr="+dr+",bt="+bt);
    	}
    	print(bt);
    	print(b);
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print(BB);
    	return BB;
    }
    public double[] set1(){
    	double A = 0;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double bt = peak(A,b);
    	double dr = 1e-1;
    	while(bt < -2){
    		double[] g = peakGradient(A,b);
    		for(double gi:g)
    			if(dr>1e-3/Math.abs(gi))
    				dr = 1e-3/Math.abs(gi);
    		for(int i=0;i<N;i++)
    			b[i] = b[i] + g[i]*dr;
    		bt = peak(A,b);
    		print(bt);
    	}
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print(BB);
    	return BB;
    }   
    public double[] set2(){
    	double A = 0;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double[] b0 = ones(N);
    	double bt = peak(A,b);
    	double dr = 0.01;
    	double a = 0.079;
    	while(bt <0.4){
    		double[] g = peakGradient(A,b);
    		for(int i=0;i<N;i++){
    			if(Math.abs(b[i]-1)>a*0.99)
    				continue;
    			double gi = Math.abs(g[i]);
    			double di = 1e-2/gi;
    			if(dr>di) dr=di;
    		}
    		for(int i=0;i<N;i++){
    			b[i] = b[i] + g[i]*dr;
    		}
    		restrict(a,b,b0);
    		bt = peak(A,b);
    		print("bt="+bt);
    	}
    	print(bt);
    	print(b);
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print(BB);
    	return BB;
    }
    public void set4(){
    	double A = 15;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double bt = peakDfp(A,b);
    	double dr = 5e-2;
    	while(bt > 0){
    		double[] g = peakDfpGradient(A,b);
			double maxg = 0;
    		for(int i=0;i<N;i++){
    			double gi = Math.abs(g[i]);
    			if(maxg < gi){
    				maxg=gi;
    			}
    		}
    		dr = 0.001/maxg;
    		for(int i=0;i<N;i++)
    			b[i] = b[i] - g[i]*dr;
    		bt = peakDfp(A,b);
    		print("bt="+bt);
    	}
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print1(BB);
    }
    private boolean oneSample1(double range){
    	Random rd = new Random();
    	double[] b = new double[N];
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1 - range + 2*range*rd.nextDouble();
    		BB[i] = b[i]*B1[i];
    	}
    	assignParameters(BB);
    	double PeM = 500;
    	//patten= dfp < 0, dfp >0, dfp <0
    	double P1 = 0.1;
    	while(dfP(P1)>=0)
    		P1 = P1*0.9;
    	boolean[] pattern = {false,true,false};
    	int ip = 0;
    	for(double P=P1;P<=PeM;P=P+0.1){
    		/*boolean re = dfP(P)>=0;
    		if(re==pattern[ip])
    			continue;
    		ip++;
    		if(ip>=pattern.length){
    			print(BB);
    			return false;
    		}*/
    		if(dfP(P)>0)
    			return false;
    	}
    	/*if(ip<pattern.length-1){
    		print(BB);
    		return false;
    	}*/
    	print(BB);
    	return true;
    }
    public LinkedList<Double> solveStationary(double A){
    	this.Ae = A;
        LinkedList<Double> solutions = new LinkedList<Double>();
        //high value, fP(h2) must be negative
        double h2 = (k_i*(i_0+V_i)/(d_I*d_i)+sgm*Ae)/(sgm*d_P*d_R*d_r/(k_RA*k_r*(r_0+V_r)))+1;
        //find all the solutions and store them in variable "solutions"
        double h1 = 0.9*h2;
        double f2 = fP(h2);
        double err = 1e-6;
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

    //old code#####################################################
    


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
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1 - range + 2*range*rd.nextDouble();
    		BB[i] = b[i]*B1[i];
    	}
    	/*String bs = "";
    	for(double bi:b)
    		bs +=  bi+"\t";
    	System.out.println(bs);*/
    	//return AeThreshold(b);
    	return 0;
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
    
    public double so(double Ae){
    	LinkedList<Double> dd = solveStationary(Ae);
    	return dd.get(0);
    }
    public double dso(double Ae){
    	double dA = 1e-2;
    	double ds = (so(Ae+dA)-so(Ae))/dA;
    	return ds;
    }
    public double ddso(double Ae){
    	double dA = 1e-2;
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
    		}
    	}
    	return Am;
    }
    double minCurAe(){
    	double Am = -1;
    	double cAm = -1;
    	for(double A=0;A<=maxA;A++){
    		double cA = curvature(A);
    		if(cAm <0 || cAm > cA){
    			Am = A;
    			cAm = cA;
    		}
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
    public void test(){
    	//test bottom(Ae)
    	/*
    	Ae = 10;
    	double ba = bottom(Ae);
    	double P = 2.621969763699587;
    	System.out.println(fP(P));
    	System.out.println(dfP(P));
    	System.out.println(ba);
    	double dp = 1e-3;
    	print(ddfP(P));
    	System.out.println((dfP(P+dp)-dfP(P))/dp);
    	*/
    	//test AeThreshold()
    	/*double range = 0.05;
    	Random rd = new Random();
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		BB[i] = (1 - range + 2*range*rd.nextDouble())*B1[i];
    	}
    	assignParameters(B1);
    	double A = AeThreshold();
    	System.out.println(bottom(A*0.95)+","+bottom(A)+","+bottom(A*1.05));
    	*/
    	//test solveStationary(Ae)
    	
    	      	//solve steady state
    	/*
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
    	k_RA = B1[0];
    	d_P = B1[1]; 
    	k_r = B1[2]; 
    	d_R = B1[3]; 
    	r_0 = B1[4]; 
    	K_r = B1[5]; 
    	V_r = B1[6]; 
    	d_r = B1[7]; 
    	V_A = B1[8]; 
    	d_A = B1[9];
    	k_i = B1[10];
    	d_I = B1[11];
    	i_0 = B1[12];
    	V_i = B1[13];
    	d_i = B1[14];
    	K_i = B1[15];
    	sgm = B1[16];
    	k_M = B1[17];
    	d_M = B1[18];
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
    	k_M = B[17];
    	d_M = B[18];
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
    	int L = d.length;
    	for(int i=0;i<L;i++){
    		if(i==0)
    			print(d[i]);
    		else
    			print(","+d[i]);
    	}
    }
    double[] ones(int K){
    	double[] x = new double[K];
    	for(int i=0;i<K;i++)
    		x[i] = 1;
    	return x;
    }
	public double precision(double f,int i){
		int d = 0;
		double f2=f;
		while(f2<1){
			f2 = f2*10;
			d++;
		}
		BigDecimal b = new BigDecimal(f); 
		double f1 = b.setScale(d+i,BigDecimal.ROUND_HALF_UP).doubleValue(); 
		return f1;
	}
	public double[] precision(double[] BB){
		double[] BB1 = new double[BB.length];
		for(int i=0;i<BB.length;i++)
			BB1[i] = precision(BB[i],3);
		return BB1;
	}
}
