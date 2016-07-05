package quorumSensing;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.Random;

public class Dimer2Backup {
	private double[] B1 = {
			0.8944
			,65.933
			,519.601
			,1.002E-4
			,0.005936
			,1.851E-4
			,6.218E-5
			,4.061E-4
			,5.474E-4
			,0.01874
			,0.07113
			,1.399
			,40.657
			,0.008229
			,0.007422
			,5.816E-4
			,7.991E-4
			,3.0E-6
			,8.467E-6
			,0.01859
			};
			double[] B2={
					0.8983
					,64.061
					,304.175
					,8.888E-5
					,0.006551
					,1.766E-4
					,6.214E-5
					,4.328E-4
					,5.553E-4
					,0.01804
					,0.07013
					,1.47
					,41.087
					,0.007297
					,0.007314
					,6.418E-4
					,8.108E-4
					,3.0E-6
					,7.974E-6
					,0.01974
			};
			double[] B3={
			0.44140500000000005
			,83.6870267955146
			,830.071322226896
			,1.449472E-4
			,0.00322218
			,2.445241568669632E-4
			,7.175642796518999E-5
			,2.32305E-4
			,4.0868119996445844E-4
			,0.024267478906827573
			,0.07640578370896509
			,0.7446
			,38.22531417212959
			,0.011900629999999999
			,0.009129179575403452
			,3.15741E-4
			,5.966885186494567E-4
			,2.9999849573220737E-6
			,1.1370596770346575E-5
			,0.0096288
			};
			double[] B4 ={
			0.46186755474213165
			,15.597197246528248
			,4229.053905359846
			,2.157057615908111E-4
			,0.004325371566803694
			,0.0020043765258987017
			,7.947908628169854E-5
			,3.126184762271783E-4
			,0.002082151972390437
			,0.0030390182892096732
			,0.022502889232934208
			,1.3261310232715051
			,277.36388960028074
			,0.01899294047137961
			,0.0030876185510485175
			,4.22641179759869E-4
			,0.0030400131757921814
			,4.215508455218751E-6
			,2.3628330758242442E-6
			,0.13959686619935044
			};

	double maxA = 200;
	double maxP = 10000;
	double pdfm = 0;
	private int N = B1.length;
	private double sgm,V_A,K_A,k_RA,d_P,r_0,i_0,d_r,d_i,V_r,V_i,K_r,K_i,k_r,k_i,d_R,d_I,d_A,k_M,d_M;
	double lowMax = 12;
	double highMin = 40;
	double highMax = 5000;
	double peak0 = -1.5;
	double[] Bfit = B1.clone();
	Dimer2Backup(){
    	assignParameters();
    }
	Dimer2Backup(double[] B){
    	this.B1 = B;
    }
    private double Ae = 0;
    public static void main(String[] arg) throws FileNotFoundException{
    	Dimer2 b = new Dimer2();
    	b.run();
    }
    public void run() throws FileNotFoundException{
    	
    	//double[] B = set2();
    	//print1(precision(B));
    	
    	//checkComponents();
    //	plotAe();
    	//set4(B1);
    	//print1(B);
    	//B1 = B2;
    	//B2 = precision(B2);
    	double[] B = fit(B1);
    	print1(precision(B1));
    	assignParameters(B);
    	plotAe();
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = (B2[i]-B1[i])/B1[i];
    	}
    	//double[] b = new double[N];
    	//for(int i=0;i<N;i++)
    		//b[i] = (B4[i]-B1[i])/B1[i];
    	//print1(b);
    	//sample("data/sample_30_june.txt",0.03,10000);
    	/*double[] b = ones(N);
    	print1(gradient(b));*/
    	Ae=30;
    	print(allComponents());
    	//sample("data/compare_0_2.txt",0.2,100000);
    	//display();
    }
    void display(){
    	String[] s = {
    			"sgm"        
    			,"V_A"  
    			,"K_A"
    			,"k_RA"
    			,"d_P"
    			,"r_0"
    			,"i_0"
    			,"d_r"
    			,"d_i"
    			,"V_r"
    			,"V_i"
    			,"K_r"
    			,"K_i"
    			,"k_r"
    			,"k_i"
    			,"d_R"
    			,"d_I"
    			,"d_A"
    			,"k_M"
    			,"d_M"
    	};
    	for(int i=0;i<N;i++)
    		print(s[i]+"\t"+B1[i]);
    }
    void checkComponents(){
    	Ae = 0;
    	print(allComponents());
    	Ae = 30;
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
    double[] fit(double[] B){
    	double tB = totalScore(B);
    	double ERR = 1e-2;
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
    	if(AT > 25)
    		return (AT-25)/25;
    	if(AT < 15)
    		return (15-AT)/15;
    	return 0;
    }
    public double totalScore(double[] B){
    	//return isASimilar(B)+isLow(B)+isHigh(B)+isAe0(B);
    	return isLow(B)+isHigh(B)+isAe(B)+isAe0(B);
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
    	for(double Ae=0.0;Ae<=80.1;Ae=Ae+0.1){
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
        return V_A*I/(K_A+I) - d_A*A + sgm*(Ae-A);
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
        double dP = 1e-3;
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
    	return dfP(p);
    }
    private double peakDfp(double Ae, double[] BB){
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
/*    private double[] peakDfpGradient(double A,double[] b){
    	double[] g = new double[N];
    	double db=1e-3;
    	double bp = peakDfp(A,b);
    	for(int i=0;i<N;i++){
        	double[] bi = (double[]) b.clone();
        	bi[i] = b[i] + db;
        	g[i] = (peakDfp(A,bi)-bp)/db;
    	}
    	return g;
    }*/
   
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
    	double a = 0.49;
    	double bt = bottom(A,b);
    	double dr = 1e-2;
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
    	double a = 0.03;
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
    public void set4(double[] B){
    	double A = 15;
    	double bt = peakDfp(A,B);
    	double db = 1e-3;
    	double dr = 1e-3;
    	double[] BB = new double[N];
    	while(bt > -0.01){
    		double[] g = new double[N];
    		for(int i=0;i<N;i++){
        		double[] Bi = B.clone();
        		Bi[i] = B[i]*(1+db);
        		double bti = peakDfp(A,Bi);
        		g[i] = (bti-bt)/db;
    		}
			double maxg = 0;
    		for(int i=0;i<N;i++){
    			double gi = Math.abs(g[i]);
    			if(maxg < gi){
    				maxg=gi;
    			}
    		}
    		dr = 0.01/maxg;
    		//print(g);
    		for(int i=0;i<N;i++)
    			BB[i] = B[i]*(1 - g[i]*dr);
    		B = BB;
    		bt = peakDfp(A,B);
    		print("bt="+bt);
    	}
    	print1(B);
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
    		double[] Ai = oneSample(range);
    		out.println(covert(Ai));
        	if(i%100==0)
        		System.out.println(i);
    	}
    	out.close();
    }
    private double[] oneSample(double range){
    	Random rd = new Random();
    	double[] b = new double[N];
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1 - range + 2*range*rd.nextDouble();
    		BB[i] = b[i]*B1[i];
    	}
    	assignParameters(BB);
    	Ae = 0;
    	double[] a1 = allComponents();
    	Ae = 80;
    	double[] a2 = allComponents();
    	double[] a = new double[a1.length*2];
    	for(int i=0;i<a1.length;i++){
    		a[i] = a1[i];
    		a[i+a1.length] = a2[i];
    	}
    	return a;

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
    	sgm=B1[0];
    	V_A=B1[1];
    	K_A=B1[2];
    	k_RA=B1[3];
    	d_P=B1[4];
    	r_0=B1[5];
    	i_0=B1[6];
    	d_r=B1[7];
    	d_i=B1[8];
    	V_r=B1[9];
    	V_i=B1[10];
    	K_r=B1[11];
    	K_i=B1[12];
    	k_r=B1[13];
    	k_i=B1[14];
    	d_R=B1[15];
    	d_I=B1[16];
    	d_A=B1[17];
    	k_M=B1[18];
    	d_M=B1[19];
    }
    private void assignParameters(double[] B){
    	sgm=B[0];
    	V_A=B[1];
    	K_A=B[2];
    	k_RA=B[3];
    	d_P=B[4];
    	r_0=B[5];
    	i_0=B[6];
    	d_r=B[7];
    	d_i=B[8];
    	V_r=B[9];
    	V_i=B[10];
    	K_r=B[11];
    	K_i=B[12];
    	k_r=B[13];
    	k_i=B[14];
    	d_R=B[15];
    	d_I=B[16];
    	d_A=B[17];
    	k_M=B[18];
    	d_M=B[19];
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
	public String covert(double[] a){
		int L = a.length;
		String s = ""+a[0];
		for(int i=1;i<L;i++){
			s=s+","+a[i];
		}
		return s;
	}
}
