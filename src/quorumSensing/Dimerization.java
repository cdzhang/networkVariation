package quorumSensing;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Random;

public class Dimerization {
	private double[] B1 = {
		1e-4,
		3e-3,
		1.6e-2,
		2e-4,
		3e-4,
		1,
		//4.8e-3,
		4e-2,
		6.0e-3,
		0.45,
		3.3e-6,
		1.6e-2,
		5e-5,
		1.5e-4,
		//2e-3,
		2e-2,
		6.0e-3,
		30,
		0.4,
		1e-5,
		1e-2
	};
	private double[] B2 ={
			1.0249848709424399E-4
			,0.002925080488180573
			,0.01639975793507904
			,1.9500536587870558E-4
			,3.0046383502596114E-4
			,0.9987745209122778
			,0.004912510778983533
			,0.005850160976361167
			,0.45976490886843013
			,3.2999993196089747E-6
			,0.016347196759766404
			,4.891792645961548E-5
			,1.5053999999999927E-4
			,0.0020361995949708103
			,0.005870151175153857
			,29.662137274453322
			,0.3900108141075393
			,1.0124909818684886E-5
			,0.009875119431918052
	};	
	private double[] B3 = {
			1.1400807680437839E-4
			,0.002514026806881732
			,0.018241292288700636
			,1.676017871254488E-4
			,3.2102564427791526E-4
			,0.9254445103704212
			,0.005156978205484437
			,0.005028053613763437
			,0.4941414257000826
			,3.299995589842403E-6
			,0.017569472913780713
			,4.459262082314368E-5
			,1.643253771241253E-4
			,0.0020056499914457636
			,0.005351114498777253
			,29.91528660397213
			,0.33520420629316183
			,1.0722311424249404E-5
			,0.009224005484558555
	};
	private double[] B4 = {
			1.1093137510577158E-4
			,0.002631139360261045
			,0.01774901999216104
			,1.754092917250148E-4
			,4.201270508373042E-4
			,0.9380602134491092
			,0.0024927612502706207
			,0.00526227869856632
			,0.3898624290619653
			,3.2999966538586362E-6
			,0.013861775286160935
			,5.588327268291611E-5
			,1.5000000050293102E-4
			,0.0017327219171294744
			,0.006705992660652299
			,33.467883459862264
			,0.35081904425742183
			,9.355873826689668E-6
			,0.01060479079371234
	};
	double[] B5 = {
			4.284865089405996E-5
			,0.006979420801638231
			,0.006855784143049592
			,4.652947201092135E-4
			,3.0365236384310616E-4
			,0.9998974382200171
			,0.019755047439057982
			,4.0035145456944027E-4
			,0.20026098282256044
			,3.3000158118205687E-6
			,0.007120390493221762
			,1.120313983093585E-4
			,1.5126151743565598E-4
			,0.013780609238675296
			,2.764800856292889E-4
			,29.927184999976156
			,0.8962469137147633
			,6.54888003964798E-6
			,0.01525871510146615
	};
	private double[] B6={
			7.061447792294671E-5
			,0.003676402657024701
			,0.01129831646767148
			,2.450935104683124E-4
			,2.990174122135705E-4
			,1.0636644600455423
			,0.028429949104901577
			,0.007352805314049413
			,0.13239242727244205
			,3.300006139397506E-6
			,0.0047072863030203135
			,6.915085992091662E-5
			,1.4829238362876199E-4
			,0.006609005665026304
			,0.00829810319051001
			,31.74092336430413
			,0.49018641388046896
			,8.655937920497144E-6
			,0.011184122288639296};
	private double[] B7={
			5.106172032139481E-5,
			0.0050741628219619726,
			0.008169875251423151,
			3.3827752146412996E-4,
			3.028797006284402E-4,
			1.0601561030473354,
			0.027457305908206706,
			5.543260330603589E-4,
			0.09084189928979636,
			3.300011710597097E-6,
			0.003229934074866397,
			1.0065372320861886E-4,
			1.5217997346444082E-4,
			0.01376940320536908,
			2.7652542181775327E-4,
			28.804877549193478,
			0.7134980826885645,
			7.3624315609793E-6,
			0.013142658612101718
	};
	double[] B8={4.024918535753244E-5
			,0.00605430751352177
			,0.006439869657205328
			,4.0362050090144466E-4
			,3.030715793593597E-4
			,1.0682459326201552
			,0.02496216091853842
			,4.281849473304441E-4
			,0.058549280134397856
			,3.300014518108974E-6
			,0.0020817521088981523
			,1.3597266533201232E-4
			,1.5195610882956718E-4
			,0.011350170482061517
			,2.0012773289851337E-4
			,31.44170092283294
			,0.8807099436450067
			,6.598682208038302E-6
			,0.01444466766721411};
	double maxA = 200;
	double maxP = 10000;
	double pdfm = 0;
	private int N = B1.length;
	private double k_RA,d_P,k_r,d_R,r_0,K_r,V_r,d_r,V_A,d_A,k_i,d_I,i_0,V_i,d_i,K_i,sgm,k_M,d_M;
    Dimerization(){
    	assignParameters();
    }
    Dimerization(double[] B){
    	this.B1 = B;
    }
    private double Ae = 23.561563490726638;
    public static void main(String[] arg) throws FileNotFoundException{
    	Dimerization b = new Dimerization();
    	b.run();
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
    double[] nextPointMin(double[] B){
    	assignParameters(B);
    	double mB = minConcentration();
    	double db = 1e-3;
    	double delta = 1e-3;
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
    	double delta = 1e-1;
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
    	while(min1<minc1||max1>maxc1||min2<minc2||max2>maxc2){
    		if(min1/minc1<min2/minc2||max1/maxc1>max2/maxc2){
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
    	}
    	return BB;
    }
    public void run() throws FileNotFoundException{
/*    	for(int i=0;i<10000000;i++){
    		if(oneSample1(0.05)){
    			print(i);
    			break;
    		}
    		if(i%1000000==0)
    			print(i);
    	}*/
    	Ae=40;
    	B1=B8;
    	assignParameters(B1);
    	//double[] B = set5(50,5000);
    	print(allComponents());
    	Ae=0;
    	print(allComponents());
    	for(double sp:solveStationary(Ae)){
    		print(sp);
    	}
    	print("#############");
    	//plotAe();
    	//tune(10);
    	//set1();
    	//plotAe();
    	//plotfP(B6,A);
/*    	Ae=40;
    	B1=B6;
    	double[] A = {0,15,40};
    	//plotfP(B6,A);
    	LinkedList<Double> ll = solveStationary(0);
    	Ae=0;
    	for(double l:ll){
    		print(l+","+fP(l));
    		print(allComponents(l));
    	}*/
    	//set1();
    	//double[] BB = set5(50,5000);
    	//print(BB);
    	//print(allComponents());
    	//assignParameters(B5);
    	//plotAe();
    }
    public void tune(int K){
    	double[] BB = new double[N];
    	for(int i=0;i<K;i++){
    		Ae = 40;
    		BB = set5(70,5000);
        	print(allComponents());
        	Ae = 0;
        	print(allComponents());
        	B1 = BB;
        	assignParameters(B1);
        	BB = set1();
        	
    	}
    	print1(BB);
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
    	assignParameters(B4);
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
    public void plotdfP(double[] B1,double[] A){
    	assignParameters(B1);
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
    		BB[i] = b[i] * B1[i];
    	print(BB);
    }
    public double[] set1(){
    	double A = 0;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double bt = peak(A,b);
    	double dr = 1e-1;
    	while(bt > -10){
    		double[] g = peakGradient(A,b);
    		for(double gi:g)
    			if(dr>1e-3/Math.abs(gi))
    				dr = 1e-3/Math.abs(gi);
    		for(int i=0;i<N;i++)
    			b[i] = b[i] - g[i]*dr;
    		bt = peak(A,b);
    		print(bt);
    	}
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print(BB);
    	return BB;
    }   
    public void set2(){
    	double A = 0;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double bt = peak(A,b);
    	double dr = 1e-3;
    	while(bt < 2){
    		double[] g = peakGradient(A,b);
    		for(int i=0;i<N;i++)
    			b[i] = b[i] + g[i]*dr;
    		bt = peak(A,b);
    	}
    	print(bt);
    	print(b);
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print(BB);
    }   
    public void set4(){
    	double A = 15;
    	double[] b = new double[N];
    	for(int i=0;i<N;i++){
    		b[i] = 1;
    	}
    	double bt = peakDfp(A,b);
    	double dr = 1e-3;
    	while(bt > -0.01){
    		double[] g = peakDfpGradient(A,b);
    		for(int i=0;i<N;i++)
    			b[i] = b[i] - g[i]*dr;
    		bt = peakDfp(A,b);
    	}
    	double[] BB = new double[N];
    	for(int i=0;i<N;i++)
    		BB[i] = b[i] * B1[i];
    	print(BB);
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
}
