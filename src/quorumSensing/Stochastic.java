package quorumSensing;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Random;

import bsh.EvalError;
import bsh.Interpreter;
public class Stochastic {
	int Nc,Nr,Np;//number of components, reactions and parameters
	double[] components;
	String[] comNames;
	double[] para;
	String[] paraNames;
	Reaction[] reactions;
	double time=0;
	Interpreter eva = new Interpreter();
	Random random = new Random();
	Stochastic(){}
	Stochastic(String input) throws IOException, EvalError{
		parseInput(input);
	}
	public static void main(String[] arg) throws EvalError, IOException{
		Interpreter i = new Interpreter();
		i.set("x",1.2);
		i.set("y",2.3);
		i.eval("x=(-x)+2");
		Object re = i.get("x");
		System.out.println(re);
		System.out.println(i.get("x").getClass());
		Stochastic s = new Stochastic();
		s.run();
	}
	public void run() throws IOException, EvalError{
		String input = "data/input.txt";
		parseInput(input);
		String outputFile = "data/output.txt";
		PrintWriter out = new PrintWriter(outputFile);
		String output="time";
		for(String com:comNames){
			output += ","+com;
		}
		out.println(output);
/*		print("Nc="+Nc);
		print("components");
		for(int i=0;i<Nc;i++){
			print(components[i]+","+comNames[i]);
		}
		print("parameters");
		for(int i=0;i<Np;i++){
			print(paraNames[i]+"="+para[i]);
		}
		print("reactions");
		for(int i=0;i<Nr;i++){
			print(reactions[i]);
		}
		print1(rate());
		print("testlog");*/
		for(int i=0;i<2000000;i++){
			react();
			if(i%100==0){
				printComponents();
				out.println(getOutput());
			}
		}
		out.close();
	}
	public void react() throws EvalError{
		double[] rates = rate();
		double t = nextReactionTime(rates);
		int i = getReaction(rates);
		performReaction(i);
		time+=t;
		if(time>=100)
			return;
	}
	private double nextReactionTime(double[] rates){
		double lambda = 0;
		for(double ri:rates)
			lambda+=ri;
		double x = random.nextDouble();
		while(x==0)
			x = random.nextDouble();
		return - Math.log(x)/lambda;
	}
	private int getReaction(double[] rates){
		double lambda = 0;
		for(double ri:rates)
			lambda+=ri;
		double x = random.nextDouble();
		double sum = 0;
		for(int i=0;i<Nr;i++){
			sum += rates[i]/lambda;
			if(sum>x)
				return i;
		}
		return Nr-1;
	}
	private void performReaction(int i) throws EvalError{
		Reaction ri = reactions[i];
		for(String expression:ri.expressions)
			eva.eval(expression);
	}
	public double[] rate() throws EvalError{
		double[] rates = new double[Nr];
		for(int i=0;i<Nr;i++){
			eva.eval("_result="+reactions[i].rateExpression);
			rates[i] = (Double) eva.get("_result");
		}
		return rates;
	}
	private void parseInput(String file) throws IOException, EvalError{
		BufferedReader in = new BufferedReader(new FileReader(file));
		String state="";
		String line="";
		LinkedList<String> com = new LinkedList<String>();
		LinkedList<String> lpa = new LinkedList<String>();
		LinkedList<Double> lpv = new LinkedList<Double>();
		LinkedList<String> rea = new LinkedList<String>();
		LinkedList<String> reaRate = new LinkedList<String>();
		LinkedList<Double> iv = new LinkedList<Double>();
		while((line=in.readLine())!=null){
			line = line.trim();
			line = line.replaceAll("//.*$","");
			if(line.length()<2)
				continue;
			if(line.length()>=2 &&line.charAt(0)=='/'&&line.charAt(1)=='/')//comment
				continue;
			if(line.charAt(0)=='#'){
				state=line.replaceAll("#", "");
				continue;
			}
			if(state.toLowerCase().contains("components")){
				String[] lc = line.split(",");
				int ln = lc.length;
				for(int i=0;i<ln;i++){
					String li = lc[i].trim();
					if(li!="")
						com.add(lc[i]);
				}
			}else if(state.toLowerCase().contains("parameters")){
				String[] lc = line.split(",");
				int ln = lc.length;
				for(int i=0;i<ln;i++){
					String li = lc[i];
					String[] lia = li.split(":");
					lpa.add(lia[0]);
					lpv.add(Double.parseDouble(lia[1]));
				}
			}else if(state.toLowerCase().contains("reactions")){
				String[] rvs = line.split(",");
				int ln = rvs.length;
				for(int i=0;i<ln;i++){
					String rvi = rvs[i];
					String[] rviv = rvi.split(":");
					rea.add(rviv[0]);
					reaRate.add(rviv[1]);
				}
			}else if(state.toLowerCase().contains("initialvalues")){
				String[] lc = line.split(",");
				for(String lci:lc){
					if(lci.trim()!="")
						iv.add(Double.parseDouble(lci));
				}
			}
		}
		in.close();
		Nc = com.size();
		components = new double[Nc];
		comNames = new String[Nc];
		for(int i=0;i<Nc;i++){
			comNames[i] = com.get(i);
			components[i] = iv.get(i);
		}
		Np = lpa.size();
		paraNames = new String[Np];
		para = new double[Np];
		for(int i=0;i<Np;i++){
			paraNames[i] = lpa.get(i);
			para[i] = lpv.get(i);
		}
		Nr = rea.size();
		reactions=new Reaction[Nr];
		for(int i=0;i<Nr;i++){
			Map<String,Integer> ri = new HashMap<String,Integer>();
			String reaction = rea.get(i);
			String[] rs = reaction.split("->");
			String[] reactants = rs[0].split("\\+");
			String[] products = rs[1].split("\\+");
			for(int j=0;j<reactants.length;j++){
				String rtj = reactants[j];
				String ms = "";
				String rt = "";
				for(int k=0;k<rtj.length();k++){
					char ck = rtj.charAt(k);
					if(ck>='0' && ck<='9'){
						ms = ms + ck;
					}else{
						rt = rtj.substring(k);
						break;
					}
				}
				if(ms=="") ms="1";
				Integer mu = Integer.parseInt(ms);
				Integer or = ri.get(rt);
				if(or==null) or = 0;
				ri.put(rt,or-mu);
			}
			for(int j=0;j<products.length;j++){
				String rtj = products[j];
				String ms = "";
				String rt = "";
				for(int k=0;k<rtj.length();k++){
					char ck = rtj.charAt(k);
					if(ck>='0' && ck<='9'){
						ms = ms + ck;
					}else{
						rt = rtj.substring(k);
						break;
					}
				}
				if(ms=="") ms="1";
				Integer mu = Integer.parseInt(ms);
				Integer or = ri.get(rt);
				if(or==null) or = 0;
				ri.put(rt,or+mu);
			}
			LinkedList<String> flag = new LinkedList<String>();
			for(String rkey:ri.keySet()){
				if(rkey.toLowerCase().equals("null")||ri.get(rkey)==0)
					flag.add(rkey);
			}
			for(String rkey:flag)
				ri.remove(rkey);
			reactions[i] = new Reaction(reaRate.get(i),ri,reaction);
			
		}
		setValues();
	}
	private void setValues() throws EvalError{
		for(int i=0;i<Nc;i++)
			eva.set(comNames[i], components[i]);
		for(int i=0;i<Np;i++)
			eva.set(paraNames[i], para[i]);
	}
	private void setParameters() throws EvalError{
		for(int i=0;i<Np;i++)
			eva.set(paraNames[i], para[i]);
	}
	private void setComponents() throws EvalError{
		for(int i=0;i<Np;i++)
			eva.set(comNames[i], components[i]);
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
    private void printComponents() throws EvalError{
    	print("time="+time+"s");
    	for(String com:comNames)
    		print(com+":"+eva.eval(com));
    }
    private String getOutput() throws EvalError{
    	String output = time +"";
    	for(String com:comNames)
    		output += ","+eva.eval(com);
    	return output;
    }
}
