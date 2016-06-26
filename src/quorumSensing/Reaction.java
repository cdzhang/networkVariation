package quorumSensing;

import java.util.Map;

public class Reaction {
	String name;
	String rateExpression;
	Map<String,Integer> reactants;
	String[] expressions;
	Reaction(String r, Map<String,Integer> re,String name){
		rateExpression = r;
		reactants = re;
		this.name = name;
		expressions = new String[re.size()];
		int i=0;
		for(String com:re.keySet()){
			String expi = com+"="+com+"+("+re.get(com)+")";
			expressions[i] = expi;
			i++;
		}
	}
	public String toString(){
		String re = "reaction:"+name + "\n";
		re = re + "rate:"+ rateExpression;
		for(String component:reactants.keySet()){
			Integer num = reactants.get(component);
			re = re + "\n"+ component + "," + num ;
		}
		return re;
	}
}
