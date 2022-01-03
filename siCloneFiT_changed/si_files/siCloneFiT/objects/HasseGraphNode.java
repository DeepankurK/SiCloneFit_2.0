package objects;

import java.util.ArrayList;

public class HasseGraphNode {
	
	public int ID;
	public String name;
	public int inDegree;
	public int outDegree;
	public int countOccurrence;
	public ArrayList<String> adjacentNextNodes;
	public ArrayList<String> nodeLeaves;
	
	public HasseGraphNode(String nodeName){
		this.name = nodeName;
		this.adjacentNextNodes = new ArrayList<>();
		this.nodeLeaves = new ArrayList<>();
		this.countOccurrence = 1;
	}
	
	public void incrCount(){
		countOccurrence++;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
