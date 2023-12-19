/**
 * 
 */
package main;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * @author Lukas Padolevicius
 *
 */
public class Algorithm3 {
	
	int[] parentList;
	int elimForestSearched = 0;
	int elimTreesSearched = 0;
	
	
	
	public int[] run() throws IOException {
		Inscriber inscriber = new Inscriber();
		BitSet[] graph = inscriber.inscribe("src/input_graphs/exact_003.gr");
		parentList = new int[graph.length];
		
		int k = 0;
		while(!elimForest(graph, k, 0)) {
			k++;
			System.out.println("new k is:"+k);
		}
		int[] res = new int[graph.length+1];
		res[0] = k;
		for(int i=1; i<res.length; i++) {
			res[i] = parentList[i-1];
		}
		return res;
	}
	
	public int shiftParent(int parent, ArrayList<Integer> removed) {
		for(int i=0; i<removed.size(); i++) {
			if(parent > removed.get(i)) {
				parent++;
			} else {break;}
		}
		return parent;
	}
	
	public BitSet[] removeVertex(BitSet[] graph, int vertex, ArrayList<Integer> removed) {
		// add vertex to removed list in order
		int ind = 0;
		if(removed.isEmpty()) {
			removed.add(vertex);
		} else {
			while(removed.get(ind) < vertex) {
				ind++;
			}
			removed.add(ind, vertex);
		}
		
		// copy bitsets before removed vertex
		BitSet[] reducedGraph = new BitSet[graph.length-1];
		for(int i = 0; i<vertex; i++) {
			reducedGraph[i] = (BitSet) graph[i].clone();
		}
		//and after removed vertex
		for(int i = vertex; i<graph.length-1; i++) {
			reducedGraph[i] = (BitSet) graph[i+1].clone();
		}
		
		//remove edges to removed vertex from its neighbours
		int next = graph[vertex].nextSetBit(0);
		while(next != -1) {
			if(next < vertex) {
				
			}
			next = graph[vertex].nextSetBit(next+1);
		}

	}
	
	public ArrayList<BitSet[]> generateComponents(BitSet[] graph, int removal) {

		return components;
	}
	
	public Boolean elimTree(BitSet[] graph, BitSet mapping, int k, int parent, ArrayList<Integer> removed) {
		if(mapping.cardinality() == 1) {
			parentList[mapping.nextSetBit(0)] = shiftParent(parent, removed);
			return true;
		} 
		int next = mapping.nextSetBit(0);
		while(next != -1) {
			parentList[next] = shiftParent(parent, removed);
			BitSet[] newGraph = removeVertex(graph, next, removed);
			
			//shift all nodes beyond the removed by 1
			mapping.clear(next);
			int next2 = mapping.nextSetBit(next);
			while(next2 != -1) {
				
				next2 = mapping.nextSetBit(next2+1);
			}
			if(elimForest(newGraph, ))

		}
		
		return false;
		
		
		// PRUNING METHOD: consider 
		else if(pruningType == "closestFOD") {
			int lim = components.size()-1;
			double threshold = k * this.param2;
			for(int i=0; i<lim; i++) {
				double closest = Double.MAX_VALUE;
				int closestInd = -1;
				for(int j=1; j<components.size(); j++) {
					int card = components.get(0)[j-1].cardinality();
					if(card != 0 && Math.abs(card - threshold) < closest) {
						closest = Math.abs(card - threshold);
						closestInd = j;
					}
				}
				//System.out.println(closestInd);
				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {return false;}
				components.get(0)[closestInd-1].clear();
			};
		}
		// PRUNING METHOD: consider 
		else if(pruningType == "closestOverFOD") {
			int lim = components.size()-1;
			double threshold = k * this.param2;
			for(int i=0; i<lim; i++) {
				double closest = Double.MAX_VALUE;
				int closestInd = -1;
				for(int j=1; j<components.size(); j++) {
					int card = components.get(0)[j-1].cardinality();
					if(card != 0 && Math.abs(card - threshold) < closest) {
						closest = Math.abs(card - threshold);
						closestInd = j;
					}
				}
				// Get the first non-empty component if threshold not reached
				if(closestInd == -1) {
					int tempInd = 1;
					while(components.get(0)[tempInd-1].cardinality() == 0) {
						tempInd++;
					}
					closestInd = tempInd;
				}
				//System.out.println(closestInd);
				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {return false;}
				components.get(0)[closestInd-1].clear();
			};
		}
		
	}

	
	public Boolean elimForest(BitSet[] graph, BitSet mapping, int k, int parent, ArrayList<Integer> removed) {
		this.elimForestSearched++;
		
		//System.out.println("Forest searched:"+elimForestSearched);
		if(k == 0 && g > 1) {return false;}
		return true;
	}
	
}
