package MCTS;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

public class NodeComponentSelection {
	
	
	
	int[] mapping;
	ArrayList<ArrayList<BitSet>> components;
	NodeVertexSelection[] compNodes;
	int score;
	int maxSize = -1;
	
	

	public NodeComponentSelection(ArrayList<BitSet> component, int vertex, int vertexIndex) {
		component.remove(vertexIndex);
		System.out.println("removing vert: "+vertex+" at index: "+vertexIndex);
		if (!component.isEmpty()) {
			this.shiftBy(vertex, component);
			this.components = this.generateComponents(component);
			compNodes = new NodeVertexSelection[components.size()];
		} else {
			this.score = 0;
			this.maxSize = 0;
		}
	}
	
	public int expand(int treedepth, int[] settings) {
		// need heur for exploration/exploitation
		if(components == null) {return treedepth;}
		int max = 0;
		for (int i = 0; i < components.size(); i++) {
			System.out.println(components);
			compNodes[i] = new NodeVertexSelection(components.get(i));
			int newMax = compNodes[i].expand(treedepth+1, settings);
			if(newMax > max) {max = newMax;}
		}
		this.score = max;
		// explored
		return max;
		
		
	}
	
	public void shiftBy(int shift, ArrayList<BitSet> graph) {
		for (int i = 0; i < graph.size(); i++) {
			graph.get(i).clear(shift);
			int temp = graph.get(i).nextSetBit(shift+1);
			while(temp > -1) {
				graph.get(i).set(temp-1);
				graph.get(i).clear(temp);
				temp = graph.get(i).nextSetBit(temp+1);
			}
		}
	}
	
	public ArrayList<ArrayList<BitSet>> generateComponents(ArrayList<BitSet> graph) {
		// setup
		ArrayList<ArrayList<BitSet>> components = new ArrayList<ArrayList<BitSet>>();
		BitSet toExamine = new BitSet();
		toExamine.set(0, graph.size());
		
		// while there are vertices unassigned to a component
		while (!toExamine.isEmpty()) {
			ArrayList<BitSet> comp = new ArrayList<BitSet>();
			BitSet remaining = new BitSet();
			BitSet tempMap = new BitSet();
			
			// pick a starting unassigned vertex
			int start = toExamine.nextSetBit(0);
			//toExamine.clear(start);
			remaining.set(start);
			//comp.add(graph.get(start));
			//System.out.println("added: "+start);
//			System.out.println("test: "+remaining.isEmpty());
			
			// while there are edges to non-taken vertices
			while(!remaining.isEmpty()) {
				System.out.println(remaining+"remains");
				int next = remaining.nextSetBit(0);
//				System.out.println("added: "+next);
				toExamine.clear(next);
				tempMap.set(next);
				remaining.clear(next);
				BitSet temp = (BitSet) graph.get(next).clone();
				temp.and(toExamine);
				remaining.or(temp);
//				System.out.println("new remaining: "+remaining);
//				System.out.println("new examining: "+toExamine);
//				
				// add the unassigned connected vertices
				
				comp.add(graph.get(next));
				}
			comp.add(0, tempMap);
			components.add(comp);
			}
		return components;
	}

	public int getMaxSize() {
		if(this.maxSize < 0) {
			for (int i = 0; i < components.size(); i++) {
				if(this.maxSize < components.get(i).size()-1) {
					this.maxSize = components.get(i).size()-1;
				}
			}
		}
		return this.maxSize;
	}

}
