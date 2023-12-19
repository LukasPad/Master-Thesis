package MCTS;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.BitSet;

public class NodeVertexSelection {
	
	
	
	
	ArrayList<BitSet> component;
	int[] mapping;
	NodeComponentSelection[] vertices;
	NodeComponentSelection parent;
	int score = Integer.MAX_VALUE;
	int bestIndex = -1;
	

	public NodeVertexSelection(ArrayList<BitSet> component) {
		BitSet mapBit = component.remove(0);
		this.component = component;
		this.parent = parent;
//		System.out.println();
//		System.out.println();
		this.mapping = new int[mapBit.cardinality()];
		this.vertices = new NodeComponentSelection[mapping.length];
		
		if(mapBit.cardinality() == 0) {
			System.out.println(" size 0 comp is: "+this.component);
		}
		
		int curr = 0;
		for (int i = 0; i < mapBit.cardinality(); i++) {
			mapping[i] = mapBit.nextSetBit(curr);
			curr = mapping[i]+1;
		}
//		System.out.println(component);
//		for (int i = 0; i < mapping.length; i++) {
//			System.out.println(mapping[i]);
//		}
	}
	
	public int expand(int treedepth, int[] settings) {
		// need heur for exploration/exploitation
		
		// random selection
		Random rand = new Random();
		System.out.println("mapping of size: "+mapping.length);
		int selectedIndex;
		boolean exploitNotExplore = true;
		// decide explore or exploit, using percentiles
		if(bestIndex > -1 && rand.nextInt(100) < settings[1]) {
			// exploit
			selectedIndex = this.bestIndex;
		} else {
			//explore
			selectedIndex = rand.nextInt(mapping.length);
		}
		int selected = mapping[selectedIndex];
		
		// unexplored
		int res;
		if(vertices[selectedIndex] == null) {
			vertices[selectedIndex] = new NodeComponentSelection(this.component, selected, selectedIndex);
			res = vertices[selectedIndex].getMaxSize();
		} else {
			res = vertices[selectedIndex].expand(treedepth, settings);
		}
		if(this.score > res) {
			this.score = res;
			this.bestIndex = selectedIndex;}
		
		return score;
		
		
	}

}
