/**
 * 
 */
package MCTS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;

import main.Inscriber;

/**
 * @author Lukas Padolevicius
 *
 */
public class MCTSPrimer {
	
	int[] parentList;
	long elimForestSearched = 0;
	long elimTreesSearched = 0;
	
	
	public long[] run(File file) throws IOException {
		// 
		Inscriber inscriber = new Inscriber();
		ArrayList<BitSet> graph = inscriber.inscribe(file, true);
		
		parentList = new int[graph.size()];
		for (int i = 0; i < graph.size(); i++) {
			System.out.println(graph.get(i));
		}
		BitSet mapping = new BitSet(graph.size());
		mapping.flip(0, graph.size());
		graph.add(0, mapping);
		System.out.println("pause");
		

		

//		for (int i = graph.size()-1; i > -1; i--) {
//			graph.get(i).clear(1);
////			if(graph.get(i).isEmpty()) { 
////				graph.remove(i);
////				this.shiftBy(i, graph);
////				}
//		}
//		graph.remove(1);
//		this.shiftBy(1, graph);
//		mapping.clear(1);
//		
//
//		for (int i = 0; i < graph.size(); i++) {
//			System.out.println(graph.get(i));
//		}
		
//		ArrayList<ArrayList<BitSet>> pr = this.generateComponents(graph);
		
		NodeVertexSelection vrt = new NodeVertexSelection(graph);
		int[] set = {1, 80};
		for (int i = 0; i < 1000; i++) {
			vrt.expand(0, set);
		}
		
//		for (int i = 0; i < pr.size(); i++) {
//			System.out.println(pr.get(i));
//		}
		
		

		long[] res = new long[graph.size()+3];
		res[0] = vrt.score;
		//res[0] = k;
		res[1] = this.elimForestSearched;
		res[2] = this.elimTreesSearched;
		for(int i=3; i<res.length; i++) {
			res[i] = parentList[i-3];
		}
		return res;
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
			System.out.println("test: "+remaining.isEmpty());
			
			// while there are edges to non-taken vertices
			while(!remaining.isEmpty()) {
				int next = remaining.nextSetBit(0);
				System.out.println("added: "+next);
				toExamine.clear(next);
				tempMap.set(next);
				remaining.clear(next);
				BitSet temp = (BitSet) graph.get(next).clone();
				temp.and(toExamine);
				remaining.or(temp);
				System.out.println("new remaining: "+remaining);
				System.out.println("new examining: "+toExamine);
				
				// add the unassigned connected vertices
				
				comp.add(graph.get(next));
				}
			comp.add(0, tempMap);
			components.add(comp);
			}
		return components;
	}
	
}
