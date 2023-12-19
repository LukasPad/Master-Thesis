/**
 * 
 */
package main;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Random;

/**
 * @author Lukas Padolevicius
 *
 */
public class GeneralPrimer {
	
	int[] parentList;
	long elimForestSearched = 0;
	long elimTreesSearched = 0;
	Random randomizer;
	
	
	public long[] run(String fileName) throws IOException {
		Inscriber inscriber = new Inscriber();
		BitSet[] graph = inscriber.inscribe(fileName);
		parentList = new int[graph.length];
	
		BitSet mapping = new BitSet(graph.length);
		mapping.flip(0, graph.length);
		
		this.randomizer = new Random();
		long[] res = new long[graph.length+3];
		//String random = list.get(randomizer.nextInt(list.size()));
		
		int k = 0;
		int runs = 0;
		long startTime = System.currentTimeMillis();
		int bestk = Integer.MAX_VALUE;
		while(runs < 100000) {
			int newBestk = elimForest(graph, mapping, k, 0);
			if(newBestk < bestk) {
				System.out.println("found k: "+newBestk+"   in: "+(System.currentTimeMillis()- startTime));
				bestk = newBestk;
				for(int i=3; i<res.length; i++) {
					res[i] = parentList[i-3];
				}
			}
			runs++;
			long time2 = System.nanoTime();
		}
		
		res[0] = bestk;
		res[1] = this.elimForestSearched;
		res[2] = this.elimTreesSearched;
		for(int i=3; i<res.length; i++) {
			res[i] = parentList[i-3];
		}
		return res;
	}
	
	public BitSet[] removeVertex(BitSet[] graph, int vertex) {

		// copy graph
		BitSet[] reducedGraph = new BitSet[graph.length];
		for(int i = 0; i<reducedGraph.length; i++) {
			if(graph[i] != null) {
				reducedGraph[i] = (BitSet) graph[i].clone();
			}
		}
		// remove vertex
		reducedGraph[vertex] = null;
		
		//remove edges to removed vertex from its neighbours
		int next = graph[vertex].nextSetBit(0);
		while(next != -1) {
			reducedGraph[next-1].clear(vertex+1);
			next = graph[vertex].nextSetBit(next+1);
		}
		return reducedGraph;
	}
	
	public ArrayList<BitSet[]> generateComponents(BitSet[] graph, BitSet mapping) {
		ArrayList<BitSet[]> components = new ArrayList<BitSet[]>();
		components.add(new BitSet[graph.length]);
		BitSet toExamine = (BitSet) mapping.clone();
		int mapIndex = 0;
		while(!toExamine.isEmpty()) {
			BitSet compMap = new BitSet();
			BitSet frontier = new BitSet();
			frontier.set(toExamine.nextSetBit(0));
			toExamine.clear(toExamine.nextSetBit(0));
			while(!frontier.isEmpty()) {
				int current = frontier.nextSetBit(0);
				toExamine.clear(current);
				compMap.set(current);
				frontier.clear(current);
				int next2 = graph[current].nextSetBit(0);
				while(next2 != -1) {
					frontier.set(next2-1);
					next2 = graph[current].nextSetBit(next2+1);
				}
				frontier.and(toExamine);
			}
			components.get(0)[mapIndex] = compMap;
			mapIndex++;
		}
		for(int i=0; i<mapIndex; i++) {
			// for every component mapping, carve out the component
			BitSet[] newComp = new BitSet[graph.length];
			int next = components.get(0)[i].nextSetBit(0);
			while(next != -1) {
				newComp[next] = graph[next];
				next = components.get(0)[i].nextSetBit(next+1);
			}
			components.add(newComp);
		}		
//		if(components.size() > 2) {
//			for(int x=0; x<components.size(); x++) {
//				for(int y=0; y<components.get(x).length; y++) {
//					System.out.println(components.get(x)[y]);
//				}
//				System.out.println();
//			}
//		}
		return components;
	}
	
	public int elimTree(BitSet[] graph, BitSet mapping, int k, int parent) {
		this.elimTreesSearched++;
		if(mapping.cardinality() == 1) {
			parentList[mapping.nextSetBit(0)] = parent;
			return k+1;
		} 
		int next = mapping.nextSetBit(0);
		boolean onlyChild = true;
		while(onlyChild) {
			next = mapping.nextSetBit(0);
			int rand = randomizer.nextInt(mapping.cardinality());
			for(int x=0;x<rand;x++) {
				next = mapping.nextSetBit(next+1);
			}
			if(graph[next].cardinality() > 1 || mapping.cardinality() < 3) {
				onlyChild = false;
			}
		}
		
		parentList[next] = parent;
		BitSet[] newGraph = removeVertex(graph, next);
		
		BitSet newMapping = (BitSet) mapping.clone();
		newMapping.clear(next);
		return elimForest(newGraph, newMapping, k+1, next+1);
	}
	
	public int elimForest(BitSet[] graph, BitSet mapping, int k, int parent) {
		this.elimForestSearched++;
						
		// make connected components
		long time1 = System.nanoTime();
		ArrayList<BitSet[]> components = generateComponents(graph, mapping);
		long time2 = System.nanoTime();
		int maxk = -1;
		for(int x=1; x<components.size(); x++) {
			maxk = Math.max(maxk, elimTree(components.get(x), components.get(0)[x-1], k, parent));
			components.get(0)[x-1].clear();
		}
		return maxk;
	}

	
}
