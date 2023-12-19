/**
 * 
 */
package main;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Random;
import java.io.File;


/**
 * @author Lukas Padolevicius
 *
 */
public class Algorithm {
	
	private int[] parentList;
	private long elimForestSearched = 0;
	private long elimTreesSearched = 0;
	private double[] threshByLayerU = new double[1000];
	private double[] threshByLayerWeightU = new double[1000];
	private int[] threshByLayerWeightoccU = new int[1000];
	private double[] threshByLayerL = new double[1000];
	private double[] threshByLayerWeightL = new double[1000];
	private int[] threshByLayerWeightoccL = new int[1000];
	private double[] threshByLayerUE = new double[1000];
	private double[] threshByLayerWeightUE = new double[1000];
	private int[] threshByLayerWeightoccUE = new int[1000];
	private double[] threshByLayerLE = new double[1000];
	private double[] threshByLayerWeightLE = new double[1000];
	private int[] threshByLayerWeightoccLE = new int[1000];
	private boolean reorder = false;
	private boolean simple_lb = false;
	private boolean path_lb = false;
	private boolean dom_rule = false;
	private int node_amount = 0;
	private int edge_amount = 0;
	private long removingTime = 0;
	private long compTime = 0;
	private int b = 0;
	private int e;
	private int v;
	private boolean reorderstart;
	private double shV;
	private double shE;
	
	
	// learning to branch techniques
	private double param;
	private String pruningType;
	
	public long[] run(File file, String pruningName, double parameter, boolean[] settings, double shiftV, double shiftE) throws IOException {
		
		this.param = parameter;
		this.reorder = settings[0];
		this.reorderstart = settings[1];
		this.simple_lb = settings[2];
		this.path_lb = settings[3];
		this.dom_rule = settings[4];
		this.pruningType = pruningName;
		this.shV = shiftV;
		this.shE = shiftE;

		Inscriber inscriber = new Inscriber();
		BitSet[] graph = inscriber.inscribe(file);
		int[] tempparams= inscriber.get_params();
		this.e = tempparams[0];
		this.v = tempparams[1];
		if(this.dom_rule && ((float) e)/((float) v) <1.275) {
			this.dom_rule = false;
		}
		parentList = new int[graph.length];
		
		long time1 = System.nanoTime();
		
		this.node_amount = graph.length;
		// find lb param b, and number of edges
		for(int i=0; i<graph.length; i++) {
			this.edge_amount += graph[i].cardinality();
			if(this.b < graph[i].cardinality()) {this.b = graph[i].cardinality();}
		}
		this.edge_amount = this.edge_amount/2;
		
		// Reorder vertices by degree
		
		if(this.reorderstart) {
			
////			for(int i=0; i< graph.length; i++) {
////				System.out.println((i)+"  is: "+graph[i]);
////			}
//			
			ArrayList<Integer> ordered = new ArrayList<Integer>();
			BitSet[] newGraph = new BitSet[graph.length];
			for(int i=0; i< graph.length; i++) {
				int ind = 0;
				while(ind < ordered.size() && graph[ordered.get(ind)].cardinality() > graph[i].cardinality()) {ind++;}
				ordered.add(ind, i);
			}
			
//			for(int i=0; i< graph.length; i++) {
//				System.out.println(ordered.get(i));
//			}
			
			// remember original indexing
			int[] oldIndexes = new int[graph.length];
			for(int i=0; i< graph.length; i++) {
				oldIndexes[ordered.get(i)] = i;
			}
			
			// substitute graph numbers
			for(int i=0; i< graph.length; i++) {
				newGraph[i] = new BitSet();
				int next = graph[ordered.get(i)].nextSetBit(0);
				while(next != -1) {
					//System.out.println("next is: "+next);
					newGraph[i].set(oldIndexes[next-1]+1);
					next = graph[ordered.get(i)].nextSetBit(next+1);
				}
			}
			
//			for(int i=0; i< graph.length; i++) {
//				System.out.println((i)+"  is: "+newGraph[i]);
//			}
			graph = newGraph;
			if(newGraph[newGraph.length-1].isEmpty()) {
				System.out.println("bad reorder");
				System.exit(0);
			}
		}

		BitSet mapping = new BitSet(graph.length);
		mapping.flip(0, graph.length);
		int k = 0;
		while(!elimForest(graph, mapping, k, 0)) {
			long time2 = System.nanoTime();
			
//			for(int x=0; x<k-1; x++) {
//				System.out.print(","+threshByLayerL[x]+","+threshByLayerWeightL[x]+","+threshByLayerU[x]+","+threshByLayerWeightU[x]);
//			}
//			System.out.println("     k= "+k);
//			for(int x=0; x<k-1; x++) {
//				System.out.print(","+threshByLayerLE[x]+","+threshByLayerWeightLE[x]+","+threshByLayerUE[x]+","+threshByLayerWeightUE[x]);
//			}
//			System.out.println("     k= "+k);
			//System.out.println("time taken for k: "+k+"  is: "+(time2-time1));
			//System.out.print((time2-time1)+",");
			time1 = System.nanoTime();
			k++;
			//System.out.println("new k is:"+k);
		}
		//System.out.println();
		//TODO: re-substitute graph numbers for param 27.0
		
		if(pruningType == "analysis") {
//			System.out.println();
//			for(int i=0; i < 5 ;i++) {
//				System.out.println(threshByLayer[i]+"  ,bob,  "+threshByLayerWeight[i]);
//			}

			//int index = 0;
//			while(threshByLayerWeight[index] != 0) {
//				index++;
//			}
			//System.out.println();
			String[][] analysis = new String[k][2];
			for(int i=0; i<k; i++) {
				analysis[i][0] = String.format("%.2f", threshByLayerL[i]/threshByLayerWeightL[i]);
				analysis[i][1] = String.format("%.2f", threshByLayerU[i]/threshByLayerWeightU[i]);
				//System.out.print(threshByLayerWeight[i]+",");
			}
			//System.out.println(threshByLayerWeight[index-1]);
			//analysis[index-1] = threshByLayer[index-1]/threshByLayerWeight[index-1];
			for(int i=1; i<analysis.length-1; i++) {
				System.out.print("["+this.threshByLayerWeightoccL[i]+"; "+this.threshByLayerWeightoccU[i]+"]"+",");
			}
			System.out.println();
			System.out.print(",vertex range,");
			for(int i=1; i<analysis.length-1; i++) {
				System.out.print("["+analysis[i][0]+"; "+analysis[i][1]+"]"+",");
			}
			
			System.out.println();
			System.out.print(",edge range,");
			String[][] analysis2 = new String[k][2];
			for(int i=0; i<k; i++) {
				analysis2[i][0] = String.format("%.2f", threshByLayerLE[i]/threshByLayerWeightLE[i]);
				analysis2[i][1] = String.format("%.2f", threshByLayerUE[i]/threshByLayerWeightUE[i]);
				//System.out.print(threshByLayerWeight[i]+",");
			}
			//System.out.println(threshByLayerWeight[index-1]);
			//analysis[index-1] = threshByLayer[index-1]/threshByLayerWeight[index-1];
			
			
			for(int i=1; i<analysis2.length-1; i++) {
				System.out.print("["+analysis2[i][0]+"; "+analysis2[i][1]+"]"+",");
			}
			//System.out.println(analysis[analysis.length-1]);
			System.out.println();

		}
		
//		System.out.println(this.removingTime);
//		System.out.println(this.compTime);
		
		long[] res = new long[graph.length+3];
		res[0] = k;
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
	
	/*
	 * mapping = vertices remaining
	 */
	public Boolean elimTree(BitSet[] graph, BitSet mapping, int k, int parent) {
		this.elimTreesSearched++;
		if(mapping.cardinality() == 1) {
			parentList[mapping.nextSetBit(0)] = parent;
			return true;
		} else {
			//OPTIMIZATION only-child rule
//			int mapS = mapping.cardinality();
//			int nex=mapping.nextSetBit(0);
//			for(int x=0; x<mapS; x++) {
//				if(graph[nex].cardinality() < 2) {
//					if(graph[nex].cardinality())
//				}
//				nex=mapping.nextSetBit(nex);
//			}
		}
		
//		for(int x=0; x<graph.length; x++) {
//			System.out.println(x+" : "+graph[x]);
//		}
//		System.out.println("mapping "+mapping.size());
//		System.out.println("mapping "+mapping);
		
		BitSet mapping2 = (BitSet) mapping.clone();
		// pruning vertices by domination rule
		if(this.dom_rule) {
			//System.out.println("dom rule ative");
			int nextx = mapping.nextSetBit(0);
			for(int x=0; x<mapping.cardinality(); x++) {
				int nexty = mapping.nextSetBit(nextx+1);
//				System.out.println(nextx+","+nexty);
//				System.out.println(graph[nextx]+","+graph[nexty]);
				for(int y=x+1; y<mapping.cardinality(); y++) {
//					if(nextx != nexty) {
					BitSet tempx = (BitSet) graph[nextx].clone();
					tempx.clear(nexty+1);
					BitSet tempy = (BitSet) graph[nexty].clone();
					tempy.clear(nextx+1);
					int maxN = Math.max(tempx.cardinality(), tempy.cardinality());
					tempx.or(tempy);
//						System.out.println(graph[nextx]+","+graph[nexty]);
					if(tempx.cardinality() == maxN) {
//						System.out.println("dom rule hapenned for:"+nextx+","+nexty);
//						System.out.println(graph[nextx]+","+graph[nexty]);
//						System.out.println(tempx+","+tempy);
						// domination rule applies
						if(tempy.cardinality() == maxN) {
							mapping.clear(nextx);
						} else {
							mapping.clear(nexty);
						}
//							System.exit(0);
					}
//					}
					nexty = mapping.nextSetBit(nexty+1);
				}
				nextx = mapping.nextSetBit(nextx+1);
			}
		}
		
		int next = mapping.nextSetBit(0);
		
		// reorder = order vertex by degree
		if(reorder && mapping.cardinality() > 2) {
			ArrayList<Integer> ord = new ArrayList<Integer>();
			while(next != -1) {
				if(graph[next].cardinality() > 1) {
					int ind = 0;
					while(ind < ord.size() && graph[ord.get(ind)].cardinality() > graph[next].cardinality()) {ind++;}
					ord.add(ind, next);
				}
				next = mapping.nextSetBit(next+1);
			}
			for(int i=0; i<ord.size(); i++) {
				BitSet[] newGraph = removeVertex(graph, ord.get(i));
				BitSet newMapping = (BitSet) mapping2.clone();
				newMapping.clear(ord.get(i));
				if(elimForest(newGraph, newMapping, k-1, ord.get(i)+1)) {return true;}
			}
		}
		else if(!reorder && mapping.cardinality() > 2) {
			int[] inds = new int[mapping.cardinality()];
			int spot = 0;
			while(next != -1) {
				inds[spot] = next;
				next = mapping.nextSetBit(next+1);
				spot +=1;
			}
			if(!this.reorderstart) {
					
				// Fisher-Yates shuffle
				int ind;
			    Random rand = new Random();
			    for (int i=inds.length-1; i>0; i--) {
			        ind = rand.nextInt(i+1);
			        if(ind != i) {
			        	inds[ind] = inds[ind]^inds[i];
			            inds[i] = inds[i]^inds[ind];
			            inds[ind] = inds[ind]^inds[i];
			        }
			    }
			}
			for(int i=0; i<inds.length; i++) {
				next = inds[i];
				if(graph[next].cardinality() > 1) {
					parentList[next] = parent;
					BitSet[] newGraph = removeVertex(graph, next);

					BitSet newMapping = (BitSet) mapping2.clone();
					newMapping.clear(next);
					if(elimForest(newGraph, newMapping, k-1, next+1)) {return true;}
				}
			}
			

			
		} else {
			while(next != -1) {
				parentList[next] = parent;
				BitSet[] newGraph = removeVertex(graph, next);
				
				BitSet newMapping = (BitSet) mapping2.clone();
				newMapping.clear(next);
				//System.out.println("after removed"+next+newGraph[next]+newMapping);
				if(elimForest(newGraph, newMapping, k-1, next+1)) {return true;}
				next = mapping.nextSetBit(next+1);
			}
		}
		return false;
	}
	
	public int simple_lb(int vertices) {
		if(vertices == 0) {return 0;}
		return (1+simple_lb((vertices-1)/b + (((vertices-1)%b == 0) ? 0 : 1)));
	}
	
	public boolean path_lb(BitSet[] comp, BitSet mapping, int k) {
		// if the component has 2 or less vertices, the path size is equal to amount of vertices
		if(mapping.size() < 3) {return (k < mapping.size());}
		
		BitSet examineable = (BitSet) mapping.clone();
		int v = examineable.nextSetBit(0);
		examineable.clear(v);
		BitSet p = new BitSet();
		p.set(v);
		
		for(int x=0; x<2; x++) {
			int u = v;
			// there's an unexplored neighbour of u
			while(examineable.intersects(comp[u])) {
				int newu = comp[u].nextSetBit(0);
				while(!examineable.get(newu)) {
					newu = comp[u].nextSetBit(newu+1);
				}
				u = newu;
				p.set(newu);
				examineable.clear(newu);
			}
		}
	    int result = (int) (Math.ceil(Math.log(p.cardinality()+1) / Math.log(2)));
	    return (result > k);
	}

	
	public Boolean elimForest(BitSet[] graph, BitSet mapping, int k, int parent) {
		this.elimForestSearched++;
		
		//System.out.println("Forest searched:"+elimForestSearched);
		if(k == 0 && mapping.cardinality() > 0) {return false;}
				
		// make connected components
		long time1 = System.nanoTime();
		ArrayList<BitSet[]> components = generateComponents(graph, mapping);
		long time2 = System.nanoTime();
		this.compTime += (time2-time1);
		//System.out.println("many components:"+components.size());
		
		
		// pruning by simple lower bound
		if(this.simple_lb) {
			for(int i=0; i<components.size()-1; i++) {
				if(simple_lb(components.get(0)[i].cardinality()) > k) {return false;}
			}
		}
		
		// pruning by path-based lower bound
		if(this.path_lb) {
			for(int i=0; i<components.size()-1; i++) {
				if(path_lb(components.get(i+1), components.get(0)[i], k)) {return false;}
			}
		}
		
		
		// PRUNING METHOD: order by size to find false components, big component first
		if(pruningType == "closestSize") {
			int lim = components.size()-1;
			for(int i=0; i<lim; i++) {
				int closest = Integer.MAX_VALUE;
				int closestInd = -1;
				for(int j=1; j<components.size(); j++) {
					int card = components.get(0)[j-1].cardinality();
					if(card != 0 && Math.abs(card - (int) param) < closest) {
						
						closest = Math.abs(card - (int) param);
						closestInd = j;
					}
				}
				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {return false;}
				components.get(0)[closestInd-1].clear();
			}
		}
		// PRUNING METHOD: consider 
		else if(pruningType == "functOfDepth") {
			int lim = components.size()-1;
			double threshold = k * this.param;
			for(int i=0; i<lim; i++) {
				double closest = Double.MAX_VALUE;
				int closestInd = -1;
				for(int j=1; j<components.size(); j++) {
					int card = components.get(0)[j-1].cardinality();
					if(card != 0 && card > threshold) {
						
						closestInd = j;
						break;
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
		// PRUNING METHOD: consider 
		else if(pruningType == "analysis") {
			int lim = components.size()-1;
			int biggestValid = -1;
			for(int i=0; i<lim; i++) {
				int closestInd = i+1;
				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {
					int val = components.get(0)[closestInd-1].cardinality();
					double weight = (double) Math.sqrt((double) val);
					this.threshByLayerU[k-1] += val*weight;
					this.threshByLayerWeightU[k-1] += weight;
					this.threshByLayerWeightoccU[k-1]++;
					
					int n=-1;
					double edg = 0;
					for(int j=0; j<components.get(closestInd).length; j++) {
						// error to fix
						if(components.get(closestInd)[j] != null) {
							edg += components.get(closestInd)[j].cardinality();
						}
					}
					edg = edg/2.0;
					
					double weightE = (double) Math.sqrt(edg);
					this.threshByLayerUE[k-1] += edg*weightE;
					this.threshByLayerWeightUE[k-1] += weightE;
					this.threshByLayerWeightoccUE[k-1]++;
					return false;
				} else {
					int val = components.get(0)[closestInd-1].cardinality();
					double weight = (((double) val)*((double) val));
					this.threshByLayerL[k-1] += val*weight;
					this.threshByLayerWeightL[k-1] += weight;
					this.threshByLayerWeightoccL[k-1]++;
					
					int n=-1;
					double edg = 0;
					for(int j=0; j<components.get(closestInd).length; j++) {
						// error to fix
						if(components.get(closestInd)[j] != null) {
							edg += components.get(closestInd)[j].cardinality();
						}
					}
					edg = edg/2.0;
					
					double weightE = (double) Math.sqrt(edg);
					this.threshByLayerLE[k-1] += edg*weightE;
					this.threshByLayerWeightLE[k-1] += weightE;
					this.threshByLayerWeightoccLE[k-1]++;
				}
				components.get(0)[closestInd-1].clear();
			};
		}
		// PRUNING METHOD: consider 
		else if(pruningType == "adaptive") {
			
			int lim = components.size()-1;
			int biggestValid = -1;
			for(int i=0; i<lim; i++) {
				int closestInd = i+1;
				if(this.threshByLayerWeightoccU[k-1] > 0 && this.threshByLayerWeightoccUE[k-1] > 0) {
					double closest = Double.MAX_VALUE;
					double LB = this.threshByLayerL[k-1]/this.threshByLayerWeightL[k-1];
					double UB = this.threshByLayerU[k-1]/this.threshByLayerWeightU[k-1];
					double LBE = this.threshByLayerLE[k-1]/this.threshByLayerWeightLE[k-1];
					double UBE = this.threshByLayerUE[k-1]/this.threshByLayerWeightUE[k-1];
					double treshV = LB + ((UB-LB)*this.shV);
					double treshE = LBE + ((UBE-LBE)*this.shE);
					double edg = 0;
					for(int j=0; j<components.get(closestInd).length; j++) {
						// error to fix
						if(components.get(closestInd)[j] != null) {
							edg += components.get(closestInd)[j].cardinality();
						}
					}
					edg = edg/2.0;
					
					for(int j=1; j<components.size(); j++) {
						int card = components.get(0)[j-1].cardinality();
						if(card != 0 && (Math.abs(card - treshV)*param +  Math.abs(edg - treshE)*(param-1))< closest) {
							
							closest = (Math.abs(card - treshV)*param +  Math.abs(edg - treshE)*(param-1));
							closestInd = j;
						}
					}
				}
				//System.out.println(closestInd);
				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {
					int val = components.get(0)[closestInd-1].cardinality();
					double weight = (double) Math.sqrt((double) val);
					this.threshByLayerU[k-1] += val*weight;
					this.threshByLayerWeightU[k-1] += weight;
					this.threshByLayerWeightoccU[k-1]++;
					
					int n=-1;
					double edg = 0;
					for(int j=0; j<components.get(closestInd).length; j++) {
						// error to fix
						if(components.get(closestInd)[j] != null) {
							edg += components.get(closestInd)[j].cardinality();
						}
					}
					edg = edg/2.0;
					
					double weightE = (double) Math.sqrt(edg);
					this.threshByLayerUE[k-1] += edg*weightE;
					this.threshByLayerWeightUE[k-1] += weightE;
					this.threshByLayerWeightoccUE[k-1]++;
					return false;
				} else {
					int val = components.get(0)[closestInd-1].cardinality();
					double weight = (((double) val)*((double) val));
					this.threshByLayerL[k-1] += val*weight;
					this.threshByLayerWeightL[k-1] += weight;
					this.threshByLayerWeightoccL[k-1]++;
					
					int n=-1;
					double edg = 0;
					for(int j=0; j<components.get(closestInd).length; j++) {
						// error to fix
						if(components.get(closestInd)[j] != null) {
							edg += components.get(closestInd)[j].cardinality();
						}
					}
					edg = edg/2.0;
					
					double weightE = (double) Math.sqrt(edg);
					this.threshByLayerLE[k-1] += edg*weightE;
					this.threshByLayerWeightLE[k-1] += weightE;
					this.threshByLayerWeightoccLE[k-1]++;
				}
				components.get(0)[closestInd-1].clear();
			};
//		}
//		// PRUNING METHOD: consider 
//		else if(pruningType == "adaptive") {
//			int lim = components.size()-1;
//			int biggestValid = -1;
//			for(int i=0; i<lim; i++) {
//				int closest = Integer.MAX_VALUE;
//				int closestInd = -1;
//				// If sufficient weight has been gathered
//				//TODO: 1.0 changed to param
//				//TODO: implement with only above threshold
//				if (threshByLayerWeight[k-1] > this.param) {
//					double thresh = threshByLayer[k-1]/threshByLayerWeight[k-1];
//					for(int j=1; j<components.size(); j++) {
//						int card = components.get(0)[j-1].cardinality();
//						if(card != 0 && Math.abs(card - thresh) < closest) {
//							closest = card;
//							closestInd = j;
//						}
//					}
//				} else {
//					//OTHERWISE continue to train
//					for(int j=1; j<components.size(); j++) {
//						int card = components.get(0)[j-1].cardinality();
//						if(card != 0 && card < closest) {
//							closest = card;
//							closestInd = j;
//						}
//					}
//				}
//				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {
//					if(biggestValid > -1) {
//						// adjust weights
//						double weight = 1.0/(1.0 + (double) ((closest - biggestValid)^2));
//						this.threshByLayer[k-1] += biggestValid*weight;
//						this.threshByLayerWeight[k-1] += weight;
//					}
//					return false;
//				} else {
//					biggestValid = components.get(0)[closestInd-1].cardinality();
//				}
//				components.get(0)[closestInd-1].clear();
//			};
//		} else if(pruningType == "adaptiveE") {
//			int lim = components.size()-1;
//			int biggestValid = -1;
//			int biggestValidE = -1;
//			for(int i=0; i<lim; i++) {
//				int closest = Integer.MAX_VALUE;
//				int closestE = Integer.MAX_VALUE;
//				double bestVal = Double.MAX_VALUE;
//				int closestInd = -1;
//				// If sufficient weight has been gathered
//				//TODO: 1.0 changed to param
//				//TODO: implement with only above threshold
//				if (this.threshByLayerWeightE[k-1] > this.param) {
//					double threshE = this.threshByLayerE[k-1]/this.threshByLayerWeightE[k-1];
//					for(int j=1; j<components.size(); j++) {
//						int card = components.get(0)[j-1].cardinality();
//						int connectivity = 0;
//						for(int x=0; x<components.get(j).length; x++) {
//							if(components.get(j)[x] != null) {
//								connectivity += components.get(j)[x].cardinality();
//							}
//						}
//						if((double) Math.abs(connectivity - threshE) < bestVal) {
//							
//							closestE = connectivity;
//							closestInd = j;
//						}
//					}
//				} else {
//					//OTHERWISE continue to train
//					for(int j=1; j<components.size(); j++) {
//						int card = components.get(0)[j-1].cardinality();
//						if(card != 0 && card < closest) {
//							closest = card;
//							closestInd = j;
//						}
//					}
//					int connectivity = 0;
//					for(int x=0; x<components.get(closestInd).length; x++) {
//						if(components.get(closestInd)[x] != null) {
//							connectivity += components.get(closestInd)[x].cardinality();
//						}
//					}
//				}
//				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {
//					if(biggestValidE > -1) {
//						// adjust weights
//						double xE = ((double) (closestE - biggestValidE))/this.edge_amount;
//						double weightE = 1.0/(1.0 + xE);
//						this.threshByLayerE[k-1] += biggestValidE*weightE;
//						this.threshByLayerWeightE[k-1] += weightE;
//					}
//					return false;	
//				} else {
//					int connectivity = 0;
//					for(int x=0; x<components.get(closestInd).length; x++) {
//						if(components.get(closestInd)[x] != null) {
//							connectivity += components.get(closestInd)[x].cardinality();
//						}
//					}
//					biggestValidE = connectivity;
//				}
//				components.get(0)[closestInd-1].clear();
//			};
//		}
//		else if(pruningType == "adaptiveMixed") {
//			int lim = components.size()-1;
//			int biggestValid = -1;
//			int biggestValidE = -1;
//			for(int i=0; i<lim; i++) {
//				int closest = Integer.MAX_VALUE;
//				int closestE = Integer.MAX_VALUE;
//				double bestVal = Double.MAX_VALUE;
//				int closestInd = -1;
//				// If sufficient weight has been gathered
//				//TODO: 1.0 changed to param
//				//TODO: implement with only above threshold
//				if (this.threshByLayerWeight[k-1] > this.param && this.threshByLayerWeightE[k-1] > this.param) {
//					double thresh = threshByLayer[k-1]/threshByLayerWeight[k-1];
//					double threshE = this.threshByLayerE[k-1]/this.threshByLayerWeightE[k-1];
//					for(int j=1; j<components.size(); j++) {
//						int card = components.get(0)[j-1].cardinality();
//						int connectivity = 0;
//						for(int x=0; x<components.get(j).length; x++) {
//							if(components.get(j)[x] != null) {
//								connectivity += components.get(j)[x].cardinality();
//							}
//						}
//						if((double) Math.abs(connectivity - threshE)/this.edge_amount + (double) Math.abs(card - thresh)/this.node_amount < bestVal) {
//							
//							closestE = connectivity;
//							closest = card;
//							closestInd = j;
//						}
//					}
//				} else {
//					//OTHERWISE continue to train
//					for(int j=1; j<components.size(); j++) {
//						int card = components.get(0)[j-1].cardinality();
//						if(card != 0 && card < closest) {
//							closest = card;
//							closestInd = j;
//						}
//					}
//					int connectivity = 0;
//					for(int x=0; x<components.get(closestInd).length; x++) {
//						if(components.get(closestInd)[x] != null) {
//							connectivity += components.get(closestInd)[x].cardinality();
//						}
//					}
//				}
//				if(!elimTree(components.get(closestInd), components.get(0)[closestInd-1], k, parent)) {
//					if(biggestValid > -1) {
//						// adjust weights
//						double x = ((double) (closest - biggestValid))/this.node_amount;
//						double weight = 1.0/(1.0 + x);
//						this.threshByLayer[k-1] += biggestValid*weight;
//						this.threshByLayerWeight[k-1] += weight;
//						double xE = ((double) (closestE - biggestValidE))/this.edge_amount;
//						double weightE = 1.0/(1.0 + xE);
//						this.threshByLayerE[k-1] += biggestValidE*weightE;
//						this.threshByLayerWeightE[k-1] += weightE;
//					}
//					return false;	
//				} else {
//					biggestValid = components.get(0)[closestInd-1].cardinality();
//					int connectivity = 0;
//					for(int x=0; x<components.get(closestInd).length; x++) {
//						if(components.get(closestInd)[x] != null) {
//							connectivity += components.get(closestInd)[x].cardinality();
//						}
//					}
//					biggestValidE = connectivity;
//				}
//				components.get(0)[closestInd-1].clear();
//			};
		} else {
			// Fisher-Yates shuffle
			int[] inds = new int[components.size()-1];
			for(int y=0; y<inds.length; y++) {
				inds[y] = y+1;
			}
			int ind;
		    Random rand = new Random();
		    for (int i=inds.length-1; i>0; i--) {
		        ind = rand.nextInt(i+1);
		        if(ind != i) {
		        	inds[ind] = inds[ind]^inds[i];
		            inds[i] = inds[i]^inds[ind];
		            inds[ind] = inds[ind]^inds[i];
		        }
		    }
		    		    
			for(int j=0; j<inds.length; j++) {
				int i = inds[j];
				//System.out.println("new component of"+i+"  with mapping:"+components.get(0)[i-1]);
				if(!elimTree(components.get(i), components.get(0)[i-1], k, parent)) {return false;}
			}
		}
//		for(int i=1; i<components.size(); i++) {
//			//System.out.println("new component of"+i+"  with mapping:"+components.get(0)[i-1]);
//			if(!elimTree(components.get(i), components.get(0)[i-1], k, parent)) {return false;}
//		}
		return true;
	}
	
	public int[] get_params() {
		int[] res = {this.e, this.v};
		return res;
	}
	
}
