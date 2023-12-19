/**
 * 
 */
package main;

import java.io.IOException;
import java.util.BitSet;
import java.io.File;

/**
 * @author Lukas Padolevicius
 *
 */
public class Main {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		int[] filenums = {11,15,5,21,9,1,13,23,4,6,8,12,22};
//		test3.set(6);
		String[] types = {"normal", "closestSize", "furthestSize","functOfDepth","closestFOD","closestOverFOD","analysis","adaptive","adaptiveE","adaptiveMixed"};
		//auto generate
		//int[] filenums = {19,25,7,3,17};
		File[] par1 = new File[filenums.length];
//		
		for(int i=0; i<filenums.length; i++) {
			par1[i] = new File("src/input_graphs/exact_"+String.format("%03d", filenums[i])+".gr");
		}
		int st = 0;

		
		File folder = new File("src/famous_graphs");
		File[] listOfFiles = folder.listFiles();
		File[] par2 = new File[listOfFiles.length];
		for(int i=0; i<listOfFiles.length; i++) {
			par2[i] = listOfFiles[i];
		}
		File[] par = new File[par1.length+par2.length];
		for(int i=0 ;i<par1.length; i++) {
			par[i] = par1[i];
		}
		for(int i=par1.length ;i<par.length; i++) {
			par[i] = par2[i-par1.length];
		}
 		
		
		
		/*
		 * Settings = 
		 * 0 -> re-order vertices at sub-problems
		 * 1 -> re-order vertices at start
		 * 2 -> prune with simple lower bound
		 * 3 -> prune with path-based lower bound
		 * 
		 * 
		 */
		int[] typeExp = {7,7,7,7,7,7,7,7,7,7,7};
		double[] paramsExp = {1.6,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0};
		boolean experiments = true;
		boolean primers = false;
		boolean analysis = false;
		int iters = 10;
	
				
		//String[] par = {"src/input_graphs/exact_005.gr"};

		for(int i=0; i<par.length; i++) {
			
			
			if(experiments) {
				Inscriber tempInc = new Inscriber();
				System.out.print(par[i]);
				for(int y=0; y<1; y++) {
					long[] tots = new long[3];
					for(int g=0; g<tots.length; g++) {
						tots[g] = 0;
					}
					BitSet[] temp = tempInc.inscribe(par[i]);
					//System.out.print(par[i]+","+types[typeExp[y]]+","+tempInc.get_params()[0]+","+tempInc.get_params()[1]);
					//System.out.print(par[i]+",occurences,");
					
					if(analysis) {
						int e = tempInc.get_params()[0];
						int v = tempInc.get_params()[1];
						
						float avgdeg = ((float) e)/ ((float) v);
						float popvar = 0;
						for(int x=0; x<v; x++) {
							//System.out.print(temp[x]);
							popvar += (temp[x].cardinality() - avgdeg)*(temp[x].cardinality() - avgdeg);
						}
						int cliqueE = (v*(v-1))/2;
						float cliqueprop = ((float) e)/((float) cliqueE);
						popvar = popvar/((float) v);
						System.out.print(","+avgdeg+","+cliqueprop+","+popvar);
						
						
					}
					
					
					
					
					for(int j=0; j<iters; j++) {
						int shift = 0;
						boolean[] settings = {true, true, true, true, true};
						//boolean[] settings = {false, false, false, false, false};
						
						long time1 = System.nanoTime();
						Algorithm alg = new Algorithm();
						long[] res = alg.run(par[i], types[7], 1, settings, paramsExp[y], paramsExp[y]);
						long time2 = System.nanoTime();
						//System.out.print(res[0]+",");
						tots[shift*3] += time2-time1;
						tots[shift*3 + 1] += res[1];
						tots[shift*3 + 2] += res[2];
						shift++;
						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
						alg = null;
						res = null;
////						settings[0] = true;
//						time1 = System.nanoTime();
//						alg = new Algorithm();
//						res = alg.run(par[i], types[7], 0, settings, 0.5, 0.5);
//						time2 = System.nanoTime();
//						tots[shift*3] += time2-time1;
//						tots[shift*3 + 1] += res[1];
//						tots[shift*3 + 2] += res[2];
//						shift++;
//						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
//						alg = null;
//						res = null;
////						settings[0] = false;
////						settings[1] = true;
//						alg = new Algorithm();
//						res = alg.run(par[i], types[7], 1, settings, 0.5, 0.5);
//						time2 = System.nanoTime();
//						tots[shift*3] += time2-time1;
//						tots[shift*3 + 1] += res[1];
//						tots[shift*3 + 2] += res[2];
//						shift++;
//						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
//						alg = null;
//						res = null;
//						settings[1] = false;
//						settings[2] = true;
//						time1 = System.nanoTime();
//						alg = new Algorithm();
//						res = alg.run(par[i], types[typeExp[y]], paramsExp[y], settings);
//						time2 = System.nanoTime();
//						tots[shift*3] += time2-time1;
//						tots[shift*3 + 1] += res[1];
//						tots[shift*3 + 2] += res[2];
//						shift++;
//						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
//						alg = null;
//						res = null;
//						settings[0] = true;
//						settings[1] = true;
//						settings[2] = false;
//						time1 = System.nanoTime();
//						alg = new Algorithm();
//						res = alg.run(par[i], types[typeExp[y]], paramsExp[y], settings);
//						time2 = System.nanoTime();
//						tots[shift*3] += time2-time1;
//						tots[shift*3 + 1] += res[1];
//						tots[shift*3 + 2] += res[2];
//						shift++;
//						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
//						alg = null;
//						res = null;
//						settings[3] = true;
//						settings[2] = false;
//						time1 = System.nanoTime();
//						alg = new Algorithm();
//						res = alg.run(par[i], types[typeExp[y]], paramsExp[y], settings);
//						time2 = System.nanoTime();
//						tots[shift*3] += time2-time1;
//						tots[shift*3 + 1] += res[1];
//						tots[shift*3 + 2] += res[2];
//						shift++;
//						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
//						alg = null;
//						res = null;
//						settings[4] = true;
//						settings[3] = false;
//						time1 = System.nanoTime();
//						alg = new Algorithm();
//						res = alg.run(par[i], types[typeExp[y]], paramsExp[y], settings);
//						time2 = System.nanoTime();
//						tots[shift*3] += time2-time1;
//						tots[shift*3 + 1] += res[1];
//						tots[shift*3 + 2] += res[2];
//						shift++;
//						//System.out.print(time2-time1+","+res[1]+","+res[2]+",");
//						alg = null;
//						res = null;
//						settings[4] = true;
					}
					
					for(int x=0; x<tots.length; x++) {
						System.out.print(","+(long) (((double) tots[x])/(double) iters));
					}
					
					//System.out.println();
				}
				System.out.println();
			}
			
			if(primers) {
				GeneralPrimer prim = new GeneralPrimer();
				//long[] res = prim.run(par[i]);
				
			}
		}
	}
	//		System.out.println("Time taken: "+(time2-time1));
	//		System.out.println("treedepth: "+res[0]);
	//		System.out.println("Elimination forest searched: "+res[1]);
	//		System.out.println("Elimination trees searched: "+res[2]);

}


//BitSet test1 = new BitSet();
//BitSet test2 = new BitSet();
//BitSet test3 = new BitSet();
//test1.set(2);
//test1.set(3);
//test1.set(5);
//test1.set(6);
//test2.set(2);
//test2.set(4);
//test2.set(5);
//test2.set(7);
////System.out.println(test2);
//test2.clear(2);
////System.out.println(test2);
////System.out.println(test2.nextSetBit(2));
////System.exit(0);
//test3.set(1);
//test3.set(5);
//test3.set(2);
//int[] filenums = {7,19,2,10,24,26};
//int[] filenums = {26,19,14,16,37,28};
//int[] filenums = {6,19};