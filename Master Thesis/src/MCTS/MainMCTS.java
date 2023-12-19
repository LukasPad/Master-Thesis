/**
 * 
 */
package MCTS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;

/**
 * @author FPC
 *
 */
public class MainMCTS {

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub

		int[] filenums = {5,15,5,21,9,1,13,23,4,6,8,12,22};
//		test3.set(6);
		String[] types = {"normal", "closestSize", "furthestSize","functOfDepth","closestFOD","closestOverFOD","analysis","adaptive","adaptiveE","adaptiveMixed"};
		//auto generate
		//int[] filenums = {19,25,7,3,17};
		File[] par1 = new File[filenums.length];
//		
		for(int i=0; i<filenums.length; i++) {
			par1[i] = new File("src/input_graphs/exact_"+String.format("%03d", filenums[i])+".gr");
		}
		MCTSPrimer mcts = new MCTSPrimer();
		System.out.println("starting on: "+par1[0]);
		long[] res = mcts.run(par1[0]);
		System.out.println("treedepth gotten: "+res[0]);

	}
	

}
