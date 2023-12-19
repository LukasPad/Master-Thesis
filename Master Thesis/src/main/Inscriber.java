/**
 * 
 */
package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;

/**
 * @author Lukas Padolevicius
 *
 */
public class Inscriber {
	
	private int e;
	private int v;
	
//	private String path = "src/input_graphs/exact_003.gr";
	
	public BitSet[] inscribe(File file) throws IOException {
//		File file = new File(path);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		// read parameters by PACE 2020 formatting
		String[] splitLine = reader.readLine().split(" ");
		while(!splitLine[0].equals("p")) {
			splitLine = reader.readLine().split(" ");
		}
		//System.out.println(splitLine[0]);
		this.v = Integer.parseInt(splitLine[2]);
		this.e = Integer.parseInt(splitLine[3]);
//		System.out.println("parameters: "+vertices+"  "+edges);
		BitSet[] result = new BitSet[this.v];
		for (int i=0; i < result.length; i++) {
			result[i] = new BitSet();
		}		
		splitLine = reader.readLine().split(" ");
		while (splitLine.length > 0) {
			if (!splitLine[0].equals("c")) {
				int shift = 0;
//				System.out.println(splitLine[0]);
//				System.out.println(splitLine[0].equals("e"));
				if(splitLine[0].equals("e")) {
					shift++;
//					System.out.println(splitLine[0]);
				}
				int vert1 = Integer.parseInt(splitLine[0+shift]);
				int vert2 = Integer.parseInt(splitLine[1+shift]);
//				System.out.println(" verts: "+vert1+"  "+vert2);
				result[vert1-1].set(vert2);
				result[vert2-1].set(vert1);
			}
			String read = reader.readLine();
			if(read == null) {
				break;
			}
			splitLine = read.split(" ");
		}
		reader.close();
//		for (int i=0; i < result.length; i++) {
//			System.out.println(result[i]);
//		}
		return result;
		
	}

	public ArrayList<BitSet> inscribe(File file, Boolean shiftTo0) throws IOException {
//		File file = new File(path);
		BufferedReader reader = new BufferedReader(new FileReader(file));
		// read parameters by PACE 2020 formatting
		String[] splitLine = reader.readLine().split(" ");
		while(!splitLine[0].equals("p")) {
			splitLine = reader.readLine().split(" ");
		}
		//System.out.println(splitLine[0]);
		this.v = Integer.parseInt(splitLine[2]);
		this.e = Integer.parseInt(splitLine[3]);
//		System.out.println("parameters: "+vertices+"  "+edges);
		BitSet[] result = new BitSet[this.v];
		for (int i=0; i < result.length; i++) {
			result[i] = new BitSet();
		}		
		splitLine = reader.readLine().split(" ");
		while (splitLine.length > 0) {
			if (!splitLine[0].equals("c")) {
				int shift = 0;
//				System.out.println(splitLine[0]);
//				System.out.println(splitLine[0].equals("e"));
				if(splitLine[0].equals("e")) {
					shift++;
//					System.out.println(splitLine[0]);
				}
				int vert1 = Integer.parseInt(splitLine[0+shift]);
				int vert2 = Integer.parseInt(splitLine[1+shift]);
//				System.out.println(" verts: "+vert1+"  "+vert2);
				result[vert1-1].set(vert2-1);
				result[vert2-1].set(vert1-1);
			}
			String read = reader.readLine();
			if(read == null) {
				break;
			}
			splitLine = read.split(" ");
		}
		reader.close();
//		for (int i=0; i < result.length; i++) {
//			System.out.println(result[i]);
//		}
		ArrayList<BitSet> convertedRes = new ArrayList<BitSet>();
		for (int i = 0; i < result.length; i++) {
			convertedRes.add(result[i]);
		}
		return convertedRes;
		
	}

	public int[] get_params() {
		int[] res = {this.e, this.v};
		return res;
	}
}
