package liemarkov;

import java.util.HashMap;
import java.util.Map;

import beast.evolution.datatype.Nucleotide;

/**
 * BasisMatrix: For use by LieMarkovModel.
 * Primarily this class provides the required 12 basis matrices as static finals:
 * BA, BA1, BB, BC, BD, BD1, BE1, BE2, BF1, BF2, BG1, BG2.
 * 
 * A Permutation enum is used to determine which nucleotide pairs are distinguished.
 * 
 * Bringing these together is the applyToRates method.
 *  
 * @author woodhams
 *
 */
public class BasisMatrix {
	private final int[][] positive; // indices of cells to add to. n x 2. Will add 1 to each cell unless isA2 is true.
	private final int[][] negative; // indices of cells to subtract from. Will subtract 1 from each cell unless isA is true
	private final boolean isA2; // for basis matrix A2, positive entries have weight 2
	private String canonicalVariable; // The variable name canonically used with this basis matrix
	
	/*
	 * For all but one case, isA2 is false, so short cut constructor for this.
	 */
	private BasisMatrix(int[][] pos, int[][] neg, String var) {
		this(pos, neg, var, false);
	}
	
	private BasisMatrix(int[][] pos, int[][] neg, String var, boolean isA2) {
		positive = pos;
		this.isA2 = isA2;
		negative = neg;
		canonicalVariable = var;
	}
	/*
	 * IMPORTANT NOTE!
	 * All the 'pictures' below are in columns-sum-to-zero order.
	 * Beast works in rows-sum-to-zero, so these will get transposed in process of 
	 * building the rate matrix.
	 */
	private static final Nucleotide nucleotideType = new Nucleotide();

	// The row/column indices for the DNA basis:
	private static int A; // 0
	private static int G; // 2
	private static int C; // 1
	private static int T; // 3
	static {
		try {
			A = nucleotideType.char2state("A");
			G = nucleotideType.char2state("G");
			C = nucleotideType.char2state("C");
			T = nucleotideType.char2state("T");
		} catch (Exception e) {
			throw new RuntimeException("Can't happen. Something here is hugely messed up, best look at the code.");
		}
	}
	/*
	 * From here on, we work in AGCT index order, using the above variables to translate to Beast's index order.
	 * We do this because the matrices are much easier to understand in this order.
	 */
	
	public enum Permutation {
		RY (new int[]{A,G,C,T}),
		WS (new int[]{A,T,C,G}),
		MK (new int[]{A,C,G,T});
		
		protected final int[] permutation;
		Permutation(int[] perm) {
			this.permutation = perm;
		}
		public static Permutation fromString(String perm) {
			if (perm == null) return RY; // default to transitions vs transversions
			/* Java 7 code:
			switch (perm) {
			case "RY":
				return RY;
			case "WS":
				return WS;
			case "MK":
				return MK;
			default:
				throw new IllegalArgumentException("Unrecognized permutation name "+perm);
			}
			*/
			/* Java 6 code: */
			if (perm == null || perm.equals("RY")) { return RY; } 
			else if (perm.equals("WS")) { return WS; }
			else if (perm.equals("MK")) { return MK; }
			else {throw new IllegalArgumentException("Unrecognized permutation name "+perm);}
		}
	}
	
	/**
	 * This is the heart of the BasisMatrix class. Increment or decrement elements of 'rateStore'
	 * by 'weight', as determined by this basis matrix permuted by 'perm'.
	 * @param rateStore
	 * @param weight
	 * @param perm
	 */
	public void applyToRates(double[][] rateStore, double weight, Permutation perm) {
		// Note transposing coordinates in building rateStore 
		for (int[] coord : negative) {
			rateStore[perm.permutation[coord[1]]][perm.permutation[coord[0]]] -= weight;
		}
		if (isA2) weight *= 2;
		for (int[] coord : positive) {
			rateStore[perm.permutation[coord[1]]][perm.permutation[coord[0]]] += weight;
		}
	}
	
	// This is just for printing out a summary of the rate matrix in symbolic form.
	public void applyToStrings(StringBuffer[][] strings) {
		// Leave in column-sum order.
		String factor = isA2 ? "2" : "";
		for (int[] coord : positive) {
			strings[coord[0]][coord[1]].append('+').append(factor).append(canonicalVariable);
		}
		for (int[] coord : negative) {
			strings[coord[0]][coord[1]].append('-').append(canonicalVariable);
		}
	}

	/*
	 * Basis matrix A2
	 * 0#--
	 * #0--
	 * --0#
	 * --#0
	 * where '#' means +2
	 */
	public static final BasisMatrix BA2 = new BasisMatrix(
			new int[][]{{0,1},{1,0},{2,3},{3,2}},
			new int[][]{{0,2},{0,3},{1,2},{1,3},{2,0},{2,1},{3,0},{3,1}}, "a2", true);

	/*
	 * Basis matrix B
	 * 00+-
	 * 00-+
	 * +-00
	 * -+00
	 */
	public static final BasisMatrix BB = new BasisMatrix(
			new int[][]{{0,2},{1,3},{2,0},{3,1}},
			new int[][]{{0,3},{1,2},{2,1},{3,0}},
			"b");
	
	/*
	 * Basis matrix C
	 * 00+-
	 * 00-+
	 * -+00
	 * +-00
	 */
	public static final BasisMatrix BC = new BasisMatrix(
			new int[][]{{0,2},{1,3},{2,1},{3,0}},
			new int[][]{{0,3},{1,2},{2,0},{3,1}},
			"c");
	
	/*
	 * Basis matrix D
	 * ++++
	 * ++++
	 * ----
	 * ----
	 */
	public static final BasisMatrix BD = new BasisMatrix(
			new int[][]{{0,0},{0,1},{0,2},{0,3},{1,0},{1,1},{1,2},{1,3}},
			new int[][]{{2,0},{2,1},{2,2},{2,3},{3,0},{3,1},{3,2},{3,3}},
			"d");
	
	/*
	 * Basis matrix D1
	 * -+00
	 * +-00
	 * 00+-
	 * 00-+
	 */
	public static final BasisMatrix BD1 = new BasisMatrix(
			new int[][]{{0,1},{1,0},{2,2},{3,3}},
			new int[][]{{0,0},{1,1},{2,3},{3,2}},
			"d1");

	/*
	 * Basis matrix E1
	 * ++++
	 * ----
	 * 0000
	 * 0000
	 */
	public static final BasisMatrix BE1 = new BasisMatrix(
			new int[][]{{0,0},{0,1},{0,2},{0,3}},
			new int[][]{{1,0},{1,1},{1,2},{1,3}},
			"e1");
	
	/*
	 * Basis matrix E2
	 * --++
	 * ++--
	 * 0000
	 * 0000
	 */
	public static final BasisMatrix BE2 = new BasisMatrix(
			new int[][]{{2,0},{2,1},{2,2},{2,3}},
			new int[][]{{3,0},{3,1},{3,2},{3,3}},
			"e2");
	
	/*
	 * Basis matrix F1
	 * ++--
	 * --++
	 * 0000
	 * 0000
	 */
	public static final BasisMatrix BF1 = new BasisMatrix(
			new int[][]{{0,0},{0,1},{1,2},{1,3}},
			new int[][]{{1,0},{1,1},{0,2},{0,3}},
			"f1");

	/*
	 * Basis matrix F2
	 * 0000
	 * 0000
	 * ++--
	 * --++
	 */
	public static final BasisMatrix BF2 = new BasisMatrix(
			new int[][]{{2,0},{2,1},{3,2},{3,3}},
			new int[][]{{3,0},{3,1},{2,2},{2,3}},
			"f2");

	/*
	 * Basis matrix G1
	 * -+00
	 * -+00
	 * +-00
	 * +-00
	 */

	public static final BasisMatrix BG1 = new BasisMatrix(
			new int[][]{{0,1},{1,1},{2,0},{3,0}},
			new int[][]{{0,0},{1,0},{2,1},{3,1}},
			"g1");
	/*
	 * Basis matrix G2
	 * 00+-
	 * 00+-
	 * 00-+
	 * 00-+
	 */
	public static final BasisMatrix BG2 = new BasisMatrix(
			new int[][]{{0,2},{1,2},{2,3},{3,3}},
			new int[][]{{0,3},{1,3},{2,2},{3,2}},
			"g2");
	
	/*
	 * All models include the Jukes Cantor matrix as a basis matrix, this matrix is not listed
	 * and is implemented in a special way by the LieMarkovModel class.
	 * 
	 *  
	 */
	public static final Map<String,BasisMatrix[]> MODEL_BASIS_MATRICES = new HashMap<String,BasisMatrix[]>();
	public static final String[] MODEL_LIST;
	static {
		MODEL_BASIS_MATRICES.put("1.1",   new BasisMatrix[]{});
		MODEL_BASIS_MATRICES.put("2.2b",  new BasisMatrix[]{BA2});
		MODEL_BASIS_MATRICES.put("3.3a",  new BasisMatrix[]{BA2,BB});
		MODEL_BASIS_MATRICES.put("3.3b",  new BasisMatrix[]{BA2,BC});
		MODEL_BASIS_MATRICES.put("3.3c",  new BasisMatrix[]{BA2,BD1});
		MODEL_BASIS_MATRICES.put("3.4",   new BasisMatrix[]{BA2,BD});
		MODEL_BASIS_MATRICES.put("4.4a",  new BasisMatrix[]{BD, BE1,BE2});
		MODEL_BASIS_MATRICES.put("4.4b",  new BasisMatrix[]{BA2,BD, BD1});
		MODEL_BASIS_MATRICES.put("4.5a",  new BasisMatrix[]{BA2,BB, BD});
		MODEL_BASIS_MATRICES.put("4.5b",  new BasisMatrix[]{BA2,BC, BD});
		MODEL_BASIS_MATRICES.put("5.6a",  new BasisMatrix[]{BA2,BB, BC, BD1});
		MODEL_BASIS_MATRICES.put("5.6b",  new BasisMatrix[]{BA2,BD, BE1,BE2});
		MODEL_BASIS_MATRICES.put("5.7a",  new BasisMatrix[]{BA2,BB, BE1,BE2});
		MODEL_BASIS_MATRICES.put("5.7b",  new BasisMatrix[]{BA2,BB, BF1,BF2});
		MODEL_BASIS_MATRICES.put("5.7c",  new BasisMatrix[]{BA2,BB, BG1,BG2});
		MODEL_BASIS_MATRICES.put("5.11a", new BasisMatrix[]{BA2,BD1,BE1,BE2});
		MODEL_BASIS_MATRICES.put("5.11b", new BasisMatrix[]{BA2,BD1,BF1,BF2});
		MODEL_BASIS_MATRICES.put("5.11c", new BasisMatrix[]{BA2,BD1,BG1,BG2});
		MODEL_BASIS_MATRICES.put("5.16",  new BasisMatrix[]{BA2,BD, BG1,BG2});
		MODEL_BASIS_MATRICES.put("6.6",   new BasisMatrix[]{BA2,BB, BC, BD, BD1});
		MODEL_BASIS_MATRICES.put("6.7a",  new BasisMatrix[]{BA2,BB, BD, BE1,BE2});
		MODEL_BASIS_MATRICES.put("6.7b",  new BasisMatrix[]{BA2,BC, BD, BE1,BE2});
		MODEL_BASIS_MATRICES.put("6.8a",  new BasisMatrix[]{BA2,BD, BD1,BE1,BE2});
		MODEL_BASIS_MATRICES.put("6.8b",  new BasisMatrix[]{BA2,BD, BD1,BG1,BG2});
		MODEL_BASIS_MATRICES.put("6.17a", new BasisMatrix[]{BA2,BB, BD, BG1,BG2});
		MODEL_BASIS_MATRICES.put("6.17b", new BasisMatrix[]{BA2,BC, BD, BG1,BG2});
		MODEL_BASIS_MATRICES.put("8.8",   new BasisMatrix[]{BA2,BD, BD1,BE1,BE2,BF1,BF2});
		MODEL_BASIS_MATRICES.put("8.10a", new BasisMatrix[]{BA2,BB, BC, BD, BD1,BE1,BE2});
		MODEL_BASIS_MATRICES.put("8.10b", new BasisMatrix[]{BA2,BB, BC, BD, BD1,BG1,BG2});
		MODEL_BASIS_MATRICES.put("8.16",  new BasisMatrix[]{BA2,BD, BD1,BE1,BE2,BG1,BG2});
		MODEL_BASIS_MATRICES.put("8.17",  new BasisMatrix[]{BA2,BB, BD, BE1,BE2,BG1,BG2});
		MODEL_BASIS_MATRICES.put("8.18",  new BasisMatrix[]{BA2,BB, BD, BE1,BE2,BF1,BF2});
		MODEL_BASIS_MATRICES.put("9.20a", new BasisMatrix[]{BA2,BB, BC, BD1,BE1,BE2,BF1,BF2});
		MODEL_BASIS_MATRICES.put("9.20b", new BasisMatrix[]{BA2,BB, BC, BD1,BF1,BF2,BG1,BG2});
		MODEL_BASIS_MATRICES.put("10.12", new BasisMatrix[]{BA2,BB, BC, BD, BD1,BE1,BE2,BF1,BF2});
		MODEL_BASIS_MATRICES.put("10.34", new BasisMatrix[]{BA2,BB, BC, BD, BD1,BE1,BE2,BG1,BG2});
		MODEL_BASIS_MATRICES.put("12.12", new BasisMatrix[]{BA2,BB, BC, BD, BD1,BE1,BE2,BF1,BF2,BG1,BG2});
		MODEL_LIST = new String[MODEL_BASIS_MATRICES.size()];
		MODEL_BASIS_MATRICES.keySet().toArray(MODEL_LIST);
	}
}
