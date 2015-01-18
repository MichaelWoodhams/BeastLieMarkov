package beast.evolution.substitutionmodel;

import java.util.HashSet;
import java.util.Set;

import liemarkov.BasisMatrix;
import beast.core.Description;
import beast.core.Citation;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.tree.Node;

/*
 * This class has quite a bit in common with, and copied from, beast.evolution.substitutionmodel.GeneralSubstitutionModel.
 * That class has assumptions of time reversibility (that rate matrix entries are the 
 * product of a relative rate and an equilibrium base frequency) which are not obeyed
 * by the Lie Markov models, which is why this isn't simply a subclass of GeneralSubstitutionModel.
 * Subclasses SubstitutionModel.Base because SiteModel requires input of this class. 
 * Then we have to change 'frequencies' input to forbidden as SubstitutionModel.Base includes it.
 * 
 * The class should be able to handle changes to 'model' mid-run, but this is untested.
 * 
 * Method setupLieMarkovRateMatrix() plus the BasisMatrix class comprise the reference implementation
 * of the distinguished pair Lie Markov models. If you're looking to reimplement them in your own
 * program, those are the bits you need.  
 */

@Citation("Woodhams, Fernández-Sánchez and Sumner, " + 
		"A New Hierarchy of Phylogenetic Models Consistent with Heterogeneous Substitution Rates, " +
        "submitted to Systematic Biology, 2015.")
@Description("Lie Markov DNA mutation models (with distinguished base pairings) " + 
		     "Woodhams, Fernández-Sánchez and Sumner, submitted to Systematic Biology, 2015.")
public class LieMarkovModel extends SubstitutionModel.Base {
    public Input<RealParameter> parametersInput =
            new Input<RealParameter>("parameters", "Parameters which define the transition rate matrix. " + 
            		"All parameters are in the range -1 to 1. The relationship between these parameters and the rate matrix is complex.", Validate.REQUIRED);
    public Input<String> modelInput = new Input<String>("model","The name of the Lie Markov model to use, e.g. '6.7a'", "3.3a", BasisMatrix.MODEL_LIST);
    // TODO: "RY" is supposed to be the default. Array is supposed to be list of valid inputs. Check that this works.
    public Input<String>  distinguishedPairInput = new Input<String>("distinguished","Which base pairing gets special treatment. Default = 'RY'" + 
    		"(i.e. transitions vs. transversions), other options are 'WS' or 'MK'.", "RY", new String[]{"RY","WS","MK"});
    public static final int STATE_COUNT = 4;

    private double[][] rateMatrix;
    private EigenSystem eigenSystem;
    private EigenDecomposition eigenDecomposition;
    private EigenDecomposition storedEigenDecomposition;
    private boolean storedUpdateMatrix = true;
    private double[] equilibriumFrequencies;
    private String lastModelInput;
    private int nDimensions;
    private BasisMatrix[] basis;

    /**
     * flag to indicate matrix needs to be updated
     */
    protected boolean updateMatrix = true;

	public LieMarkovModel() {
		super();
        frequenciesInput.setRule(Validate.FORBIDDEN); // can't take frequencies as input
	}

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        updateMatrix = true;
        eigenSystem = new DefaultEigenSystem(STATE_COUNT);
        rateMatrix = new double[STATE_COUNT][STATE_COUNT];
        equilibriumFrequencies = new double[STATE_COUNT];
        lastModelInput = modelInput.get();
        basis = BasisMatrix.MODEL_BASIS_MATRICES.get(lastModelInput);
        nDimensions = basis.length;
        parametersInput.get().setDimension(nDimensions);
    } // initAndValidate

    /*
     * The setupRateMatrix() method is where all the Lie Markov magic happens.
     * Everything aside from this method is Beast boilerplate code.
     */
    /*
     * As written, keeps no memory of the model or permutation, and checks them every time
     * this method is called. In the usual circumstance where the model and permutation are
     * hard-coded in the xml, this is unnecessary work. Perhaps there should be an alternative.
     */
    /*
     *  TODO: does updateMatrix get set false (so this gets called) in all necessary circumstances?
     *  I.e. changes in any of parametersInput, modelInput, distinguishedPairInput?
     */
    /*
     * Updates rateMatrix, eigenDecomposition, equilibriumFrequencies 
     */
    public void setupRateMatrix() {
    	if (!lastModelInput.equals(modelInput.get())) {
			try {
				initAndValidate();
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
    	}
    	BasisMatrix.Permutation perm = BasisMatrix.Permutation.fromString(distinguishedPairInput.get());
    	
    	double[] params = new double[nDimensions];
		for (int i=0; i<nDimensions; i++) {
			params[i] = parametersInput.get().getArrayValue(i);
		}
		
		setupLieMarkovRateMatrix(rateMatrix, params, perm);
		// rateMatrix is complete. Now we update eigenDecomposition and equilibriumFrequencies.
		
        eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix); // overwrites rateMatrix
        // Update equilibriumFrequencies
        double[] eValues = eigenDecomposition.getEigenValues();
        // One of the eigenvalues should be zero (up to truncation errors). It is not always #0, so need to search for it
        double[] ieValues = eigenDecomposition.getImEigenValues();
        double minVal = Double.MAX_VALUE;
        int zeroIndex = -1;
        if (ieValues==null) {
	        for (int i=0; i<STATE_COUNT; i++) {
	        	double mag = Math.abs(eValues[i]);
	        	if (mag < minVal) {
	        		minVal = mag;
	        		zeroIndex = i;
	        	}
	        }
        } else {
	        for (int i=0; i<STATE_COUNT; i++) {
	        	double mag = eValues[i]*eValues[i]+ieValues[i]+ieValues[i];
	        	if (mag < minVal) {
	        		minVal = mag;
	        		zeroIndex = i;
	        	}
	        }
        }
        // Paranoia sanity check:
        if (Math.abs(eValues[zeroIndex])>1e-12) {
        	throw new RuntimeException("Can't happen: rate matrix has no zero eigenvalue!");
        }
        double sum=0;
        minVal = Double.MAX_VALUE;
        double[] eVectors = eigenDecomposition.getEigenVectors();
        for (int i=0; i<STATE_COUNT; i++) {
        	equilibriumFrequencies[i]=eVectors[STATE_COUNT*i+zeroIndex];
        	sum += equilibriumFrequencies[i];
        }
        // Normalize:
        for (int i=0; i<STATE_COUNT; i++) {
        	equilibriumFrequencies[i] /= sum;
        	if (equilibriumFrequencies[i]<0) {
        		throw new RuntimeException("Can't happen: negative equilibrium frequency");
        	}
        }
	}
    
    /*
     * This method plus the BasisMatrix class comprise the reference implementation of
     * the Lie Markov models. 
     */
    private void setupLieMarkovRateMatrix(double[][] rateMatrix, double[] params, BasisMatrix.Permutation perm) {
		// Zero rateMatrix
		for (int i=0; i<STATE_COUNT; i++) {
			for (int j=0; j<STATE_COUNT; j++) {
				rateMatrix[i][j] = 0;
			}
		}
		
		// Add weighted basis matrices and find 'saturation' s, which is how close we are to boundary
		// of hypercube.
		double s=0; 
		for (int i=0; i<nDimensions; i++) {
			basis[i].applyToRates(rateMatrix, params[i], perm);
			s = Math.max(s, Math.abs(params[i]));
		}

		// find min off-diagonal (will always be negative, unless all off-diagonals are zero.)
		double min = Double.MAX_VALUE;
		for (int row=0; row<STATE_COUNT; row++) {
			for (int col=0; col<STATE_COUNT; col++) {
				if (col != row) min = Math.min(min,rateMatrix[row][col]);
			}
		}

		/*
		 *  If we call current matrix X, we want perturbation matrix P
		 *  P = k X 
		 *  with k chosen such that minOffDiagonal(P) = -1,
		 *  then final rate matrix is
		 *  Q = (A + s P)/3
		 *  where s = saturation value (found above) and A = Jukes Cantor rate matrix 
		 *  (and always-present basis matrix)
		 *  A = (-3  1  1  1)
		 *      ( 1 -3  1  1)
		 *      ( 1  1 -3  1)
		 *      ( 1  1  1 -3)
		 *     
		 */
		 // 'min' is always non-positive, r is always non-negative, hence scale is non-negative
		double scale = (min==0) ? 0 : -s/min/3; // min==0 iff r==0.
		for (int row=0; row<STATE_COUNT; row++) {
			for (int col=0; col<STATE_COUNT; col++) {				
				if (col == row) {
					// diagonal: A/3 contribution = -1
					rateMatrix[row][col]=rateMatrix[row][col]*scale-1;
				} else {
					// off-diagonal: A/3 contribuiton = 1/3
					rateMatrix[row][col]=rateMatrix[row][col]*scale+1./3;
				}
			}
		}
   }
    
    @Override
    public int getStateCount() {
    	return STATE_COUNT;
    }
 
    /*
     * I couldn't find an example in the standard Beast code which does other than this.
     * Lacking an example on how to flatten the rate matrix, I'll just stick with this.
     */
    @Override
    public double[] getRateMatrix(Node node) {
        return null;
    }
    
    @Override
    public double[] getFrequencies() {
    	return equilibriumFrequencies;
    }

	/*
	 * Based on beast.evolution.substitutionmodel.GeneralSubstitutionModel#getTransitionProbabilities(...)
	 */
    @Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        double distance = (fStartTime - fEndTime) * fRate;

        int i, j, k;
        double temp;

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupRateMatrix();
                updateMatrix = false;
            }
        }

        double[] iexp = new double[STATE_COUNT * STATE_COUNT];
        // Eigen vectors
        double[] Evec = eigenDecomposition.getEigenVectors();
        // inverse Eigen vectors
        double[] Ievc = eigenDecomposition.getInverseEigenVectors();
        // Eigen values
        double[] Eval = eigenDecomposition.getEigenValues();
        for (i = 0; i < STATE_COUNT; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (j = 0; j < STATE_COUNT; j++) {
                iexp[i * STATE_COUNT + j] = Ievc[i * STATE_COUNT + j] * temp;
            }
        }

        int u = 0;
        for (i = 0; i < STATE_COUNT; i++) {
            for (j = 0; j < STATE_COUNT; j++) {
                temp = 0.0;
                for (k = 0; k < STATE_COUNT; k++) {
                    temp += Evec[i * STATE_COUNT + k] * iexp[k * STATE_COUNT + j];
                }

                matrix[u] = Math.abs(temp);
                u++;
            }
        }
    } // getTransitionProbabilities


    /* 
     * For non-Nucleotide data, returns 'false' (correct Beast 2.2 behaviour)
     * rather than throw an exception (Beast 2.1 behaviour.) 
     */
    @Override
    public boolean canHandleDataType(DataType dataType) {
    	return (dataType instanceof Nucleotide);
    }
    
    /*
     * The Lie Markov models which can return complex diagonalization are
     * precisely those which include basis matrix "C".
     */
	public static final Set<String> COMPLEX_MODELS = new HashSet<String>();
	static {
		COMPLEX_MODELS.add("3.3b");
		COMPLEX_MODELS.add("4.5b");
		COMPLEX_MODELS.add("5.6a");
		COMPLEX_MODELS.add("6.6");
		COMPLEX_MODELS.add("6.7b");
		COMPLEX_MODELS.add("6.17b");
		COMPLEX_MODELS.add("8.10a");
		COMPLEX_MODELS.add("8.10b");
		COMPLEX_MODELS.add("9.20a");
		COMPLEX_MODELS.add("10.12");
		COMPLEX_MODELS.add("10.34");
		COMPLEX_MODELS.add("12.12");
	}
    @Override
    public boolean canReturnComplexDiagonalization() {
    	return COMPLEX_MODELS.contains(modelInput.get());
    }

    
   
    /*
     * Code from here on is largely copied from GeneralSubstitutionModel with little understanding.
     */
	
	/**
     * CalculationNode implementation follows *
     */
    @Override
    public void store() {
        storedUpdateMatrix = updateMatrix;
        storedEigenDecomposition = eigenDecomposition.copy();
        super.store();
    }

    /**
     * Restore the additional stored state
     */
    @Override
    public void restore() {
        updateMatrix = storedUpdateMatrix;
        // To restore all this stuff just swap the pointers...
//        double[] tmp1 = storedRelativeRates;
//        storedRelativeRates = relativeRates;
//        relativeRates = tmp1;
        EigenDecomposition tmp = storedEigenDecomposition;
        storedEigenDecomposition = eigenDecomposition;
        eigenDecomposition = tmp;
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        updateMatrix = true;
        return true;
    }


    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */
    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        synchronized (this) {
            if (updateMatrix) {
                setupRateMatrix();
                eigenDecomposition = eigenSystem.decomposeMatrix(rateMatrix);
                updateMatrix = false;
            }
        }
        return eigenDecomposition;
    }

}
