/* LAPACKMatrixOperationsWrapper.java: Interface to numeric linear algebra computations like
 *                  matrix opns in CLAPACK.
 * Author: Mani Narayanan
 * Adapted to Java: Barbara Engelhardt (9/06)
 */

package sifter.components;

//import sifter_components.LAPACKVectorOperationsWrapper;
import org.netlib.util.*;
import org.netlib.lapack.Dgetrf;
import org.netlib.lapack.Dgetri;
//import org.netlib.lapack.Dsyevd;
//import org.netlib.lapack.Dgees;
import org.netlib.lapack.Dgeev;
import org.netlib.blas.Dgemm;


public class LAPACKMatrixOperationsWrapper {
	public double[] matrix;
	public int rows;
	public int cols;
	public boolean isDiagonal;
	private LAPACKMatrixOperationsWrapper ZStore;
	private LAPACKMatrixOperationsWrapper DStore;
	private LAPACKMatrixOperationsWrapper ZinvStore;
	private LAPACKVectorOperationsWrapper dStore;

	////////////////////////////////////////////
	// Constructors, setters
	////////////////////////////////////////////

	// Initializes diagonal matrix with vector elements
	// on diagonal
	public LAPACKMatrixOperationsWrapper(LAPACKVectorOperationsWrapper v) {
		matrix = new double[v.length()*v.length()];

		for (int i = 0; i < v.length(); i++) {
			for (int j = 0; j < v.length(); j++) {
				matrix[(i*v.length())+j] = 0.0;
			}

			matrix[(i*v.length())+i] = v.vector[i];
		}

		rows = v.length();
		cols = v.length();
		isDiagonal = true;
		ZStore = null;
		DStore = null;
		dStore = null;
	}

	// assumes square matrix
	public LAPACKMatrixOperationsWrapper(double[] m) {
		matrix = new double[m.length];
		rows = (int)Math.sqrt(m.length);
		cols = rows;

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[index(i, j)] = m[index(i, j)];
			}
		}

		ZStore = null;
		dStore = null;
	}

	// assumes square matrix
	public LAPACKMatrixOperationsWrapper(double[][] m) {
		rows = m.length;

		for (int i = 0; i < rows; i++) {
			matrix = new double[m[0].length*m.length];
			cols = m[0].length;

			for (int j = 0; j < cols; j++) {
				matrix[index(i, j)] = m[i][j];
			}
		}

		ZStore = null;
		DStore = null;
		dStore = null;
	}

	// Initializes matrix (r rows, c cols) with all zeros
	public LAPACKMatrixOperationsWrapper(int r, int c) {
		matrix = new double[r*c];

		for (int i = 0; i < r*c; i++)
			matrix[i] = 0.0;

		rows = r;
		cols = c;
		isDiagonal = true;
		ZStore = null;
		DStore = null;
		dStore = null;
	}

	public void set(int r, int c, double s) {
		if (r < 0 || c < 0) {
			System.out.print("Error: trying to access negative ");
			System.out.print("indexed rows/columns in LAPACKMatrixOperationsWrapper (" + r + "," + c + ")");
		}

		if (r >= rows || c >= cols) {
			System.out.print("Error: trying to set out of ");
			System.out.print("bounds rows/columns in LAPACKMatrixOperationsWrapper (" + r + "," + c + ")");
			System.out.println("when there are only (" + rows + "," + cols + ")");
			return;
		}

		matrix[index(r, c)] = s;

		// maintain diagonal setting
		//if(isDiagonal && r != c) {
		//    if(s != 0.0) isDiagonal = false;
		//} else if (r != c) {
		//    checkDiagonal();
		//	}
		ZStore = null;
		DStore = null;
		dStore = null;
	}


	public void add(int r, int c, double s) {
		if (r < 0 || c < 0) {
			System.out.print("Error: trying to access negative ");
			System.out.print("indexed rows/columns in LAPACKMatrixOperationsWrapper (" + r + "," + c + ")");
		}

		if (r >= rows || c >= cols) {
			System.out.print("Error: trying to add to out of ");
			System.out.print("bounds rows/columns in LAPACKMatrixOperationsWrapper (" + r + "," + c + ")");
			System.out.println("when there are only (" + rows + "," + cols + ")");
			return;
		}

		matrix[index(r, c)] += s;

		// maintain diagonal setting
		//if(isDiagonal && r != c) {
		//    if(s != 0.0) isDiagonal = false;
		//} else if (r != c) {
		//    checkDiagonal();
		//}
		ZStore = null;
		DStore = null;
		dStore = null;
	}

	public double get(int r, int c) {
		if (r < 0 || c < 0) {
			System.out.print("Error: trying to access negative ");
			System.out.print("indexed rows/columns in LAPACKMatrixOperationsWrapper (" + r + "," + c + ")");
		}

		if (r >= rows || c >= cols) {
			System.out.print("Error: trying to access out of ");
			System.out.print("bounds rows/columns in LAPACKMatrixOperationsWrapper (" + r + "," + c + ")");
			System.out.println("when there are only (" + rows + "," + cols + ")");
			return 0;
		}

		return matrix[index(r, c)];
	}

	/*
	private void checkDiagonal()
	{
	for(int i = 0; i < rows; i++) {
	  for(int j = 0; j < cols; j++) {
	if(i != j &&  get(i,j) != 0.0) {
	    isDiagonal = false;
	    return;
	}
	  }
	}
	isDiagonal = true;
	}*/

	// Compute eigen decomposition: S = Z D Zinv.
	// But S is symmetric, so Zinv = Ztransp.
	// The diagD is just the diagonal elems of D (i.e., eigvals of S).
	// Note that we assume that S is symmetric and
	// only read the upper-triangular part of S.
	// (verified using octave on certain R matrices; maybe write a
	// testNumericalgb() fn later to test it
	// on simple matrices like diagonal or block-diagonal R..)
	public void diagonalizeLAPACKMatrixOperationsWrapper(LAPACKMatrixOperationsWrapper Z, LAPACKVectorOperationsWrapper diagD) {
		// Check that matrix is square
		if (rows != cols) {
			System.out.print("Error: Cannot diagonalize ");
			System.out.println("non-square matrix.");
			return;
		}

		//interface with jlapack
		jlapackDsyevr(Z, diagD);
	}

	/* Supposedly one of the robust and fastest eigen decomposer
	 * (based on the sth called the RRR method).
	 * To interface with this clapack routine below,
	 * I prepare all arguments and
	 * interface the matrices of clapack and ublas properly.
	 *
	 (subroutine) int dsyevr_(char *jobz, char *range, char *uplo, integer *n,
	 doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	 il, integer *iu, doublereal *abstol, integer *m, doublereal *w,
	 doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
	 integer *lwork, integer *iwork, integer *liwork, integer *info); */

	public void jlapackDsyevr(LAPACKMatrixOperationsWrapper Z, LAPACKVectorOperationsWrapper outW) {
		if (rows != cols) {
			System.out.print("Error: Trying to get eigenvalues of");
			System.out.println(" non-square matrix.");
			return;
		}

		String jobvl = "N"; // left eigenvectors not computed
		String jobvr = "V"; // right eigenvectors are computed
		int N = rows;
		int LDA = rows;
		//String range = "A"; // compute all eigvals
		//String sort = "S"; // compute all eigvals
		//String select = "N"; // compute all eigvals
		//String uplo = "U"; // upper triangle has input and is destroyed on exit

		//double vl, vu;  int il, iu; //not referenced when range='A'
		//vl = vu = 0.0; il = iu = 0;
		//double abstol = 0.0; //use default tol.

		/* mostly outputs */
		//org.netlib.util.intW M = new intW(0); //num. of eigvals found
		int LDZ = N;
		//low-prior.: not sure what this support for Z is for?
		//int[] isuppz = new int[2*N]; // integer vector

		if (!Assert(outW.length == N, "Error: Eigenvectors not right length: "
		            + outW.length))
			return;

		if (!Assert(Z.rows == N, "Error: Eigenvalues not right rows: "
		            + Z.rows + ", " + N))
			return;

		if (!Assert(Z.cols == N, "Error: Eigenvalues not right cols: " + Z.cols))
			return;

		/* do a workspace query */
		double[] work = new double[1];
		int lwork = -1; //workspace query
		// This should not be necessary
		//LAPACKMatrixOperationsWrapper WR = new LAPACKMatrixOperationsWrapper(N,N);
		LAPACKVectorOperationsWrapper WR = new LAPACKVectorOperationsWrapper(N);
		LAPACKMatrixOperationsWrapper WI = new LAPACKMatrixOperationsWrapper(N, N);
		LAPACKMatrixOperationsWrapper VS = new LAPACKMatrixOperationsWrapper(N, N);
		org.netlib.util.intW info = new intW(0);

		Dgeev.dgeev(jobvl, jobvr, N, Z.matrix, 0, N,
		            WR.vector, 0,
		            WI.matrix, 0,
		            null, 0, N,
		            VS.matrix, 0, N, work, 0,
		            lwork, info);

		//Dsyevd.dsyevd(jobz, uplo, N, Z.matrix, 0,
		//	      LDA, outW.vector, 0,
		//	      work, 0, lwork,
		//	      iwork, 0, liwork, info);
		//	Dsyevr.dsyevr(jobz, range, uplo, N, matrix, matrix.length,
		//      LDA, vl, vu, il, iu, abstol,
		//      M, W, W.length, Z.matrix, Z.matrix.length, LDZ,
		//      work, work.length-1, lwork,
		//      iwork, iwork.length-1, liwork, info);
		if (!Assert(info.val == 0,
		            "dsyevr_ workspace query returned with error " + info))
			return;

		lwork = (int)work[0];
		work = new double[lwork];

		/* do the actual call */
		//System.out.println("Performing the actual call...");
		Dgeev.dgeev(jobvl, jobvr, N, Z.matrix, 0, N,
		            WR.vector, 0,
		            WI.matrix, 0,
		            null, 0, N,
		            VS.matrix, 0, N,
		            work, 0,
		            lwork, info);
		//Dsyevd.dsyevd(jobz, uplo, N, Z.matrix, 0,
		//	      LDA, outW.vector, 0,
		//	      work, 0, lwork,
		//	      iwork, 0, liwork, info);
		//Dsyevr.dsyevr(jobz, range, uplo, N, A, LDA, vl, vu, il, iu,
		//	abstol,
		//	M, W, Z, LDZ, isuppz, work, lwork, iwork, liwork, info);
		//dsyevd_(jobz, uplo, N, A, LDA, W, work, lwork, iwork, liwork, info);
		//dsyev_(jobz, uplo, N, A, LDA, W, work, lwork, info);
		//System.out.println("Finished deev");
		//WR.print();
		//VS.print();
		outW.vector = WR.vector;
		Z.matrix = VS.matrix;

		if (!Assert(info.val == 0,
		            "dsyevd_ actual call returned with error " + info))
			return;

		if (!Assert(LDZ == N, "LDZ == N") ||
		    !Assert(LDA == N, "LDA == N") || !Assert(N == rows, "N == rows"))
			return;
	}

	public LAPACKMatrixOperationsWrapper matrixExponential() {
		LAPACKMatrixOperationsWrapper mex;

		// Get eigenvalues (d) and eigenvector matrix (Z)
		if (ZStore == null) {
			System.out.println("Recomputing.." + matrix.length);
			ZStore = new LAPACKMatrixOperationsWrapper(matrix);
			dStore = new LAPACKVectorOperationsWrapper(rows);
			jlapackDsyevr(ZStore, dStore);
			ZinvStore = ZStore.matrixInverse();
		}
		else {
			System.out.println("NOT Recomputing.." + matrix.length);

		}

		// Compute Z exp(D) Zinv using jlapack functions
		LAPACKVectorOperationsWrapper d = new LAPACKVectorOperationsWrapper(dStore.vector);
		d.exponential();
		LAPACKMatrixOperationsWrapper D = new LAPACKMatrixOperationsWrapper(d);
		mex = matrixMultiplication(ZStore, D);
		mex = matrixMultiplication(mex, ZinvStore);

		return mex;
	}


	public LAPACKMatrixOperationsWrapper matrixExponential(double t) {
		LAPACKMatrixOperationsWrapper mex;

		// Get eigenvalues (d) and eigenvector matrix (Z)
		// If these are already computed, don't bother
		if (ZStore == null) {
			//System.out.println("Recomputing ZStore "+matrix.length);
			ZStore = new LAPACKMatrixOperationsWrapper(matrix);
			DStore = new LAPACKMatrixOperationsWrapper(ZStore.rows, ZStore.cols);
			dStore = new LAPACKVectorOperationsWrapper(rows);
			jlapackDsyevr(ZStore, dStore);
			ZinvStore = ZStore.matrixInverse();
		}

		//System.out.println("Printing Z");
		//ZStore.print();
		//System.out.println("Printing Zinv");
		//ZinvStore.print();

		// Compute Z exp(D) Zinv using jlapack functions
		DStore.expDiagonal(dStore.vector, t);
		mex = matrixMultiplication(ZStore, DStore);
		mex = matrixMultiplication(mex, ZinvStore);

		return mex;
	}

	public void expDiagonal(double[] v, double scale) {
		for (int i = 0; i < v.length; i++) {
			set(i, i, Math.exp(v[i]*scale));
		}
	}

	// Compute Zinv using jlapack
	public LAPACKMatrixOperationsWrapper matrixInverse() {
		// N is the order of the matrix
		int N = rows;
		// LDA = leading dimension of matrix
		int LDA = rows;
		// IPIV = integer array
		int[] IPIV = new int[N];
		// LWORK = integer
		int LWORK = N * 2;
		// WORK = double array
		double[] WORK = new double[LWORK];
		// INFO = integer
		org.netlib.util.intW INFO = new intW(0);
		LAPACKMatrixOperationsWrapper Zinv = new LAPACKMatrixOperationsWrapper(matrix);
		Dgetrf.dgetrf(rows, cols, Zinv.matrix,
		              0, LDA, IPIV, 0, INFO);

		if (!Assert(INFO.val == 0,
		            "dgetrf actual call returned with error " + INFO.val))
			return null;

		Dgetri.dgetri(N, Zinv.matrix, 0, LDA,
		              IPIV, 0, WORK, 0, LWORK, INFO);

		if (!Assert(INFO.val == 0,
		            "dgetri actual call returned with error " + INFO.val))
			return null;

		return Zinv;
	}

	public void scale(double s) {
		for (int i = 0; i < matrix.length; i++) {
			matrix[i] = matrix[i] * s;
		}
	}

	/* some helpers to deal with interface between clapack and lapack style
	 * matrices/vectors. may be move these to defutils.cpp later if it's worth
	 * it. */
	//matrix also created like an array for use with clapack.
	/*public T* arrayCreate(int N)
	{
	assert(N > 0);
	T *arr = new T[N];
	if (arr == null) {
	  cerr << "Error: memory allocation failed.\n";
	  assert(0);
	}
	return arr;
	} */

	/*    public void arrayDestroy(T *arr)
	{
	delete [] arr;
	//or delete arr;
	}*/

	/*    public void arrayCopy(const T *srcArr, int N, LAPACKVectorOperationsWrapper tgt)
	{
	int i;
	assert(N >= 0 && N <= (int)tgt.size());
	for (i = 0; i < N; i++)
	  {
	tgt[i] = srcArr[i];
	  }
	} */

	//should int in arguments be integer??
	/*public void matrixCopy(LAPACKMatrixOperationsWrapper src, int LDA, int N)
	{
	int i, j;
	if(!Assert(LDA >= (int)src.size1()) || !Assert(N == (int)src.size2()))
	  return;
	for (i = 0; i < src.size1(); i++) {
	  for (j = 0; j < src.size2(); j++) {
	//tgtArr in col-major order, so (i,j)th elem in posn j*LDA+i.
	tgtArr[j*LDA + i] = src(i,j);
	  }
	}
	} */

	/*    public void matrixCopy()
	{
	//actually i should also input M, the actual # of rows in srcArr holds the
	//data. but its ok, i know that LDA=M for my particular usage above.
	assert(LDA == (int)tgt.size1() && N == (int)tgt.size2());
	//assert(M == tgt.size1() && N == tgt.size2());
	for (int i = 0; i < LDA; i++) {
	  for (int j = 0; j < N; j++) {
	tgt(i,j) = srcArr[j*LDA + i];
	  }
	}
	} */

	private boolean Assert(boolean b, String s) {
		if (!b) {
			System.out.println(s);
		}

		return(b);
	}

	public LAPACKMatrixOperationsWrapper matrixMultiplication(LAPACKMatrixOperationsWrapper A, LAPACKMatrixOperationsWrapper B) {
		String transa = "n";
		String transb = "n";
		int M = A.rows;
		int N = B.cols;
		int K = A.cols;
		double alpha = 1;
		int LDA = M;
		int LDB = K;
		double beta = 0.0;
		LAPACKMatrixOperationsWrapper C = new LAPACKMatrixOperationsWrapper(M, N);
		int LDC = M;
		Dgemm.dgemm(transa, transb, M, N, K,
		            alpha, A.matrix, 0, LDA,
		            B.matrix, 0, LDB, beta,
		            C.matrix, 0, LDC);
		return C;
	}

	public void print() {
		/*System.out.println("Printing matrix");

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				System.out.print(matrix[index(i, j)] + " ");
			}

			System.out.println();
		}*/

		//double[][] twoDMatrix = new double[rows][cols];
		//MatConv.copyOneDintoTwoD(twoDMatrix, matrix);
		//System.out.println("Printing two d:");
		//for(int i = 0; i < rows; i++) {
		//    for(int j = 0; j < cols; j++) {
		//	System.out.print(twoDMatrix[i][j]+" ");
		//   }
		//   System.out.println();
		//}
	}


	public void sumRows() {
		double sum = 0.0;
		System.out.println("Printing row sums");

		for (int i = 0; i < rows; i++) {
			sum = 0.0;

			for (int j = 0; j < cols; j++) {
				sum += matrix[index(i, j)];
			}

			System.out.print(sum + " ");
		}

		System.out.println();
	}

	public void sumCols() {
		System.out.println("Printing column sums");
		double sum = 0.0;

		for (int i = 0; i < cols; i++) {
			sum = 0.0;

			for (int j = 0; j < rows; j++) {
				sum += matrix[index(j, i)];
			}

			System.out.print(sum + " ");
		}

		System.out.println();
	}

	// Normalizes by row
	public void normalize() {
		double sum = 0.0;

		for (int i = 0; i < rows; i++) {
			sum = 0.0;

			for (int j = 0; j < cols; j++) {
				sum += matrix[index(i, j)];
			}

			for (int j = 0; j < cols; j++) {
				matrix[index(i, j)] = matrix[index(i, j)] / sum;
			}
		}
	}

	// Normalizes by row
	public void zero() {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				matrix[index(i, j)] = 0.0;
			}
		}
	}

	// Normalizes by row
	public void normalizeNotFirstCol() {
		double sum = 0.0;

		for (int i = 0; i < rows; i++) {
			sum = 0.0;

			for (int j = 1; j < cols; j++) {
				sum += matrix[index(i, j)];
			}

			for (int j = 1; j < cols; j++) {
				matrix[index(i, j)] = matrix[index(i, j)] / sum;
			}
		}
	}

	public void matrixScale(double a) {
		for (int i = 0; i < matrix.length; i++) {
			matrix[i] = matrix[i] * a;
		}
	}

	public void elementwiseMultiply(double[][] m) {
		if (!Assert(m.length == rows, "Error in elementwiseMultiply: "
		            + "rows do not match"))
			return;

		for (int i = 0; i < rows; i++) {
			if (!Assert(m[0].length == cols, "Error in elementwiseMultiply: "
			            + "columns do not match"))
				return;

			for (int j = 0; j < cols; j++) {
				matrix[index(i, j)] = matrix[index(i, j)] * m[i][j];
			}
		}
	}

	public void elementwiseDivide(double[][] m) {
		if (!Assert(m.length == rows, "Error in elementwiseDivide: "
		            + "rows do not match"))
			return;

		for (int i = 0; i < rows; i++) {
			if (!Assert(m[0].length == cols, "Error in elementwiseDivide: "
			            + "columns do not match"))
				return;

			for (int j = 0; j < cols; j++) {
				matrix[index(i, j)] = matrix[index(i, j)] / m[i][j];
			}
		}
	}

	private int index(int r, int c) {
		return ((c*rows) + r);
	}

	// Make sure rows sum to 0
	public void rowsSumToZero() {
		double sum = 0.0;

		for (int i = 0; i < rows; i++) {
			sum = 0.0;

			for (int j = 0; j < cols; j++) {
				if (i != j)
					sum += get(i, j);
			}

			set(i, i, -sum);
		}
	}
}
