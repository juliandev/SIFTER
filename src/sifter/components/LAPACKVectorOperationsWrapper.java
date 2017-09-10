/* LAPACKVectorOperationsWrapper.java: Interface to numeric linear algebra computations like
 *                  matrix opns in CLAPACK.
 * Author: Barbara Engelhardt (9/06)
 */

package sifter.components;

public class LAPACKVectorOperationsWrapper {

	public double[] vector;
	public int length;


	public LAPACKVectorOperationsWrapper() {
		length = 0;
		vector = null;
	}

	public LAPACKVectorOperationsWrapper(int len) {
		vector = new double[len];
		length = len;
	}

	// puts exp(v[i]*s) into the new vector
	public LAPACKVectorOperationsWrapper(double[] v, double s) {
		vector = new double[v.length];
		length = v.length;

		for (int i = 0; i < length; i++) {
			vector[i] = Math.exp(v[i] * s);
		}
	}

	public LAPACKVectorOperationsWrapper(double[] v) {
		vector = new double[v.length];
		length = v.length;

		for (int i = 0; i < length; i++) {
			vector[i] = v[i];
		}
	}

	public int length() {
		return length;
	}

	public void exponential() {
		for (int i = 0; i < length; i++) {
			vector[i] = Math.exp(vector[i]);
		}
	}

	// this version includes a scaling parameter s
	public void exponential(double s) {
		for (int i = 0; i < length; i++) {
			vector[i] = Math.exp(vector[i] * s);
		}
	}

	public void getDiag(double[] mat) {
		int len = (int)(Math.log(mat.length) / Math.log(2));

		if (len == vector.length) {
			for (int i = 0; i < len; i++) {
				vector[i] = mat[i*len+i];
			}
		}
	}

	public void print() {
		System.out.println("Printing vector");

		for (int i = 0; i < length; i++) {
			System.out.print(vector[i] + " ");
		}

		System.out.println();
	}
}
