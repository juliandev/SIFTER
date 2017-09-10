/*
 * Created on Jun 18, 2005
 *
 *
 * @author Barbara Engelhardt (bee@cs.berkeley.edu)
 */
package sifter.components;


public class ProteinFunctionMarkovState {

	private int[] meter;
	private int maxSum;
	private int currentSum;
	private int currentIndex;
	private int[] pps;
	private int[] ppsIndices;
	private int currentPositiveSum;
	private int numChildren;
	private int numFunctions;
	private boolean ppsStart;
	private boolean start;
	private boolean emptyPowerSet;
	private int[] childTemp;
	private int[] childMarker;

	public ProteinFunctionMarkovState(int length, int maxsum) {
		maxSum = maxsum;
		meter = new int[length];
		//System.out.println("Power set of size "+(number*length)+", number "+number+" length "+length);
		Reset();
		numChildren = 1;
		numFunctions = length;
		childTemp = null;
		childMarker = null;
	}

	public ProteinFunctionMarkovState(int number, int length, int maxsum) {
		meter = new int[number*length];
		maxSum = maxsum;
		//System.out.println("Power set of size "+(number*length)+", number "+number+" length "+length);
		Reset();
		numChildren = number;
		numFunctions = length;

		if (number == 1) {
			childTemp = null;
			childMarker = null;
		}
		else {
			childTemp = new int[length];
			childMarker = new int[number];
		}
	}

	public ProteinFunctionMarkovState(int length) {
		meter = new int[length];
		maxSum = length;
		//System.out.println("Power set of size "+(number*length)+", number "+number+" length "+length);
		Reset();
		numChildren = 1;
		numFunctions = length;
		childTemp = null;
		childMarker = null;
	}

	public void Reset() {
		for (int i = 0; i < meter.length; i++)
			meter[i] = 0;

		currentSum = 0;
		pps = null;
		ppsIndices = null;
		currentPositiveSum = 0;
		ppsStart = true;
		start = true;

		if (childTemp != null) {
			for (int i = 0; i < numChildren; i++)
				childMarker[i] = 0;

			for (int i = 0; i < numFunctions; i++)
				childTemp[i] = 0;
		}

		currentIndex = 0;
	}

	public int length() {
		return meter.length;
	}

	public int elementAt(int i) {
		if (i >= meter.length) {
			System.out.println("WARNING: accessing an element in power set that doesn't exist: " + i);
			return -1;
		}

		return meter[i];
	}

	public boolean hasNext() {
		if (numChildren > 1 && maxSum < meter.length)
			return hasNextChildren(meter);

		return hasNext(meter);
	}

	public boolean hasNext(int[] v) {
		//System.out.println("in has next with ");
		//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
		//System.out.println();
		if (v == null)
			return false;

		for (int i = v.length - 1; i >= v.length - maxSum && i >= 0; i--) {
			if (v[i] == 0)
				return true;
		}

		return false;
	}

	public boolean hasNextChildren(int[] v) {
		//System.out.println("in has next with ");
		//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
		//System.out.println();
		if (v == null)
			return false;

		for (int i = 0; i < childMarker.length; i++) {
			for (int j = 0; j < childTemp.length; j++) {
				childTemp[j] = v[i*childTemp.length+j];

				if (hasNext(childTemp))
					return true;
			}
		}

		return false;
	}

	public void getNext() {
		if (meter == null)
			return;

		if (start) {
			start = false;
			return;
		}

		if (maxSum == meter.length)
			currentSum = getNext(meter, currentSum);
		else
			if (numChildren == 1)
				currentSum = getNextApproximate(meter, currentSum, maxSum);
			else
				currentSum = getNextApproximateChildren(meter, currentSum, maxSum);

		currentIndex++;
	}

	/* do not include the case when all meter is 0.
	 * Doesnt work for approximation yet. */
	public void getNextNonZeros() {
		if (meter == null)
			return;

		if (maxSum == meter.length)
			currentSum = getNext(meter, currentSum);
		else
			if (numChildren == 1)
				currentSum = getNextApproximate(meter,
				                                currentSum, maxSum);
			else
				currentSum = getNextApproximateChildren(meter,
				                                        currentSum,
				                                        maxSum);

		currentIndex++;
	}


	public int getNext(int[] v, int sum) {
		//System.out.println("IN get next with ");
		//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
		//System.out.println();
		if (v == null)
			return 0;

		if (v.length == 1) {
			v[0] = 1;
			return 1;
		}
		else {
			int total = 0;

			for (int j = v.length - 1; j >= 0 ; j--) {
				boolean nextRotate = true;

				for (int innerK = 0; innerK < j; innerK++) {
					if (v[innerK] == 0)
						nextRotate = false;
				}

				if (nextRotate)
					v[j] = (v[j] + 1) % 2;

				if (v[j] == 1)
					total++;
			}

			//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
			//System.out.println();
			return total;
		}
	}

	public int getCurrentSum() {
		return currentSum;
	}

	// This is functional, although it doesn't really
	// work in the way it should for children yet.
	public int getNextApproximate(int[] v, int sum, int limit) {
		//System.out.println("In get next approximate");
		//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
		//System.out.println();
		if (v == null)
			return 0;

		if (v.length == 1) {
			v[0] = 1;
			return 1;
		}
		else {
			int total = 0;

			for (int j = v.length - 1; j >= 0 ; j--) {
				boolean nextRotate = true;
				int innerTotal = 0;

				for (int innerK = j - 1; innerK >= 0; innerK--) {
					if (v[innerK] == 0 &&
					    total + innerTotal < limit)
						nextRotate = false;

					if (v[innerK] == 1)
						innerTotal++;
				}

				if (nextRotate) {
					v[j] = 1;

					for (int k = 0; k < j; k++) {
						v[k] = 0;
					}

					total++;
					j = 0; // break
				}

				if (v[j] == 1)
					total++;
			}

			//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
			//System.out.println();
			return total;
		}
	}

	public int getNextApproximateChildren(int[] v, int sum, int limit) {
		//System.out.println("In get next approximate children");
		//for(int i = 0; i < v.length; i++) System.out.print(v[i]);
		//System.out.println();
		if (v == null)
			return 0;

		if (v.length == 1) {
			v[0] = 1;
			return 1;
		}
		else
			if (numChildren > 1) {
				int total = sum;

				for (int i = 0; i < numChildren; i++) {
					int tempTotal = 0;

					for (int j = 0; j < childTemp.length; j++) {
						childTemp[j] = v[i*childTemp.length+j];
						tempTotal += childTemp[j];
					}

					if (hasNext(childTemp)) {
						// get next on this child
						total -= tempTotal;
						total += getNextApproximate(childTemp, tempTotal, maxSum);

						for (int j = 0; j < childTemp.length; j++) {
							v[i*childTemp.length+j] = childTemp[j];
						}

						i = numChildren; //break
					}
					else {
						// reset this child
						total -= tempTotal;

						for (int j = 0; j < childTemp.length; j++) {
							v[i*childTemp.length+j] = 0;
						}
					}
				}

				return total;
			}

		return 0;
	}
	// 01011 -> 01100 -> 01101 -> 01110 -> 01111

	// This makes a power set out of only the positive variables.
	public void positivePowerSet() {
		int len = 0;

		for (int i = 0; i < meter.length; i++) {
			if (meter[i] == 1)
				len++;
		}

		if (len == 0) {
			// empty power set
			emptyPowerSet = true;
		}

		if (len > 0) {
			pps = new int[len];
			ppsIndices = new int[len];
			emptyPowerSet = false;
		}

		len = 0;

		for (int i = 0; i < meter.length; i++) {
			if (meter[i] == 1) {
				pps[len] = 0;
				ppsIndices[len++] = i;
			}
		}

		ppsStart = true;
		currentPositiveSum = 0;
	}

	public int positiveLength() {
		//if(pps == null) return 0;
		//return pps.length;
		return currentPositiveSum;
	}

	public boolean positiveHasNext() {
		if (emptyPowerSet == true) {
			emptyPowerSet = false;
			return true;
		}

		if (pps == null)
			return false;

		//System.out.println("In positive has next");
		return hasNext(pps);
	}

	public void positiveNext() {
		if (pps == null)
			return;

		if (ppsStart) {
			ppsStart = false;
			//for(int i = 0; i < pps.length; i++)
			//	 System.out.print(pps[i]);
			//System.out.println();
			return;
		}

		currentPositiveSum = getNext(pps, currentPositiveSum);
		//for(int i = 0; i < pps.length; i++)
		//System.out.print(pps[i]);
		//System.out.println();
	}

	public int positivePSAndNegativeLength() {
		// Also the case when pps is all negative
		if (pps == null)
			return meter.length;
		else
			return (meter.length - currentSum + currentPositiveSum);
	}

	public int positivePSAndNegativeIndex(int j) {
		// Also the case when pps is all negative
		if (pps == null)
			return j;
		else {
			for (int i = 0; i < meter.length; i++) {
				if (meter[i] == 0 && j == 0)
					return i;

				if (meter[i] == 0)
					j--;
			}

			for (int i = 0; i < pps.length; i++) {
				if (j == 0 && pps[i] == 1)
					return ppsIndices[i];

				if (pps[i] == 1)
					j--;
			}

			return 0;
		}
	}

	// TODO: make sure sums are in working order.
	// this function returns the molecular function index
	// of the integer j
	public int positivePSAndNegativeFunctionIndex(int j) {
		// Also the case when pps is all negative
		if (pps == null)
			return (j % numFunctions);
		else {
			for (int i = 0; i < meter.length; i++) {
				if (meter[i] == 0 && j == 0)
					return i % numFunctions;

				if (meter[i] == 0)
					j--;
			}

			for (int i = 0; i < pps.length; i++) {
				if (j == 0 && pps[i] == 1)
					return ppsIndices[i] % numFunctions;

				if (pps[i] == 1)
					j--;
			}

			return 0;
		}
	}

	// TODO: make sure sums are in working order.
	// Returns the index of the current child
	public int positivePSAndNegativeChildIndex(int j) {
		// also the case when all negative
		if (pps == null)
			return (j / numFunctions);
		else {
			//int len = meter.length - pps.length + currentPositiveSum;
			for (int i = 0; i < meter.length; i++) {
				if (meter[i] == 0 && j == 0)
					return (i / numFunctions);

				if (meter[i] == 0)
					j--;
			}

			for (int i = 0; i < pps.length; i++) {
				if (j == 0 && pps[i] == 1)
					return (ppsIndices[i] / numFunctions);

				if (pps[i] == 1)
					j--;
			}

			return 0;
		}
	}

	public int functionIndex(int j) {
		return j % (numFunctions);
	}

	public int index(int j) {
		return (j / (numFunctions));
	}

	public int setIndex() {
		return currentIndex;
		/*int ind = 0;
		int pow = 1;
		if(meter == null) return 0;
		for(int i = 0; i < meter.length; i++) {
		ind += pow*meter[i];
		pow *= 2;
		}
		return ind;*/
	}

	// Assumes inclusion of empty set
	public int powerSetSize() {
		int sum = 1;

		for (int i = 1; i <= maxSum; i++) {
			sum += (int)GenericMathFunctions.choose(meter.length, i);
		}

		return sum;
	}

	// count the number of elements in power set
	public int powerSetSum() {
		int sum = 1;

		for (int i = 1; i <= maxSum; i++) {
			sum += meter[i];
		}

		return sum;
	}

	public void printPowerSet() {
		if (meter == null)
			return;

		for (int i = 0; i < meter.length; i++)
			System.out.print(meter[i]);

		System.out.println();
	}

	public void printPositive() {
		if (meter == null || pps == null)
			return;

		for (int i = 0; i < meter.length; i++)
			System.out.print(meter[i]);

		System.out.print(" ");

		for (int i = 0; i < pps.length; i++)
			System.out.print(pps[i]);

		System.out.println();
	}

	// Returns -1 when not off by one
	// Returns length when equal
	// otherwise returns index of off-by-one index
	public int offByOne(ProteinFunctionMarkovState ps2) {
		if (ps2.length() != meter.length)
			return -1;

		int index = -1;

		for (int i = 0; i < meter.length; i++) {
			if (meter[i] != ps2.elementAt(i)) {
				if (index >= 0)
					return -1;
				else
					index = i;
			}
		}

		if (index == -1)
			return meter.length;

		return index;
	}

}
