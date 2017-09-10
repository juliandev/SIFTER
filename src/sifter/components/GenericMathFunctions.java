/**
 * Unions, utility stuff, choose, file IO.
 *
 * @author Barbara Engelhardt
 */

package sifter.components;

import java.util.*;
import java.io.*;


public class GenericMathFunctions {

	public static double choose(int n, int k) {
		double total = 1;
		double divisor = 1;
		k = Math.max(k, n - k);

		for (int i = n; i > k; i--) {
			total *= (double)i;
		}

		for (int i = 1; i <= (n - k); i++) {
			divisor *= (double)i;
		}

		total /= divisor;
		return total;
	}

	// Computes the size of the power set of n >= 2
	public static int powerSetSize(int n) {
		int total = 0;

		for (int i = 1; i < (n + 1) / 2; i++) {
			total += (2 * choose(n, i));
		}

		if (n % 2 == 0) {
			total += (choose(n, n / 2));
		}

		return total + 2;
	}

	// REFACTORME: use templated return types!
	public static <T extends Object> Vector<T> intersection(Vector<T> a, Vector<T> b) {
		Vector<T> inter = new Vector<T>();

		for (int i = 0; i < a.size(); i++) {
			if (b.contains(a.elementAt(i))) {
				inter.add(a.elementAt(i));
			}
		}

		return inter;
	}

	public static <T extends Object> Vector<T> union(Vector<T> a, Vector<T> b) {
		Vector<T> union = b;

		for (int i = 0; i < a.size(); i++) {
			if (!b.contains(a.elementAt(i))) {
				union.add(a.elementAt(i));
			}
		}

		return union;
	}

	public static double logSafe(double x) {
		if (x <= 0)
			return -100000000;
		else
			return Math.log(x);
	}

	public static double getLogSum(double log_a, double log_b) {
		double v;

		if (log_a < log_b) {
			v = log_b + Math.log(1 + Math.exp(log_a - log_b));
		}
		else {
			v = log_a + Math.log(1 + Math.exp(log_b - log_a));
		}

		if (v >= 0.0)
			return (Math.log(Math.exp(log_a) + Math.exp(log_b)));

		return(v);
	}

	// Are these two doubles equal up to 5 significant digits?
	public static boolean areEqual(double a, double b) {
		if ((int)Math.round(a) != (int)Math.round(b))
			return false;

		a = a - (int)Math.floor(a);
		b = b - (int)Math.floor(b);

		if ((int)Math.round(a*10000) != (int)Math.round(b*10000))
			return false;

		return true;
	}

	@SuppressWarnings({ "unused", "resource" })
	public static Hashtable<String, Double> readInHashtable(String filename) {
		int lineno = 0;
		String str = null;
		Hashtable<String, Double> params = new Hashtable<String, Double>();

		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));

			while ((str = in.readLine()) != null) {
				String type = null;
				double val = -1.0;
				StringTokenizer st = new StringTokenizer(str);

				// type of scale param
				if (st.hasMoreTokens()) {
					type = st.nextToken();
				}

				// value of scale param
				if (st.hasMoreTokens()) {
					val = Double.parseDouble(st.nextToken());
				}

				if (type != null && val > -0.05) {
					params.put(type, new Double(val));
				}

				lineno++;
			}
		}
		catch (IOException ioe) {
			System.err.println("Couldn't get data file "
			                   + filename + "from the directory");
			System.exit(1);
		}

		return params;
	}


	// Assumes a vector of doubles
	@SuppressWarnings({ "unused", "resource" })
	public static Vector<Double> readInVector(String filename) {
		int lineno = 0;
		int i = 0;
		String str = null;
		Vector<Double> params = new Vector<Double>();

		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));

			while ((str = in.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(str);

				while (st.hasMoreTokens()) {
					params.add(i++, new Double(st.nextToken()));
				}

				lineno++;
			}
		}
		catch (IOException ioe) {
			System.err.println("Couldn't get data file "
			                   + filename + "from the directory");
			System.exit(1);
		}

		return params;
	}

	public static <T extends Object> void randomizeVector(Vector<T> v) {
		Random rand = new Random();

		for (int i = 0; i < v.size(); i++) {
			swapVectorElements(v, i, rand.nextInt(v.size()));
		}
	}

	public static <T extends Object> void swapVectorElements(Vector<T> v, int i, int j) {
		if (i == j)
			return;

		if (i > j) {
			int k = j;
			j = i;
			i = k;
		}

		v.insertElementAt(v.elementAt(i), j);
		v.setElementAt(v.elementAt(j + 1), i);
		v.removeElementAt(j + 1);
	}

}
