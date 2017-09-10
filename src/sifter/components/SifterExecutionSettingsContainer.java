// $Id: SifterExecutionSettingsContainer.java,v 1.2 2005/05/27 06:56:02 sprite Exp $

package sifter.components;
import java.util.HashMap;

/** Class that stores options, constants, and other preferences
 * that can be modified by the user.
 *
 * I created this to store nearly any variable setting, for the classes that
 * Barbara uses. Using the commandlinearguments class is not acceptable because
 * we may want to extend this to a GUI in the future.
 *
 * TODO: Handle any type of value? So far it only handles objects.
 *
 * @author Steven Chan <steven@berkeley.edu>
 * @author Barbara Engelhardt (bee@cs.berkeley.edu)
 **/
public class SifterExecutionSettingsContainer {
	public SifterExecutionSettingsContainer() {
		this.h = new HashMap<String, Object>();

		// THESE ARE NOT DEFAULTS. DEFAULTS EXPLICITLY SET
		// IN Sifter.java
		setSetting("output_filename", "./output.rdata");
		setSetting("bestguess-min", new Float(0.2));
		setSetting("bestguess-max", new Float(10.0));
		setSetting("bestguess-incr", new Float(0.2));
		setSetting("verbose", new Boolean(false));
		setSetting("em", new Boolean(false));
		setSetting("folds", new Integer(0));
		setSetting("truncation", new Integer(0));
		setSetting("cutoff", new Double(0.01));
	}

	/** Set option "key" to "value".
	 * <br /><b>Example: </b> setOption("output_file", "./output.rdata");
	 * @param key Key value string.
	 * @param value Value of the key
	 **/
	public void setSetting(String key, Object value) {
		h.put(key, value);
	}

	/** Get the value of option "key".
	 * @param key Key value object. (Most of the time, it's a string.)
	 *
	 **/
	public Object getSetting(Object key) {
		if (h.containsKey(key))
			return h.get(key);

		return null;
	}

	/** Stores all the key-value pairs.
	 **/
	private HashMap<String, Object> h;

	/** Set option "key" to "value".
	 * <br /><b>Example: </b> setOption("output_file", true);
	 * @param key Key value object. (Most of the time, you'd want a string.)
	 * @param value Value of the key
	 */
	/*    public void setSetting(String string, boolean b) {
	        h.put(string, "b");
	    }*/

	public String toString() {
		String retval = "";

		retval += h.entrySet().toString();

		return retval;
	}
}

