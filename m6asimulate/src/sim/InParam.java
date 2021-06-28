package sim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

public class InParam {
	
	private static InParam args = new InParam();
	private HashMap<String, Object> params = null;
	
	private InParam() {
		params = new HashMap<>();
		params.put("-bn", "10");
		params.put("-en", "16");
		params.put("-sl", "300");
		params.put("-rl", "150");
		params.put("-mrl", "5");
		params.put("-mp", "0.5");
		params.put("-rs", "100");
		params.put("-q", "J");
	}
	
	public static InParam getParams() {
		return args;
	}
	
	public void putParams(String[] args) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].charAt(0) == '-') {
				params.put(args[i], i + 1 < args.length ? args[i + 1] : ""); 
			}
		}
		try {
			putBackgoundNum((String) params.get("-bn"));
			putEnrichNum((String) params.get("-en"));
			putSegmentLen((String) params.get("-sl"));
			putReadLen((String) params.get("-rl"));
			putMinReadLen((String) params.get("-mrl"));
			putMutationProportion((String) params.get("-mp"));
			putRandSegSize((String) params.get("-rs"));
		}
		catch (Exception e) {
			System.err.println(e.getClass().getSimpleName() + ": " + e.getMessage());
		}
	}
	
	public boolean check() {
		boolean out = true;
		out = !params.containsKey("-h");
		try {
			out = out && getGenomeFileName() != null;
			out = out && getReferenceFileName() != null;
			out = out && getPeakFileName() != null;
			out = out && getMutationFileName() != null;
			out = out && getIPFileName() != null;
			out = out && getInputFileName() != null;
			if (!out) {
				throw new FileNotFoundException("Lack parameter(s)");
			}
			FileReader fr = new FileReader(new File(getGenomeFileName()));
			fr.close();
			fr = new FileReader(new File(getReferenceFileName()));
			fr.close();
			
			out = out && getBackgoundNum() > 0;
			out = out && getEnrichNum() >= 0.0;
			out = out && getReadLen() > 0;
			out = out && getMutationProportion() >= 0.0 && getMutationProportion() <= 1.0;
			if (!out) {
				throw new NumberFormatException("Exsits negative number(s) or mutation proportion greater than 1.0");
			}
			if (getSegmentLen() <= getReadLen()) {
				params.put("-sl", params.get("-rl"));
				System.out.println("Warning: Segment length is no greater than read length, simulating single-end fastq using read length");
			}
			if (getSegmentLen() - getReadLen() <= getRandSegSize()) {
				params.put("-rs", getSegmentLen() - getReadLen());
				System.out.println("Warning: Segment length is no greater than the sum of random segment length and read length, random segment length is cut");
			}
			return out;
		} catch (Exception e) {
			System.err.println(e.getClass().getSimpleName() + ": " + e.getMessage());
			help();
			return false;
		}
	}
	
	private void help() {
		InputStream is = Object.class.getResourceAsStream("/Help_Document");
		try {
			InputStreamReader isr = new InputStreamReader(is, "UTF-8");
			BufferedReader reader = new BufferedReader(isr);
			String line = null;
			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public String getGenomeFileName() {
		return (String) params.get("-g");
	}
	
	public String getReferenceFileName() {
		return (String) params.get("-r");
	}
	
	public String getPeakFileName() {
		return (String) params.get("-p");
	}
	
	public String getMutationFileName() {
		return (String) params.get("-m");
	}
	
	public String getIPFileName() {
		return (String) params.get("-ip");
	}
	
	public String getInputFileName() {
		return (String) params.get("-in");
	}
	
	public String getCircFileName() {
		return (String) params.get("-circ");
	}
	
	public String getQuality() {
		return (String) params.get("-q");
	}
	
	private void putBackgoundNum(String s) {
		params.put("-bn", Integer.parseInt(s));
	}
	
	public int getBackgoundNum() {
		return (int) params.get("-bn");
	}
	
	private void putEnrichNum(String s) {
		params.put("-en", Double.parseDouble(s));
	}
	
	public double getEnrichNum() {
		return (double) params.get("-en");
	}
	
	private void putSegmentLen(String s) {
		params.put("-sl", Integer.parseInt(s));
	}
	
	public int getSegmentLen() {
		return (int) params.get("-sl");
	}
	
	private void putMinReadLen(String s) {
		params.put("-mrl", Integer.parseInt(s));
	}
	
	public int getMinReadLen() {
		return (int) params.get("-mrl");
	}
	
	private void putRandSegSize(String s) {
		params.put("-rs", Integer.parseInt(s));
	}
	
	public int getRandSegSize() {
		return (int) params.get("-rs");
	}
	
	private void putReadLen(String s) {
		params.put("-rl", Integer.parseInt(s));
	}
	
	public int getReadLen() {
		return (int) params.get("-rl");
	}
	
	private void putMutationProportion(String s) {
		params.put("-mp", Double.parseDouble(s));
	}
	
	public double getMutationProportion() {
		return (double) params.get("-mp");
	}
	
	public boolean isPaired() {
		return getReadLen() != getSegmentLen();
	}
	
	public boolean isControled() {
		return params.get("-uc") == null;
	}
	
	public boolean isAnnoted() {
		return params.get("-ua") == null;
	}
}
