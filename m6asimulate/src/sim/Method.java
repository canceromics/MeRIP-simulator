package sim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.AccessibleObject;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.*;

import genome.*;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import note.*;

public class Method {
	
	private final static NormalDistribution mut_per_dis = new NormalDistribution(0.0, 0.05);
	private final static NormalDistribution OR_dis = new NormalDistribution(0.0, 0.9);
	
	private final static char[] base_units = {'A', 'G', 'C', 'T'};
	private static String FixMapQ = null;
	private static int id_count = 0;
	private static HashMap<String, ArrayList<Mutation>> muts = null;
	private static InParam args = null;
	private final static int[] real_read = {1,36,222,548,1187,1982,3055,4231,5815,7547,9517,11760,14159,17008,20023,23516,27328,31099,35594,
			40294,46261,52970,62635,74133,86894,101858,119138,136700,155572,176409,401071,404055,428248,452741,478697,505872,533741,561985,
			592168,623776,655799,687947,722339,757025,792746,828189,865288,904037,944077,985317,1027454,1084572,1239489,1418217,1660947,
			1938452,2245296,2568393,2910301,3263897,3623931,3980022,4390208,4825790,5276488,5755329,6242831,6747273,7282460,7819999,8379786,
			8954017,9549646,10184293,10818876,11473150,12152280,12861224,13580622,14330465,15061394,15811365,16617938,17445991,18290369,
			19166564,20072868,20988116,21950236,22918837,23910047,24902197,25936977,26985149,28084077,29298832,30928220,33417147,37314527,
			77863330,95825613,109005185,154826218}; 
	private final static UniformIntegerDistribution real_read_dis = new UniformIntegerDistribution(1,real_read[real_read.length - 1]);
	
	public static void run(InParam in_args) {
		args = in_args;
		HashMap<String, ArrayList<Gene>> gene_info = loadGeneInfo(args.getReferenceFileName());
		HashMap<String, ArrayList<Peak>> peak_info = Files.exists(Paths.get(args.getPeakFileName())) ? loadPeakInfo(args.getPeakFileName()) : new HashMap<>();
		for (Entry<String, ArrayList<Peak>> entry : peak_info.entrySet()) {
			Collections.sort(entry.getValue());
		}
		loadGenomeInfo(args, gene_info, peak_info);
		muts = loadMutationInfo(args.getMutationFileName());
//		HashMap<String, IntervalTree<CircleClip>> circ_info = loadCircs(args.getCircFileName());
		HashMap<String, IntervalTree<CircleClip>> circ_info = loadCircs(null);
		for (Entry<String, ArrayList<Mutation>> entry : muts.entrySet()) {
			Collections.sort(entry.getValue());
		}
		printNow("Load complete!");
		FixMapQ = buildString(args.getReadLen(), "J");
		ensureMutDes(peak_info);
		makeMutTimes();
		
		writeFile(args.getMutationFileName() + "_used", muts);
		writeFile(args.getInputFileName().substring(0, args.getInputFileName().lastIndexOf("/") + 1) + "sim_annotation.txt", new ArrayList<>());
//		appendSimFile(args, gene_info, peak_info, muts, circ_info);
		printNow("Add unvisited peaks!");
		appendPeak(peak_info, muts, circ_info, gene_info);
//		writeFile(args.getMutationFileName() + "_used", muts);
		writeTreeFile(args.getCircFileName() + "_used", circ_info);
	}

	public static void test() {
		ArrayList<Double> region = new ArrayList<>();
		int[] count = new int[10];
		for (double i = -3.95; i < 4.03; i += 0.1) {
			System.out.println(sampleMutPerOR(0.5));
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(new File("D:\\work\\workspace\\data\\OR_value")))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				++count[binarySearch(Math.log(Double.parseDouble(line)), region, 0, region.size())];
			}
			System.out.println(count);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static HashMap<String, ArrayList<Gene>> loadGeneInfo(String file_name){
		return readFile(file_name, 'e');
	}
	
	public static HashMap<String, ArrayList<Site>> loadSiteInfo(String file_name){
		return readFile(file_name, 's');
	}
	
	public static HashMap<String, ArrayList<Peak>> loadPeakInfo(String file_name){
		return readFile(file_name, 'p');
	}
	
	public static HashMap<String, ArrayList<Mutation>> loadMutationInfo(String file_name){
		return readFile(file_name, 'm');
	}
	
	public static void loadGenomeInfo(InParam args, HashMap<String, ArrayList<Gene>> gene_info, HashMap<String, ArrayList<Peak>> peak_info){
		HashMap<String, ArrayList<Mutation>> mut = new HashMap<>();
		HashMap<String, ArrayList<CircleClip>> circ = new HashMap<>();
		boolean mut_flag = Files.exists(Paths.get(args.getMutationFileName()));
		boolean peak_flag = Files.exists(Paths.get(args.getPeakFileName()));
		boolean circ_flag = args.getCircFileName() == null || Files.exists(Paths.get(args.getCircFileName()));
		File fi = new File(args.getGenomeFileName());
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(fi));
			String line = null;
			String chr = null;
			StringBuffer bases = new StringBuffer();
			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>'){
					IntervalTree<Gene> gene_region = new IntervalTree<>();
					IntervalTree<Exon> exon_region = new IntervalTree<>();
					if (gene_info.containsKey(chr)) {
						for (int i = 0; i < gene_info.get(chr).size(); i++) {
							Gene gene = gene_info.get(chr).get(i);
							gene.checkPos();
							gene.setBases(bases);
							gene_region.put(gene.getStart() + 1, gene.getEnd(), gene);
							putInExonTree(exon_region, gene);
						}
					}
					if (chr != null && (!mut_flag || !peak_flag)) {
						ArrayList<Peak> peaks = peak_flag ? peak_info.get(chr) : new ArrayList<>(); 
						ArrayList<CircleClip> circs = new ArrayList<>();
						ArrayList<Mutation> muts = simMutPeak(args, peaks, circs, gene_region, exon_region, bases);
						if (!mut_flag) {
							mut.put(chr, muts);
						}
						if (!circ_flag) {
							circ.put(chr, circs);
						}
						peak_info.put(chr, peaks);
					}
					String[] temp_chr = line.split("\\s+");
					chr = temp_chr[0].substring(1);
					bases.setLength(0);
				}
				else {
					bases.append(line);
				}
			}
			reader.close();
			IntervalTree<Gene> gene_region = new IntervalTree<>();
			IntervalTree<Exon> exon_region = new IntervalTree<>();
			if (gene_info.containsKey(chr)) {
				for (int i = 0; i < gene_info.get(chr).size(); i++) {
					Gene gene = gene_info.get(chr).get(i);
					gene.checkPos();
					gene.setBases(bases);
					gene_region.put(gene.getStart(), gene.getEnd() - 1, gene);
					putInExonTree(exon_region, gene);
				}
			}
			if (chr != null && (!mut_flag || !peak_flag)) {
				ArrayList<Peak> peaks = Files.exists(Paths.get(args.getPeakFileName())) ? peak_info.get(chr) : new ArrayList<>(); 
				ArrayList<CircleClip> circs = new ArrayList<>();
				ArrayList<Mutation> muts = simMutPeak(args, peaks, circs, gene_region, exon_region, bases); 
				if (!mut_flag) {
					mut.put(chr, muts);
					writeFile(args.getMutationFileName(), mut);
				}
				if (!circ_flag) {
					circ.put(chr, circs);
					writeFile(args.getCircFileName(), circ);
				}
				if (!peak_flag) {
					peak_info.put(chr, peaks);
					writeFile(args.getPeakFileName(), peak_info);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}finally {
			if (reader != null){
				try{
					reader.close();
				}
				catch(IOException e1){
				}
			}
		}
	}
	
	private static ArrayList<Mutation> simMutPeak(InParam args, ArrayList<Peak> peaks, ArrayList<CircleClip> circ, IntervalTree<Gene> gene_region, IntervalTree<Exon> exon_region, StringBuffer bases) {
		String[] des = {"Peak", "Gene", "None", "Exon"};
		int read_len = args.getReadLen();
		UniformIntegerDistribution uid = new UniformIntegerDistribution(1, bases.length());
		ArrayList<Integer> peak_region = new ArrayList<>();
		ArrayList<Mutation> out =  new ArrayList<>();
		ArrayList<Integer> list = new ArrayList<>();
		if (peaks != null) {
			if (peaks.size() == 0) {
				int count = bases.length() / 50000;
				int none_count = count / 50;
				while (count > 0) {
					int s = uid.sample();
					Iterator<Node<Gene>> nodes = gene_region.overlappers(s, s);
					if (nodes.hasNext()) {
						if (--count >= 0) {
							Iterator<Node<Exon>> exon_nodes = exon_region.overlappers(s, s);
							if (exon_nodes.hasNext()) {
								peaks.add(new Peak(s, des[3]));
							}
							else {
								peaks.add(new Peak(s, des[1]));
							}
						}
					}
					else {
						if (--none_count >= 0) {
							peaks.add(new Peak(s, des[2]));
						}
					}
				}
				Collections.sort(peaks);
			}
			for (Peak peak : peaks) {
				list.addAll(randomNoRepeat(peak.getEnd() - read_len, peak.getEnd() + read_len, 4));
				mergerAddRegion(peak_region, new Exon(null, peak.getStart() - read_len, peak.getEnd() + read_len));
				if (Math.random() < 0.5) {
					CircleClip cc = simCircPeak(gene_region, exon_region, peak, Math.random() < 0.625, Math.random() < 0.75, bases);
					if (cc != null) {
						circ.add(cc);
					}
				}
			}
		}
		UniformIntegerDistribution buid = new UniformIntegerDistribution(0, base_units.length - 1);
		for (int pos : list) {
			char c = base_units[buid.sample()];
			while (c == bases.charAt(pos - 1)) {
				c = base_units[buid.sample()];
			}
			out.add(new Mutation(pos, bases.charAt(pos), c, des[0]));
		}
		HashSet<Integer> gene_set = new HashSet<>();
		HashSet<Integer> none_set = new HashSet<>();
		int num = bases.length() / 5000 - list.size();
		int none_num = num / 6;
		num -= none_num;
		while (gene_set.size() < num || none_set.size() < none_num) {
			int s = uid.sample();
			if ((binarySearch(s, peak_region, 0, peak_region.size()) & 1) == 0) {
				if (gene_region.overlappers(s, s).hasNext()) {
					if (gene_set.size() < num && bases.charAt(s - 1) != 'N') {
						gene_set.add(s);
						if (gene_set.size() % 10 == 0) {
							CircleClip cc = simCircPeak(gene_region, exon_region, null, true, Math.random() < 0.75, bases);
							if (cc != null) {
								circ.add(cc);
							}
						}
					}
				}
				else {
					if (none_set.size() < none_num && bases.charAt(s - 1) != 'N') {
						none_set.add(s);
						if (none_set.size() % 10 == 0) {
							CircleClip cc = simCircPeak(gene_region, exon_region, null, true, Math.random() < 0.75, bases);
							if (cc != null) {
								circ.add(cc);
							}
						}
					}
				}
			}
		}
		for (int pos : gene_set) {
			char c = base_units[buid.sample()];
			while (c == bases.charAt(pos - 1)) {
				c = base_units[buid.sample()];
			}
			out.add(new Mutation(pos, bases.charAt(pos - 1), c, des[1]));
		}
		for (int pos : none_set) {
			char c = base_units[buid.sample()];
			while (c == bases.charAt(pos - 1)) {
				c = base_units[buid.sample()];
			}
			out.add(new Mutation(pos, bases.charAt(pos - 1), c, des[2]));
		}
		Collections.sort(out);
		return out;
	}
	
	private static CircleClip simCircPeak(IntervalTree<Gene> gene_region, IntervalTree<Exon> exon_region, Peak peak, boolean in_gene, boolean in_exon, StringBuffer bases){
		StringBuilder descript = new StringBuilder();
		int start = 0;
		int end = 0;
		if (peak != null) {
			descript.append("Peak:");
		}
		else {
			int i = (int) (Math.random() * gene_region.size());
			Node<Gene> node = gene_region.findByIndex(i);
			i = (int) (Math.random() * (node.getEnd() - node.getStart() + 1)) + node.getStart();
			peak = new Peak(i, null);
		}
		start = Math.max(2, peak.getStart() - 2 * args.getReadLen() + 1 + (int) (Math.random() * 2 * args.getReadLen()));
		end = Math.min(start + 101 + (int) (Math.random() * 50000), bases.length() - 2);
		Iterator<Node<Gene>> gene_nodes = gene_region.overlappers(peak.getEnd(), peak.getEnd());
		if (in_gene && gene_nodes.hasNext()) {
			Iterator<Node<Exon>> exon_nodes = exon_region.overlappers(peak.getEnd(), peak.getEnd());
			if (in_exon && exon_nodes.hasNext()) {
				return simCircExon(exon_nodes.next().getValue(), bases, descript);
			}
			descript.append("Gene");
			Gene gene = gene_nodes.next().getValue();
			start = Math.max(start, gene.getStart());
			end = Math.min(end, gene.getEnd());
		}
		return new CircleClip(start, end, descript.toString(), isGTAG(bases.substring(end, end + 2), bases.substring(start - 2, start)));
	}
	
	private static CircleClip simCircExon(Exon exon, StringBuffer bases, StringBuilder descript){
		Transcript script = exon.getScript();
		descript.append("Exon");
		int index = script.getExons().indexOf(exon);
		int left_index = (int) (Math.random() * (index + 1));
		int right_index = (int) (Math.random() * (script.getExons().size() - index)) + index;
		while (!withIn(script.getExons().get(right_index).getEnd() - script.getExons().get(left_index).getStart(), 100, 200000)) {
			left_index = (int) (Math.random() * (index + 1));
			right_index = (int) (Math.random() * (script.getExons().size() - index)) + index;
			if (script.getExonLength() < 400 || script.getExonLength() > 200000) {
				return null;
			}
		}
		return new CircleClip(script.getExons().get(left_index).getStart(), script.getExons().get(right_index).getEnd(), descript.toString(), 
				isGTAG(bases.substring(script.getExons().get(right_index).getEnd(), script.getExons().get(right_index).getEnd() + 2),
						bases.substring(script.getExons().get(left_index).getStart() - 2, script.getExons().get(left_index).getStart())));
	}
	
	private static boolean isGTAG(String gt, String ag) {
		String[] AG = {"AG","AC","AG","AC","AT","GC"};
		String[] GT = {"GT","AT","GC","CT","GT","CT"};
		for (int i = 0; i < AG.length; i++) {
			if (AG[i].equals(ag) && GT[i].equals(gt)) {
				return true;
			}
		}
		return false;
	}
	
	private static <T> HashMap<String, ArrayList<T>> readFile(String file_name, char mode){
		HashMap<String, ArrayList<T>> out = null;
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))){
			switch(mode) {
			case 'e':
				out = readExons(reader);
				break;
			case 's':
				out = readSites(reader);
				break;
			case 'm':
				out = readMutations(reader);
				break;
			case 'p':
				out = readPeaks(reader);
				break;
			default:
				System.err.println("Error Mode Reading file: " + file_name);
			}
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}
	
	@SuppressWarnings("unchecked")
	private static <T> HashMap<String, ArrayList<T>> readExons(BufferedReader reader) throws NumberFormatException, IOException {
		HashMap<String, ArrayList<T>> out = new HashMap<>();
		String line = null;
		Gene gene = null;
		Transcript script = null;
		String gene_id = "gene_id";
		String gene_symbol = "gene_name";
		String script_id = "script_id";
		char sep = '"';
		while ((line = reader.readLine()) != null) {
			String[] cols = line.split("\t");
			if (cols.length < 7) {
				continue;
			}
			String key = cols[2].toLowerCase();
			String chr = cols[0];
			int start = Integer.parseInt(cols[3]) - 1;
			int end = Integer.parseInt(cols[4]);
			char strand = cols[6].charAt(0);
			if ("gene".equals(key)) {
				if (!out.containsKey(chr)) {
					out.put(chr, new ArrayList<>());
				}
				gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, null, start, end, new ArrayList<>());
				out.get(chr).add((T) gene);
			}
			else if ("transcript".equals(key)) {
				if (!out.containsKey(chr)) {
					out.put(chr, new ArrayList<>());
					gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, null, 0, 0, new ArrayList<>());
					out.get(chr).add((T) gene);
				}
				else if (!getIDSym(line, gene_id, sep).equals(gene.getGene_id())){
					gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, null, 0, 0, new ArrayList<>());
					out.get(chr).add((T) gene);
				}
				script = new Transcript(gene, getIDSym(line, script_id, sep), start, end, new ArrayList<>());
				gene.getScripts().add(script);
			}
			else if (!"CDS".equals(key)) {
				if (!out.containsKey(chr)) {
					out.put(chr, new ArrayList<>());
					gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, null, 0, 0, new ArrayList<>());
					out.get(chr).add((T) gene);
					script = new Transcript(gene, getIDSym(line, script_id, sep), 0, 0, new ArrayList<>());
					gene.getScripts().add(script);
				}
				else if (!getIDSym(line, gene_id, sep).equals(gene.getGene_id())){
					gene = new Gene(chr, getIDSym(line, gene_id, sep), getIDSym(line, gene_symbol, sep), strand, null, 0, 0, new ArrayList<>());
					out.get(chr).add((T) gene);
					script = new Transcript(gene, getIDSym(line, script_id, sep), 0, 0, new ArrayList<>());
					gene.getScripts().add(script);
				}
				else if (!getIDSym(line, script_id, sep).equals(script.getScript_id())) {
					script = new Transcript(gene, getIDSym(line, script_id, sep), 0, 0, new ArrayList<>());
					gene.getScripts().add(script);
				}
				Exon exon = new Exon(script, start, end);
				script.getExons().add(exon);
			}
		}
		return out;
	}
	
	private static String getIDSym(String line, String key, char spl_sym) {
		int index = line.indexOf(key);
		if (index != -1) {
			index += key.length() + 2;
			return line.substring(index, line.indexOf(spl_sym, index));
		}
		return "";
	}
	
	@SuppressWarnings("unchecked")
	private static <T> HashMap<String, ArrayList<T>> readSites(BufferedReader reader) throws NumberFormatException, IOException {
		HashMap<String, ArrayList<T>> out = new HashMap<>();
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] cols = line.split("\t");
			String chr = cols[0];
			Site site = new Site(Integer.parseInt(cols[1]));
			if (!out.containsKey(chr)) {
				out.put(chr, new ArrayList<>());
			}
			out.get(chr).add((T) site);
		}
		return out;
	}
	
	@SuppressWarnings("unchecked")
	private static <T> HashMap<String, ArrayList<T>> readPeaks(BufferedReader reader) throws NumberFormatException, IOException {
		HashMap<String, ArrayList<T>> out = new HashMap<>();
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] cols = line.split("\t");
			String chr = cols[0];
			String des = cols.length > 3 ? cols[3] : null;
			Peak peak = new Peak(Integer.parseInt(cols[1]), des);
			if (!out.containsKey(chr)) {
				out.put(chr, new ArrayList<>());
			}
			out.get(chr).add((T) peak);
		}
		return out;
	}
	
	@SuppressWarnings("unchecked")
	private static <T> HashMap<String, IntervalTree<T>> loadCircs(String name) {
		HashMap<String, IntervalTree<T>> out = new HashMap<>();
		if (name == null) {
			return out;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(name)))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				String chr = cols[0];
				CircleClip c = null;
				if (cols.length > 4) {
					c = new CircleClip(Integer.parseInt(cols[1]) - 1, Integer.parseInt(cols[2]), cols[3], Boolean.parseBoolean(cols[4]));
				}
				else {
					c = new CircleClip(Integer.parseInt(cols[1]) - 1, Integer.parseInt(cols[2]), null, false);
				}
				if (!out.containsKey(chr)) {
					out.put(chr, new IntervalTree<>());
				}
				out.get(chr).put(c.getStart() + 1, c.getEnd(), (T) c);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return out;
	}
	
	@SuppressWarnings("unchecked")
	private static <T> HashMap<String, ArrayList<T>> readMutations(BufferedReader reader) throws NumberFormatException, IOException {
		HashMap<String, ArrayList<T>> out = new HashMap<>();
		HashMap<String, String> stat = new HashMap<>();
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] cols = line.split("\t");
			String chr = cols[0];
			char alt = cols.length > 4 ? cols[3].charAt(0) : 'N';
			char snp = cols.length > 4 ? cols[4].charAt(0) : 'N';
			String des = null;
			if (cols.length > 6) {
				if (!stat.containsKey(cols[5])) {
					stat.put(cols[5], cols[5]);
					des = cols[5];
				}
				else {
					des = stat.get(cols[5]);
				}
			}
			Mutation mut = new Mutation(Integer.parseInt(cols[1]), alt, snp, des);
			if (!out.containsKey(chr)) {
				out.put(chr, new ArrayList<>());
			}
			out.get(chr).add((T) mut);
		}
		return out;
	}
	
	private static void appendPeak(HashMap<String, ArrayList<Peak>> peak_info, HashMap<String, ArrayList<Mutation>> muts, HashMap<String, IntervalTree<CircleClip>> circ_info, HashMap<String, ArrayList<Gene>> gene_info){
		BufferedWriter writer_annote = null;
		int dual_len = args.getSegmentLen();
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(args.getGenomeFileName())));
				BufferedWriter writer_input = new BufferedWriter(new FileWriter(new File (args.getInputFileName() + "_L1.fastq"), true));
				BufferedWriter writer_ip = new BufferedWriter(new FileWriter(new File (args.getIPFileName() + "_L1.fastq"), true));
				BufferedWriter writer_input2 = new BufferedWriter(new FileWriter(new File (args.getInputFileName() + "_L2.fastq"), true));
				BufferedWriter writer_ip2 = new BufferedWriter(new FileWriter(new File (args.getIPFileName() + "_L2.fastq"), true))){
			if (args.isAnnoted()) {
				writer_annote = new BufferedWriter(new FileWriter(new File(args.getInputFileName()
						+ "sim_annotation.txt"), true));
			}
			String line = null;
			String chr = null;
			StringBuffer bases = new StringBuffer();
			while ((line = reader.readLine()) != null) {
				if (line.charAt(0) == '>'){
					printNow("appending " + chr);
					int ip_lines = 0;
					int input_lines = 0;
//					if (peak_info.containsKey(chr)) {
//						ArrayList<Peak> peak_list = peak_info.get(chr);
//						for (Peak peak : peak_list) {
//							if (peak.isVisited()) {
//								continue;
//							}
//							String base_seq = bases.substring(peak.getStart() - dual_len - args.getRandSegSize(), peak.getStart() + dual_len + args.getRandSegSize());
//							ArrayList<Integer> mut_points = new ArrayList<>();
//							String mut_seq = getMutBases(base_seq, muts.get(chr), peak.getStart() - dual_len - args.getRandSegSize(), mut_points, 0);
//							if (args.isControled()) {
//								double[] mut_per = {sampleMutPer(args.getMutationProportion())};
//								writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildControlFastqs(base_seq, mut_seq, chr, null, -1,
//										true, peak.getStart() - dual_len - args.getRandSegSize(), mut_per,
//										getMutRegionPercent(0, base_seq.length() - dual_len, mut_points)));
//								writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildControlFastqs(base_seq, mut_seq, chr, null, 
//										dual_len + args.getRandSegSize(), true, peak.getStart() - dual_len - args.getRandSegSize(), mut_per, 
//										getMutRegionPercent(args.getRandSegSize(), dual_len + args.getRandSegSize(), mut_points)));
//								mut_per[0] = sampleMutPerOR(mut_per[0]);
//								writeWithAnnote(writer_input, writer_input2, writer_annote, buildControlFastqs(base_seq, mut_seq, chr, null,
//										-1, false, peak.getStart() - dual_len - args.getRandSegSize(), mut_per,
//										getMutRegionPercent(0, base_seq.length() - dual_len, mut_points)));
//							}
//							else {
//								writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildFastqs(base_seq, mut_seq, chr, null, -1, true, 
//										peak.getStart() - dual_len - args.getRandSegSize(), false));
//								writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildFastqs(base_seq, mut_seq, chr, null,
//										dual_len + args.getRandSegSize(), true, peak.getStart() - dual_len - args.getRandSegSize(), false));
//								writeWithAnnote(writer_input, writer_input2, writer_annote, buildFastqs(base_seq, mut_seq, chr, null, -1,
//										false, peak.getStart() - dual_len - args.getRandSegSize(), false));
//							}
//						}
//					}					
					if (muts.containsKey(chr)) {
						ArrayList<Mutation> mut_list = muts.get(chr);
						for (Mutation mut : mut_list) {
							String base_seq = bases.substring(mut.getStart() - dual_len - args.getRandSegSize(), mut.getStart() + dual_len + args.getRandSegSize());
							if (args.isControled()) {
								ip_lines += writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildMutControlFastqs(base_seq, chr, mut, -1,
										true, mut.getStart() - dual_len - args.getRandSegSize()));
								input_lines += writeWithAnnote(writer_input, writer_input2, writer_annote, buildMutControlFastqs(base_seq, chr, mut, -1,
										false, mut.getStart() - dual_len - args.getRandSegSize()));
							}
						}
					}
//					if (gene_info.containsKey(chr)) {
//						int add_ip = (ip_lines - input_lines) / 2 / 3;
//						if (add_ip > 0) {
//							writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildBackControlFastqs(chr, bases, gene_info.get(chr), muts.get(chr), add_ip, true));
//							writeWithAnnote(writer_input, writer_input2, writer_annote, buildBackControlFastqs(chr, bases, gene_info.get(chr), muts.get(chr), add_ip * 4, false));
//						}
//						ip_lines = 0;
//						input_lines = 0;
//					}
					if (circ_info.containsKey(chr)) {
						IntervalTree<CircleClip> circ_tree = circ_info.get(chr);
						Iterator<Node<CircleClip>> nodes = circ_tree.iterator();
						while (nodes.hasNext()) {
							Node<CircleClip> node = nodes.next();
							if (node.getValue().isVisited()) {
								continue;
							}
							String base_seq = bases.substring(node.getValue().getEnd() - dual_len - args.getRandSegSize(), node.getValue().getEnd()) + 
									bases.substring(node.getValue().getStart(), node.getValue().getStart() + dual_len + args.getRandSegSize());
							ArrayList<Exon> exons = new ArrayList<>(2);
							Transcript script = new Transcript(null, null, node.getValue().getStart(), node.getValue().getEnd(), exons);
							exons.add(new Exon(script, node.getValue().getEnd() - dual_len - args.getRandSegSize(), node.getValue().getEnd()));
							exons.add(new Exon(script, node.getValue().getStart(), node.getValue().getStart() + dual_len + args.getRandSegSize()));
							writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildFastqs(base_seq, base_seq, chr, script,
									dual_len + args.getRandSegSize(), true, 0, true));
							writeWithAnnote(writer_input, writer_input2, writer_annote, buildFastqs(base_seq, base_seq, chr, script, -1,
									false, 0, true));
						}
					}
					String[] temp_chr = line.split("\\s+");
					chr = temp_chr[0].substring(1);
					bases.setLength(0);
				}
				else {
					bases.append(line);
				}
			}
//			if (peak_info.containsKey(chr)) {
//				ArrayList<Peak> peak_list = peak_info.get(chr);
//				for (Peak peak : peak_list) {
//					if (peak.isVisited()) {
//						continue;
//					}
//					String base_seq = bases.substring(peak.getStart() - dual_len - args.getRandSegSize(), peak.getStart() + dual_len + args.getRandSegSize());
//					ArrayList<Integer> mut_points = new ArrayList<>();
//					String mut_seq = getMutBases(base_seq, muts.get(chr), peak.getStart() - dual_len + args.getRandSegSize(), mut_points, 0);
//					if (args.isControled()) {
//						double[] mut_per = {sampleMutPer(args.getMutationProportion())};
//						writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildControlFastqs(base_seq, mut_seq, chr, null, -1,
//								true, peak.getStart() - dual_len - args.getRandSegSize(), mut_per,
//								getMutRegionPercent(0, base_seq.length() - dual_len, mut_points)));
//						writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildControlFastqs(base_seq, mut_seq, chr, null, 
//								dual_len + args.getRandSegSize(), true, peak.getStart() - dual_len - args.getRandSegSize(), mut_per, 
//								getMutRegionPercent(args.getRandSegSize(), dual_len + args.getRandSegSize(), mut_points)));
//						mut_per[0] = sampleMutPerOR(mut_per[0]);
//						writeWithAnnote(writer_input, writer_input2, writer_annote, buildControlFastqs(base_seq, mut_seq, chr, null,
//								-1, false, peak.getStart() - dual_len - args.getRandSegSize(), mut_per,
//								getMutRegionPercent(0, base_seq.length() - dual_len, mut_points)));
//					}
//					else {
//						writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildFastqs(base_seq, mut_seq, chr, null, -1, true, 
//								peak.getStart() - dual_len - args.getRandSegSize(), false));
//						writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildFastqs(base_seq, mut_seq, chr, null,
//								dual_len + args.getRandSegSize(), true, peak.getStart() - dual_len - args.getRandSegSize(), false));
//						writeWithAnnote(writer_input, writer_input2, writer_annote, buildFastqs(base_seq, mut_seq, chr, null, -1,
//								false, peak.getStart() - dual_len - args.getRandSegSize(), false));
//					}
//					
//				}
//			}
			if (muts.containsKey(chr)) {
				ArrayList<Mutation> mut_list = muts.get(chr);
				for (Mutation mut : mut_list) {
					String base_seq = bases.substring(mut.getStart() - dual_len - args.getRandSegSize(), mut.getStart() + dual_len + args.getRandSegSize());
					if (args.isControled()) {
						writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildMutControlFastqs(base_seq, chr, mut, -1,
								true, mut.getStart() - dual_len - args.getRandSegSize()));
						writeWithAnnote(writer_input, writer_input2, writer_annote, buildMutControlFastqs(base_seq, chr, mut, -1,
								false, mut.getStart() - dual_len - args.getRandSegSize()));
					}
				}
			}
			if (circ_info.containsKey(chr)) {
				IntervalTree<CircleClip> circ_tree = circ_info.get(chr);
				Iterator<Node<CircleClip>> nodes = circ_tree.iterator();
				while (nodes.hasNext()) {
					Node<CircleClip> node = nodes.next();
					if (node.getValue().isVisited()) {
						continue;
					}
					String base_seq = bases.substring(node.getValue().getEnd() - dual_len - args.getRandSegSize(), node.getValue().getEnd()) + 
							bases.substring(node.getValue().getStart(), node.getValue().getStart() + dual_len + args.getRandSegSize());
					ArrayList<Exon> exons = new ArrayList<>(2);
					Transcript script = new Transcript(null, null, node.getValue().getStart(), node.getValue().getEnd(), exons);
					exons.add(new Exon(script, node.getValue().getEnd() - dual_len - args.getRandSegSize(), node.getValue().getEnd()));
					exons.add(new Exon(script, node.getValue().getStart(), node.getValue().getStart() + dual_len + args.getRandSegSize()));
					writeWithAnnote(writer_ip, writer_ip2, writer_annote, buildFastqs(base_seq, base_seq, chr, script,
							dual_len + args.getRandSegSize(), true, node.getValue().getEnd() - dual_len - args.getRandSegSize(), true));
					writeWithAnnote(writer_input, writer_input2, writer_annote, buildFastqs(base_seq, base_seq, chr, script, -1,
							false, node.getValue().getEnd() - dual_len - args.getRandSegSize(), true));
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		finally {
			if (writer_annote != null) {
				try {
					writer_annote.close();
				} catch (IOException e) {
				}
			}
		}
	}
	
	public static ArrayList<Fastq> buildFastqs() {
		return null;
	}
	
	public static ArrayList<Fastq> buildFastqs(String bases, String mut_bases, String chr, Transcript script, int peak, boolean ip_flag, int offset, boolean circ_flag) {
		ArrayList<Fastq> out = new ArrayList<>();
		if (bases.length() < args.getSegmentLen()) {
			return out;
		}
		int max_seg = bases.length() > args.getSegmentLen() + args.getRandSegSize() ? args.getSegmentLen() + args.getRandSegSize() : bases.length();
		BinomialDistribution rld = new BinomialDistribution((int) (args.getReadLen() / 0.95 + 0.99), 0.95);
		BinomialDistribution rsld = new BinomialDistribution(args.getSegmentLen() * 2, 0.5);
		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, bases.length() - args.getSegmentLen() + args.getRandSegSize());
		if (peak < 0) {
			int base_num = args.getBackgoundNum() / 2;
			base_num = ip_flag ? base_num / 4 : base_num;
			for (int i = 0; i < base_num; i++) {
				int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
				int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
				int start = sampleInRange(uid, 0, bases.length() - seg_len);
				String id = buildFastqID(chr, script, circ_flag);
				String annote = null;
				if (args.isAnnoted()) {
					StringBuffer anno = new StringBuffer();
					ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, script, start + offset, read_len, ip_flag);
					if (args.isPaired()) {
						mutations.addAll(buildFastqAnno(anno, chr, script, start + offset + seg_len - read_len, read_len, ip_flag));
					}
					anno.append(ensureMut(mutations, ip_flag, false));
					annote = anno.toString();
				}
				out.add(new Fastq(id, bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
				if (args.isPaired()) {
					out.add(new Fastq(id, reverseFliq(bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
				}
			}
			for (int i = 0; i < base_num; i++) {
				int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
				int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
				int start = sampleInRange(uid, 0, bases.length() - seg_len);
				String id = buildFastqID(chr, script, circ_flag);
				String annote = null;
				if (args.isAnnoted()) {
					StringBuffer anno = new StringBuffer();
					ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, script, start + offset, read_len, ip_flag);
					if (args.isPaired()) {
						mutations.addAll(buildFastqAnno(anno, chr, script, start + offset + seg_len - read_len, read_len, ip_flag));
					}
					anno.append(ensureMut(mutations, ip_flag, true));
					annote = anno.toString();
				}
				out.add(new Fastq(id, mut_bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
				if (args.isPaired()) {
					out.add(new Fastq(id, reverseFliq(mut_bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
				}
			}
		}
		else {
			int base_num = args.getBackgoundNum() / 2;
			base_num = ip_flag ? (int) (base_num * args.getEnrichNum()) : base_num;
			for (int i = 0; i < base_num * (1.0 - args.getMutationProportion()); i++) {
				int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
				int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
				int start = sampleInRange(uid, Math.max(0, peak - seg_len + 1), Math.min(peak, bases.length() - seg_len));
				String id = buildFastqID(chr, script, circ_flag);
				String annote = null;
				if (args.isAnnoted()) {
					StringBuffer anno = new StringBuffer();
					ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, script, start + offset, read_len, ip_flag);
					if (args.isPaired()) {
						mutations.addAll(buildFastqAnno(anno, chr, script, start + offset + seg_len - read_len, read_len, ip_flag));
					}
					anno.append(ensureMut(mutations, ip_flag, false));
					annote = anno.toString();
				}
				out.add(new Fastq(id, bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
				if (args.isPaired()) {
					out.add(new Fastq(id, reverseFliq(bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
				}
			}
			for (int i = 0; i < base_num * (1.0 - args.getMutationProportion()); i++) {
				int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
				int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
				int start = sampleInRange(uid, Math.max(0, peak - seg_len + 1), Math.min(peak, bases.length() - seg_len));
				String id = buildFastqID(chr, script, circ_flag);
				String annote = null;
				if (args.isAnnoted()) {
					StringBuffer anno = new StringBuffer();
					ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, script, start + offset, read_len, ip_flag);
					if (args.isPaired()) {
						mutations.addAll(buildFastqAnno(anno, chr, script, start + offset + seg_len - read_len, read_len, ip_flag));
					}
					anno.append(ensureMut(mutations, ip_flag, true));
					annote = anno.toString();
				}
				out.add(new Fastq(id, mut_bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
				if (args.isPaired()) {
					out.add(new Fastq(id, reverseFliq(mut_bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
				}
			}
		}
		return out;
	}
	
	public static ArrayList<Fastq> buildControlFastqs(String bases, String mut_bases, String chr, Transcript script, int peak, boolean ip_flag, int offset, double[] mut_per, double mut_region_per) {
		ArrayList<Fastq> out = new ArrayList<>();
		if (bases.length() < args.getSegmentLen()) {
			return out;
		}
		int max_seg = bases.length() > args.getSegmentLen() + args.getRandSegSize() ? args.getSegmentLen() + args.getRandSegSize() : bases.length();
		BinomialDistribution rld = new BinomialDistribution((int) (args.getReadLen() / 0.95 + 0.99), 0.95);
		BinomialDistribution rsld = new BinomialDistribution(args.getSegmentLen() * 2, 0.5);
		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, bases.length() - args.getSegmentLen() + args.getRandSegSize());
		
		int base_num = args.getBackgoundNum() / 2;
		base_num = ip_flag ? (peak < 0 ? base_num / 4 : (int) (base_num * args.getEnrichNum())) : base_num;
		int mut_num = (int) Math.round((base_num * mut_region_per));
		int no_mut_num = (int) Math.round(mut_num * (1.0 - mut_per[0]));
		base_num = base_num - mut_num;
		mut_num = mut_num - no_mut_num;
		int mut_count = 0;
		int no_mut_count = 0;
		
		int try_time = 1;
		while(base_num > 0 || mut_num > 0 || no_mut_num > 0) {
			int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
			int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
			int start = peak < 0 ? sampleInRange(uid, 0, bases.length() - seg_len) :
				sampleInRange(uid, Math.max(0, peak - seg_len + 1), Math.min(peak, bases.length() - seg_len));
			StringBuffer anno = new StringBuffer();
			ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, script, start + offset, read_len, ip_flag);
			if (args.isPaired()) {
				mutations.addAll(buildFastqAnno(anno, chr, script, start + offset + seg_len - read_len, read_len, ip_flag));
			}
			if (mutations.size() > 0) {
				if (mut_num > 0) {
					--mut_num;
					++mut_count;
					String id = buildFastqID(chr, script, false);
					String annote = null;
					anno.append(ensureMut(mutations, ip_flag, true));
					annote = anno.toString();
					out.add(new Fastq(id, mut_bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
					if (args.isPaired()) {
						out.add(new Fastq(id, reverseFliq(mut_bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
					}
				}
				else if (no_mut_num > 0) {
					--no_mut_num;
					++no_mut_count;
					String id = buildFastqID(chr, script, false);
					String annote = null;
					anno.append(ensureMut(mutations, ip_flag, false));
					annote = anno.toString();
					out.add(new Fastq(id, bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
					if (args.isPaired()) {
						out.add(new Fastq(id, reverseFliq(bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
					}
				}
				else {
					++try_time;
				}
			}
			else if (base_num > 0) {
				--base_num;
				String id = buildFastqID(chr, script, false);
				String annote = null;
				anno.append(ensureMut(mutations, ip_flag, false));
				annote = anno.toString();
				out.add(new Fastq(id, bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
				if (args.isPaired()) {
					out.add(new Fastq(id, reverseFliq(bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
				}
			}
			else {
				++try_time;
			}
			if (try_time % 10000 == 0) {
//				System.out.printf("Warning: try over %d times in %s\n", try_time, script.getScript_id());
				int t = base_num;
				base_num = no_mut_num;
				no_mut_num = t;
				++try_time;
			}
		}
		if (mut_count != 0 || no_mut_count != 0) {
			mut_per[0] = (double) mut_count / (double) (mut_count + no_mut_count);
		}
		return out;
	}
	
	public static ArrayList<Fastq> buildMutControlFastqs(String bases, String chr, Mutation mut, int peak, boolean ip_flag, int offset) {
		ArrayList<Fastq> out = new ArrayList<>();
		if (bases.length() < args.getSegmentLen()) {
			return out;
		}
		int max_seg = bases.length() > args.getSegmentLen() + args.getRandSegSize() ? args.getSegmentLen() + args.getRandSegSize() : bases.length();
		BinomialDistribution rld = new BinomialDistribution((int) (args.getReadLen() / 0.95 + 0.99), 0.95);
		BinomialDistribution rsld = new BinomialDistribution(args.getSegmentLen() * 2, 0.5);
		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, bases.length() - args.getSegmentLen() + args.getRandSegSize());
		
		int try_time = 1;
		while(mut.ensureUsedTime(ip_flag)) {
			int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
			int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
			int start = peak < 0 ? sampleInRange(uid, 0, bases.length() - seg_len) :
				sampleInRange(uid, Math.max(0, peak - seg_len + 1), Math.min(peak, bases.length() - seg_len));
			StringBuffer anno = new StringBuffer();
			ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, null, start + offset, read_len, ip_flag);
			if (args.isPaired()) {
				mutations.addAll(buildFastqAnno(anno, chr, null, start + offset + seg_len - read_len, read_len, ip_flag));
			}
			if (mutations.size() > 0 && mutations.contains(mut)) {
				if (mut.ensureUsedTime(ip_flag, true)) {
					String mut_string = ensureMut(mutations, ip_flag, true);
					if (mut_string == null) {
						++try_time;
						continue;
					}
					String id = buildFastqID(chr, null, false);
					String annote = null;
					anno.append(mut_string);
					annote = anno.toString();
					String mut_bases = getMutBases(bases, mutations, offset);
					out.add(new Fastq(id, mut_bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
					if (args.isPaired()) {
						out.add(new Fastq(id, reverseFliq(mut_bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
					}
				}
				else if (mut.ensureUsedTime(ip_flag, false)) {
					String mut_string = ensureMut(mutations, ip_flag, false);
					if (mut_string == null) {
						++try_time;
						continue;
					}
					String id = buildFastqID(chr, null, false);
					String annote = null;
					anno.append(mut_string);
					annote = anno.toString();
					out.add(new Fastq(id, bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
					if (args.isPaired()) {
						out.add(new Fastq(id, reverseFliq(bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
					}
				}
				else {
					++try_time;
				}
			}
			else {
				++try_time;
			}
			if (try_time % 50000 == 0) {
				System.out.println("debug");
			}
		}
		return out;
	}
	
	public static ArrayList<Fastq> buildBackControlFastqs(String chr, StringBuffer bases, ArrayList<Gene> gene_info, ArrayList<Mutation> mut_info, int num, boolean ip_flag) {
		ArrayList<Fastq> out = new ArrayList<>();
		if (bases.length() < args.getSegmentLen()) {
			return out;
		}
		ArrayList<IntRegion> gene_region = mergerRegion(gene_info);
		ArrayList<Integer> region_length = getRegionSizes(gene_region);
		int effect_length = region_length.get(0);
		for (int i = 1; i < region_length.size(); i++) {
			effect_length += region_length.get(i);
			region_length.set(i, effect_length);
		}
		if (effect_length < args.getSegmentLen()) {
			return out;
		}
		int max_seg = args.getSegmentLen() + args.getRandSegSize();
		BinomialDistribution rld = new BinomialDistribution((int) (args.getReadLen() / 0.95 + 0.99), 0.95);
		BinomialDistribution rsld = new BinomialDistribution(args.getSegmentLen() * 2, 0.5);
		UniformIntegerDistribution uid = new UniformIntegerDistribution(1, effect_length - args.getSegmentLen() + args.getRandSegSize());
		
		int try_time = 1;
		while(num > 0) {
			int seg_len = sampleCutInRange(rsld, args.getSegmentLen() - args.getRandSegSize(), max_seg);
			int read_len = sampleReadLength(rld, args.getReadLen() - args.getRandSegSize(), args.getReadLen());
			int start = sampleInRange(uid, 1, effect_length - seg_len);
			int start_index = binarySearch(start, region_length, 0, region_length.size() - 1);
			if (start_index > 0) {
				start = start - region_length.get(start_index - 1);
			}
			start = start + gene_region.get(start_index).getStart() - 1;
			start_index = binarySearch(start, mut_info, 0, mut_info.size());
			if (start_index < mut_info.size() && start >= mut_info.get(start_index).getStart() - seg_len) {
				if (++try_time % 10000 == 0) {
					System.out.println("debug");
				}
				continue;
			}
			StringBuffer anno = new StringBuffer();
			if (num > 0) {
				ArrayList<Mutation> mutations = buildFastqAnno(anno, chr, null, start, read_len, ip_flag);
				if (args.isPaired()) {
					mutations.addAll(buildFastqAnno(anno, chr, null, start + seg_len - read_len, read_len, ip_flag));
				}
				String id = buildFastqID(chr, null, false);
				String annote = null;
				annote = anno.toString();
				out.add(new Fastq(id, bases.substring(start, start + read_len), buildString(read_len, args.getQuality()), annote));
				if (args.isPaired()) {
					out.add(new Fastq(id, reverseFliq(bases.substring(start + seg_len - read_len, start + seg_len)), buildString(read_len, args.getQuality()), annote));
				}
				--num;
			}
		}
		return out;
	}
	
	public static String ensureMut(ArrayList<Mutation> mutations, boolean ip_flag, boolean mut_flag) {
		if (mutations == null || mutations.size() == 0) {
			return "";
		}
		StringBuilder sb = new StringBuilder();
		sb.append('\t');
		for (int i = 0; i < mutations.size(); i++) {
			if (mutations.get(i).ensureUsedTime(ip_flag, mut_flag)) {
				mutations.get(i).reduceUsedTime(ip_flag, mut_flag);
				if (mut_flag) {
					sb.append(mutations.get(i).getStart() + 1);
					sb.append(mutations.get(i).getMut());
					sb.append(':');
				}
			}
			else {
				return null;
			}
		}
		sb.setLength(sb.length() - 1);
		return sb.toString();
	}
	
	public static double getMutRegionPercent(int start, int end, ArrayList<Integer> points) {
		int last_end = start;
		int sum = 0;
		int region = Math.min(args.getSegmentLen(), 2 * args.getReadLen());
		for (int i = 0; i < points.size(); i++) {
			if (points.get(i) > start) {
				if (points.get(i) >= end) {
					sum += points.get(i) - end >= region ? 0 : points.get(i) - last_end > region ? region - points.get(i) + end : end - last_end;
					break;
				}
				sum += Math.min(points.get(i) - last_end, region);
			}
			last_end = Math.max(last_end, points.get(i));
		}
		return (double) sum / (double) (end - start);
	}
	
	private static void makeMutTimes() {
		for (Entry<String, ArrayList<Mutation>> entry : muts.entrySet()) {
			Mutation last_mut = null;
			for (Mutation mut : entry.getValue()) {
				if (last_mut != null && mut.getEnd() <= last_mut.getEnd() + args.getReadLen()) {
					mut.cloneTimes(last_mut);
				}
				else {
					mut.makeMutTimes(args);
				}
				last_mut = mut;
			}
		}
	}
	
	public static void ensureMutDes(HashMap<String, ArrayList<Peak>> peak_info) {
		for (Entry<String, ArrayList<Mutation>> entry : muts.entrySet()) {
			for (Mutation mut : entry.getValue()) {
				if (peak_info.containsKey(entry.getKey())) {
					ArrayList<Peak> peak_list = peak_info.get(entry.getKey());
					int index = binarySearch(mut.getStart(), peak_list, 0, peak_list.size());
					if ((index > 0 && peak_list.get(index - 1).getStart() + args.getReadLen() >= mut.getStart())
							|| (index < peak_list.size() && peak_list.get(index).getStart() - args.getReadLen() <= mut.getStart())) {
						if (!"Peak".equals(mut.getDescription())){
							mut.setDescription("Peak");
						}
					}
					else if ("Peak".equals(mut.getDescription())){
						mut.setDescription("None");
					}
				}
				else if ("Peak".equals(mut.getDescription())){
					mut.setDescription("None");
				}
			}
		}
	}
	
	public static String reverseFliq(String bases) {
		if (bases == null || bases.length() == 0) {
			return bases;
		}
		HashMap<Character, Integer> unit_map = new HashMap<>();
		for (int i = 0; i < base_units.length; i++) {
			unit_map.put(base_units[i], i + 1);
		}
		StringBuilder sb = new StringBuilder();
		for (int i = bases.length() - 1; i >= 0; i--) {
			sb.append(unit_map.containsKey(bases.charAt(i)) ? base_units[base_units.length - unit_map.get(bases.charAt(i))] : bases.charAt(i));
		}
		return sb.toString();
	}
	
	public static String getMutBases(String bases, ArrayList<Mutation> muts, int start) {
		StringBuffer sb = new StringBuffer();
		sb.append(bases);
		for (int i = 0; i < muts.size(); i++) {
			muts.get(i).visited();
			sb.setCharAt(muts.get(i).getStart() - start, muts.get(i).getMut());
		}
		return sb.toString();
	}
	
	public static String getMutBases(String bases, ArrayList<Mutation> muts, int start, ArrayList<Integer> mut_points, int point_start) {
		if (bases == null || muts == null) {
			return bases;
		}
		StringBuffer sb = new StringBuffer();
		sb.append(bases);
		int index = binarySearch(start, muts, 0, muts.size());
		for (int i = index; i < muts.size() && muts.get(i).getStart() < start + bases.length(); i++) {
			muts.get(i).visited();
			sb.setCharAt(muts.get(i).getStart() - start, muts.get(i).getMut());
			if (mut_points != null) {
				mut_points.add(muts.get(i).getStart() - start + point_start);
			}
		}
		return sb.toString();
	}
	
	private static String buildFastqID(String chr, Transcript script, boolean circ) {
		StringBuffer out = new StringBuffer("@");
		out.append(chr);
//		out.append(':');
//		out.append(script == null ? "None" : script.getScript_id());
		if (circ) {
			out.append(":Circ");
		}
		makeUniqueID(out);
		return out.toString();
	}
	
	private static ArrayList<Mutation> buildFastqAnno(StringBuffer anno, String chr, Transcript script, int start, int len, boolean ip_flag) {
		ArrayList<Mutation> out = new ArrayList<>();
		char spl = '\t';
		if (anno.length() <= 0) {
			anno.append(id_count + 1);
			anno.append(spl);
			anno.append(chr);
			anno.append(spl);
			anno.append(script == null ? "None" : script.getScript_id());
		}
		if (script == null) {
			anno.append(spl);
			anno.append(start);
			anno.append('-');
			anno.append(start + len);
			int left_index = binarySearch(start, muts.get(chr), 0, muts.get(chr).size());
			int right_index = binarySearch(start + len, muts.get(chr), 0, muts.get(chr).size());
			for (int i = left_index; i < right_index; i++) {
				if (muts.get(chr).get(i).ensureUsedTime(ip_flag)) {
					out.add(muts.get(chr).get(i));
				}
			}
		}
		else {
			int index = 0;
			while (start >= script.getExons().get(index).getLength()) {
				start -= script.getExons().get(index).getLength();
				index++;
				index = index % script.getExons().size();
			}
			while (len > 0) {
				anno.append(spl);
				anno.append(start + script.getExons().get(index).getStart() + 1);
				anno.append('-');
				int end = len >= script.getExons().get(index).getLength() - start ? script.getExons().get(index).getEnd() : script.getExons().get(index).getStart() + start + len;
				anno.append(end);
				len -= script.getExons().get(index).getLength() - start;
				int left_index = binarySearch(start + script.getExons().get(index).getStart(), muts.get(chr), 0, muts.get(chr).size());
				int right_index = binarySearch(end, muts.get(chr), 0, muts.get(chr).size());
				for (int i = left_index; i < right_index; i++) {
					if (muts.get(chr).get(i).ensureUsedTime(ip_flag)) {
						out.add(muts.get(chr).get(i));
					}
				}
				start = 0;
				++index;
				index = index % script.getExons().size();
				spl = ':';
			}
		}
		return out;
	}
	
	private static void makeUniqueID(StringBuffer id) {
		if (id == null) {
			return;
		}
		id.append(':');
		id.append(++id_count);
	}
	
	public static <T> void writeTreeFile(String file_name, HashMap<String, IntervalTree<T>> tree) {
		File fo = new File (file_name);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo))) {
			for (Entry<String, IntervalTree<T>> entry : tree.entrySet()){
				Iterator<Node<T>> nodes = entry.getValue().overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
				while (nodes.hasNext()) {
					Node<T> node = nodes.next();
					writer.write(entry.getKey());
					writer.write('\t');
					writer.write(node.getValue().toString());
					writer.newLine();
				}
				writer.flush();
			}
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void writeFile(String file_name, ArrayList<T> list) {
		File fo = new File (file_name);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo))) {
			for (T t : list){
				writer.write(t.toString());
				writer.newLine();
			}
			writer.flush();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void writeUsedFile(String file_name, HashMap<String, ArrayList<Mutation>> map, int time) {
		File fo = new File (file_name);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo))) {
			for (Entry<String, ArrayList<Mutation>> entry : map.entrySet()){
				for (int i = 0; i < entry.getValue().size(); i++) {
					if (entry.getValue().get(i).getLeftTime(false, true) >= time) {
						writer.write(entry.getKey());
						writer.write('\t');
						writer.write(entry.getValue().get(i).toString());
						writer.write('\n');
					}
				}			
			}
			writer.flush();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void writeUsedFile(String file_name, String ip_file, HashMap<String, ArrayList<Mutation>> map, int time) {
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(new File (file_name)));
				BufferedWriter writer2 = new BufferedWriter(new FileWriter(new File (ip_file)))) {
			for (Entry<String, ArrayList<Mutation>> entry : map.entrySet()){
				for (int i = 0; i < entry.getValue().size(); i++) {
					if (entry.getValue().get(i).getLeftTime(false, true) >= time) {
						writer.write(entry.getKey());
						writer.write('\t');
						writer.write(entry.getValue().get(i).toString());
						writer.write('\n');
					}
					if (entry.getValue().get(i).getLeftTime(true, true) >= time) {
						writer2.write(entry.getKey());
						writer2.write('\t');
						writer2.write(entry.getValue().get(i).toString());
						writer2.write('\n');
					}
				}			
			}
			writer.flush();
			writer2.flush();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T extends Site> void writeVisitedFile(String file_name, HashMap<String, ArrayList<T>> map) {
		File fo = new File (file_name);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo))) {
			for (Entry<String, ArrayList<T>> entry : map.entrySet()){
				for (int i = 0; i < entry.getValue().size(); i++) {
					if (entry.getValue().get(i).isVisited()) {
						writer.write(entry.getKey());
						writer.write('\t');
						writer.write(entry.getValue().get(i).toString());
						writer.write('\n');
					}
				}			
			}
			writer.flush();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void writeFile(String file_name, HashMap<String, ArrayList<T>> map) {
		File fo = new File (file_name);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo))) {
			for (Entry<String, ArrayList<T>> entry : map.entrySet()){
				for (int i = 0; i < entry.getValue().size(); i++) {
					writer.write(entry.getKey());
					writer.write('\t');
					writer.write(entry.getValue().get(i).toString());
					writer.write('\n');
				}			
			}
			writer.flush();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void appendFile(String file_name, ArrayList<T> list) {
		File fo = new File (file_name);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo, true))){
			for (T t : list){
				writer.write(t.toString());
				writer.write('\n');
			}
			writer.flush();
			writer.close();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void appendFile(String file_name, String file_name2, ArrayList<T> list) {
		File fo = new File (file_name);
		File fo2 = new File (file_name2);
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(fo, true));
				BufferedWriter writer2 = new BufferedWriter(new FileWriter(fo2, true))){
			for (int i = 0; i + 1 < list.size(); ++i){
				writer.write(list.get(i).toString());
				writer.write('\n');
				writer2.write(list.get(++i).toString());
				writer2.write('\n');
			}
			writer.flush();
			writer2.flush();
		}
		catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static <T> void writeInTwo(BufferedWriter writer1, BufferedWriter writer2, ArrayList<T> list) throws IOException {
		for (int i = 0; i + 1 < list.size(); ++i){
			writer1.write(list.get(i).toString());
			writer1.newLine();
			writer2.write(list.get(++i).toString());
			writer2.newLine();
		}
		writer1.flush();
		writer2.flush();
	}
	
	public static <T> void write(BufferedWriter writer, ArrayList<T> list) throws IOException {
		for (int i = 0; i + 1 < list.size(); ++i){
			writer.write(list.get(i).toString());
			writer.newLine();
		}
		writer.flush();
	}
	
	public static void writeAnnote(BufferedWriter writer, ArrayList<Fastq> list) throws IOException {
		for (int i = 0; i + 1 < list.size(); ++i){
			writer.write(list.get(i).getAnnotation());
			writer.newLine();
			if (args.isPaired()) {
				++i;
			}
		}
		writer.flush();
	}
	
	public static int writeWithAnnote(BufferedWriter writer1, BufferedWriter writer2, BufferedWriter writer, ArrayList<Fastq> list) throws IOException {
		writeInTwo(writer1, writer2, list);
		if (writer != null) {
			writeAnnote(writer, list);
		}
		writer.flush();
		return list.size();
	}
	
	private static void appendSimFile(InParam args, HashMap<String, ArrayList<Gene>> gene_info, HashMap<String, ArrayList<Peak>> peak_info, HashMap<String, ArrayList<Mutation>> mut_info, HashMap<String, IntervalTree<CircleClip>> circ_info) {
		BufferedWriter writer_annote = null;
		try (BufferedWriter writer_input = new BufferedWriter(new FileWriter(new File (args.getInputFileName() + "_L1.fastq")));
				BufferedWriter writer_ip = new BufferedWriter(new FileWriter(new File (args.getIPFileName() + "_L1.fastq")));
				BufferedWriter writer_input2 = new BufferedWriter(new FileWriter(new File (args.getInputFileName() + "_L2.fastq")));
				BufferedWriter writer_ip2 = new BufferedWriter(new FileWriter(new File (args.getIPFileName() + "_L2.fastq")))){
			if (args.isAnnoted()) {
				writer_annote = new BufferedWriter(new FileWriter(new File(args.getInputFileName()
						+ "sim_annotation.txt")));
			}
			for (Entry<String, ArrayList<Gene>> entry : gene_info.entrySet()) {
				printNow("Simulate in " + entry.getKey());
				for (int i = 0; i < entry.getValue().size(); i++) {
					for (int j = 0; j < entry.getValue().get(i).getScripts().size(); j++) {
						Transcript script = entry.getValue().get(i).getScripts().get(j);
						double[] mut_per = {sampleMutPer(args.getMutationProportion())};
						writeWithAnnote(writer_input, writer_input2, writer_annote, entry.getValue().get(i).getScripts().get(j).
								simulate(args, peak_info.get(entry.getKey()), mut_info.get(entry.getKey()), false, mut_per));
						mut_per[0] = sampleMutPerOR(mut_per[0]);
						writeWithAnnote(writer_ip, writer_ip2, writer_annote, entry.getValue().get(i).getScripts().get(j).
								simulate(args, peak_info.get(entry.getKey()), mut_info.get(entry.getKey()), true, mut_per));
						if (!circ_info.containsKey(entry.getKey())) {
							continue;
						}
						Iterator<Node<CircleClip>> nodes = circ_info.get(entry.getKey()).overlappers(script.getStart(), script.getEnd());
						while (nodes.hasNext()) {
							Node<CircleClip> node = nodes.next();
							writeWithAnnote(writer_input, writer_input2, writer_annote, script.
									simulateCirc(args, peak_info.get(entry.getKey()), node.getValue(), false));
						}
					}
				}
			}
		}
		catch(IOException e) {
			e.printStackTrace();
		}
		finally {
			if (writer_annote != null) {
				try {
					writer_annote.close();
				} catch (IOException e) {
				}
			}
		}
	}
	
	public static HashSet<Integer> randomNoRepeat(int start, int end, int num) {
		if (num > end - start) {
			return null;
		}
		HashSet<Integer> out = new HashSet<>();
		UniformIntegerDistribution uid = new UniformIntegerDistribution(start, end);
		while (out.size() < num) {
			out.add(uid.sample());
		}
		return out;
	}
	
	public static int binarySearch(int target, int[] list, int start, int end) {
		if (start >= end) {
			return start;
		}
		int mid = (start + end) >> 1;
		return target > list[mid] ? binarySearch(target, list, mid + 1, end) : binarySearch(target, list, start, mid);
	}

	public static int binarySearchToMiddle(double target, double[] list) {
		int index = binarySearch(target, list, 0, list.length);
		if (index >= list.length / 2) {
			--index;
		}
		return index;
	}
	
	public static int binarySearch(double target, double[] list, int start, int end) {
		if (start >= end) {
			return start;
		}
		int mid = (start + end) >> 1;
		return target > list[mid] ? binarySearch(target, list, mid + 1, end) : binarySearch(target, list, start, mid);
	}
	
	public static int binarySearch(int target, ArrayList<? extends IntRegion> list, int start, int end) {
		if (start >= end) {
			return start;
		}
		int mid = (start + end) >> 1;
		return target > list.get(mid).getStart() ? binarySearch(target, list, mid + 1, end) : binarySearch(target, list, start, mid);
	}

	public static <T extends Comparable<? super T>> int binarySearch(T target, ArrayList<T> list, int start, int end) {
		if (start >= end) {
			return start;
		}
		int mid = (start + end) >> 1;
		return target.compareTo(list.get(mid)) > 0 ? binarySearch(target, list, mid + 1, end) : binarySearch(target, list, start, mid);
	}
	
	public static <T extends IntRegion> ArrayList<IntRegion> mergerRegion(ArrayList<T> list){
		if (list == null || list.size() == 0) {
			return null;
		}
		ArrayList<IntRegion> out = new ArrayList<>();
		Collections.sort(list);
		out.add(new IntRegion(list.get(0).getStart(), list.get(0).getEnd()));
		IntRegion last_region = out.get(0);
		for (int i = 1; i < list.size(); i++) {
			IntRegion ir = list.get(i);
			if (ir.getStart() <= last_region.getEnd()) {
				if (ir.getEnd() > last_region.getEnd()) {
					last_region.resetStartAndEnd(last_region.getStart(), ir.getEnd());
				}
			}
			else {
				out.add(last_region = new IntRegion(ir.getStart(), ir.getEnd()));
			}
		}
		return out;
	}
	
	public static <T extends IntRegion> ArrayList<Integer> getRegionSizes(ArrayList<T> list){
		if (list == null || list.size() == 0) {
			return null;
		}
		ArrayList<Integer> out = new ArrayList<>();
		Collections.sort(list);
		for (int i = 0; i < list.size(); i++) {
			out.add(list.get(i).getLength());
		}
		return out;
	}
	
	public static void putInExonTree(IntervalTree<Exon> exon_tree, Gene gene) {
		if (gene == null) {
			return;
		}
		for (Transcript script : gene.getScripts()) {
			for (Exon exon : script.getExons()) {
				exon_tree.put(exon.getStart(), exon.getEnd() - 1, exon);
			}
		}
	}
	
	public static <T extends IntRegion> void mergerAddRegion(ArrayList<Integer> out, T region){
		int index_left = binarySearch(region.getStart(), out, 0, out.size());
		int index_right = binarySearch(region.getEnd(), out, 0, out.size());
		if (index_left == index_right) {
			if ((index_left & 1) == 0) {
				if (index_right < out.size() && out.get(index_right) == out.get(index_right + 1)) {
					out.set(index_left, region.getStart());
				}
				else {
					out.add(index_left, region.getEnd());
					out.add(index_left, region.getStart());
				}
			}
		}
		else {
			if ((index_left & 1) == 0) {
				out.set(index_left, region.getStart());
				if ((index_right & 1) == 0) {
					if (index_right < out.size() && out.get(index_right) == region.getEnd()) {
						out.remove(index_right);
						out.remove(--index_right);
					}
					else {
						out.set(--index_right, region.getEnd());
					}
				}
				for (int i = index_right - 1; i > index_left; i--) {
					out.remove(i);
				}
			}
			else {
				if ((index_right & 1) == 0) {
					if (index_right < out.size() && out.get(index_right) == region.getEnd()) {
						out.remove(index_right);
						out.remove(--index_right);
					}
					else {
						out.set(--index_right, region.getEnd());
					}
				}
				for (int i = index_right - 1; i >= index_left; i--) {
					out.remove(i);
				}
			}
		}
	}
	
	public static boolean checkAllGene(HashMap<String, ArrayList<Gene>> gene_info) {
		boolean out = true;
		for (Entry<String, ArrayList<Gene>> genes : gene_info.entrySet()) {
			for (Gene gene : genes.getValue()) {
				out &= gene.checkPos();
			}
		}
		return out;
	}
	
	public static String buildString(int length, String chars) {
		if (chars == null || chars.length() < 1) {
			return null;
		}
		if (chars.length() == 1) {
			if (FixMapQ != null){
				if (length <= FixMapQ.length()) {
					return FixMapQ.substring(0, length);
				}
				StringBuilder sb = new StringBuilder();
				length -= FixMapQ.length();
				sb.append(FixMapQ);
				while (--length >= 0) {
					sb.append(chars);
				}
				return FixMapQ = sb.toString();
			}
			StringBuilder sb = new StringBuilder();
			while (--length >= 0) {
				sb.append(chars);
			}
			return FixMapQ = sb.toString();
		}
		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, chars.length() - 1);
		StringBuffer sb = new StringBuffer();
		while(--length >= 0) {
			sb.append(chars.charAt(uid.sample()));
		}
		return sb.toString();
	}
	
	public static int sampleRealReadLength() {
		return binarySearch(real_read_dis.sample(), real_read, 0, real_read.length);
	}
	
	public static <T extends AbstractIntegerDistribution> int sampleReadLength(T aid, int min_len, int max_len) {
		int out = aid.sample();
		while (out < min_len) {
			out = aid.sample();
		}
		out = out > max_len ? max_len : out;
		return out;
	}
	
	public static <T extends AbstractIntegerDistribution> int sampleCutInRange(T aid, int lower, int upper) {
		int s = aid.sample();
		while (s < lower || s > upper) {
			s = aid.sample();
		}
		return s;
	}
	
	public static <T extends AbstractIntegerDistribution>  int sampleInRange(T aid, int lower, int upper) {
		int length = upper - lower + 1;
		int multi = (aid.getSupportUpperBound() - aid.getSupportLowerBound() + 1) / length;	
		if (multi == 0) {
			System.err.println("Cannot sample!");
			return 0;
		}
		int s = aid.sample();
		while (s > length * multi) {
			s = aid.sample();
		}
		s = lower + (s - aid.getSupportLowerBound()) % length;
		return s;
	}
	
	public static double sampleMutPer(double mut_per) {
		double s = mut_per_dis.sample();
		while(s < - mut_per || s > 1.0 - mut_per) {
			s = mut_per_dis.sample();
		}
		return s + mut_per;
	}
	
	public static double sampleMutPerOR(double mut_per) {
		return 1.0 - 1.0 / ( 1.0 + sampleORValue() * mut_per / (1.0 - mut_per));
	}
	
	private static double sampleORValue() {
		return Math.exp(Math.floor((OR_dis.sample() + 0.15) / 0.3) * 0.3);
	}
	
	public static <T extends AbstractRealDistribution> double sampleCutInRange(T ard, double lower, double upper) {
		double s = ard.sample();
		while (s < lower || s > upper) {
			s = ard.sample();
		}
		return s;
	}
	
	public static boolean withIn(int target, int lower, int upper) {
		if (lower > upper) {
			lower ^= upper;
			upper ^= lower;
			lower ^= upper;
		}
		return target >= lower && target <= upper;
	}
	
	public static InParam getParams() {
		return args;
	}
	
	public static void printNow(String prefix) {
		Date time = new Date();
		System.out.printf("%s at %tF %tT\n", prefix, time, time);
	}
	
	public static String toString(Object obj) {
		if (obj == null) {
			return null;
		}
		StringBuilder sb = new StringBuilder();
		HashSet<Object> visited = new HashSet<>();
		visited.add(null);
		toString(sb, obj, visited);
		return sb.toString();
	}
	
	private static void toString(StringBuilder sb, Object obj, HashSet<Object> visited) {
		if (obj == null) {
			sb.append("null");
		}
		if (visited.add(obj)) {
			Class<?> class1 = obj.getClass();
			if (class1 == String.class) {
				sb.append(obj);
			}
			else if (class1.isArray()) {
				sb.append(class1.getComponentType());
				sb.append("[]{");
				if (Array.getLength(obj) > 0) {
					if (class1.getComponentType().isPrimitive()) {
						sb.append(Array.get(obj, 0));
					}
					else {
						toString(sb, Array.get(obj, 0), visited);
					}
				}
				for (int i = 1; i < Array.getLength(obj); i++) {
					sb.append(',');
					if (class1.getComponentType().isPrimitive()) {
						sb.append(Array.get(obj, i));
					}
					else {
						toString(sb, Array.get(obj, i), visited);
					}
				}
				sb.append("}\n");
			}
			else {
				sb.append(class1.getName());
				do {
					sb.append('[');
					Field[] fields = class1.getDeclaredFields();
					AccessibleObject.setAccessible(fields, true);
					for (Field field : fields) {
						if (!Modifier.isStatic(field.getModifiers())) {
							if (sb.charAt(sb.length() - 1) != '[') {
								sb.append(',');
							}
							sb.append(field.getName());
							sb.append('=');
							try {
								Object val = field.get(obj);
								if (field.getType().isPrimitive()) {
									sb.append(val);
								}
								else {
									toString(sb, val, visited);
								}
							} catch (Exception e) {
								e.printStackTrace();
							}
							
						}
					}
					sb.append(']');
					class1 = class1.getSuperclass();
				}while(class1 != null);
			}
		}
		else {
			sb.append("...");
		}
	}
}
