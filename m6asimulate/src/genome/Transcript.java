package genome;

import java.util.ArrayList;
import java.util.Collections;

import note.CircleClip;
import note.Fastq;
import note.Mutation;
import note.Site;
import sim.*;

public class Transcript extends IntRegion{

	private Gene gene = null;
	private String script_id = null;
	private ArrayList<Exon> exons = null;
	private int exon_length = 0;
	private String exon_bases = null;
	
	public Transcript(Gene gene, String script_id, int start, int end, ArrayList<Exon> exons) {
		super(start, end);
		this.gene = gene;
		this.script_id = script_id;
		this.exons = exons;
	}

	public String getChr() {
		return gene==null? null : gene.getChr_symbol();
	}
	
	public Gene getGene() {
		return gene;
	}
	
	public int getExonLength() {
		return exon_length;
	}
	
	public String getScript_id() {
		return script_id;
	}

	public ArrayList<Exon> getExons() {
		return exons;
	}
	
	public String getBases() {
		return gene.getBases() == null ? null : gene.getBases().substring(getStart() - gene.getStart(), getEnd() - gene.getStart());
	}
	
	public String getExonBases(String bases) {
		if (exon_bases != null) {
			return exon_bases;
		}
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < exons.size(); i++) {
			sb.append(bases.substring(exons.get(i).getStart(), exons.get(i).getEnd()));
		}
		return exon_bases = sb.toString();
	}
	
	public int getExonIndex(int base_index) {
		int exon_index = Method.binarySearch(base_index + 1, exons, 0, exons.size());
		if (exon_index > 0 && base_index < exons.get(exon_index - 1).getEnd()) {
			return exon_index - 1;
		}
		return -1;
	}
	
	public int getExonSeqIndex(int base_index, int seq_index) {
		int exon_index = getExonIndex(base_index);
		if (exon_index < 0) {
			return Integer.MIN_VALUE;
		}
		while (seq_index > base_index - exons.get(exon_index).getStart()) {
			seq_index -= base_index - exons.get(exon_index).getStart();
			if (--exon_index >= 0) {
				base_index = exons.get(exon_index).getEnd();
			}
			else {
				break;
			}
		}
		if (exon_index >= 0) {
			base_index = base_index - seq_index;
			seq_index = 0;
			for (int i = 0; i < exon_index; i++) {
				seq_index += exons.get(i).getLength();
			}
			seq_index += base_index - exons.get(exon_index).getStart();
			return seq_index;
		}
		return -seq_index;
	}
	
	public String getExonBases(ArrayList<Mutation> muts, ArrayList<Integer> mut_points) {
		if (gene.getBases() == null) {
			return null;
		}
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < exons.size(); i++) {
			sb.append(Method.getMutBases(gene.getBases().substring(exons.get(i).getStart() - gene.getStart(), exons.get(i).getEnd() - gene.getStart()), muts, exons.get(i).getStart(), mut_points, sb.length()));
		}
		return sb.toString();
	}
	
	public int inExons(int pos) {
		int out = 0;
		if (exons != null) {
			for (int i = 0; i < exons.size(); i++) {
				if (pos < exons.get(i).getEnd() && pos >= exons.get(i).getStart()) {
					out += pos - exons.get(i).getStart();
					return out;
				}
				out += exons.get(i).getEnd() - exons.get(i).getStart();
			}
		}
		return -2;
	}
	
	public int inExons(Site site) {
		int out = 0;
		int pos = site.getStart();
		if (exons != null) {
			for (int i = 0; i < exons.size(); i++) {
				if (pos < exons.get(i).getEnd() && pos >= exons.get(i).getStart()) {
					out += pos - exons.get(i).getStart();
					site.visited();
					return out;
				}
				out += exons.get(i).getEnd() - exons.get(i).getStart();
			}
		}
		return -1;
	}
	
	public ArrayList<Fastq> simulateCirc(InParam args, ArrayList<? extends Site> peaks, CircleClip circ, boolean ip_flag){
		ArrayList<Fastq> out = new ArrayList<>();
		if (circ != null) {
			ArrayList<Exon> old_exons = exons;
			int left_index = Method.binarySearch(circ.getStart(), exons, 0, exons.size());
			int right_index = Method.binarySearch(circ.getEnd(), exons, 0, exons.size());
			exons = new ArrayList<>(old_exons.subList(left_index, right_index));
			if (exons.size() > 0 && circ.getStart() == exons.get(0).getStart() && circ.getEnd() == exons.get(exons.size() - 1).getEnd()) {
				String base_seq = getExonBases(null, null);
				if (base_seq != null) {
					int junction = base_seq.length();
					base_seq += base_seq;
					circ.visited();
					out = Method.buildFastqs(base_seq, base_seq, getChr(), this, junction, ip_flag, 0, true);
				}
			}
			exons = old_exons;
		}
		return out;
	}
	
	public ArrayList<Fastq> simulate(InParam args, ArrayList<? extends Site> peaks, ArrayList<Mutation> muts, boolean ip_flag, double[] mut_per){
		ArrayList<Fastq> out = new ArrayList<>();
		String base_seq = getExonBases(null, null);
		if (base_seq == null) {
			return out;
		}
		ArrayList<Integer> mut_points = new ArrayList<>();
		String mut_seq = getExonBases(muts, mut_points);
		if (Method.getParams().isControled()) {
			out = Method.buildControlFastqs(base_seq, mut_seq, getChr(), this, -1, ip_flag, 0, 
					mut_per, Method.getMutRegionPercent(0, base_seq.length() - args.getSegmentLen() + 1, mut_points));
		}
		else {
			out = Method.buildFastqs(base_seq, mut_seq, getChr(), this, -1, ip_flag, 0, false);
		}
		if (peaks == null) {
			return out;
		}
		int peak_left = Method.binarySearch(getStart(), peaks, 0, peaks.size());
		int peak_right = Method.binarySearch(getEnd(), peaks, 0, peaks.size());
		for (int i = peak_left; ip_flag && i < peak_right; i++) {
			int offset = inExons(peaks.get(i));
			if (offset < 0) {
				continue;
			}
			int left = Math.max(0, offset - args.getSegmentLen() - args.getRandSegSize() + 1);
			int peak = Math.min(offset, args.getSegmentLen() + args.getRandSegSize() - 1);
			int right = Math.min(base_seq.length(), offset + args.getSegmentLen() + args.getRandSegSize());
			if (Method.getParams().isControled()) {
				out.addAll(Method.buildControlFastqs(base_seq.substring(left, right), mut_seq.substring(left, right), getChr(), this, peak, ip_flag, left, 
							mut_per, Method.getMutRegionPercent(Math.max(0, offset - args.getSegmentLen()), Math.min(base_seq.length(), offset + args.getSegmentLen()) - args.getSegmentLen() + 1, mut_points)));
			}
			else {
				out.addAll(Method.buildFastqs(base_seq.substring(left, right), mut_seq.substring(left, right), getChr(), this, peak, ip_flag, left, false));
			}
		}
		return out;
	}

	public boolean checkPos() {
		boolean out = true;
		if (exons != null) {
			exon_length = 0;
			Collections.sort(exons);
			for (int i = exons.size() - 1; i >= 0; i--) {
				out &= getStart() <= exons.get(i).getStart()
						&& getEnd() >= exons.get(i).getEnd();
				exon_length += exons.get(i).getLength();
				if (i > 0 && exons.get(i).getStart() <= exons.get(i - 1).getEnd()) {
					exons.get(i - 1).resetStartAndEnd(exons.get(i - 1).getStart(), Math.max(exons.get(i).getEnd(), exons.get(i - 1).getEnd()));
					exons.remove(i);
				}
			}
			resetStartAndEnd(getStart() > exons.get(0).getStart() ? exons.get(0).getStart() : getStart(), 
					getEnd() < exons.get(exons.size() - 1).getEnd() ? exons.get(exons.size() - 1).getEnd() : getEnd());
		}
		return out;
	}
}
