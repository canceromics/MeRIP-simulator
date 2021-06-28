package note;

import java.util.ArrayList;

import genome.Exon;
import genome.Transcript;
import sim.InParam;
import sim.Method;

public class Fastq {
	
	private String id = "@";
	private String seq = null;
	private String add_info = "+";
	private String quali = null;
	private String annote = null;
	private int len = 0;
	private Boolean mut_flag = null;
	private Fastq pair = null;
	
	public static void main(String[] args) {
		InParam.getParams().putParams(args);
		Fastq fq = new Fastq("@t", 15);
		String s = "123456789012345678901234567890";
		
		int i = fq.buildSeq(s, 0, 0, null);
		i = fq.buildSeq(s, 0 + i, 10, null);
		
		Transcript script = new Transcript(null, null, 0, 0, new ArrayList<>());
		script.getExons().add(new Exon(script, 0, 2));
		script.getExons().add(new Exon(script, 4, 7));
		script.getExons().add(new Exon(script, 10, 13));
		script.getExons().add(new Exon(script, 21, 27));
		script.checkPos();
		
		i = script.getExonIndex(-1);
		i = script.getExonIndex(Integer.MAX_VALUE);
		i = script.getExonIndex(5);
		i = script.getExonIndex(7);
		i = script.getExonIndex(4);
		
		i = script.getExonSeqIndex(-1, 0);
		i = script.getExonSeqIndex(6, -1);
		i =	script.getExonSeqIndex(Integer.MAX_VALUE, 40);
		i = script.getExonSeqIndex(24, -4);
		i = script.getExonSeqIndex(7, 5);
		i = script.getExonSeqIndex(6, 5);
		i = script.getExonSeqIndex(10, 4);
		
		fq.seq = null;
		fq.len = 6;
		i = fq.buildSeq(s, 10, 4, script, null);
		
		System.out.println(fq);
	}
	
	public Fastq(String id, int len) {
		super();
		this.id = id;
		this.len = len;
	}
	
	public Fastq(String id, String seq, String add_info, String quali, String annote) {
		super();
		this.id = id;
		this.seq = seq;
		this.add_info = add_info;
		this.quali = quali;
		this.annote = annote;
	}

	public Fastq(String id, String seq, String quali, String annote) {
		super();
		this.id = id;
		this.seq = seq;
		this.quali = quali;
		this.annote = annote;
	}
	
	public boolean buildMutSeq(ArrayList<Mutation> mut_list, boolean mut_flag) {
		if (seq == null || annote == null) {
			return false;
		}
		if (this.mut_flag != null) {
			return this.mut_flag == mut_flag;
		}
		if (mut_flag) {
			StringBuilder mut_seq = new StringBuilder();
			StringBuilder mut_anno = new StringBuilder();
			mut_anno.append(annote);
			mut_anno.append('\t');
			mut_seq.append(seq);
			
			String[] cols = annote.split("\t");
			String[] regions = cols[2].split(":");
			for (String region : regions) {
				String[] bounder = region.split("-");
				int start = Integer.parseInt(bounder[0]) - 1;
				int end = Integer.parseInt(bounder[1]);
				start = Method.binarySearch(start, mut_list, 0, mut_list.size());
				end = Method.binarySearch(end, mut_list, 0, mut_list.size());
				for (;start < end; ++start) {
					
				}
			}
		}
		else {
			
		}
		this.mut_flag = mut_flag;
		return true;
	}
	
	public boolean isDefinedMut() {
		return mut_flag != null;
	}
	
	public Boolean getMutFlag() {
		return mut_flag;
	}
	
	public int buildSeq(String bases, int base_index, int seq_index, Transcript script, ArrayList<Mutation> mut_list) {
		if (seq != null) {
			return 0;
		}
		if (script == null || script.getExons() == null || script.getExons().size() == 0) {
			return -1;
		}
		int exon_start = script.getExonSeqIndex(base_index, seq_index);
		if (exon_start + len < InParam.getParams().getMinReadLen()) {
			return -1;
		}
		if (exon_start < 0) {
			len += exon_start;
			exon_start = 0;
		}
		if (exon_start + InParam.getParams().getMinReadLen() > script.getExonLength()) {
			return -1;
		}
		if (exon_start + len > script.getExonLength()) {
			len = script.getExonLength() - exon_start;
		}
		
		
		seq = script.getExonBases(bases).substring(exon_start, exon_start + len);
		StringBuilder sa = new StringBuilder();
		sa.append(id);
		sa.append('\t');
		sa.append(script.getScript_id());
		sa.append('\t');
		StringBuilder sq = null;
		StringBuilder smut = null;
		boolean mut_flag = false;
		for (Exon exon : script.getExons()) {
			if (exon.getLength() <= exon_start) {
				exon_start -= exon.getLength();
			}
			else {
				if (mut_list != null && mut_list.size() > 0) {
					sa.append(exon_start + exon.getStart() + 1);
					sa.append('-');
					sa.append(Math.min(exon.getEnd(), exon_start + exon.getStart() + len));
					sa.append(':');
					int left_index = Method.binarySearch(exon_start + exon.getStart(), mut_list, 0, mut_list.size());
					int right_index = Method.binarySearch(Math.min(exon.getEnd(), exon_start + exon.getStart() + len), mut_list, 0, mut_list.size());
					for (;left_index < right_index; ++left_index) {
						if (!mut_flag) {
							sq = new StringBuilder();
							sq.append(seq);
							smut = new StringBuilder();
							smut.append('\t');
							mut_flag = true;
						}
						sq.setCharAt(mut_list.get(left_index).getStart() - base_index, mut_list.get(left_index).getMut());
						smut.append(String.valueOf(mut_list.get(left_index).getEnd()) + mut_list.get(left_index).getMut());
						smut.append(':');
					}
				}
				if (exon.getEnd() >= exon_start + len) {
					break;
				}
				exon_start = 0;
				len -= exon.getLength();
			}
		}
		sa.setLength(sa.length() - 1);
		if (mut_flag) {
			seq = sq.toString();
			smut.setLength(smut.length() - 1);
			sa.append(smut);
		}
		annote = sa.toString();
		len = seq.length();
		quali = Method.buildString(len, InParam.getParams().getQuality());
		return 0;
	}
	
	public int buildSeq(String bases, int base_index, int seq_index, ArrayList<Mutation> mut_list) {
		if (seq != null) {
			return 0;
		}
		if (seq_index > base_index) {
			if (seq_index + InParam.getParams().getMinReadLen() > base_index + len) {
				return seq_index - base_index;
			}
			else {
				len += base_index - seq_index;
				seq_index = base_index;
			}
		}
		if (bases.length() < len - seq_index + base_index) {
			if (bases.length() < base_index - seq_index + InParam.getParams().getMinReadLen()) {
				return bases.length() - len + seq_index - base_index;
			}
			else {
				len = bases.length() + seq_index - base_index;
			}
		}
		base_index = base_index - seq_index;
		if (mut_list == null || mut_list.size() == 0) {
			seq = bases.substring(base_index, base_index + len);
		}
		else {
			StringBuilder sb = new StringBuilder();
			sb.append(bases.substring(base_index, base_index + len));
			int left_index = Method.binarySearch(base_index, mut_list, 0, mut_list.size());
			int right_index = Method.binarySearch(base_index + len, mut_list, 0, mut_list.size());
			for (;left_index < right_index; ++left_index) {
				sb.setCharAt(mut_list.get(left_index).getStart() - base_index, mut_list.get(left_index).getMut());
			}
			seq = sb.toString();
		}
		quali = Method.buildString(len, InParam.getParams().getQuality());
		return 0;
	}
	
	public String getID() {
		return id;
	}
	
	public String getAnnotation() {
		return annote;
	}
	
	public int seqLen() {
		return this.seq.length();
	}
	
	private boolean checkID() {
		return id.length() > 0 && id.charAt(0) == '@';
	}
	
	private boolean checkSeq() {
		for (int i = 0; i < seq.length(); i++) {
			switch(seq.charAt(i)) {
			case 'A':
			case 'T':
			case 'U':
			case 'G':
			case 'C':
			case 'N':
				break;
			default:
				return false;
			}
		}
		return true;
	}
	
	private boolean checkAdd_Info() {
		return add_info.length() > 0 && add_info.charAt(0) == '+';
	}
	
	private boolean checkQuality() {
		if (quali.length() != seq.length()) {
			return false;
		}
		for (int i = 0; i < quali.length(); i++) {
			if (quali.charAt(i) < 33 || quali.charAt(i) > 126) {
				return false;
			}
		}
		return true;
	}
	
	public boolean check() {
		return checkID() && checkSeq() && checkAdd_Info() && checkQuality();
	}

	@Override
	public String toString() {
		StringBuffer tmp = new StringBuffer();
		tmp.append(id);
		tmp.append('\n');
		tmp.append(seq);
		tmp.append('\n');
		tmp.append(add_info);
		tmp.append('\n');
		tmp.append(quali);
		return tmp.toString();
	}
	
	public String toStringWithCheck() {
		if (!check()) {
			System.err.println("Warning: Wrong read:\n" + toString());
			return null;
		}
		return toString();
	}
}
