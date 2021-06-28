package genome;

import java.util.ArrayList;
import java.util.Collections;

public class Gene extends IntRegion{
	
	private String chr_symbol = null;
	private String gene_id = null;
	private String gene_symbol = null;
	private String bases = null;
	private char strand = '.';
	private ArrayList<Transcript> scripts = null;
	
	public Gene(String chr_symbol, String gene_id, String gene_symbol, char strand, String bases, int start, int end, ArrayList<Transcript> scripts) {
		super(start, end);
		this.chr_symbol = chr_symbol;
		this.gene_symbol = gene_symbol;
		this.gene_id = gene_id;
		this.strand = strand;
		this.bases = bases;
		this.scripts = scripts;
	}

	public String getChr_symbol() {
		return chr_symbol;
	}
	
	public String getGene_id() {
		return gene_id;
	}

	public String getGene_symbol() {
		return gene_symbol;
	}
	
	public String getBases() {
		return bases;
	}

	public char getStrand() {
		return strand;
	}
	
	public void mergeGene(Gene gene) {
		this.gene_id = this.gene_id + "," + gene.gene_id;
		this.gene_symbol = this.gene_symbol + ',' + gene_symbol;
		this.bases = mergeBases(gene.getBases(), gene.getStart(), gene.getEnd());
		this.strand = this.strand == gene.strand ? this.strand : '.';
		this.scripts.addAll(gene.scripts);
		checkPos();
		resetStartAndEnd(Math.min(getStart(), gene.getStart()), Math.max(getEnd(), gene.getEnd()));
	}
	
	private String mergeBases(String base2, int start2, int end2) {
		if (base2 == null) {
			return bases;
		}
		if (bases == null) {
			return base2;
		}
		StringBuilder sb = new StringBuilder();
		if (getStart() <= start2) {
			sb.append(bases);
			if (end2 > getEnd()) {
				sb.append(base2.substring(getEnd() - start2, end2 - start2));
			}
		}
		else {
			sb.append(base2);
			if (getEnd() > end2) {
				sb.append(bases.substring(end2 - getStart(), getEnd() - getStart()));
			}
		}
		return sb.toString();
	}
	
	public void setBases(StringBuffer bases) {
		this.bases = bases.length() >= getEnd() ? bases.substring(getStart(), getEnd()) : null;
	}
	
	public boolean checkPos() {
		boolean out = true;
		if (scripts != null) {
			for (int i = 0; i < scripts.size(); i++) {
				out &= scripts.get(i).checkPos()
						&& getStart() <= scripts.get(i).getStart()
						&& getEnd() >= scripts.get(i).getEnd();
			}
			Collections.sort(scripts);
			resetStartAndEnd(Math.min(getStart(), scripts.get(0).getStart()), Math.max(getEnd(), scripts.get(scripts.size() - 1).getEnd()));
		}
		return out;
	}

	public ArrayList<Transcript> getScripts() {
		return scripts;
	}
	
}
