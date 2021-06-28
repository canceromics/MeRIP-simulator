package genome;

public class Exon extends IntRegion{

	private Transcript script = null;
	private Exon same_position = null;
	private Exon next_exon = null;
	
	public Exon(Transcript script, int start, int end) {
		super(start, end);
		this.script = script;
	}

	Exon addDownStream(Exon exon) {
		Exon head = new Exon(null, 0, 0);
		Exon tail = head;
		Exon p = this;
		Exon q = exon;
		while (p != null && q != null) {
			if (p.getStart() < q.getStart()) {
				tail.next_exon = p;
				tail = p;
				p = p.next_exon;
			}
			else {
				tail.next_exon = q;
				tail = q;
				q = q.next_exon;
			}
		}
		if (p == null) {
			tail.next_exon = q;
		}
		else {
			tail.next_exon = p;
		}
		return head.next_exon;
	}
	
	public Exon getDownStream() {
		return next_exon;
	}
	
	public void addAnotherExon(Exon exon) {
		if (exon != null) {
			exon.same_position = this.same_position;
		}
		this.same_position = exon;
	}
	
	public Exon getAnotherExon() {
		return same_position;
	}
	
	public String getChr() {
		return script==null? null : script.getChr();
	}
	
	public Gene getGene() {
		return script==null? null : script.getGene();
	}
	
	public Transcript getScript() {
		return script;
	}
	
}
