package note;

import genome.IntRegion;

public class Site extends IntRegion {
	
	private boolean visited = false;
	
	public Site(int pos) {
		super(pos - 1, pos);
	}
	
	public void visited() {
		visited = true;
	}
	
	public boolean isVisited() {
		return visited;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getStart() + 1);
		sb.append('\t');
		sb.append(getEnd());
		return sb.toString();
	}
}
