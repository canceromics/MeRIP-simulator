package note;

import genome.IntRegion;

public class CircleClip extends IntRegion {

	private String description = null;
	private boolean gt_ag = false;
	private boolean visited = false;
	
	public CircleClip(int start, int end, String description, boolean gt_ag) {
		super(start, end);
		this.description = description;
		this.gt_ag = gt_ag;
	}

	public boolean isGTAG() {
		return gt_ag;
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
		sb.append('\t');
		sb.append(description);
		sb.append('\t');
		sb.append(gt_ag);
		sb.append('\t');
		sb.append(visited);
		return sb.toString();
	}
}
