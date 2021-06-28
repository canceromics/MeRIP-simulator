package genome;

public class IntRegion implements Comparable<IntRegion>{
	
	private int start;
	private int end;
	
	public IntRegion(int start, int end){
		super();
		this.start = start > end ? end : start;
		this.end = start ^ end ^ this.start;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}
	
	public int getLength() {
		return end - start;
	}
	
	public void resetStartAndEnd(int start, int end) {
		this.start = start > end ? end : start;
		this.end = start ^ end ^ this.start;
	}
	
	@Override
	public int compareTo(IntRegion anotherRegion) {
		return start < anotherRegion.start ? -1 : 
			start > anotherRegion.start ? 1 : 
			end < anotherRegion.end ? -1 : 
			end > anotherRegion.end ? 1 : 0;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getStart());
		sb.append('\t');
		sb.append(getEnd());
		return sb.toString();
	}
	
}
