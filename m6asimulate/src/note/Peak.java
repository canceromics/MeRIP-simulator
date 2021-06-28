package note;

public class Peak extends Site {
	
	private String description = null;

	public Peak(int pos, String description) {
		super(pos);
		this.description = description;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		if (description != null) {
			sb.append('\t');
			sb.append(description);
		}
		return sb.toString();
	}

}
