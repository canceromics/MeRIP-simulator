package note;

import sim.InParam;
import sim.Method;

public class Mutation extends Site{

	private char alt = 'N';
	private char mut = 'N';
	private String description = null;
	private int no_mut_time = 0;
	private int used_time = 0;
	private int ip_used_time = 0;
	private int ip_no_mut_time = 0;
	private double ip_mut_per = 0.0;
	private double mut_per = 0.0;
	
	public static void main(String[] args) {
		Mutation mutation = new Mutation(0, 'T', 'C', null);
		mutation.mut_per = 0.6;
		mutation.ip_mut_per = Method.sampleMutPerOR(mutation.mut_per);
		mutation.no_mut_time = 5;
		mutation.used_time = 6;
		mutation.ip_used_time = 53;
		int[] i = mutation.ensureMutPer(false);
		i = mutation.ensureMutPer(true);
		System.out.println(i);
	}
	
	public Mutation(int pos, char alt, char mut, String description) {
		super(pos);
		this.alt = alt;
		this.mut = mut;
		this.description = description;
	}
	
	public char getMut() {
		return mut;
	}
	
	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}

	public void cloneTimes(Mutation mut) {
		this.used_time = mut.used_time;
		this.ip_used_time = mut.ip_used_time;
		this.no_mut_time = mut.no_mut_time;
		this.ip_no_mut_time = mut.ip_no_mut_time;
	}
	
	public int[] ensureMutPer(boolean ip_flag) {
		int[] out = {0, 0};
		double mut_time = ip_flag ? ip_used_time : used_time;
		double back_time = ip_flag ? ip_no_mut_time : no_mut_time;
		double per = ip_flag ? ip_mut_per : mut_per;
		double bias = mut_time / (mut_time + back_time) - per;
		if (Math.abs(bias) < 1e-6) {
			return out;
		}
		if (bias < 0) {
			out[0] = (int) (per * back_time / (1.0 - per) - mut_time);
			bias = Math.abs((mut_time + out[0]) / (mut_time + back_time + out[0]) - per);
			if (Math.abs((mut_time + out[0] + 1) / (mut_time + back_time + out[0] + out[1] + 1) - per) < bias) {
				out[0] += 1;
				bias = Math.abs((mut_time + out[0]) / (mut_time + back_time + out[0] + out[1]) - per);
			}
			if (Math.abs((mut_time + out[0]) / (mut_time + back_time + out[0] + out[1] + 1) - per) < bias) {
				out[1] += 1;
			}
			return out;
		}
		if (bias > 0) {
			out[1] = (int) ((1.0 - per) * mut_time / per - back_time);
			bias = Math.abs((mut_time + out[1]) / (mut_time + back_time + out[1]) - per);
			if (Math.abs((mut_time + out[0]) / (mut_time + back_time + out[0] + out[1] + 1) - per) < bias) {
				out[1] += 1;
				bias = Math.abs((mut_time + out[0]) / (mut_time + back_time + out[0] + out[1]) - per);
			}
			if (Math.abs((mut_time + out[0] + 1) / (mut_time + back_time + out[0] + out[1] + 1) - per) < bias) {
				out[0] += 1;
			}
			return out;
		}
		return out;
	}
	
	public int incUsedTime(boolean ip_flag, boolean mut_flag) {
		return mut_flag ? (ip_flag ? ++ip_used_time : ++used_time) : (ip_flag ? ++ip_no_mut_time : ++no_mut_time);
	}
	
	public boolean ensureUsedTime(boolean ip_flag) {
		return 0 < (ip_flag ? ip_used_time + ip_no_mut_time : used_time + no_mut_time);
	}
	
	public boolean ensureUsedTime(boolean ip_flag, boolean mut_flag) {
		return 0 < (mut_flag ? (ip_flag ? ip_used_time : used_time) : (ip_flag ? ip_no_mut_time : no_mut_time));
	}
	
	public boolean reduceUsedTime(boolean ip_flag, boolean mut_flag) {
		return 0 <= (mut_flag ? (ip_flag ? --ip_used_time : --used_time) : (ip_flag ? --ip_no_mut_time : --no_mut_time));
	}
	
	public int getLeftTime(boolean ip_flag, boolean mut_flag) {
		return mut_flag ? (ip_flag ? ip_used_time : used_time) : (ip_flag ? ip_no_mut_time : no_mut_time);
	}
	
	public void makeMutTimes(InParam args) {
		int base_num = args.getBackgoundNum();
		int ip_base_num = "Peak".equals(description) ? (int) (base_num / 4 + base_num * args.getEnrichNum()) : base_num / 4;
		used_time = (int) (base_num * Method.sampleMutPer(args.getMutationProportion()));
		no_mut_time = base_num - used_time;
		ip_used_time = (int) (ip_base_num * Method.sampleMutPerOR((double) used_time / (double) (used_time + no_mut_time)));
		ip_no_mut_time = ip_base_num - ip_used_time;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		sb.append('\t');
		sb.append(alt);
		sb.append('\t');
		sb.append(mut);
		if (description != null) {
			sb.append('\t');
			sb.append(description);
		}
		sb.append('\t');
		sb.append(isVisited());
		sb.append('\t');
		sb.append(ip_used_time);
		sb.append('\t');
		sb.append(ip_no_mut_time);
		sb.append('\t');
		sb.append(used_time);
		sb.append('\t');
		sb.append(no_mut_time);
		sb.append('\t');
		sb.append((double) (ip_used_time * no_mut_time) / (double) (ip_no_mut_time * used_time));
		return sb.toString();
	}
}
