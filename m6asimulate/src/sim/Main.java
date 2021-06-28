package sim;

public class Main {

	public static void main(String[] args) {
		InParam in_args = InParam.getParams();
		in_args.putParams(args);
		if (in_args.check()) {
			Method.printNow("Start program");
			Method.run(in_args);
		}
//		Method.test();
		Method.printNow("Finished!");
	}

}
