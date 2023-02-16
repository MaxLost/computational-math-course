import org.math.computational.functions.*;

import org.junit.Test;

public class FunctionTest {

	@Test
	public void test(){

		double x = 2;

		Function eq = t -> Math.pow(t, 2) + t;
		System.out.println(eq.evaluate(x));
	}

}
