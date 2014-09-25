package cluster;

public class Vector {
	double x;
	double y;
	
	public Vector(double x1, double y1) {
		x=x1;
		y=y1;
	}

	public Vector mult(double d) {
		x *= d;
		y *= d;
		return this;
	}
	
	public static double distance(Vector v1, Vector v2) {
		double dx = v1.x - v2.x;
		double dy = v1.y - v2.y;
		return Math.sqrt(dx * dx + dy * dy);
	}
	
	public static double sqDistance(Vector v1, Vector v2) {
		double dx = v1.x - v2.x;
		double dy = v1.y - v2.y;
		return dx * dx + dy * dy;
	}
	
	public double magnitude() {
		return Math.sqrt(x * x + y * y );
	}
	
	public double sqMagnitude() {
		return x * x + y * y ;
	}
	
	public static Vector subtract(Vector v1, Vector v2) {
		Vector v = new Vector(v1.x - v2.x, v1.y - v2.y);
		return v;
	}
	
	public static Vector add(Vector v1, Vector v2) {
		Vector v = new Vector(v1.x + v2.x, v1.y + v2.y);
		return v;
	}

	public void subtract(Vector d) {
		x -=d.x;
		y -=d.y;
	}

	public void add(Vector d) {
		x +=d.x;
		y +=d.y;
	}
}
