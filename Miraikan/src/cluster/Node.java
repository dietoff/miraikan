package cluster;

public class Node {
	Vector pos;
	int col;
	Vector goal;
	Vector origin;
	double diameter;
	boolean target = false;
	private double size;
	private Cluster parent;
	
	public Node() {
		col = 0xffffffff;
		pos = new Vector(0,0);
		goal = new Vector(0,0);
	}
	
	public Node(int color, int x, int y) {
		col = color;
		pos = new Vector(x,y);
		goal = new Vector(x,y);
		origin = new Vector(x,y);
	}
	public Node(int color, double x, double y) {
		col = color;
		pos = new Vector(x,y);
		goal = new Vector(x,y);
		origin = new Vector(x,y);
	}

	public double getSize() {
		return size;
	}

	public void setSize(double size) {
		this.size = size;
	}

	public Cluster getParent() {
		return parent;
	}

	public void setParent(Cluster parent) {
		this.parent = parent;
	}
	
}
