package cluster;

import java.util.ArrayList;

public class Cluster {
	String name = "name";
	ArrayList<Node> nodes;
	int color;
	Vector center;
	double radius;
	float factor;

	public Cluster() {
		nodes = new ArrayList<Node>();
		name = "cluster";
		color = 0xffffffff;
		center = new Vector(0,0);
		radius = 3; //not used now
		factor = 25;
	}

	public Cluster(String n1, Integer col1) {
		nodes = new ArrayList<Node>();
		name = n1;
		color = col1;
		center = new Vector(0,0);
		radius = 3; //not used now
		factor = 1;
	}

	public Cluster(String n1, Integer col1, float factors) {
		nodes = new ArrayList<Node>();
		name = n1;
		color = col1;
		center = new Vector(0,0);
		radius = 3; //not used now
		factor = factors;
	}

	public Node addNode(int x, int y) {
		Node n = new Node(color,x,y);
		nodes.add(n);
		n.setParent(this);
		return n;
	}
	public Node addNode(double x, double y) {
		Node n = new Node(color,x,y);
		nodes.add(n);
		n.setParent(this);
		return n;
	}

	public void updateCenter() {
		double xavg=0;
		double yavg=0;

		// plain old average
		for (Node n:nodes) {
			xavg +=n.pos.x;
			yavg +=n.pos.y;
		}
		xavg = xavg/nodes.size();
		yavg = yavg/nodes.size();
		center.x = xavg;
		center.y = yavg;
	}

	public void attractNodes(double attractFactor, boolean rad) {
		double t; //treshhold before attract is active
		if (rad) {
			t = 0.01;
		} else
			t = 0.5;
		
		for (Node n:nodes) {
			Vector d = Vector.subtract(center, n.pos);
			double magnitude = d.sqMagnitude();
			double sqThreshold = t*t;
			if(magnitude>sqThreshold) {
				d.mult(attractFactor);
				n.pos.add(d);
			}
		}
	}

	public void setCenter(int i, int h) {
		center.x=i;
		center.y=h;
	}
	
	

	public void repelNodes(double repelFactor, double maxDistance, double diameter, double mindistance) {

		for (Node n1 : nodes) {
			for (Node n2 : nodes) {
				if (n1 != n2) {
					double distance = Math.max(
							Vector.sqDistance(n1.pos, n2.pos), mindistance);
					if (distance < maxDistance) {
						Vector d = Vector.subtract(n2.pos, n1.pos);
						double dn1 = diameter*n1.getSize()/2;
						double dn2 = diameter*n2.getSize()/2;
						d.mult(repelFactor * dn1 * dn2 / distance);
						n2.pos.add(d);
						n1.pos.subtract(d);
					}
				}
			}
		}
	}
}
