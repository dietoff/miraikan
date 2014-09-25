package cluster;

import java.util.ArrayList;
import java.util.List;

public class Cell {
	List<Node> nodes;
	int x; 
	int y;
	
	public Cell( int a, int b) {
		x=a;
		y=b;
		nodes = new ArrayList<Node>();
	}
}
