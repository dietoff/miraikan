package cluster;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

public class Grid {
	HashMap<Integer,Cell> c;
	public Grid() {
		c= new HashMap<Integer,Cell>();
	}
	
	public void addNode(int x, int y, Node n) {
		
		int id = id(x, y);
		
		if (!c.containsKey(id)) {
			Cell cell = new Cell(x, y);
			cell.nodes.add(n);
			c.put(id, cell);
			
		} else {
			Cell cell = c.get(id);
			cell.nodes.add(n);
		}
	}

	private int id(int x, int y) {
		int id = (x << 16) + y;
		return id;
	}
	public Cell getCell(int x, int y) {
		 return c.get(id(x,y));
	}
	
	public Collection<Cell> getCells() {
		return c.values();
	}
}
