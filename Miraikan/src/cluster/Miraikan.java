package cluster;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import gov.nasa.worldwind.geom.Angle;
import gov.nasa.worldwind.geom.LatLon;
import processing.core.PApplet;
import processing.core.PFont;
import processing.core.PGraphics;
import processing.core.PImage;
import processing.data.Table;
import processing.data.TableRow;
import processing.pdf.*;

public class Miraikan extends PApplet {
//	private static final double cutoff_lat = 85.0511; // dont consider pixels north of this latitude
	private static double cutoff_lat = 88.97; // dont consider pixels north of this latitude
	private static final String exportfile = "/Users/offenhuber_d/Downloads/capture/frame-####.pdf";
	private static String inputimg = "data/250_night.png";
	boolean rec = false;
	int wScreen = 1000;
	int hScreen = (int)(wScreen*1.5);

	double wCoord,hCoord;
	
	float diameter; // pixel size determined once image is loaded
	
	// physics properties undistorted
	double attractFactor = 0.1; // attract strength between nodes
	double repelFactor = 0.45 * attractFactor; // repel strength between nodes
	double maxDistance = wScreen*wScreen/36; // don't repel beyond this squared distance
	double maxDistanceM = PI*PI/16; // don't repel beyond this squared distance
	
	private Cluster[] clusters;
	private Node[] nodes;
	//private PImage image;
	
	float minX,maxX;
	float minY,maxY;
	
//	private boolean rad = true; // use radians (mercator mode)
	private boolean cluster = false;
	private int cols;
	private String tblfile = "data/1deg_grid.csv";
	private char mode ='m';
	private int counter = 0;
	private int start=0;


	Comparator<Node> cmpY = new Comparator<Node>() {
	      public int compare(Node o1, Node o2) {
	        int comp = Double.valueOf(Math.abs(o1.pos.y)).compareTo(Double.valueOf(Math.abs(o2.pos.y)));
	        if (comp == 0)  comp = o1.getParent().color-o2.getParent().color;
	        if (comp == 0)	comp = -Double.valueOf(o1.getSize()).compareTo(Double.valueOf(o2.getSize()));
	        return comp;
	      }
	};


	Comparator<Node> cmpX = new Comparator<Node>() {
	      public int compare(Node o1, Node o2) {
	        int comp = Double.valueOf(o1.pos.x).compareTo(Double.valueOf(o2.pos.x)); // dirty hack TODO
	        if (comp == 0)  comp = o1.getParent().color-o2.getParent().color;
	        if (comp == 0)	comp = -Double.valueOf(o1.getSize()).compareTo(Double.valueOf(o2.getSize()));
	        return comp;
	      }
	};

	Comparator<Node> cmpColor = new Comparator<Node>() {
		public int compare(Node o1, Node o2) {
			int c1 = o1.getParent().color;
			int c2 = o2.getParent().color;
			int b1 = (int) (p.red(c1)+p.green(c1)+p.blue(c1));
			int b2 = (int) (p.red(c2)+p.green(c2)+p.blue(c2));
			int comp = b2-b1;
			if (comp == 0)	comp = c2-c1;
			if (comp == 0)  comp = Double.valueOf(Math.abs(o1.origin.x)).compareTo(Double.valueOf(Math.abs(o2.origin.x)));
			if (comp == 0)	comp = -Double.valueOf(o1.getSize()).compareTo(Double.valueOf(o2.getSize()));
			return comp;
		}
	};
	Comparator<Node> cmpColorSize = new Comparator<Node>() {
		public int compare(Node o1, Node o2) {
			int c1 = o1.getParent().color;
			int c2 = o2.getParent().color;
			int b1 = (int) (p.red(c1)+p.green(c1)+p.blue(c1));
			int b2 = (int) (p.red(c2)+p.green(c2)+p.blue(c2));
			int comp = b2-b1;
			if (comp == 0)	comp = c2-c1;
			if (comp == 0)	comp = -Double.valueOf(o1.getSize()).compareTo(Double.valueOf(o2.getSize()));
			if (comp == 0)  comp = Double.valueOf(Math.abs(o1.origin.x)).compareTo(Double.valueOf(Math.abs(o2.origin.x)));
			return comp;
		}
	};
	private Miraikan p;
	private boolean stack = true;

	public void setup() {
		p = this;
		size(wScreen, hScreen,P3D);
		 loadImgRad(inputimg, false);
		
		for (int i=0;i<clusters.length; i++) println("color:" + clusters[i].color + ", " +clusters[i].nodes.size());
		println("tot. n. nodes:"+ nodes.length);
	}

	private void loadTblRad(String file, String name, int color, double cutoff, int mode, float factor) { // mode: 1 - linear, 2- sqrt, 3-log10
		Table table = loadTable(file,"header");
		Cluster c = new Cluster(color + "", color);
		double radFactor = PI/180.0; // for conversion
		
		double[] ctmp = table.getDoubleColumn("x");
		HashMap<Double, Double> tmp = new HashMap<Double,Double>();
		for (Double dd:ctmp) tmp.put(dd, dd);
		cols = tmp.keySet().size();
		
		minX = minY = Float.MAX_VALUE;
		maxX = maxY = Float.MIN_VALUE;
		
		// determine max / min
		double[] col = table.getDoubleColumn(name);
		ArrayList<Double> t = new ArrayList<Double>();
		for (double d:col) t.add(d);
		Collections.sort(t);
		int percentile = (int) (t.size()*.999); // to get rid of extreme outliers on the max end
		 Object[] array = t.toArray();
		Double dd = t.get(percentile);
		
		float min = 0;
		float max = (float) (dd*1.0f);
		
		for (TableRow r:table.rows()) {
			double x = r.getDouble("x")*radFactor+PI;
			double y = r.getDouble("y")*radFactor;
			
			if (x < minX) minX = (float) x;
			if (y < minY) minY = (float) y;
			if (x > maxX) maxX = (float) x;
			if (y > maxY) maxY = (float) y;
			
			float val = r.getFloat(name);
			if (val>cutoff){
			Node n = c.addNode(x, y);
			n.goal = new Vector(x, y);
			n.origin = new Vector(x, y);
			
			
			float valN;
			switch (mode){
			case 2:
			valN = map(val, min, max, 0f, 1f);
			n.setSize(Math.sqrt(valN)*factor);
			case 3:
			valN = map(val, min, max, 1f, 10f);
			n.setSize(Math.log10(valN)*factor);
			break;
			default:
			valN = map(val, min, max, 0f, 1f);
			n.setSize(valN*factor);
			break;
			}
			}
		}
		diameter = maxX / (float)cols;
		HashMap<Integer,Cluster> map = new HashMap<Integer, Cluster>();
		for (int i=0; i< clusters.length; i++)  map.put(i, clusters[i]);
		map.put(clusters.length, c);
		registerLists(map);
		//Cluster[] cl = {c};
		//clusters = cl;
		//nodes = c.nodes.toArray(new Node[0]);
	}
	
	private void loadImgRad(String inputimg, boolean brightness) {
		PImage image = loadImage(inputimg);
		minX=0;
		maxX=2*PI;
		minY=PI/2;
		maxY=-PI/2;
		diameter = maxX / (float)image.width;
		cols = image.width;
	
		HashMap<Integer, Cluster> map = new HashMap<Integer, Cluster>();
		int ignore = this.color(0, 0, 0);
		
		double e = 90-cutoff_lat;
		double d = image.height*e/180f;
		
		double i = d;
		while (i < image.height-d) { // offset at the poles
			
			double degY = 90-i*180/(double)image.height;
			double radY = Math.toRadians(degY); // calc y in radians
			double dotWidth=0;
			double j = 0;
			while (j< image.width - dotWidth) {
				 dotWidth = 1/Math.cos(radY); // dot width through mercator local distortion factor
				 j += dotWidth/2; // make sure edge of the dot is inside screen (left edge)
				
				double rest = image.width%dotWidth; // how much space left after fitting an int nr of dots?
				double nrDots = (int) (image.width/dotWidth);
				double addSpacing = 0;
				if (nrDots>0) addSpacing = rest/nrDots; // we will add that later to have envely spaced dots betw. the edges
				
				double degX = j*360/(double)image.width; // we go from 0 to 360 instead of -180 to 180 to make calculations easier
				double radX = Math.toRadians(degX); // calc x in radians
	
				int color=0;
				float size=1;
				if (brightness) {
					size = brightness(image.get((int)j, (int)i))/255.0f;
					color = this.color(255);
					
				} else {
					int r = (int) red(image.get((int)j, (int)i));
					int g = (int) green(image.get((int)j, (int)i));
					int b = (int) blue(image.get((int)j, (int)i));
					color = this.color(r, g, b);
				}
	
				if (color != ignore) {
					Cluster c;
					if (map.containsKey(color)) {
						c = map.get(color);
					} else {
						c = new Cluster(color + "", color);
						map.put(color, c);
					}
					double x = radX;
					double y = radY;
					
					Node n = c.addNode(x, y);
					//n.setSize(Math.random());
					n.setSize(size);
					
				}
				j +=(dotWidth/2+addSpacing);
			}
			i++;
		}
		
		registerLists(map);
	}

	private void setTransparent(){
		for (Node n:nodes) {
			float r = this.red(n.col);
			float g = this.green(n.col);
			float b = this.blue(n.col);
			n.col = this.color(r,g,b,0);
		}
	}
	
	private void registerLists(HashMap<Integer, Cluster> map) {
		clusters = map.values().toArray(new Cluster[0]);
		ArrayList<Node> ntmp = new ArrayList<Node>();
		for (int i = 0; i < clusters.length; i++) {
			ntmp.addAll(clusters[i].nodes);
		}
		Collections.sort(ntmp, cmpColor);
		nodes = ntmp.toArray(new Node[0]);
		setTransparent();
	}
	
	private void sortLists(Comparator cmp) {
		ArrayList<Node> ntmp = new ArrayList<Node>();
		for (int i = 0; i < clusters.length; i++) {
			ntmp.addAll(clusters[i].nodes);
		}
		Collections.sort(ntmp, cmp);
		nodes = ntmp.toArray(new Node[0]);
	}

	
	public void draw() {
		blendMode(ADD);
		if (rec)
			beginRecord(PDF, exportfile);

		switch (mode) {
		case 'v': // sort horizontal
			histogramRad(nodes);
			break;
		case 'h':
			horizontalRad(nodes); 
			break;
		case 'm':
			map(nodes); //map
			break;
		case 'f':
			force(nodes); // free forces
			break;
		}

		if (rec) {
			 endRecord();
			//saveFrame(exportfile);
			//rec = false;
		}
	}

	/**
	 * set up node placements in the bar-chart - Mercator version
	 * TODO lots of them
	 * @param nodes
	 */
	private void horizontalRad(Node[] nodes) {
		if (nodes.length ==0) return;
		List<Node> l = new ArrayList<Node>();
		for (Node nn:nodes) l.add(nn);

		//Collections.sort(l, cmpColor);
		
			double y = Math.toRadians(25); //starting near the top
			Iterator<Node> n = l.iterator();
			Node node = l.get(0);
				float x = 0;
				while (n.hasNext()) {
					node = n.next();
					
					double mercDia = mercatorDiameter(y, diameter);
					node.goal = new Vector(x * mercDia+mercDia/2, y);
					double d = maxX/mercDia;
					double g = d/(int)d;
					x+=g;
					if (x > d-0.9 ) {
						x = 0;
						y -= diameter;
					}
				}

		moveToGoal(nodes);
		renderNodesRad(nodes);
	}
	
	
/**
	 * set up node placements horizontally
	 * 
	 * @param nodes
	 */
	private void histogramRad(Node[] nodes) {
		if (nodes.length ==0) return;
		int slots = (int) (cols*.75);
		HashMap<Integer,ArrayList<Node>> hist = new HashMap<Integer,ArrayList<Node>>();
		
		// init histogram
		for (int i = 0; i<slots;i++) {
			hist.put(i, new ArrayList<Node>()); 
		}
		// file in nodes
		for (Node nn:nodes) {
			int slot = (int) (nn.origin.x*slots/2/PI);
			hist.get(slot).add(nn);
		}
		
		// sort nodes & display
		for (int j = 0; j<slots;j++) {
			ArrayList<Node> ls = hist.get(j);
			Collections.sort(ls,cmpColorSize);
			
			double column= 0,vertN=0, vertS=0;
			Iterator<Node> i = ls.iterator();
			int lastCol = 0;
			while (i.hasNext()) {
				Node next = i.next();
				
				if ( !stack  && next.col != lastCol ) {
					vertN = vertS = 0;
					lastCol = next.col;
				}
				
				next.goal.x = j*PI/slots*2;
					if (next.origin.y>0) {
						vertN += diameter*next.getSize()/2;
						next.goal.y = vertN;
						vertN += diameter*next.getSize()/2;
					} else {
						vertS -= diameter*next.getSize()/2;
						next.goal.y = vertS;
						vertS -= diameter*next.getSize()/2;
					}
			}
			
		}
		
		moveToGoal(nodes);
		renderNodesRad(nodes);
	}
	

	/**
	 * set up node placements in the map
	 * works for both rad
	 * @param nodes
	 */
	private void map(Node[]  nodes) {
		{
			for (Node n : nodes)
				n.goal = new Vector(n.origin.x, n.origin.y);

			moveToGoal(nodes);
			renderNodesRad(nodes);
		}
	}

	/**
	 * this calculates the position of the nodes in the force model
	 * 
	 * @param nodes
	 */
	private void force(Node[]  nodes) {

		// repel nodes
		if(cluster) repPerCluster(clusters); else repelOpt(nodes);
		
		// attract
		for (int i = 0; i < clusters.length; i++) {
			Cluster c = clusters[i];
			c.updateCenter(); // calculate new cluster center
			c.attractNodes(attractFactor);// make nodes move tow the cluster center
		}

		boundaryRad(nodes);  // make borders of the screen repelling
		
		renderNodesRad(nodes);  // render our nodes and text
	}

	private void moveToGoal(Node[] nodes) {
		if (nodes.length==0) return;
		// make nodes move towards their goals, over multiple frames
		int range = 300;  // larger is slower
		int df = nodes.length/range; // divides the number of nodes in [range] units, defining the chunk of nodes the loop goes through
		int f = ((frameCount-start)*df)%nodes.length; // determine start frame based on current frame 
		int j = 150; // loop through j chunks at the time - larger is faster, less frames
		for (int i = Math.max(0,f-df*j); i < f+df*j; i++) { // actually start earlier, to allow the first nodes enough time to settle 
			Node n = nodes[i%nodes.length];
			Vector d = Vector.subtract(n.goal, n.pos);
			d.mult(.1f);
			n.pos.add(d);
		}
		fadeIn();
	}
	

	private void fadeIn() {
		int range = 50; // larger is slower, more granular
		int df = nodes.length/range; 
		int f = (frameCount-start)*df%nodes.length;
		int j = 5; // larger is faster, less frames
		for (int i = Math.max(0,f-df*j); i < f+df*j; i++) {
			Node n = nodes[i%nodes.length];
			float r = this.red(n.col);
			float g = this.green(n.col);
			float b = this.blue(n.col);
			float a = Math.min(255, this.alpha(n.col)+50); // increase alpha by this
			n.col = this.color(r,g,b,a);
		}
	}
	
	private void repPerCluster(Cluster[] clusters2) {
		double minDistance;
//		if (rad){
			minDistance = 0.00001;
			maxDistance = wScreen/3.0;
//		}
//		else{
//			minDistance = 0.01;
//			maxDistance = PI/6.0;
//		}
		
		for (int i = 0; i < clusters.length; i++) {
			Cluster c = clusters[i];
			 c.repelNodes(repelFactor, maxDistanceM, diameter, minDistance);// make nodes repel from cluster center
		}
		
		for (int i = 0; i < clusters.length; i++) {
			for (int j = 0; j < clusters.length; j++) {
				if (i != j) {
					Cluster c1 = clusters[i];
					Cluster c2 = clusters[j];
					c1.updateCenter();
					c2.updateCenter();

					double distance = Math.max(Vector.sqDistance(c1.center, c2.center),minDistance);
					
					float dia = diameter / 40f;
				

					if (distance < maxDistance) {
						Vector d = Vector.subtract(c2.center, c1.center);
						
						d.mult(repelFactor * dia / distance);
						for (Node n : c2.nodes) {
							n.pos.add(d);
						}
						for (Node n : c1.nodes) {
							n.pos.subtract(d);
						}
					}
				}
			}
		}

	}

	/**
	 * method for calculating repelling forces between nodes, optimized
	 * 
	 * @param nodes
	 */
	private void repelOpt(Node[] nodes) {
		int gridx = 1, gridy = 1;

		float cellLength = diameter * 10;
		int ox = (int) (maxX / cellLength);
		int oy = (int) (PI / cellLength);

		double d = Math.random();
		gridx = ox - (int) (ox / 2 * d);
		gridy = oy - (int) (oy / 2 * d);

		// long startTime = System.nanoTime();

		Grid grid;
		grid = makeGridRad(nodes, gridx, gridy);
		iterGrid(grid);

		// long endTime = System.nanoTime();
		// System.out.println("2 execution time: " + (endTime-startTime) +
		// "ms");
	}


	private void iterGrid(Grid grid) {
		int cutoffdist = 9; // max distance for rough calc (squared)
		int precdist = 4; // max distance for precision calc (squared)

		for (Cell cell:grid.getCells()) {
			for (Node n : cell.nodes) {
				for (Cell c:grid.getCells()){
					double dist = distanceSq(new Vector(cell.x,cell.y), new Vector(c.x,c.y)); //topological distance between cells
					if (dist <= cutoffdist) {
						if (dist > precdist) repelCellRad(n, c.nodes);
						else {
							for (Node n2 : c.nodes) repelNodesRad(n, n2);
						}
					}
				}
			}
		}
	}

	private Grid makeGridRad(Node[]  nodes, int gridx, int gridy) {
		Grid grid = new Grid();

		for (Node n : nodes) {
			int nx = (int) (Math.max(0,Math.min(maxX, n.pos.x) * gridx / (float) maxX));
			int ny = (int) (Math.max(0,Math.min(PI, n.pos.y) * gridy / (float) PI));
			grid.addNode(nx,ny, n);
		}
//		print();
		return grid;
	}

	
	
	private void repelNodesRad(Node n1, Node n2) {
		if (n1 != n2) {
			//float dia1 = (float) (diameter / Math.max(Math.cos((n1.pos.y+n2.pos.y)/2), 0.5)); // don't increase rep. force beyond 60 deg latitude
			//float dia2 = (float) (diameterRad / Math.max(Math.cos(n2.pos.y), 0.5));

			double mindistance = 0.00001;
			double distance = Math.max(Vector.sqDistance(n1.pos, n2.pos),mindistance);
			
			if (distance < maxDistanceM) {
				Vector d = Vector.subtract(n2.pos, n1.pos);
				
				float d1 = (float) (diameter*n1.getSize()*0.5);
				float d2 = (float) (diameter*n2.getSize()*0.5);
				
				d.mult(repelFactor * d1 * d2 / distance); 
				
				//d.mult(repelFactor * dia1 * dia1 / 4 / distance);
				n2.pos.add(d);
				n1.pos.subtract(d);
			}
		}
	}

	private void repelCellRad(Node n1, List<Node> n2) {
		if (n2 == null)
			return;

		double meanx = 0;
		double meany = 0;
		int count = 0;
		for (Node n : n2) {
			meanx += n.pos.x;
			meany += n.pos.y;
			count++;
		}
		meanx = meanx / count;
		meany = meany / count;
		
		float dia = (float) (diameter / Math.max(Math.cos(meany), 0.5));
		
		Vector p = new Vector(meanx, meany);

		double minDistance = 0.0001;
		double distance = Math.max(Vector.sqDistance(n1.pos, p), minDistance);
		if (distance < maxDistanceM) {
			Vector d = Vector.subtract(p, n1.pos);
			d.mult(repelFactor * dia * dia / 4 / distance);
			for (Node n : n2)
				n.pos.add(d);
			n1.pos.subtract(d);
		}
	}
	/**
	 * make edges of the screen repel nodes (for force-based mode)
	 * here we are in lat/long space, PI/2 > y > -PI/2
	 * @param nodes
	 */
	private void boundaryRad(Node[]  nodes) {
		float r = PI/2;
		double rep = repelFactor/2f;
		
		for (Node n : nodes) {
			if (n.pos.x > maxX - r)
				n.pos.x -= rep * 0.5 * (n.pos.x - (maxX - r));
			if (n.pos.x < r)
				n.pos.x += rep * 0.5 *  (r - n.pos.x);
			
			n.pos.y = Math.min(n.pos.y, PI/2.01);
			n.pos.y = Math.max(n.pos.y, -PI/2.01);
			
			float y1 = minY - r;
			if (n.pos.y > y1)
				n.pos.y -= rep * (n.pos.y - y1);
			float y2 = maxY + r;
			if (n.pos.y < y2) //+r because maxY is negative in radians
				n.pos.y += rep * (-n.pos.y + y2);
				
		}
	}

	/**
	 * render our nodes to the buffer
	 * 
	 * @param nodes
	 */
	private void renderNodes(Node[]  nodes) {
		//background(0, 20, 45);
		background(0);
		smooth();
		noStroke();
		
		for (Node n : nodes) {
			int col = n.col;
			int a = 0xddffffff;
			float size = (float) ((float) diameter * n.getSize());
			fill(col & a);
			ellipse((float) n.pos.x, (float) n.pos.y, size, size);
		}
	}
	private void renderNodesRad(Node[]  nodes) {
		//transform to mercator first!
		//background(0, 20, 45);
		background(0);
		smooth();
		noStroke();
		for (Node n : nodes) {
			//int col = n.col;
			//int a = 0xddffffff;
			//fill(col & a);
			fill(n.col);
			float newx = (float) n.pos.x; // starts at 0, ends at 2pi
			newx = (float) (newx*wScreen/(Math.PI*2f)); // convert back to screen units
			float size = (float) n.getSize();
			float newy = (float) mercLat2Y(n.pos.y);
			float newDia = (float) mercatorDiameter(n.pos.y, wScreen * size / (float)cols);
			
			ellipse(newx,newy,newDia,newDia);
		}
	}

	/*
	 * convenience functions for controlling the simulation
	 */
	public void keyReleased() {
		switch (key) {
		case '1':
			mode = 'v';
			start = frameCount;
			sortLists(cmpColor);
			break;
		case '2':
			mode = 'h';
			start = frameCount;
			sortLists(cmpColorSize);
			break;
		case '3':
			mode = 'f';
			break;
		case '4':
			mode = 'm';
			start = frameCount;
			break;
		case 'q':
			start = frameCount;
			loadTblRad(tblfile, "pop_sum", 0xffffffff, 250000,2,1.3f);
//			clusters = new Cluster[]{clusters[0]};
//
//			ArrayList<Node> ntmp = new ArrayList<Node>();
//			for (int i = 0; i < clusters.length; i++) {
//				ntmp.addAll(clusters[i].nodes);
//			}
//			nodes = ntmp.toArray(new Node[0]);
			break;
		case 'w':
			start = frameCount;
			loadTblRad(tblfile, "water_sum", 0xff33bbff, 1000,2,1);
			break;
		case 'e':
			start = frameCount;
			loadTblRad(tblfile, "mineral_sum", 0xff666666, 1,2,5f);
			break;
		case 'r':
			start = frameCount;
			loadTblRad(tblfile, "acc_mean", 0xffbbff00, 1,2,1.4f);
			break;
		case 't':
			start = frameCount;
			loadTblRad(tblfile, "inv_acc_mean", 0xffbbff00, 0,2,1.4f);
			break;
		case 'y':
			counter ++;
			counter = 1+counter%12;
			clusters=new Cluster[0];
			nodes = new Node[0];
			start = frameCount;
			loadImgRad("data/rain/"+counter+".png", true); 
			break;
		case 'z':
			break;
		case '0':
			start = frameCount;
			inputimg = "data/veg90.png";
			loadImgRad(inputimg, false); 
			break;
		case '6':
			start = frameCount;
			inputimg = "data/63prc_500.png";
			loadImgRad(inputimg, false); 
			break;
		case '8':
			start = frameCount;
			inputimg = "data/250_night.png";
			loadImgRad(inputimg, false); 
			break;
		case '9':
			start = frameCount;
			inputimg = "data/500_night.png";
			stack = false;
			loadImgRad(inputimg, false); 
			break;
		case '5':
			start = frameCount;
			inputimg = "data/veg_cities.png";
			loadImgRad(inputimg, false); 
			break;
		case '7':
			start = frameCount;
			inputimg = "data/access2.png";
			stack = true;
			loadImgRad(inputimg, false);
			break;
		case '-':
			clusters=new Cluster[0];
			nodes = new Node[0];
			break;
		case 'c':
			cluster = !cluster;
			break;
		case 'x':
			rec = !rec;
			break;
		case 's':
			//repelFactor = 0.45f* attractFactor;
			for (Node n:nodes) n.setSize(Math.random());
			break;
		}
	}

	private void randomize() {
		for(Cluster c:clusters) 
			for (Node n:c.nodes) {
				n.pos.x += 10*(Math.random()-0.5);
				n.pos.y += 10*(Math.random()-0.5);
			}
	}

	private double distanceSq(Vector a, Vector b) {
		double dx = a.x - b.x;
		double dy = a.y - b.y;
		return dx*dx+dy*dy;
	}
	
	private double mercLat2Y(double rad) { // in radians
		double mercY = Math.log(Math.tan(PI/4+rad/2));
		double mercMax = Math.log(Math.tan(PI/4+Math.toRadians(cutoff_lat)/2)); //start at 85.0511 deg for a square map
		double mercMin = Math.log(Math.tan(PI/4+Math.toRadians(-cutoff_lat)/2));
		double y     = wScreen*(mercMax-mercY)/(2*PI);
		return y;
	}

	private double mercatorDiameter(double y, float diameter) {
		return Math.abs((1/Math.cos(y)*diameter));
	}
	
	private LatLon greatCircle(double prc, double x1, double y1, double x2, double y2) {
		LatLon pos1 = new LatLon(Angle.fromRadians(x1), Angle.fromRadians(y1));
		LatLon pos2 = new LatLon(Angle.fromRadians(x2), Angle.fromRadians(y2));
		LatLon interPoint = LatLon.interpolateGreatCircle(prc, pos1, pos2);
		return interPoint;
	}
	
	public static void main(String args[]) {
	    PApplet.main(new String[] { "cluster.Miraikan" });
	    System.out.println("Working Directory = " +
	              System.getProperty("user.dir"));
	  }

}
	
