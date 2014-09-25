package viewer;


import java.io.File;

import processing.core.PApplet;
import processing.core.PImage;
import processing.core.PShape;
import processing.core.PApplet;

public class Viewer extends PApplet {

	private PImage img;
	private PShape globe;
	private File file;// = "/Volumes/Dropbox/testfiles_dietmar/migration streams/migration_bundle_cities-01.png";
	private Float roty = 0f;
	private int pointer=0;
	private File[] list;
	private boolean rec=false;
	
	public void setup(){
		size(800, 600, P3D);
		perspective(PI/11, (float)(width)/(float)(height), 
	            1, 2000);
		noStroke();
		
		file = new File("/Users/offenhuber_d/Downloads/capture/out");
		list = file.listFiles();
		getFile(list);
		
		sphereDetail(64);
		smooth(8);
		noStroke();
		globe = createShape(SPHERE, 200); 
		globe.setTexture(img);
	}
	
	public void getFile(File[] list2) {
		pointer++;
		pointer = pointer%list2.length;
		if (!list2[pointer].getName().endsWith(".png")) {
			pointer++;
			getFile(list2);
		} else {
			String absolutePath = list2[pointer].getAbsolutePath();
			img = loadImage(absolutePath);
		}
	}

	public void draw(){
		background(100);
		translate(width/2,height/2,-1000);
		rotate(roty,0 ,1 , 0);
		//lights();
		globe.setTexture(img);
		shape(globe);
		if (rec) {
			saveFrame("/Users/offenhuber_d/Downloads/capture/screenshots/shot###.png");
			rec = false;
		}
	}
	
	public void keyReleased() {
		switch (key){
		
		case ',':
			roty +=0.1f;
			break;
		case '.':
			roty -=0.1f;
			break;
		case 'z':
			getFile(list);
			break;
		case 'x':
			rec=true;
			break;
		}
	}
	
	
	public static void main(String args[]) {
	    PApplet.main(new String[] { "viewer.Viewer" });
	  }

	
}
