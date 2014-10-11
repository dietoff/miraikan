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
	private float roty = 0f;
	private int pointer=0;
	private File[] list;
	private boolean rec=false;
	private float rotx = 0f;
	private int oldX=0;
	private int oldY=0;

	public void setup(){
		size(800, 800, P3D);
		perspective(PI/12, (float)(width)/(float)(height), 
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
		pointer = pointer%list2.length+1;
		if (!list2[pointer].getName().endsWith(".png")) {
			pointer++;
			getFile(list2);
		} else {
			String absolutePath = list2[pointer].getAbsolutePath();
			img = loadImage(absolutePath);
		}
		pointer++;
	}

	public void draw(){
		background(100);
		translate(width/2,height/2,-1000);
		rotate(rotx,1 ,0 , 0);
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
		case 'a':
			rotx +=0.1f;
			break;
		case 'z':
			rotx -=0.1f;
			break;
		case '1':
			getFile(list);
			break;
		case 'x':
			rec=true;
			break;
		}
	}
	public void mousePressed() {
		oldX = mouseX;
		oldY = mouseY;
	}
	
	public void mouseDragged() {
		float deltaX = oldX - mouseX;
		float deltaY = oldY - mouseY;
		rotx = Math.max(-PI/2f, Math.min(PI/2f, rotx + deltaY/200f));
		roty -= deltaX/180f;
		oldX = mouseX;
		oldY = mouseY;
	}

	public static void main(String args[]) {
		PApplet.main(new String[] { "viewer.Viewer" });
	}


}
