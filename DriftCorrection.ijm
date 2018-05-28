//3D Drift Corrects series of ZTC stacked cells
// Author, David Podorefsky - 1/7/18
dir = getDirectory("Choose Directory");
outputfol = dir + "Output/";
File.makeDirectory(outputfol); //Make output directory
phaselist = getFileList(dir + "Stacks/"); //Vector of all the folders of phases
for (i=0; i<phaselist.length; i++) {
	phase = phaselist[i];
	File.makeDirectory(outputfol + phase); //Make output phase folder
	celllist = getFileList(dir + "Stacks/" + phase);
	for (j=0; j<celllist.length; j++) {
		cell = celllist[j];
		File.makeDirectory(outputfol + phase + cell); //Make output cell folder
		pl = lengthOf(phase)-1;
		cl = lengthOf(cell)-1;
		p = substring(phase,0,pl);
		s = substring(cell,0,cl);
		print("Cell: " + p + "," + s);
		slist = getFileList(dir + "Stacks/" + phase + cell);
		t = slist.length / 2 / 96; //Total images /2 channels /96 layers		
		cellimage = dir + "Stacks/" + phase + cell + slist[1]; //first cell image
		run("Image Sequence...", "open=cellimage use convert");
		run("Stack to Hyperstack...", "order=xyczt(default) channels=2 slices=96 frames=t display=Grayscale");
		run("Correct 3D drift", "channel=1 only=0 lowest=1 highest=96");
		saveAs("Tiff", outputfol + phase + cell + p + "," + s + ",t=" + t);
		close(); 
		close();
	}
}