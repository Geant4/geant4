// Created by
// - Ph. Barberet, J. Bordes
//   Bordeaux U., France
//   E-mail: barberet@lp2ib.in2p3.fr
// - L. Morelli
//   Politecnico di Milano, Italy

// Show progress
showProgress(0);

// The phantom file will be saved in the directory chosen by the user
dir = getDirectory("Choose the output directory");

// Get voxel size and image dimensions
getVoxelSize(voxelWidth, voxelHeight, depth, unit);
getDimensions(imgWidth, imgHeight, channels, slices, frames);

// User settings dialog
title = "Phantom settings";
threshold1 = 0;
threshold2 = 0;
threshold3 = 0;
Dialog.createNonBlocking(title);
Dialog.addString("Output file name (.dat):", "phantom");
Dialog.addNumber("Threshold red [0:255]:", 30);
Dialog.addNumber("Threshold green [0:255]:", 30);
Dialog.addNumber("Threshold blue [0:255]:", 30);

numberSlices = 0;
items = newArray("RGB", "RBG", "BRG", "BGR", "GRB", "GBR");  //Definition of color priority order (1st color priority, 2nd color priority, 3rd color priority)
Dialog.addChoice("Priority", items);
Dialog.show();

// Read dialog parameters
filename = Dialog.getString();
threshold1 = Dialog.getNumber();
threshold2 = Dialog.getNumber();
threshold3 = Dialog.getNumber();
priority = Dialog.getChoice();

// Print output directory
print(dir);

// Generate file path
path2file = dir + filename + ".dat";

for (num = 0; File.exists(path2file); num++) {
    newfilename = filename + "_" + num;
    path2file = dir + newfilename + ".dat";
}

// Open temporary file for writing
tempF = File.open(dir + "_temp.dat");

// Display file parameters
W = getWidth(); // Image width in voxels
H = getHeight(); // Image height in voxels

print("Voxel size : ", voxelWidth, " ", voxelHeight, " ", depth, " ", unit);
print("Number of slices : ", slices);
print("Definition : ", W, "*", H);
print("Thresholds : ", threshold1, threshold2, threshold3);

// Display number of voxels
showStatus("Voxels count");


// Initialize voxel counters
numberVoxels1 = 0;
numberVoxels2 = 0;
numberVoxels3 = 0;


// Initialize a string to store lines of data
linesToWrite = "";
linesArray = newArray("");

// Loop through the image to write voxel coordinates and material in the phantom file
for(k=0; k< nSlices; k++)
{
  showProgress(k/(nSlices));
  setSlice(k+1);

   for(j=0; j<H; j++)
    {
       for(i=0; i< W; i++)
      {

        v=getPixel(i,j);
        red = (v>>16)&0xff;     //Extracting red color data - bits 23-16
        green = (v>>8)&0xff;    //Extracting green color data - bits 15-8
        blue = v&0xff;          //Extracting blue color data - bits 7-0

        //voxel coordinates (real units)
	    x=i*voxelWidth;
	    y=j*voxelWidth;
	    z=k*depth;

		material = 0;
        if (priority=="RGB")  //Red has priority over blue, which has priority over green, if 2 or 3 of these colors are greater than their threshold.
        {
          if (red>=threshold1)  {
            numberVoxels1 +=1;
            material = 1;}
          else if (green>=threshold2) {
            numberVoxels2 +=1;
            material = 2;}
          else if (blue>=threshold3) {
            numberVoxels3 +=1;
            material = 3;}


        }

        else if (priority=="RBG")
        {
          if (red>=threshold1)  {
            numberVoxels1 +=1;
            material = 1;}
          else if (blue>=threshold3) {
            numberVoxels3 +=1;
            material = 3}
          else if (green>=threshold2) {
            numberVoxels2 +=1;
            material = 2;}
        }

        else if (priority=="BRG")
        {
          if (blue>=threshold3) {
            numberVoxels3 +=1;
            material = 3;}
          else if (red>=threshold1)  {
            numberVoxels1 +=1;
            material = 1;}
          else if (green>=threshold2) {
            numberVoxels2 +=1;
            material = 2;}
        }

        else if (priority=="BGR")
        {
          if (blue>=threshold3) {
            numberVoxels3 +=1;
            material = 3;}
          else if (green>=threshold2) {
            numberVoxels2 +=1;
            material = 2;}
          else if (red>=threshold1)  {
            numberVoxels1 +=1;
            material = 1;}
        }

        else if (priority=="GBR")
        {
          if (green>=threshold2) {
            numberVoxels2 +=1;
            material = 2;}
          else if (blue>=threshold3) {
            numberVoxels3 +=1;
            material = 3;}
          else if (red>=threshold1)  {
            numberVoxels1 +=1;
            material = 1;}
        }

        else if (priority=="GRB")
        {
          if (green>=threshold2) {
          numberVoxels2 +=1;
          material = 2;}
          else if (red>=threshold1)  {
          numberVoxels1 +=1;
          material = 1;}
          else if (blue>=threshold3) {
          numberVoxels3 +=1;
          material = 3;}
        }

        // Append the line to the list of lines to write
        if (material != 0){
        	print(tempF, d2s(x,4) + "  \t" + d2s(y,4) + "  \t" + d2s(z,4) + "  \t" + material + "\n");
        }

        }
    }
}

numberVoxels=numberVoxels1+numberVoxels2+numberVoxels3;

// Close temporary file
File.close(tempF);

// Open main file for writing
F = File.open(path2file);

// Write header in main file
print(F, numberVoxels + "\t" + numberVoxels1 + "\t" + numberVoxels2 + "\t" + numberVoxels3 + "\n");
print(F, imgWidth * voxelWidth + "\t" + imgHeight * voxelWidth + "\t" + slices * depth + "\t" + unit + "\n");
print(F, voxelWidth + "\t" + voxelWidth + "\t" + depth + "\t" + unit + "\n");

// Read data from temporary file and write to main file
data = File.openAsString(dir + "_temp.dat");
print(F, data);

// Close main file
File.close(F);

// Delete temporary file
File.delete(dir + "_temp.dat");

// Show completion messages
showProgress(1)

if (num > 0) {
    showMessage("WARNING: '" + filename + ".dat' file already exists.\nNew file: '" + newfilename + ".dat'");
}
showStatus("Completed");
showMessage("Completed");
