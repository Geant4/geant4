#include "G4RTJpeg.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTJpegCoder.hh"
#include <fstream.h>

G4RTJpegMaker::G4RTJpegMaker()
{;}

G4RTJpegMaker::~G4RTJpegMaker()
{;}

void G4RTJpegMaker::CreateFigureFile(G4String fileName,
              int nColumn, int nRow,
              u_char* colorR,
              u_char* colorG,
	      u_char* colorB)
{
        G4JpegCoder aFigure(colorR,colorG,colorB);
        G4JpegProperty aProperty;
        aProperty.nColumn = nColumn;
        aProperty.nRow = nRow;
        aProperty.Units = 0;
        aProperty.HDensity = 1;
        aProperty.VDensity = 1;
        aProperty.ExtensionCode = 0;
				aProperty.Comment = "Geant4 Ray Tracer Version 1.0 by M.Asai K.Minamimoto C.Kishinaga";

        aFigure.SetJpegProperty(aProperty);
        aFigure.DoCoding();

        char* jpegAddress;
        int jpegSize;

        aFigure.GetJpegData(&jpegAddress,jpegSize);

        ofstream ofs;
        ofs.open(fileName,ios::out);
        ofs.write(jpegAddress,jpegSize);
        ofs.close();

}


