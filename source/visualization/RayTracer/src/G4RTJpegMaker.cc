// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTJpegMaker.cc,v 1.4 2000-03-09 15:36:37 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//

#include "G4RTJpeg.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTJpegCoder.hh"
#include <g4std/fstream>

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

        G4std::ofstream ofs;
        ofs.open(fileName,G4std::ios::out);
        ofs.write(jpegAddress,jpegSize);
        ofs.close();

}
