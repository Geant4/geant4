//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4RTJpegMaker.cc,v 1.7 2001-11-22 17:29:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//

#include "G4RTJpeg.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTJpegCoder.hh"
#include "g4std/fstream"

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
#ifdef G4USE_STD_NAMESPACE
        ofs.open(fileName,G4std::ios::out|G4std::ios::trunc|G4std::ios::binary);
#else
	ofs.open(fileName,ios::out|ios::trunc);
#endif
        ofs.write(jpegAddress,jpegSize);
        ofs.close();

}









