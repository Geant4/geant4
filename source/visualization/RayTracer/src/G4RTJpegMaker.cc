//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4RTJpegMaker.cc 77479 2013-11-25 10:01:22Z gcosmo $
//
//
//

#include "G4RTJpeg.hh"
#include "G4RTJpegMaker.hh"
#include "G4RTJpegCoder.hh"
#include <fstream>

G4RTJpegMaker::G4RTJpegMaker()
{;}

G4RTJpegMaker::~G4RTJpegMaker()
{;}

void G4RTJpegMaker::CreateFigureFile(const G4String& fileName,
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

        std::ofstream ofs;
        ofs.open(fileName,std::ios::out|std::ios::trunc|std::ios::binary);
        ofs.write(jpegAddress,jpegSize);
        ofs.close();
}
