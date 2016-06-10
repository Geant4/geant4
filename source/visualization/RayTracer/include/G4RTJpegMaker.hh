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
// $Id: G4RTJpegMaker.hh 77479 2013-11-25 10:01:22Z gcosmo $
//
//

// class description:
//
//  This is a concrete class of G4VFigureFileMaker.
//  This class converts 8 bits unsighned integer arrays which represent RGB of
// each pixel to JPEG code and stores it to a file. The only one public method
// of this class will be invoked by G4RayTracer.
//

#ifndef G4RTJpegMaker_H
#define G4RTJpegMaker_H 1

#include "G4VFigureFileMaker.hh"
#include "G4RTJpeg.hh"

class G4RTJpegMaker : public G4VFigureFileMaker
{
  public:
    G4RTJpegMaker();
	virtual ~G4RTJpegMaker();

  public:
    virtual void CreateFigureFile(const G4String& fileName,
	          int nColumn,int nRow,
		  u_char* colorR, u_char* colorG,
		  u_char* colorB);
};

#endif

