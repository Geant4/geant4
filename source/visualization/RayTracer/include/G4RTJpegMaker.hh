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
// $Id: G4RTJpegMaker.hh,v 1.5 2001-07-11 10:09:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    virtual void CreateFigureFile(G4String fileName,
	          int nColumn,int nRow,
		  u_char* colorR, u_char* colorG,
		  u_char* colorB);
};

#endif

