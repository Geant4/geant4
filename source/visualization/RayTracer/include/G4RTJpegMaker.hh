// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTJpegMaker.hh,v 1.3 2000-03-09 15:36:34 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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

