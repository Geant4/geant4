// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VFigureFileMaker.hh,v 1.2 2000-03-09 15:36:35 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//


#ifndef G4VFigureFileMaker_H
#define G4VFigureFileMaker_H 1

#include "globals.hh"

class G4VFigureFileMaker
{
  public:
    G4VFigureFileMaker() {;}
    virtual ~G4VFigureFileMaker() {;}

  public:
    virtual void CreateFigureFile(G4String fileName,
              int nColumn,int nRow,
              unsigned char* colorR, unsigned char* colorG, 
              unsigned char* colorB) = 0;
};

#endif

