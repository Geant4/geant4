
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

