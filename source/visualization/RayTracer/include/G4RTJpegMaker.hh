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

