///////////////////////////////////////////////////////////////////////////////
// File: CMSRotationMatrixFactory.hh
// Date: 03/98 I. Gonzalez
// Modification: 16/12/99 I.G. Updated for STL.  Needs revision!!! 
//               27/03/00 S.B. In OSCAR
// Description: CMSRotationFactory is a singleton class to get from a file
//              the information to build CMS Rotation matrices.
///////////////////////////////////////////////////////////////////////////////
#ifndef CMSRotationMatrixFactory_h
#define CMSRotationMatrixFactory_h 1

#include "g4std/map"
#include "G4RotationMatrix.hh"

typedef G4RotationMatrix* G4RotationMatrixPtr;
typedef G4std::map<G4String, G4RotationMatrixPtr, less<G4String> > G4RotationMatrixTable;
typedef G4std::map<G4String, G4RotationMatrixPtr, less<G4String> >::iterator G4RotationMatrixTableIterator;


//typedef RWTPtrOrderedVector<G4RotationMatrix> G4RotationMatrixTable;
//  Is an instantiation of the template class RWTPtrOrderedVector.

class CMSRotationMatrixFactory {
public:
  ~CMSRotationMatrixFactory();

  static CMSRotationMatrixFactory* getInstance();
  static CMSRotationMatrixFactory* getInstance(const G4String& rotfile);
  static void setFileName(const G4String& rotfile);

  G4RotationMatrix* findMatrix(const G4String&);
  G4RotationMatrix* AddMatrix(const G4String& name, 
			      G4double th1, G4double phi1,  //Axis angles
			      G4double th2, G4double phi2,  //in rads
			      G4double th3, G4double phi3); //

private:
  CMSRotationMatrixFactory();

private:
  static CMSRotationMatrixFactory* instance;
  static G4String file;

  G4RotationMatrixTable theMatrices; //Where the matrices are stored.
};

ostream& operator<<(ostream&, const G4RotationMatrix &);
#endif
