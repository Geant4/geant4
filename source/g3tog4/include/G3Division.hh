// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3Division.hh,v 1.2 1999-12-05 17:50:01 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, V.Berejnoi, 27 Sep 99

#ifndef G3_DIVISION_H
#define G3_DIVISION_H

#include "globals.hh"
#include "geomdefs.hh"

class G3VolTableEntry;
class G4VPhysicalVolume;
class G4LogicalVolume;

enum G3DivType { kDvn, kDvn2, kDvt, kDvt2 };

class G3Division 
{
  public:
    G3Division(G3DivType type, G3VolTableEntry* vte, G3VolTableEntry* mvte, 
               G4int nofDivision, G4int iaxis, G4int nmed, G4double c0, 
	       G4double step);
    G3Division(G3VolTableEntry* vte, G3VolTableEntry* mvte,
               const G3Division& division);
    virtual ~G3Division();
    
    // methods
    void UpdateVTE();
    G4VPhysicalVolume* CreatePVReplica();   
    
  private:
    // methods
    void SetRangeAndAxis();
    void CreateSolid(G4String shape, G4double par[], G4int npar);
    G3VolTableEntry* CreateEnvelope(G4String shape, G4double hi, G4double lo, 
                       G4double par[], G4int npar);
    G3VolTableEntry* Dvn ();
    G3VolTableEntry* Dvn2();
    G3VolTableEntry* Dvt ();
    G3VolTableEntry* Dvt2();
    void Exception(G4String where, G4String what);
    
    // data members 
    G3DivType         fType;
    G3VolTableEntry*  fVTE;    
    G3VolTableEntry*  fMVTE;    
    G4int             fNofDivisions;  // ndiv/ndvmx    
    G4int             fIAxis;
    G4int             fNmed;
    G4double          fC0;
    G4double          fStep;
    G4double          fLowRange;  
    G4double          fHighRange;
    G4double          fWidth;
    G4double          fOffset;
    EAxis             fAxis;
};

#endif //G3_DIVISION_H
