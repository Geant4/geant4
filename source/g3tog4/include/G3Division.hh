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
// $Id: G3Division.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class description:
//
// This class tranforms G3 divided volumes to G4 replicated volumes. 
// UpdateVTE() method checks parameters of mother volume
// and in case they are complete the solid that will be replicated 
// is created. In case of division with offset an additinal envelope 
// VTE (G3VolTableEntry instance) is created.
// CreatePVReplica() methods creates the G4PVReplica instance.

// ----------------------
//
// by I.Hrivnacova, V.Berejnoi, 27 Sep 99

#ifndef G3DIVISION_HH
#define G3DIVISION_HH 1

#include "globals.hh"
#include "geomdefs.hh"

class G3VolTableEntry;
class G4VPhysicalVolume;
class G4LogicalVolume;

enum G3DivType { kDvn, kDvn2, kDvt, kDvt2 };

class G3Division 
{
  public: // with description

    G3Division(G3DivType type, G3VolTableEntry* vte, G3VolTableEntry* mvte, 
               G4int nofDivision, G4int iaxis, G4int nmed, G4double c0, 
	       G4double step);
    G3Division(G3VolTableEntry* vte, G3VolTableEntry* mvte,
               const G3Division& division);
    virtual ~G3Division();
    
    // methods
    void UpdateVTE();
    void CreatePVReplica();   
    
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
