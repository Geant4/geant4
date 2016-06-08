// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPVReplica.ddl,v 1.5 2000/11/02 12:42:11 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// class G4PPVReplica
//
//  Persistent-capable class of G4PVReplica
//
// History:
// 05.11.99  Y.Morita  First non-stub version

#ifndef G4PPVREPLICA_DDL
#define G4PPVREPLICA_DDL

#include "G4Pglobals.hh"

#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"
#include "G4VPhysicalVolume.hh"
// #include "G4RotationMatrix.hh"

#include "G4PVPhysicalVolume.hh"

class G4PPVReplica
 : public G4PVPhysicalVolume
{
public:
  G4PPVReplica( G4VPhysicalVolume *PhysVol,
                HepRef(G4PLogicalVolume) persLogVol);

  ~G4PPVReplica();

  G4VPhysicalVolume* MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother );

  virtual G4bool IsMany() const;
  virtual G4int GetCopyNo() const;

protected:
  virtual void  SetCopyNo(G4int CopyNo);

protected:
  EAxis faxis;
  G4Pint fnReplicas;
  G4Pdouble fwidth,foffset;
    
  G4Pint    fcopyNo;
};

#endif



