// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPVParameterised.ddl,v 1.4 1999/12/15 14:51:23 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// class G4PPVParameterised
//
//  Persistent-capable class of G4PVParameterised
//
// History:
// 05.11.99  Y.Morita  First non-stub version

#ifndef G4PPVPARAMETERISED_DDL
#define G4PPVPARAMETERISED_DDL 1

#include "G4PersistentSchema.hh"
#include "G4PPVReplica.hh"

class G4PPVParameterised
 : public G4PPVReplica
{
public:
  G4PPVParameterised(G4VPhysicalVolume *PhysVol,
                       HepRef(G4PLogicalVolume) persLogVol);

  ~G4PPVParameterised();

  G4VPhysicalVolume* MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother );

private:
//    d_Ref<G4VPVParameterisation> fparam; 
};

#endif

