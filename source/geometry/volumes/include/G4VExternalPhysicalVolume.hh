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
// G4VExternalPhysicalVolume
//
// Class description:
//
// Base class to represent a physical volume managed by an external 
// sub-navigator.
// 
// Intial assumptions:
//   * volume type is similar to G4PVPlacement -- not replicated
//   * external navigator may provide 'many'/Boolean operation

// Author: J.Apostolakis, October 2019
// ----------------------------------------------------------------------
#ifndef G4VEXTERNALPHYSICSVOLUME_HH
#define G4VEXTERNALPHYSICSVOLUME_HH

#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"

class G4VExternalPhysicalVolume : public G4VPhysicalVolume
{
  public:  // with description

    G4VExternalPhysicalVolume( G4RotationMatrix* pRot,
                               const G4ThreeVector& tlate,
                               G4LogicalVolume* pCurrentLogical,                            
                               const G4String& pName,
                               G4VPhysicalVolume* pMother );
       
    G4VExternalPhysicalVolume(const G4VExternalPhysicalVolume&) = delete;
    G4VExternalPhysicalVolume& operator=(const G4VExternalPhysicalVolume&) = delete;
      // Forbidden copy constructor and assignment operator.

    ~G4VExternalPhysicalVolume() override;
      // Default destructor.

    G4bool CheckOverlaps(G4int res=1000, G4double tol=0.,
                         G4bool verbose=true, G4int maxErr=1) override = 0;
      // Verifies if the placed volume is overlapping with existing
      // daughters or with the mother volume. Provides default resolution
      // for the number of points to be generated and verified.
      // A tolerance for the precision of the overlap check can be specified,
      // by default it is set to maximum precision.
      // Reports a maximum of overlaps errors according to parameter in input.
      // Returns true if the volume is overlapping.

  public:  // without description

    G4VExternalPhysicalVolume(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    EVolume VolumeType() const final;
    G4bool IsMany() const final;
   
    G4bool IsReplicated() const final;
    G4bool IsParameterised() const final;
    G4VPVParameterisation* GetParameterisation() const final;
    void GetReplicationData(EAxis& axis,
                            G4int& nReplicas,
                            G4double& width,
                            G4double& offset,
                            G4bool& consuming) const final;

    // Methods for use with Regular Navigation
    //
    G4bool  IsRegularStructure() const final; 
    G4int   GetRegularStructureId() const final; 

    void SetMany(G4bool overlap);
      // Set the 'many' flag
   
  private:

    G4bool fMany = false;
};

#endif

