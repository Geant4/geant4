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
// $Id: G4ReplicatedSlice.hh 66356 2012-12-18 09:02:32Z gcosmo $
// 
// class G4ReplicatedSlice
//
// Class description:
//
// G4ReplicatedSlice represents many touchable detector elements differing
// only in their positioning. The elements' positions are calculated by means
// of a simple linear formula.
// 
// Division may occur along:
//
// o Cartesian axes (kXAxis,kYAxis,kZAxis)
//
//   The divisions, of specified width have coordinates of
//   form (-width*(nReplicas-1)*0.5+n*width,0,0) where n=0.. nReplicas-1
//   for the case of kXAxis, and are unrotated.
//
// o Radial axis (cylindrical polar) (kRho)
//
//   The divisions are cons/tubs sections, centred on the origin
//   and are unrotated.
//   They have radii of width*n+offset to width*(n+1)+offset
//                      where n=0..nReplicas-1
//
// o Phi axis (cylindrical polar) (kPhi)
//   The divisions are `phi sections' or wedges, and of cons/tubs form
//   They have phi of offset+n*width to offset+(n+1)*width where
//   n=0..nReplicas-1

// History:
// -------
// Apr/20/2010 - Initial version extended from G4PVDivision
// ----------------------------------------------------------------------
#ifndef G4ReplicatedSlice_HH 
#define G4ReplicatedSlice_HH 1

#include "geomdefs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VDivisionParameterisation.hh"

class G4LogicalVolume;
class G4VSolid;

class G4ReplicatedSlice : public G4VPhysicalVolume
{
  public:  // with description
    
    G4ReplicatedSlice(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4LogicalVolume* pMotherLogical,
                      const EAxis pAxis,
                      const G4int nReplicas,
                      const G4double width,
                      const G4double half_gap,
                      const G4double offset );
      // Constructor with number of divisions and width

    G4ReplicatedSlice(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4LogicalVolume* pMotherLogical,
                      const EAxis pAxis,
                      const G4int nReplicas,
                      const G4double half_gap,
                      const G4double offset );
      // Constructor with number of divisions 

    G4ReplicatedSlice(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4LogicalVolume* pMotherLogical,
                      const EAxis pAxis,
                      const G4double width,
                      const G4double half_gap,
                      const G4double offset );
      // Constructor with width

  public:  // without description

    G4ReplicatedSlice(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4VPhysicalVolume* pMotherPhysical,
                      const EAxis pAxis,
                      const G4int nReplicas,
                      const G4double width,
                      const G4double half_gap,
                      const G4double offset);
      // Constructor in mother physical volume

    G4ReplicatedSlice(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4VPhysicalVolume* pMotherPhysical,
                      const EAxis pAxis,
                      const G4int nReplicas,
                      const G4double half_gap,
                      const G4double offset );
      // Constructor with number of divisions 

    G4ReplicatedSlice(const G4String& pName,
                            G4LogicalVolume* pLogical,
                            G4VPhysicalVolume* pMotherPhysical,
                      const EAxis pAxis,
                      const G4double width,
                      const G4double half_gap,
                      const G4double offset );
      // Constructor with width

  public:  // with description

    virtual ~G4ReplicatedSlice();

    virtual G4bool IsMany() const;
    virtual G4int GetCopyNo() const;
    virtual void  SetCopyNo(G4int CopyNo);
    virtual G4bool IsReplicated() const;
    virtual G4VPVParameterisation* GetParameterisation() const;
    virtual void GetReplicationData( EAxis& axis,
                                     G4int& nReplicas,
                                     G4double& width,
                                     G4double& offset,
                                     G4bool& consuming ) const;
    EAxis  GetDivisionAxis() const;
    G4bool IsParameterised() const;

  public:  // without description

    G4bool IsRegularStructure() const; 
    G4int  GetRegularStructureId() const; 
      // Methods to identify volume that can apply 'regular' navigation.
      // Currently divisions do not qualify for this.

  private:

    void CheckAndSetParameters( const EAxis pAxis,
                                const G4int nDivs,
                                const G4double width,
                                const G4double half_gap,
                                const G4double offset, 
                                      DivisionType divType,
                                      G4LogicalVolume* pMotherLogical,
                                const G4LogicalVolume* pLogical );

    void SetParameterisation(       G4LogicalVolume* motherLogical,
                              const EAxis pAxis,
                              const G4int nReplicas,
                              const G4double width,
                              const G4double half_gap,
                              const G4double offset, 
                                    DivisionType divType );

    void ErrorInAxis( EAxis axis, G4VSolid* solid );

  private:

    G4ReplicatedSlice(const G4ReplicatedSlice&);
    const G4ReplicatedSlice& operator=(const G4ReplicatedSlice&);
      // Private copy constructor and assignment operator.

  protected:

    EAxis faxis;             // axis of optimisation
    EAxis fdivAxis;          // axis of division
    G4int fnReplicas;
    G4double fwidth,foffset;
    G4int    fcopyNo;
    G4VDivisionParameterisation *fparam; 
};

#endif
