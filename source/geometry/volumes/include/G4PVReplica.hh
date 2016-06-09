//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PVReplica.hh,v 1.4 2005/11/11 22:39:00 japost Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// class G4PVReplica
//
// Class description:
//
// Represents many touchable detector elements differing only in their
// positioning. The elements' positions are calculated by means of a simple
// linear formula, and the elements completely fill the containing mother
// volume.
// 
// G4PVReplica(const G4String& pName,
//                   G4LogicalVolume *pLogical,
//                   G4LogicalVolume *pMother,
//             const EAxis pAxis,
//             const G4int nReplicas,
//             const G4double width,
//             const G4double offset=0);
//
// Replication may occur along:
//
// o Cartesian axes (kXAxis,kYAxis,kZAxis)
//
//   The replications, of specified width have coordinates of
//   form (-width*(nReplicas-1)*0.5+n*width,0,0) where n=0.. nReplicas-1
//   for the case of kXAxis, and are unrotated.
//
// o Radial axis (cylindrical polar) (kRho)
//
//   The replications are cons/tubs sections, centred on the origin
//   and are unrotated.
//   They have radii of width*n+offset to width*(n+1)+offset
//                      where n=0..nReplicas-1
//
// o Phi axis (cylindrical polar) (kPhi)
//   The replications are `phi sections' or wedges, and of cons/tubs form
//   They have phi of offset+n*width to offset+(n+1)*width where
//   n=0..nReplicas-1

// History:
// 29.07.95 P.Kent         First non-stub version
// 26.10.97 J.Apostolakis  Added constructor that takes mother logical volume
// 16.02.98 J.Apostolakis  Added copy number
// ----------------------------------------------------------------------
#ifndef G4PVREPLICA_HH
#define G4PVREPLICA_HH

#include "G4VPhysicalVolume.hh"

class G4PVReplica : public G4VPhysicalVolume
{
  public:  // with description

    G4PVReplica(const G4String& pName,
                      G4LogicalVolume* pLogical,
                      G4LogicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
                const G4double width,
                const G4double offset=0);

  public:  // without description

    G4PVReplica(const G4String& pName,
                      G4LogicalVolume* pLogical,
                      G4VPhysicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
                const G4double width,
                const G4double offset=0);

    G4PVReplica(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

  public:  // with description

    virtual ~G4PVReplica();

    G4bool IsMany() const;
    G4bool IsReplicated() const;

    virtual G4int GetCopyNo() const;
    virtual void  SetCopyNo(G4int CopyNo);
    virtual G4bool IsParameterised() const;
    virtual G4VPVParameterisation* GetParameterisation() const;
    virtual G4int GetMultiplicity() const;
    virtual void GetReplicationData(EAxis& axis,
                                    G4int& nReplicas,
                                    G4double& width,
                                    G4double& offset,
                                    G4bool& consuming) const;

    // Methods for specialised geometries
    G4bool  IsRegularStructure()   const; 
    G4int  GetRegularStructureId() const;
      // return values
    virtual void SetRegularStructureId( G4int Code ); 
      // This method must set a unique code for each type
      //  of regular structure.
      //  - It must be called only during detector construction.
      //  - It can also be used to prepare a/any corresponding 
      //    special navigation 'conditions'.
 private:

    void CheckAndSetParameters(
                         const EAxis pAxis,
                         const G4int nReplicas,
                         const G4double width,
                         const G4double offset);

    G4PVReplica(const G4PVReplica&);
    const G4PVReplica& operator=(const G4PVReplica&);

 private:
    G4int fRegularStructureCode; 

 protected:

    EAxis faxis;
    G4int fnReplicas;
    G4double fwidth,foffset;
    
    G4int    fcopyNo;

  private:
    G4int fRegularVolsId; 

 };

#endif
