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
// $Id: G4PVReplica.hh,v 1.6 2010-07-05 13:29:12 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

    virtual void SetRegularStructureId( G4int Code ); 
      // This method must set a unique code for each type of regular structure.
      // - It must be called only during detector construction.
      // - It can also be used to prepare any corresponding special
      //  navigation 'conditions'.

    G4bool IsRegularStructure() const; 
    G4int GetRegularStructureId() const;
      // Accessors for specialised geometries

  protected:

    EAxis faxis;
    G4int fnReplicas;
    G4double fwidth,foffset;
    G4int    fcopyNo;

  private:

    void CheckAndSetParameters(const EAxis pAxis, const G4int nReplicas,
                               const G4double width, const G4double offset);
    G4PVReplica(const G4PVReplica&);
    const G4PVReplica& operator=(const G4PVReplica&);

  private:

    G4int fRegularStructureCode; 
    G4int fRegularVolsId;
};

#endif
