// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVParameterised.hh,v 1.1 1999-01-07 16:07:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVParameterised
//
// Represents many touchable detector elements differing in their
// positioning and dimensions. Both are calculated by means
// of a G4VParameterisation object. The positioning is assumed to
// be dominant along a cartesian axis (specified).
// 
// G4PVParameterised(const G4String& pName,
//	       G4LogicalVolume *pLogical,
//	       G4VPhysicalVolume *pMother,
//             const EAxis pAxis,
//             const G4int nReplicas,
//             G4VPVParameteriastion *pParam)
//
// Replicate the volume nReplicas Times using the paramaterisation pParam,
// withing the mother volume pMother. The positioning of the replicas
// is dominant along the specified axis 
//
// G4PVParameterised(const G4String& pName,
//	       G4LogicalVolume *pLogical,
//	       G4LogicalVolume *pMotherLogical,
//             const EAxis pAxis,
//             const G4int nReplicas,
//             G4VPVParameteriastion *pParam)
//
//  Almost exactly similar to first constructor, changing only mother 
//  pointer's type to LogicalVolume.
//
// History:
// 29.07.95 P.Kent First non-stub version

#ifndef G4PVPARAMETERISED_HH
#define G4PVPARAMETERISED_HH

#include "G4PVReplica.hh"

class G4PVParameterised : public G4PVReplica
{
public:
    G4PVParameterised(const G4String& pName,
		G4LogicalVolume* pLogical,
		G4VPhysicalVolume* pMother,
                const EAxis pAxis,
                const G4int nReplicas,
		G4VPVParameterisation *pParam);

    G4PVParameterised(const G4String& pName,
		G4LogicalVolume* pLogical,
		G4LogicalVolume* pMotherLogical,
                const EAxis pAxis,
                const G4int nReplicas,
		G4VPVParameterisation *pParam);

    virtual G4VPVParameterisation* GetParameterisation() const;
    virtual void GetReplicationData(EAxis& axis,
                                   G4int& nReplicas,
				   G4double& width,
                                   G4double& offset,
                                   G4bool& consuming) const;
private:
    G4VPVParameterisation *fparam; 
};

#endif




