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
// $Id: G4PVParameterised.cc,v 1.7 2002-10-23 16:11:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4PVParameterised
//
// Implementation
//
// ********************************************************************

#include "G4PVParameterised.hh"

G4PVParameterised::G4PVParameterised( const G4String& pName,
                                            G4LogicalVolume* pLogical,
                                            G4VPhysicalVolume* pMother,
                                      const EAxis pAxis,
                                      const G4int nReplicas,
                                            G4VPVParameterisation *pParam)
  : G4PVReplica(pName, pLogical, pMother, pAxis, nReplicas, 0, 0),
    fparam(pParam)
{
#ifdef G4VERBOSE  
  if ((pMother) && (pMother->IsParameterised()))
  {
    G4cout << "WARNING - A parameterised volume is being placed" << G4endl
           << "          inside another parameterised volume !" << G4endl
           << "          To make sure that no overlaps are generated," << G4endl
           << "          you should verify the mother replicated shapes" << G4endl
           << "          are of the same type and dimensions." << G4endl
           << "             Mother physical volume: " << pMother->GetName() << G4endl
           << "             Parameterised volume: " << pName << G4endl
           << "  (To switch this warning off, compile with G4_NO_VERBOSE)" << G4endl;
  }
#endif
}

G4PVParameterised::G4PVParameterised( const G4String& pName,
                                            G4LogicalVolume* pLogical,
                                            G4LogicalVolume* pMotherLogical,
                                      const EAxis pAxis,
                                      const G4int nReplicas,
                                            G4VPVParameterisation *pParam)
  : G4PVReplica(pName, pLogical, pMotherLogical, pAxis, nReplicas, 0, 0),
    fparam(pParam)
{
}

G4PVParameterised::~G4PVParameterised()
{
}

G4VPVParameterisation* G4PVParameterised::GetParameterisation() const
{
  return fparam;
}

G4bool G4PVParameterised::IsParameterised() const
{
  return true;
}

void G4PVParameterised::GetReplicationData( EAxis& axis,
                                            G4int& nReplicas,
                                            G4double& width,
                                            G4double& offset,
                                            G4bool& consuming) const
{
  axis = faxis;
  nReplicas = fnReplicas;
  width = fwidth;
  offset = foffset;
  consuming = false;
}
