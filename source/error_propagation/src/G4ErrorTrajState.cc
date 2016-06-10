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
// $Id: G4ErrorTrajState.cc 78318 2013-12-11 15:02:40Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorTrajState.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ErrorPropagatorData.hh"

#include <iomanip>

//--------------------------------------------------------------------------
G4ErrorTrajState::G4ErrorTrajState( const G4String& partType,
                                    const G4Point3D& pos,
                                    const G4Vector3D& mom,
                                    const G4ErrorTrajErr& errmat)
  : fParticleType(partType), fPosition(pos), fMomentum(mom), fCharge(0.),
    fError(errmat), theTSType(G4eTS_FREE), theG4Track(0)
{
  iverbose = G4ErrorPropagatorData::verbose();
}


//--------------------------------------------------------------------------
G4int G4ErrorTrajState::PropagateError( const G4Track* )
{ 
   std::ostringstream message;
   message << "Wrong trajectory state type !" << G4endl
           << "Called for trajectory state type: " << G4int(GetTSType());
   G4Exception("G4ErrorTrajState::PropagateError()", "GEANT4e-Error",
               FatalException, message);
   return -1; 
}


//--------------------------------------------------------------------------
void G4ErrorTrajState::UpdatePosMom( const G4Point3D& pos,
                                     const G4Vector3D& mom )
{
  fPosition = pos;
  fMomentum = mom;
}


//--------------------------------------------------------------------------
void G4ErrorTrajState::SetData( const G4String& partType,
                                const G4Point3D& pos, const G4Vector3D& mom )
{
  fParticleType = partType;
  BuildCharge();
  fPosition = pos;
  fMomentum = mom;
}


//--------------------------------------------------------------------------
void G4ErrorTrajState::BuildCharge()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(fParticleType); 
  if( particle == 0)
  {
    std::ostringstream message;
    message << "Particle type not defined: " << fParticleType;
    G4Exception( "G4ErrorTrajState::BuildCharge()", "GEANT4e-error",
                  FatalException, message);
  }
  else
  {
    fCharge = particle->GetPDGCharge();
  }
}


//------------------------------------------------------------------------
void G4ErrorTrajState::DumpPosMomError( std::ostream& out ) const
{
  out << *this;
}


//--------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const G4ErrorTrajState& ts)
{
  //  long mode = out.setf(std::ios::fixed,std::ios::floatfield);
  out  
    << " G4ErrorTrajState of type " << ts.theTSType << " : partycle: "
    << ts.fParticleType << "  position: " << std::setw(6) << ts.fPosition
    << "              momentum: " << ts.fMomentum
    << "   error matrix ";
  G4cout << ts.fError << G4endl;

  return out;
}

