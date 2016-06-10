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
// $Id: G4ErrorTrajState.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
// Class Description:
//
// Base class for the trajectory state

// History:
// - Created:   P. Arce
// --------------------------------------------------------------------

#ifndef G4ErrorTrajState_hh
#define G4ErrorTrajState_hh

#include "globals.hh"
#include "G4Track.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4ErrorTrajErr.hh"

enum G4eTSType{ G4eTS_FREE, G4eTS_OS };

class G4ErrorTrajState
{
 public:  // with description

  G4ErrorTrajState() : fCharge(0.), theG4Track(0), iverbose(0) {}

  G4ErrorTrajState( const G4String& partType, const G4Point3D& pos,
                    const G4Vector3D& mom,
                    const G4ErrorTrajErr& errmat = G4ErrorTrajErr(5,0) );
    // Constructor by providing particle, position and momentum

  virtual ~G4ErrorTrajState(){}

  void SetData( const G4String& partType, const G4Point3D& pos,
                const G4Vector3D& mom );
    // Set particle, position and momentum

  void BuildCharge();
    // Build charge based on particle type

  friend
  std::ostream& operator<<(std::ostream&, const G4ErrorTrajState& ts);

  virtual G4int PropagateError( const G4Track* ); 
    // Propagate the error along the step

  virtual G4int Update( const G4Track* ){ return -1; }
    // Update parameters from G4Track
 
  void UpdatePosMom( const G4Point3D& pos, const G4Vector3D& mom );
    // Update position and momentum

  void DumpPosMomError( std::ostream& out = G4cout ) const;
    // Dump position, momentum and error

  virtual void Dump( std::ostream& out = G4cout ) const = 0;
    // Abstract method to dump all TrajState parameters
  
  // Set and Get methods 

  const G4String& GetParticleType() const
    { return fParticleType;}
  void SetParticleType( const G4String& partType )
    { fParticleType = partType;}

  G4Point3D GetPosition() const
    { return fPosition; }
  virtual void SetPosition( const G4Point3D pos )
    { fPosition = pos; }

  G4Vector3D GetMomentum() const
    { return fMomentum; }
  virtual void SetMomentum( const G4Vector3D& mom )
    { fMomentum = mom; }

  G4ErrorTrajErr GetError() const
    { return fError; }
  virtual void SetError( G4ErrorTrajErr em )
    { fError = em; }

  G4Track* GetG4Track() const
    { return theG4Track; }
  void SetG4Track( G4Track* trk )
    { theG4Track = trk; }

  G4double GetCharge() const
    { return fCharge; }
  void SetCharge( G4double ch )
    { fCharge = ch; }

  virtual G4eTSType GetTSType() const
    { return theTSType; }

 protected:

  G4String fParticleType;
  G4Point3D fPosition;
  G4Vector3D fMomentum;
  G4double fCharge;
  G4ErrorTrajErr fError;

  G4eTSType theTSType;

  G4Track* theG4Track;

  G4int iverbose;
};

#endif
