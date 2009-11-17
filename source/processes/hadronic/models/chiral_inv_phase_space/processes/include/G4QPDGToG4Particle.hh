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
// $Id: G4QPDGToG4Particle.hh,v 1.1 2009-11-17 10:36:54 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QPDGToG4Particle header ----------------
//                 by Mikhail Kossov, December 2003.
//  Header of the singletone class of the CHIPS Simulation Branch in GEANT4
// ------------------------------------------------------------------------
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************
// Short description: This is a helper class, which converts the PDG-defined
// G4QHadrons of the CHIPS model to the G4 particles, defined by the singetones.
// -----------------------------------------------------------------------------

#ifndef G4QPDGToG4Particle_hh
#define G4QPDGToG4Particle_hh

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

class G4QPDGToG4Particle
{
  // Constructor/Destructor
protected:
  G4QPDGToG4Particle();        // the Default Construction is protected - Singelton
public:
  ~G4QPDGToG4Particle();       // Destructor is public because of Windows compilation error

  // Member Functions

public:
  // Pointers to Particles of the Singeltone of the CHIPS World
  static G4QPDGToG4Particle* Get();
  G4ParticleDefinition* GetParticleDefinition(G4int PDGCode);
  void DefineAllParticles();

// Body
private:
};
#endif
