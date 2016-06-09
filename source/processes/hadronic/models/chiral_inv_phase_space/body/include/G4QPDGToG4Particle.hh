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
// $Id: G4QPDGToG4Particle.hh,v 1.3 2005/08/30 07:15:14 mkossov Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//      ---------------- G4QPDGToG4Particle header ----------------
//                 by Mikhail Kossov, December 2003.
//  Header of the singletone class of the CHIPS Simulation Branch in GEANT4
// ------------------------------------------------------------------------
// ****************************************************************************************
// ********* This HEADER is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

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
