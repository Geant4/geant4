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
// $Id: G4LowEnergyPolarizedRayleigh.hh,v 1.5 2006-06-29 19:36:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//
// File name:     G4LowEnergyPolarizedRayleigh.hh
//
// Author:        Capra Riccardo
//
// Creation date: May 2005
//
// History:
// -----------
// 02 May 2005  R. Capra         1st implementation
//
//----------------------------------------------------------------
        

#ifndef   G4LowEnergyPolarizedRayleigh_hh
 #define  G4LowEnergyPolarizedRayleigh_hh
 
 #include "G4VLowEnergyDiscretePhotonProcess.hh"

// G4LowEnergyPolarizedRayleigh
// The class implements the polarized Rayleigh process
//
//          Specialization of G4VLowEnergyDiscretePhotonProcess class.
//          This is a preliminary implementation of the Rayleigh polarized process.
//          Reference articles are:
// 
// Rayleigh process:                     The Quantum Theory of Radiation
//                                       W. Heitler,       Oxford at the Clarendon Press, Oxford (1954)                                                 
// Scattering function:                   A simple model of photon transport
//                                        D.E. Cullen,      Nucl. Instr. Meth. in Phys. Res. B 101 (1995) 499-510                                       
// Polarization of the outcoming photon:  Beam test of a prototype detector array for the PoGO astronomical hard X-ray/soft gamma-ray polarimeter
//                                        T. Mizuno et al., Nucl. Instr. Meth. in Phys. Res. A 540 (2005) 158-168                                        
// 
class G4LowEnergyPolarizedRayleigh : public G4VLowEnergyDiscretePhotonProcess
{
public:
  //   Class constructor
  //   processName The name of the process
  G4LowEnergyPolarizedRayleigh(const G4String &processName = "polarLowEnRayleigh");

  //   Class destructor
  inline virtual ~G4LowEnergyPolarizedRayleigh(void) {}

  //   Calculates the new photon direction and polarization
  //
  //   In the code the following definitions are used
  //   aTrack The track
  //   aStep  The step
  //   The proposed changes for the photon
  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
   
private:

  //   Hides copy constructor
  G4LowEnergyPolarizedRayleigh(const G4LowEnergyPolarizedRayleigh&);

  //   Hides assignment operator
  G4LowEnergyPolarizedRayleigh & operator=(const G4LowEnergyPolarizedRayleigh&);

  //   Generates \f$cos \left ( \theta\right )\f$ of the scattered photon
  //   incomingPhotonEnergy The energy of the incoming photon
  //   zAtom Atomic number
  //   \f$cos \left ( \theta\right )\f$
  G4double GenerateCosTheta(G4double incomingPhotonEnergy, G4int zAtom) const;
   
  //   Generates \f$\phi\f$ of the scattered photon
  //   cosTheta \f$cos \left ( \theta\right )\f$ of the scattered photon
  //    \f$\phi\f$
  G4double GeneratePhi(G4double cosTheta) const;

  //   Generates the polarization direction \f$\beta\f$ in the plane x, y relative to the x direction
  //   \f$\beta\f$
  G4double GeneratePolarizationAngle(void) const;
   
  //   Lower energy limit due to the available data in the physics tables
  G4double intrinsicLowEnergyLimit;

  //   Higher energy limit due to the available data in the physics tables
  G4double intrinsicHighEnergyLimit;
};

#endif /* G4LowEnergyPolarizedRayleigh_hh */
