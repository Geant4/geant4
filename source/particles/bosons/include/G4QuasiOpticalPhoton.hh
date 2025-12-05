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
// Class Description:
// G4QuasiOpticalPhoton is a non-physical particle definition designed to
// encapsulate the conditions for generating a burst of optical photons in
// processes such as scintillation or Cerenkov radiation. Its role is to serve
// as a compact carrier of optical emission metadata—stored via
// G4VUserTrackInformation—including the number of photons to generate, parent
// particle charge, information about the step and track state, and associated
// material properties. This enables deferred or offloaded optical photon
// production and propagation, particularly for hybrid or heterogeneous
// workflows. The term *quasi* (suggested by Daren Sawkey) refers to the
// particle’s representation as a collective distribution of optical photons-
// a synthetic construct used to manage secondary particle generation in the
// Geant4 tracking and optical processes.
#ifndef G4QuasiOpticalPhoton_h
#define G4QuasiOpticalPhoton_h 1

#include "G4ParticleDefinition.hh"

class G4QuasiOpticalPhoton : public G4ParticleDefinition
{
  public:
    static G4QuasiOpticalPhoton* Definition();
    static G4QuasiOpticalPhoton* QuasiOpticalPhotonDefinition();
    static G4QuasiOpticalPhoton* QuasiOpticalPhoton();

  private:
    G4QuasiOpticalPhoton() {}
    ~G4QuasiOpticalPhoton() override = default;

    static G4QuasiOpticalPhoton* theInstance;
};

#endif
