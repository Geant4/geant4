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
// G4PhysicsLogVector
//
// Class description:
//
// A physics vector which has values of energy-loss, cross-section,
// and other physics values of a particle in matter in a given
// range of energy, momentum, etc. The scale of energy/momentum
// bins is logarithmic.

// Authors:
// - 02 Dec. 1995, G.Cosmo: Structure created based on object model
// - 03 Mar. 1996, K.Amako: Implemented the 1st version
// Revisions:
// - 11 Nov. 2000, H.Kurashige : Use STL vector for dataVector and binVector
// --------------------------------------------------------------------
#ifndef G4PhysicsLogVector_hh
#define G4PhysicsLogVector_hh 1

#include "G4PhysicsVector.hh"
#include "globals.hh"

class G4PhysicsLogVector : public G4PhysicsVector
{
public:
  // The vector will be filled from external file using Retrieve() method
  explicit G4PhysicsLogVector(G4bool spline = false);

  // Energies will be computed and filled at construction, values will be 
  // filled with zeros. Required Nbin > 1 and Emax > Emin > 0.
  // Use PutValue(..) to fill the data vector
  explicit G4PhysicsLogVector(G4double Emin, G4double Emax, std::size_t Nbin,
                              G4bool spline = false);

  ~G4PhysicsLogVector() override = default;

protected:

  void Initialise() final;
};

#endif
