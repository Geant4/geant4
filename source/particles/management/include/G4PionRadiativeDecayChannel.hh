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
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History:
//               30 July 2007 P.Gumplinger
//               Samples Radiative Pion Decay
//               Reference: 
//                    TRIUMF/PIENU Technote: 
//                    "Inclusion of pi->enug in the Monte Carlo" 
//                    by M. Blecher
//
// ------------------------------------------------------------
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4PionRadiativeDecayChannel_h
#define G4PionRadiativeDecayChannel_h 1

#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4VDecayChannel.hh"

class G4PionRadiativeDecayChannel : public G4VDecayChannel
{
  // Class Decription
  // This class describes radiative pion decay kinemtics, but
  // gives incorrect energy spectrum for neutrino

  public:  // With Description

    //Constructors 
    G4PionRadiativeDecayChannel(const G4String& theParentName,
                                G4double        theBR);
    //  Destructor
    virtual ~G4PionRadiativeDecayChannel();

  public:  // With Description

    virtual G4DecayProducts *DecayIt(G4double);

  private:

    G4double beta;

    G4double cib;
    G4double csdp;
    G4double csdm;
    G4double cif;
    G4double cig;

    G4double xl, yl, xu, yu, d2wmax;

    G4double D2W(const G4double x, const G4double y);

};

inline G4double G4PionRadiativeDecayChannel::D2W(const G4double x, 
                                                 const G4double y)
{
  G4double d2w = cib*(1.-y)*(1.+((1.-x)*(1.-x)))/((x*x)*(x+y-1.)) +
                 csdp*(1.-x)*((x+y-1.)*(x+y-1.)) +
                 csdm*(1.-x)*((1.-y)*(1.-y)) +
                 cif*(x-1.)*(1.-y)/x +
                 cig*(1.-y)*(1.-x+(x*x)/(x+y-1.))/x;
  return d2w;
}

#endif
