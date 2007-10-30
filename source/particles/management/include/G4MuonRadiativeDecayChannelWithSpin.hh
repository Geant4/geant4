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
//               25 July 2007 P.Gumplinger
//               Samples Radiative Muon Decay
//               Reference: 
//                    TRIUMF/TWIST Technote: 
//                    "Radiative muon decay" by P. Depommier and A. Vacheret
//
// ------------------------------------------------------------
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4MuonRadiativeDecayChannelWithSpin_h
#define G4MuonRadiativeDecayChannelWithSpin_h 1

#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4VDecayChannel.hh"

class G4MuonRadiativeDecayChannelWithSpin : public G4VDecayChannel
{
  // Class Decription
  // This class describes radiative muon decay kinemtics, but
  // gives incorrect energy spectrum for neutrinos

  public:  // With Description

    //Constructors
    G4MuonRadiativeDecayChannelWithSpin(const G4String& theParentName,
                                        G4double        theBR);
    //  Destructor
    virtual ~G4MuonRadiativeDecayChannelWithSpin();

  public:  // With Description

    virtual G4DecayProducts *DecayIt(G4double);

    void SetPolarization(G4ThreeVector);
    const G4ThreeVector& GetPolarization() const;

  private:

    G4ThreeVector parent_polarization;

    G4double fron(G4double Pmu, G4double x, G4double y,
                  G4double cthetaE, G4double cthetaG, G4double cthetaEG);

    // rn3dim generates random vectors, uniformly distributed over the surface 
    // of a sphere of given radius.
    // See http://wwwinfo.cern.ch/asdoc/shortwrupsdir/v131/top.html

    void rn3dim(G4double& x, G4double& y, G4double& z, G4double xlong);

    G4double atan4(G4double x, G4double y);

    G4double EMMU;
    G4double EMASS;

};

inline void G4MuonRadiativeDecayChannelWithSpin::
                                        SetPolarization(G4ThreeVector polar)
{
  parent_polarization = polar;
}

inline const G4ThreeVector& G4MuonRadiativeDecayChannelWithSpin::
                                        GetPolarization() const
{
  return parent_polarization;
}

inline void G4MuonRadiativeDecayChannelWithSpin::rn3dim(G4double& x,
                                                        G4double& y,
                                                        G4double& z,
                                                        G4double xlong)
{
      G4double a = 0.; G4double b = 0.; G4double c = 0.; G4double r = 0.;

      do {
         a = G4UniformRand() - 0.5;
         b = G4UniformRand() - 0.5;
         c = G4UniformRand() - 0.5;
         r = a*a + b*b + c*c;
      } while (r > 0.25);

      G4double rinv = xlong/(std::sqrt(r)); 
      x = a * rinv;
      y = b * rinv;
      z = c * rinv;

      return;
}

inline G4double G4MuonRadiativeDecayChannelWithSpin::atan4(G4double x,
                                                           G4double y)
{
      G4double phi = 0.;

      G4double pi = 4.*std::atan(1.);

      if        (x==0. && y>0.){
         phi = 0.5*pi;
      } else if (x==0. && y<0.){
         phi = 1.5*pi;
      } else if (y==0. && x>0.){
         phi = 0.;
      } else if (y==0. && x<0.){
         phi = pi;
      } else if (x>0.         ){
         phi = std::atan(y/x);
      } else if (x<0.         ){
         phi = std::atan(y/x) + pi;
      }

      return phi;
}

#endif
