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
//  Author: Dennis Wright (SLAC)
//  Date:   5 August 2011
//

#ifndef G4BetaDecayCorrections_h
#define G4BetaDecayCorrections_h 1

#include "G4BetaDecayType.hh"


class G4BetaDecayCorrections
{
  public:
    G4BetaDecayCorrections(const G4int Z, const G4int A);
    
    ~G4BetaDecayCorrections() {};

    G4double FermiFunction(const G4double& W);

    G4double ShapeFactor(const G4BetaDecayType&,
                         const G4double& p_e, const G4double& e_nu);

  private:
    G4double ModSquared(const G4double& x, const G4double& y);
    G4double Gamma(const G4double& arg);

    const G4int Z;
    const G4int A;
    G4double alphaZ;
    G4double Rnuc;    // Nuclear radius  
    G4double V0;      // Electron screening potential  
    G4double gamma0;

    G4double gc[6];   // Real gamma function polynomial coefficients 
};

#endif

