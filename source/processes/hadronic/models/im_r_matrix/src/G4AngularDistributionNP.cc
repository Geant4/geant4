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
#include "G4AngularDistributionNP.hh"
#include "Randomize.hh"
#include "G4ios.hh"

// Initialization of static data arrays:
#include "G4AngularDistributionNPData.hh"
#include "Randomize.hh"


G4double G4AngularDistributionNP::CosTheta(G4double s, G4double m1, G4double m2) const
{
    G4int verboseLevel=1;

    G4double ek= ((s - sqr(m1) -sqr(m2) )/(2*m1) - m1 )/GeV   ;    // kinetic energy in GeV

    // Find energy bin

    G4int je1 = 0;
    G4int je2 = NENERGY - 1;
    do {
      G4int midBin = (je1 + je2)/2;
      if (ek < elab[midBin])
        je2 = midBin;
      else
        je1 = midBin;
    } while (je2 - je1 > 1);
    //    G4int j;
    //abs(ek-elab[je1]) < abs(ek-elab[je2]) ? j = je1 : j = je2;
    G4double delab = elab[je2] - elab[je1];

    // Sample the angle

    G4float sample = G4UniformRand();
    G4int ke1 = 0;
    G4int ke2 = NANGLE - 1;
    G4double dsig = sig[je2][0] - sig[je1][0];
    G4double rc = dsig/delab;
    G4double b = sig[je1][0] - rc*elab[je1];
    G4double sigint1 = rc*ek + b;
    G4double sigint2 = 0.;

    if (verboseLevel > 1) G4cout << "sample=" << sample << G4endl
                                 << ek << " " << ke1 << " " << ke2 << " "
                                 << sigint1 << " " << sigint2 << G4endl;

    do {
      G4int midBin = (ke1 + ke2)/2;
      dsig = sig[je2][midBin] - sig[je1][midBin];
      rc = dsig/delab;
      b = sig[je1][midBin] - rc*elab[je1];
      G4double sigint = rc*ek + b;
      if (sample < sigint) {
        ke2 = midBin;
        sigint2 = sigint;
      }
      else {
        ke1 = midBin;
        sigint1 = sigint;
      }
      if (verboseLevel > 1)G4cout << ke1 << " " << ke2 << " "
                                  << sigint1 << " " << sigint2 << G4endl;
    } while (ke2 - ke1 > 1);

    // sigint1 and sigint2 should be recoverable from above loop

    //    G4double dsig = sig[je2][ke1] - sig[je1][ke1];
    //    G4double rc = dsig/delab;
    //    G4double b = sig[je1][ke1] - rc*elab[je1];
    //    G4double sigint1 = rc*ek + b;

    //    G4double dsig = sig[je2][ke2] - sig[je1][ke2];
    //    G4double rc = dsig/delab;
    //    G4double b = sig[je1][ke2] - rc*elab[je1];
    //    G4double sigint2 = rc*ek + b;

    dsig = sigint2 - sigint1;
    rc = 1./dsig;
    b = ke1 - rc*sigint1;
    G4double kint = rc*sample + b;
    G4double theta = (0.5 + kint)*pi/180.;

    //    G4int k;
    //abs(sample-sig[j][ke1]) < abs(sample-sig[j][ke2]) ? k = ke1 : k = ke2;
    //    G4double theta = (0.5 + k)*pi/180.;

    if (verboseLevel > 1) {
      G4cout << "   energy bin " << je1 << " energy=" << elab[je1] << G4endl;
      G4cout << "   angle bin " << kint << " angle=" << theta/degree << G4endl;
    }
    G4double costh= cos(theta);
    return costh;
}

G4double G4AngularDistributionNP::Phi() const
{
    return twopi * G4UniformRand();
}
