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
 // G4 Low energy model: n-n or p-p scattering
 // F.W. Jones, L.G. Greeniaus, H.P. Wellisch
 //
 // For further comments see G4LEnpData.hh and G4LEnp.cc
 //

#ifndef G4AngularDistributionNP_h
#define G4AngularDistributionNP_h 1

#include "globals.hh"

#include "G4VAngularDistribution.hh"

class G4AngularDistributionNP : public G4VAngularDistribution
{
 private:

//   enum { NENERGY=22, NANGLE=180 };
   enum { NENERGY=39, NANGLE=180 };

 public:

  G4AngularDistributionNP() { }

  virtual ~G4AngularDistributionNP() { }

  virtual G4double CosTheta(G4double s, G4double m1, G4double m2) const;

  virtual G4double Phi() const;

 private:

   static const G4float sig[NENERGY][NANGLE];
   static const G4float pcm[NENERGY], elab[NENERGY],
     dsigmax[NENERGY], sigtot[NENERGY];
};

#endif

