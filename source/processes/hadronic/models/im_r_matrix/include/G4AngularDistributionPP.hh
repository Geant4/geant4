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
 //
 // G4 Low energy model: n-n or p-p scattering
 // F.W. Jones, L.G. Greeniaus, H.P. Wellisch
 //
 // For further comments see G4LEnpData.hh and G4LEnp.cc
 //

#ifndef G4AngularDistributionPP_h
#define G4AngularDistributionPP_h 1

#include "globals.hh"

#include "G4VAngularDistribution.hh"

class G4AngularDistributionPP : public G4VAngularDistribution
{
 private:

   enum { NENERGY=22, NANGLE=180 };

 public:

  G4AngularDistributionPP() { }

  virtual ~G4AngularDistributionPP() { }

  virtual G4double CosTheta(G4double s, G4double m1, G4double m2) const;

  virtual G4double Phi() const;

 private:

   static G4float sig[NENERGY][NANGLE];
   static G4float pcm[NENERGY], elab[NENERGY],
     dsigmax[NENERGY], sigtot[NENERGY];
};

#endif

