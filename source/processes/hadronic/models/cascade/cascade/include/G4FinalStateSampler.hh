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
// $Id: G4FinalStateSampler.hh,v 1.2 2009/09/17 18:11:57 dennis Exp $
// GEANT4 tag $Name: geant4-09-03 $
//
// Author: D. H. Wright
// Date:   26 March 2009
//
 
#ifndef G4FinalStateSampler_h
#define G4FinalStateSampler_h 1

// Class Description:
//
// Implementation base class for sampling partial cross sections and final
// states for inelastic hadron-nucleon interactions

#include "globals.hh"
#include <vector>
#include "G4FastVector.hh"
#include "G4ReactionProduct.hh"


class G4FinalStateSampler
{
 public:  // with description
    
   G4FinalStateSampler()
   { }
     
   virtual ~G4FinalStateSampler()
   { }
    
   enum {pi0=7, pip=3, pim=5, kp=11, km=13, k0=15, k0b=17, pro=1, neu=2, 
         lam=21, sp=23, s0=25, sm=27, xi0=29, xim=31, om=33, ap=35, an=37};

 protected:  // with description
    
   std::pair<G4int, G4double> interpolateEnergy(G4double ke) const;

   G4int sampleFlat(std::vector<G4double> sigma) const;

   void CheckQnums(G4FastVector<G4ReactionProduct,256> &vec,
                   G4int &vecLen,
                   G4ReactionProduct &currentParticle,
                   G4ReactionProduct &targetParticle,
                   G4double Q, G4double B, G4double S);

 private:
   
   static const G4double energyScale[30];

};
 
#endif
