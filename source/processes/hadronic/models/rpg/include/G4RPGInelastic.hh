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
// $Id: G4RPGInelastic.hh 81944 2014-06-06 15:57:38Z gcosmo $
//
// Author: D. H. Wright
// Date:   26 May 2007
//
 
#ifndef G4RPGInelastic_h
#define G4RPGInelastic_h 1

// Class Description:
//
// Base class for re-parameterized Gheisha-style models.

#include <vector>

#include "globals.hh"
#include "G4FastVector.hh"
#include "G4HadronicInteraction.hh"
#include "G4ReactionProduct.hh"
#include "Randomize.hh"
#include "G4RPGFragmentation.hh"
#include "G4RPGTwoCluster.hh"
#include "G4RPGTwoBody.hh"
#include "G4RPGStrangeProduction.hh"
#include "G4RPGPionSuppression.hh"

enum{ GHADLISTSIZE=256};


class G4RPGInelastic : public G4HadronicInteraction
 {
 public:  // with description
    
   G4RPGInelastic(const G4String& modelName = "RPGInelastic"); 
    
   virtual ~G4RPGInelastic()
   { }
    
 protected:  // with description
    
   G4double Pmltpc(G4int np, G4int nm, G4int nz, G4int n,
                   G4double b, G4double c);

   G4int Factorial(G4int n);
    
   G4bool MarkLeadingStrangeParticle(const G4ReactionProduct& currentParticle,
                                     const G4ReactionProduct& targetParticle,
                                     G4ReactionProduct& leadParticle);
    
   void SetUpPions(const G4int np, const G4int nm, const G4int nz,
                   G4FastVector<G4ReactionProduct,256> &vec,
                   G4int &vecLen);
    
   //   void Rotate(G4FastVector<G4ReactionProduct,256> &vec, G4int &vecLen);

   void GetNormalizationConstant(const G4double availableEnergy,
                                 G4double &n,
                                 G4double &anpn);
    
   void CalculateMomenta(G4FastVector<G4ReactionProduct,256> &vec,
                         G4int &vecLen,
                         const G4HadProjectile *originalIncident,
                         const G4DynamicParticle *originalTarget,
                         G4ReactionProduct &modifiedOriginal,
                         G4Nucleus &targetNucleus,
                         G4ReactionProduct &currentParticle,
                         G4ReactionProduct &targetParticle,
                         G4bool &incidentHasChanged,
                         G4bool &targetHasChanged,
                         G4bool quasiElastic);
    
   void SetUpChange(G4FastVector<G4ReactionProduct,256> &vec,
                    G4int &vecLen,
                    G4ReactionProduct &currentParticle,
                    G4ReactionProduct &targetParticle,
                    G4bool &incidentHasChanged);

   G4RPGFragmentation fragmentation;

   G4RPGTwoCluster twoCluster;

   G4RPGPionSuppression pionSuppression;

   G4RPGStrangeProduction strangeProduction;

   G4RPGTwoBody twoBody;

   std::pair<G4int, G4double> interpolateEnergy(G4double ke) const;

   G4int sampleFlat(std::vector<G4double> sigma) const;

   void CheckQnums(G4FastVector<G4ReactionProduct,256> &vec,
                   G4int &vecLen,
                   G4ReactionProduct &currentParticle,
                   G4ReactionProduct &targetParticle,
                   G4double Q, G4double B, G4double S);

   enum {pi0, pip, pim, kp, km, k0, k0b, pro, neu, 
         lam, sp, s0, sm, xi0, xim, om, ap, an};

 protected:

   G4ParticleDefinition* particleDef[18];

 private:
   
   G4double cache;
   G4ThreeVector what;

   static const G4double energyScale[30];

 };
 
#endif
