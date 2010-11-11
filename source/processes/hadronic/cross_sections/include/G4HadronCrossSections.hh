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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 Hadron physics class -- header file
// F.W. Jones, TRIUMF, 03-DEC-96
//  
// This class encapsulates cross section data and interpolations 
// from the Geant3/Gheisha routine GHESIG.
// For further comments see G4HadronCrossSections.cc.
//
// Note: this is implemented as a SINGLETON class
//
// 27-MAR-97 FWJ: first version for Alpha release
// 14-APR-97 FWJ: class name changed from G4LCrossSectionData
//    to G4HadronicCrossSections
// 14-APR-98 FWJ: rewritten as class G4HadronCrossSections
//    and adapted to G4CrossSectionDataSet/DataStore class design.
// 26-JUN-98 FWJ: added elastic/inelastic caching to improve performance.
// 06-NOV-98 FWJ: added first-order correction for low-energy
//    inelastic cross sections
//


#ifndef G4HadronCrossSections_h
#define G4HadronCrossSections_h 1
 
#include "globals.hh"
#include "G4Element.hh"
#include "G4VProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4HadTmpUtil.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4AntiNeutron.hh"
#include "G4Lambda.hh"
#include "G4AntiLambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaZero.hh" 
#include "G4SigmaMinus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4XiZero.hh"
#include "G4XiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiXiMinus.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiOmegaMinus.hh"


enum { TSIZE=41, NPARTS=35, NELAB=17, NCNLW=15, NFISS=21 };

class G4HadronCrossSections
{
public:

   G4HadronCrossSections()
     : prevParticleDefinition(0), prevElement(0), lastEkx(0.),
       lastEkxPower(0.), verboseLevel(0)
   {}

   ~G4HadronCrossSections()
   {}

   static G4HadronCrossSections* Instance()
   {
      if (!theInstance) theInstance = new G4HadronCrossSections();
      return theInstance;
   }

   G4bool IsApplicable(const G4DynamicParticle* aParticle,
                       const G4Element* )
   {
     return (GetParticleCode(aParticle) > 0);
   }

   G4bool IsApplicable(const G4DynamicParticle* aParticle,
                         G4double /*ZZ*/, G4double /*AA*/)
   {
     return (GetParticleCode(aParticle) > 0);
   }

   G4double GetElasticCrossSection(const G4DynamicParticle*, const G4Element*);

   G4double GetElasticCrossSection(const G4DynamicParticle*,
                                   G4int /*ZZ*/, G4int /*AA*/);

   G4double GetInelasticCrossSection(const G4DynamicParticle*, 
                                     const G4Element*);

   G4double GetInelasticCrossSection(const G4DynamicParticle*,
                                     G4int /*ZZ*/, G4int /*AA*/);

   G4double GetCaptureCrossSection(const G4DynamicParticle* aParticle,
                                   const G4Element* anElement)
   {
     G4int Z = G4lrint(anElement->GetZ());
     G4int A = G4lrint(anElement->GetN());
     return GetCaptureCrossSection(aParticle, Z, A);
   }

   G4double GetCaptureCrossSection(const G4DynamicParticle*,
                                   G4int /*ZZ*/, G4int /*AA*/);

   G4double GetFissionCrossSection(const G4DynamicParticle* aParticle,
                                   const G4Element* anElement)
   {
     G4int Z = G4lrint(anElement->GetZ());
     G4int A = G4lrint(anElement->GetN());
     return GetFissionCrossSection(aParticle, Z, A);
   }

   G4double GetFissionCrossSection(const G4DynamicParticle*,
                                   G4int /*ZZ*/, G4int /*AA*/);


   static void SetCorrectInelasticNearZero(G4bool value)
   {
      correctInelasticNearZero = value;
   }

   static G4bool GetCorrectInelasticNearZero()
   {
      return correctInelasticNearZero;
   }

   void SetVerboseLevel(G4int value)
   {
      verboseLevel = value;
   }

   G4int GetVerboseLevel()
   {
      return verboseLevel;
   }

private:

   G4int GetParticleCode(const G4DynamicParticle*);

   void CalcScatteringCrossSections(const G4DynamicParticle*, 
                                    G4int /*ZZ*/, G4int /*AA*/);

   static G4HadronCrossSections* theInstance;

   G4double sigelastic, siginelastic;
   G4ParticleDefinition* prevParticleDefinition;
   G4Element* prevElement;
   G4double prevKineticEnergy;
   G4double lastEkx, lastEkxPower;

   static G4bool correctInelasticNearZero;

   G4int verboseLevel;

// The following arrays are declared static to allow the use of initializers.  
// They are initialized in G4HadronCrossSections.cc, thus providing some 
// data hiding.

   static G4float plab[TSIZE];
   static G4float csel[NPARTS][TSIZE];
   static G4float csin[NPARTS][TSIZE];

   static G4float cspiel[3][TSIZE];
   static G4float cspiin[3][TSIZE];

   static G4float cspnel[3][TSIZE];
   static G4float cspnin[3][TSIZE];

   static G4float elab[NELAB];
   static G4float cnlwat[NCNLW], cnlwel[NCNLW][NELAB], cnlwin[NCNLW][NELAB];

   static G4float cscap[100];

   static G4float ekfiss[NFISS], csfiss[4][NFISS];

   static G4float alpha[NPARTS], alphac[TSIZE];

   static G4float partel[35], partin[35];
   static G4int   icorr[35], intrc[35];

   static G4float csa[4];
   static G4int ipart2[7];
};
#endif
