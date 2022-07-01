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
// 20110907  M. Kelsey -- Fix constness for accessors

#ifndef G4HadSecondary_hh
#define G4HadSecondary_hh

#include "globals.hh"

class G4DynamicParticle;
class G4ParticleDefinition;

class G4HadSecondary
{
public:

  G4HadSecondary(G4DynamicParticle * aT, G4double aWeight = 1, G4int mod=-1);
  ~G4HadSecondary();

  inline G4DynamicParticle * GetParticle()   {return theP;}
  inline const G4DynamicParticle* GetParticle() const {return theP;}
  inline G4double GetWeight() const          {return theWeight;}
  inline void SetWeight(G4double aW)         {theWeight= aW;}
  inline void SetTime(G4double aT)           {theTime = aT;}
  inline G4double GetTime() const            {return theTime;}
  inline void SetCreatorModelID(G4int id)    {theCreatorModel = id;}
  inline G4int GetCreatorModelID() const     {return theCreatorModel;}
  inline const G4ParticleDefinition* GetParentResonanceDef() const {return theParentResonanceDef;}
  inline void SetParentResonanceDef(const G4ParticleDefinition* parentDef) {theParentResonanceDef = parentDef;}
  inline G4int GetParentResonanceID() const {return theParentResonanceID;}
  inline void SetParentResonanceID(const G4int parentID) {theParentResonanceID = parentID;}
   
private:
  
   G4HadSecondary(){};
   
   G4DynamicParticle * theP; 
   G4double theWeight;
   G4double theTime;
   G4int theCreatorModel;
   const G4ParticleDefinition* theParentResonanceDef = nullptr;
   G4int theParentResonanceID;
};

#endif
