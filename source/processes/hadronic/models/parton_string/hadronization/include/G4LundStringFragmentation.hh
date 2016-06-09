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
// $Id: G4LundStringFragmentation.hh,v 1.4 2007/04/24 14:55:23 gunter Exp $
// GEANT4 tag $Name: geant4-08-03 $ Maxim Komogorov
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------

#ifndef G4LundStringFragmentation_h
#define G4LundStringFragmentation_h 1

#include "G4VLongitudinalStringDecay.hh"

//**************************************************************************************************************

class G4LundStringFragmentation: public G4VLongitudinalStringDecay
    {
public:
    G4LundStringFragmentation();
//    G4LundStringFragmentation(G4double sigmaPt);
     G4LundStringFragmentation(const G4LundStringFragmentation &right);
     virtual ~G4LundStringFragmentation();
     virtual G4KineticTrackVector* FragmentString(const G4ExcitedString& theString);

public:
    const G4LundStringFragmentation & operator=(const G4LundStringFragmentation &right);
    int operator==(const G4LundStringFragmentation &right) const;
    int operator!=(const G4LundStringFragmentation &right) const;


private:
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py);      

   virtual void Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom, G4double AntiMass, G4double InitialMass); 
   virtual G4bool StopFragmenting(const G4FragmentingString  * const string);
   virtual G4bool IsFragmentable(const G4FragmentingString * const string);
   virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron, G4FragmentingString * string);
   virtual G4bool SplitLast(G4FragmentingString * string, 
		    G4KineticTrackVector * LeftVector,
		    G4KineticTrackVector * RightVector);
   void SetMinimalStringMass(const G4FragmentingString  * const string);		    
   void SetMinimalStringMass2(const G4double aValue);		    

private:
   G4double MinimalStringMass;
   G4double MinimalStringMass2;
   G4double WminLUND;		    
};

//**************************************************************************************************************
// Class G4LundStringFragmentation 
#endif


