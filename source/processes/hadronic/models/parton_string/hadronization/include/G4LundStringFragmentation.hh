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
// $Id: G4LundStringFragmentation.hh,v 1.7 2010-09-20 12:46:23 vuzhinsk Exp $
// GEANT4 tag $Name: not supported by cvs2svn $ Maxim Komogorov
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
    G4LundStringFragmentation(const G4LundStringFragmentation &right);
    virtual ~G4LundStringFragmentation();

    const G4LundStringFragmentation & operator=(const G4LundStringFragmentation &right);
    int operator==(const G4LundStringFragmentation &right) const;
    int operator!=(const G4LundStringFragmentation &right) const;

    virtual G4KineticTrackVector* FragmentString(const G4ExcitedString& theString);

private:
   void SetMinimalStringMass(const G4FragmentingString  * const string);		    
   void SetMinimalStringMass2(const G4double aValue);	

   virtual G4bool StopFragmenting(const G4FragmentingString  * const string);
   virtual G4bool IsFragmentable(const G4FragmentingString * const string);

   virtual G4bool SplitLast(G4FragmentingString * string, 
		            G4KineticTrackVector * LeftVector,
		            G4KineticTrackVector * RightVector);

   virtual void Sample4Momentum(G4LorentzVector* Mom,     G4double Mass, 
                                G4LorentzVector* AntiMom, G4double AntiMass, 
                                G4double InitialMass); 

   virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron, 
                                        G4FragmentingString * string,
                                        G4FragmentingString * newString); // Uzhi

   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, 
                                  G4int PartonEncoding,  
                                  G4ParticleDefinition* pHadron,
                                  G4double Px, G4double Py);      

   G4double lambda(G4double s, G4double m1_Sqr, G4double m2_Sqr);

private:
// ------ For estimation of a minimal string mass ---------------
   G4double Mass_of_light_quark;
   G4double Mass_of_heavy_quark;
   G4double Mass_of_string_junction;
// ------ An estimated minimal string mass ----------------------
   G4double MinimalStringMass;
   G4double MinimalStringMass2;
// ------ Minimal invariant mass used at a string fragmentation -
   G4double WminLUND;		    

   G4int          Meson[3][3][6];
   G4double MesonWeight[3][3][6];

   G4int          Baryon[3][3][3][4];
   G4double BaryonWeight[3][3][3][4];

   G4double Prob_QQbar[3];
};

//**************************************************************************************************************
// Class G4LundStringFragmentation 
#endif


