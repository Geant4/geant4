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
// $Id: G4LundStringFragmentation.hh 102717 2017-02-20 10:37:13Z gcosmo $
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
    virtual ~G4LundStringFragmentation();

    virtual G4KineticTrackVector* FragmentString(const G4ExcitedString& theString);

private:
    // not implemented to protect/forbid use
    G4LundStringFragmentation(const G4LundStringFragmentation &right);
    const G4LundStringFragmentation & operator=(const G4LundStringFragmentation &right);
    int operator==(const G4LundStringFragmentation &right) const;
    int operator!=(const G4LundStringFragmentation &right) const;

private:
   void SetMinMasses();   // Uzhi 23 Dec. 2016
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

   virtual G4KineticTrack * Splitup(G4FragmentingString *string,
                            G4FragmentingString *&newString);

   virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron, 
                                        G4FragmentingString * string,
                                        G4FragmentingString * newString);

   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, 
                                  G4int PartonEncoding,  
                                  G4ParticleDefinition* pHadron,
                                  G4double Px, G4double Py);      

   G4double lambda(G4double s, G4double m1_Sqr, G4double m2_Sqr);

   virtual G4ParticleDefinition * DiQuarkSplitup(G4ParticleDefinition* decay, 
                                         G4ParticleDefinition *&created);

private:
   // Internal methods introduced to improve the code structure (AR Nov 2011)

   G4bool Loop_toFragmentString(const G4ExcitedString & theStringInCMS, // * &
                                G4KineticTrackVector * & LeftVector, 
                                G4KineticTrackVector * & RightVector);

   G4bool Diquark_AntiDiquark_belowThreshold_lastSplitting(G4FragmentingString * & string,
                                                           G4ParticleDefinition * & LeftHadron,
                                                           G4ParticleDefinition * & RightHadron);

   G4bool Diquark_AntiDiquark_aboveThreshold_lastSplitting(G4FragmentingString * & string,
                                                           G4ParticleDefinition * & LeftHadron,
                                                           G4ParticleDefinition * & RightHadron);

   G4bool Quark_AntiQuark_lastSplitting(G4FragmentingString * & string,
                                        G4ParticleDefinition * & LeftHadron,
                                        G4ParticleDefinition * & RightHadron);

   G4bool Quark_Diquark_lastSplitting(G4FragmentingString * & string,
                                      G4ParticleDefinition * & LeftHadron,
                                      G4ParticleDefinition * & RightHadron );

   G4int SampleState(void); 

private:
// ------ For estimation of a minimal string mass ---------------
   G4double Mass_of_light_quark;
   G4double Mass_of_heavy_quark;
   G4double Mass_of_string_junction;

   G4double minMassQQbarStr[3][3];
   G4double minMassQDiQStr[3][3][3];

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

// ------ To improve the code structure
   G4ParticleDefinition * FS_LeftHadron[35], * FS_RightHadron[35];
   G4double FS_Weight[35];
   G4int NumberOf_FS;

};

//**************************************************************************************************************
// Class G4LundStringFragmentation 
#endif


