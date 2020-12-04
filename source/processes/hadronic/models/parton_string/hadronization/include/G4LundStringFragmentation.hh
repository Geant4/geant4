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
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------

#ifndef G4LundStringFragmentation_h
#define G4LundStringFragmentation_h 1

#include "G4VLongitudinalStringDecay.hh"

//******************************************************************************

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
    G4bool operator==(const G4LundStringFragmentation &right) const;
    G4bool operator!=(const G4LundStringFragmentation &right) const;

  private:
    virtual G4bool StopFragmenting(const G4FragmentingString  * const string);
    virtual G4bool IsItFragmentable(const G4FragmentingString * const string);

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

    // virtual G4ParticleDefinition * QuarkSplitup(G4ParticleDefinition* decay, 
    //                                             G4ParticleDefinition *&created);

    virtual G4ParticleDefinition * DiQuarkSplitup(G4ParticleDefinition* decay, 
                                                  G4ParticleDefinition *&created);

  private:
    G4bool Loop_toFragmentString(const G4ExcitedString & theStringInCMS,
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

    G4double Tmt;  // "Temperature" for sampling Pt of hadrons

};

//******************************************************************************
// Class G4LundStringFragmentation 

#endif

