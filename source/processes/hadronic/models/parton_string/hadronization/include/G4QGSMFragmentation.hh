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
#ifndef G4QGSMFragmentation_h
#define G4QGSMFragmentation_h 1

#include "G4VLongitudinalStringDecay.hh"

//******************************************************************************
class G4QGSMFragmentation:public G4VLongitudinalStringDecay
{
  public:
    G4QGSMFragmentation();
    ~G4QGSMFragmentation();
    virtual G4KineticTrackVector* FragmentString(const G4ExcitedString& theString);

  private:
    // not implemented to protect/forbid use
    G4QGSMFragmentation(const G4QGSMFragmentation &right);
    const G4QGSMFragmentation & operator=(const G4QGSMFragmentation &right);
    G4bool operator==(const G4QGSMFragmentation &right) const;
    G4bool operator!=(const G4QGSMFragmentation &right) const;

  private:
    virtual G4bool StopFragmenting(const G4FragmentingString  * const string);
    virtual G4bool IsItFragmentable(const G4FragmentingString * const string);

    virtual G4bool SplitLast(G4FragmentingString * string, 
	 	             G4KineticTrackVector * LeftVector, 
                             G4KineticTrackVector * RightVector);

    virtual void Sample4Momentum(G4LorentzVector* Mom, G4double Mass, G4LorentzVector* AntiMom,
                                 G4double AntiMass, G4double InitialMass); 

    virtual G4KineticTrack * Splitup(G4FragmentingString *string, 
                                     G4FragmentingString *&newString);

    // The hadron can be producet at QuarkSplitup or DiQuarkSplitup
    // virtual G4ParticleDefinition * QuarkSplitup(G4ParticleDefinition* decay,
    //	  	   		                   G4ParticleDefinition *&created);

    virtual G4ParticleDefinition * DiQuarkSplitup(G4ParticleDefinition* decay,
                                                  G4ParticleDefinition *&created);

    virtual G4LorentzVector * SplitEandP(G4ParticleDefinition * pHadron, 
                                         G4FragmentingString * string, 
                                         G4FragmentingString * newString);

    virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding, 
                                   G4ParticleDefinition* pHadron, G4double Px, G4double Py);

  private:
    // model parameters
    G4double arho;
    G4double aphi;
    G4double aJPs;        // alpha_J/Psi
    G4double aUps;        // alpha_Y
    G4double an;
    G4double ala;

    G4double alaC;
    G4double alaB;
    G4double aXi;
    G4double aXiC;
    G4double aXiB;
    G4double aXiCC;
    G4double aXiCB;
    G4double aXiBB;

    G4double aksi;
    G4double alft;

    G4double FFq2q[5][5][2];
    G4double FFq2qq[5][15][2];
    G4double FFqq2q[15][5][2];
    G4double FFqq2qq[15][5][2];

  private:
    void SetFFq2q();
    void SetFFq2qq();
    void SetFFqq2qq();
    void SetFFqq2q();

    G4int IndexDiQ[5][5];
  };

// Class G4QGSMFragmentation 

#endif

