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
//$Id: G4GeneralPhaseSpaceDecay.hh,v 1.1 1997/05/21
// ----------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 20th May 1998
//
//      Note: this class is a generalization of the 
//            G4PhaseSpaceDecayChannel one
// ----------------------------------------------------------------
#ifndef G4GeneralPhaseSpaceDecay_h
#define G4GeneralPhaseSpaceDecay_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"
#include "G4HadronicException.hh"

class G4GeneralPhaseSpaceDecay : public G4VDecayChannel
{
  public:
    //Constructors 
      G4GeneralPhaseSpaceDecay(G4int Verbose = 1);

      G4GeneralPhaseSpaceDecay(const G4String& theParentName,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 = "",
			       const G4String& theDaughterName3 = "");

      G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                               G4double        theParentMass,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 = "",
			       const G4String& theDaughterName3 = "");

      G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                               G4double        theParentMass,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 ,
			       const G4String& theDaughterName3 ,
			       const G4double * masses);

      G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                               G4double        theParentMass,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2 ,
			       const G4String& theDaughterName3 ,
			       const G4String& theDaughterName4 ,
			       const G4double * masses);

    //  Destructor
      virtual ~G4GeneralPhaseSpaceDecay();

  public:
     G4double GetParentMass() const;
     void SetParentMass(const G4double aParentMass);
     virtual G4DecayProducts* DecayIt(G4double mass=0.0);   
     static G4double Pmx(G4double e, G4double p1, G4double p2);

  protected:
     G4DecayProducts* OneBodyDecayIt();
     G4DecayProducts* TwoBodyDecayIt();
     G4DecayProducts* ThreeBodyDecayIt();
     G4DecayProducts* ManyBodyDecayIt();
     
  private:
     G4double parentmass;
     const G4double * theDaughterMasses;
    
};  



inline G4double G4GeneralPhaseSpaceDecay::GetParentMass() const
{
  return parentmass;
}

inline void G4GeneralPhaseSpaceDecay::SetParentMass(const G4double aParentMass)
{
  parentmass = aParentMass;
}



inline
 G4double G4GeneralPhaseSpaceDecay::Pmx(G4double e, G4double p1, G4double p2)
{
   // calculate momentum of daughter particles in two-body decay
   if (e-p1-p2 < 0 )
   {  
        throw G4HadronicException(__FILE__, __LINE__, "G4GeneralPhaseSpaceDecay::Pmx energy in cms > mass1+mass2");	
   }
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return std::sqrt(ppp);
   else       return -1.;
}

#endif
