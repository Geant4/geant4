// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KL3DecayChannel.hh,v 1.2 1999-10-28 23:24:11 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//       7 June 1997 H.Kurashige
// ------------------------------------------------------------
#ifndef G4KL3DecayChannel_h
#define G4KL3DecayChannel_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4KL3DecayChannel :public G4VDecayChannel
{
  public:  // With Description
    //Constructors 
      G4KL3DecayChannel(const G4String& theParentName,
			G4double        theBR,
			const G4String& thePionName,
			const G4String& theLeptonName,
			const G4String& theNutrinoName);
    //  Destructor
      virtual ~G4KL3DecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);     

  protected:
     // assignment of daughter particles for arrays of daughters[] etc.
     enum{idPi=0, idLepton=1, idNutrino=2}; 
     G4double daughterM[3];

  protected:
     // calcurate momentum of daughters
     //     results will be stored in Edaughter[3] and Pdaughter[3]
     void PhaseSpace(G4double Mparent,
		     const G4double* Mdaughter,
		     G4double*       Edaughter,
		     G4double*       Pdaughter);

  protected:
     // Dalitz Plot Density
     // KL3 decay   Dalitz Plot Density
     //               see Chounet et al Phys. Rep. 4, 201
     //  arguments
     //    Epi: kinetic enregy of pion
     //    El:  kinetic enregy of lepton (e or mu)
     //    Enu: kinetic energy of nutrino
     //  constants
     //    pLambda : linear energy dependence of f+
     //    pXi0    : = f+(0)/f-
     //    pNorm   : normalization factor
     G4double   DalitzDensity(G4double Epi, G4double El,  G4double Enu);  
  private:
     // constants used in DalitzDensity() 
     //   Kon mass
     G4double massK;
     //   coefficients
     G4double   pLambda;
     G4double   pXi0;
     G4double   pNormalization;

  public:
     void SetDalitzParameter(G4double aLambda, G4double aXi );
     G4double GetDalitzParameterLambda() const;
     G4double GetDalitzParameterXi() const;
};  

inline 
 void G4KL3DecayChannel::SetDalitzParameter(G4double aLambda, G4double aXi)
{
   pLambda  = aLambda;
   pXi0      = aXi;
}

inline 
 G4double G4KL3DecayChannel::GetDalitzParameterLambda() const
{
  return  pLambda;
}

inline 
 G4double G4KL3DecayChannel::GetDalitzParameterXi() const
{
  return  pXi0;
}


#endif




