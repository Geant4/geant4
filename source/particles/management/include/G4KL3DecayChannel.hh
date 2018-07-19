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
// $Id: G4KL3DecayChannel.hh 95906 2016-03-02 10:56:50Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
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
 
  protected:
    // Copy constructor and assignment operator
      G4KL3DecayChannel(const G4KL3DecayChannel &);
      G4KL3DecayChannel & operator=(const G4KL3DecayChannel &);

  private:
      G4KL3DecayChannel();

  public:  // With Description
     virtual G4DecayProducts *DecayIt(G4double);     

  protected:
     // assignment of daughter particles for arrays of daughters[] etc.
     enum{idPi=0, idLepton=1, idNutrino=2}; 
     //G4double daughterM[3];

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
     G4double   DalitzDensity(G4double parentmass, G4double Epi, G4double El,  G4double Enu,
                              G4double massPi, G4double massL , G4double massNu );
  private:
     // constants used in DalitzDensity() 
     //   coefficients
     G4double   pLambda;
     G4double   pXi0;

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




