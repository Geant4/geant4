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
// G4KL3DecayChannel

// Author: H.Kurashige, 30 May 1997 
// --------------------------------------------------------------------
#ifndef G4KL3DecayChannel_hh
#define G4KL3DecayChannel_hh 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VDecayChannel.hh"

class G4KL3DecayChannel : public G4VDecayChannel
{
  public:

    G4KL3DecayChannel(const G4String& theParentName,
                            G4double  theBR,
                      const G4String& thePionName,
                      const G4String& theLeptonName,
                      const G4String& theNutrinoName);
    virtual ~G4KL3DecayChannel();
      // Constructor & destructor

    virtual G4DecayProducts* DecayIt(G4double);     

    inline void SetDalitzParameter(G4double aLambda, G4double aXi );
    inline G4double GetDalitzParameterLambda() const;
    inline G4double GetDalitzParameterXi() const;

  protected:

    G4KL3DecayChannel(const G4KL3DecayChannel&);
    G4KL3DecayChannel& operator=(const G4KL3DecayChannel&);
      // Copy constructor and assignment operator

    enum { idPi=0, idLepton=1, idNutrino=2 }; 
      // Assignment of daughter particles for arrays of daughters[] etc.

    void PhaseSpace(G4double Mparent,
                    const G4double* Mdaughter,
                    G4double*       Edaughter,
                    G4double*       Pdaughter);
      // Calculate momentum of daughters

    G4double DalitzDensity(G4double parentmass, G4double Epi, G4double El,
                           G4double Enu, G4double massPi, G4double massL,
                           G4double massNu );
      // Dalitz Plot Density
      // KL3 decay   Dalitz Plot Density, see Chounet et al Phys. Rep. 4, 201
      //  Arguments
      //    Epi: kinetic enregy of pion
      //    El:  kinetic enregy of lepton (e or mu)
      //    Enu: kinetic energy of nutrino
      //  Constants
      //    pLambda : linear energy dependence of f+
      //    pXi0    : = f+(0)/f-
      //    pNorm   : normalization factor

  private:

    G4KL3DecayChannel();

  private:

    G4double pLambda = 0.0;
    G4double pXi0 = 0.0;
      // Used in DalitzDensity() coefficients
};  

// ------------------------
// Inline methods
// ------------------------

inline 
void G4KL3DecayChannel::SetDalitzParameter(G4double aLambda, G4double aXi)
{
  pLambda = aLambda;
  pXi0    = aXi;
}

inline 
G4double G4KL3DecayChannel::GetDalitzParameterLambda() const
{
  return pLambda;
}

inline 
G4double G4KL3DecayChannel::GetDalitzParameterXi() const
{
  return pXi0;
}

#endif
