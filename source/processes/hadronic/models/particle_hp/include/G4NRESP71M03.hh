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
#ifndef G4NRESP71M03_HH
#define G4NRESP71M03_HH

//#include "G4PhysicalConstants.hh"
//#include "G4SystemOfUnits.hh"
//#include "Randomize.hh"

#include "G4ReactionProduct.hh"
#include "globals.hh"

class G4NRESP71M03 
{
   public:
     
      G4NRESP71M03(){;}; 
      ~G4NRESP71M03(){;}; 

      void DKINMA(G4ReactionProduct *p1, G4ReactionProduct *p2, G4ReactionProduct *p3, G4ReactionProduct *p4, const G4double Q, const G4double costhcm3);

      G4int ApplyMechanismI_NBeA2A(G4ReactionProduct &neut, G4ReactionProduct &carb, G4ReactionProduct *theProds, const G4double QI);
      G4int ApplyMechanismII_ACN2A(G4ReactionProduct &neut, G4ReactionProduct &carb, G4ReactionProduct *theProds, const G4double QI);

      G4int ApplyMechanismABE(G4ReactionProduct &neut, G4ReactionProduct &carb, G4ReactionProduct *theProds);

   private:
      // Defining the arrays with the angular distribution data 
      static const G4int ndist = 32; // Number of angular distributions.
      static const G4int nrhos = 51; // Number of Rho values.
      // Energies for which angular distributions are given.
      static const G4double BEN2[ndist];
      // Angular distribution of alpha particles for each energy value in BEN2.
      static const G4double B2[ndist][nrhos];

};

#endif
