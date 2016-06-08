//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4DiffractiveStringBuilder.hh"

//***************************************************************************************************

G4DiffractiveStringBuilder::G4DiffractiveStringBuilder()
   {
   }

G4DiffractiveStringBuilder::G4DiffractiveStringBuilder(const G4DiffractiveStringBuilder &right)
   {
   }

G4DiffractiveStringBuilder::~G4DiffractiveStringBuilder()
   {
   }

//***************************************************************************************************

G4ExcitedString* G4DiffractiveStringBuilder::BuildString(G4PartonPair * aPair)       
    {
    return  new G4ExcitedString(aPair->GetParton1(), aPair->GetParton2(), aPair->GetDirection());
    }
