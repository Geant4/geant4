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
// -------------------------------------------------------------------
// $Id: G4DNAFinalStatesPolicyScreenedRutherford.hh,v 1.1 2007-10-08 09:18:43 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNAFinalStatesPolicyScreenedRutherford_HH
#define G4DNAFinalStatesPolicyScreenedRutherford_HH 1

#include "G4DNACrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNAFinalStatesPolicyScreenedRutherford  
{
 protected:
   G4DNAFinalStatesPolicyScreenedRutherford() {}
   ~G4DNAFinalStatesPolicyScreenedRutherford() {}

   G4double EmfietzoglouRandomizeCosTheta(G4double k, G4int z);
   G4double BrennerRandomizeCosTheta(G4double k, G4int z);

  private:

   G4double EmfietzoglouScreeningFactor(G4double k, G4int z);
   G4double BrennerCalculatePolynomial(G4double k, const G4double *vector, G4int size);

   // Hides default constructor and assignment operator as private
   G4DNAFinalStatesPolicyScreenedRutherford(const G4DNAFinalStatesPolicyScreenedRutherford & copy);
   G4DNAFinalStatesPolicyScreenedRutherford & operator=(const G4DNAFinalStatesPolicyScreenedRutherford & right);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DNAFinalStatesPolicyScreenedRutherford.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif 
