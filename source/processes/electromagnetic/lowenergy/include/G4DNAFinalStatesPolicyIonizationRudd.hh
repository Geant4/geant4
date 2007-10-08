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
// $Id: G4DNAFinalStatesPolicyIonizationRudd.hh,v 1.1 2007-10-08 09:18:43 sincerti Exp $
// -------------------------------------------------------------------
//

#ifndef G4DNAFinalStatesPolicyIonizationRudd_HH
#define G4DNAFinalStatesPolicyIonizationRudd_HH 1

#include "G4DNACrossSectionDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4DNAFinalStatesPolicyIonizationRudd  
{
 protected:
   G4DNAFinalStatesPolicyIonizationRudd() {}
   ~G4DNAFinalStatesPolicyIonizationRudd() {}

   G4double RandomizeEjectedElectronEnergy(const G4Track& track, G4double incomingParticleEnergy, G4int shell) ;
   void RandomizeEjectedElectronDirection(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4double
                                           outgoingParticleEnergy, G4double & cosTheta, G4double & phi );
   G4double EnergyConstant(G4int ionizationLevel);

  private:

   double DifferentialCrossSection(const G4Track& track, G4double k, G4double energyTransfer, G4int shell);
   G4double CorrectionFactor(G4ParticleDefinition * aParticleDefinition, G4double k);

   G4double S_1s(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber);

   G4double S_2s(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber);


   G4double S_2p(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber);

   G4double R(G4double t, G4double energyTransferred, G4double slaterEffectiveCharge, G4double shellNumber) ;

   G4double slaterEffectiveCharge[3];
   G4double sCoefficient[3];

   typedef std::map<double, std::map<double, double> > TriDimensionMap;
   TriDimensionMap DiffCrossSectionData[6];

   std::vector<double> TdummyVec;
 
   typedef std::map<double, std::vector<double> > VecMap;
   VecMap vecm;

   // Hides default constructor and assignment operator as private
   G4DNAFinalStatesPolicyIonizationRudd(const G4DNAFinalStatesPolicyIonizationRudd & copy);
   G4DNAFinalStatesPolicyIonizationRudd & operator=(const G4DNAFinalStatesPolicyIonizationRudd & right);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4DNAFinalStatesPolicyIonizationRudd.icc"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif 
