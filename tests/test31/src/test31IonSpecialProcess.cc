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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: test31IonSpecialProcess.cc,v 1.2 2003-06-19 14:46:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//------------ test31IonSpecialProcess physics process -----------------------------
//                   by Michel Maire, April 1996
//
// 28-05-96, DoIt() small change in ElecDirection, by M.Maire
// 10-06-96, simplification in ComputeMicroscopicCrossSection(), by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 13-09-96, small changes in DoIt for better efficiency. Thanks to P.Urban
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 05-03-97, new Physics scheme, M.Maire
// 28-03-97, protection in BuildPhysicsTable, M.Maire
// 07-04-98, remove 'tracking cut' of the scattered gamma, MMa
// 04-06-98, in DoIt, secondary production condition:
//                                     range>std::min(threshold,safety)
// 13-08-98, new methods SetBining()  PrintInfo()
// 15-12-98, cross section=0 below 10 keV
// 28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
// 13-07-01, DoIt: suppression of production cut for the electron (mma)
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 17-09-01, migration of Materials to pure STL (mma)
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 17-04-02, LowestEnergyLimit = 1*keV     
// -----------------------------------------------------------------------------

#include "test31IonSpecialProcess.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// constructor
 
test31IonSpecialProcess::test31IonSpecialProcess(const G4String& processName)
  : G4VDiscreteProcess (processName)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31IonSpecialProcess::~test31IonSpecialProcess()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31IonSpecialProcess::BuildPhysicsTable(const G4ParticleDefinition&)
// Build cross section and mean free path tables
{
  PrintInfoDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double test31IonSpecialProcess::GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition)
{
  *condition = NotForced;
   
   
   /*
   const G4ElementTable* theElementTable = G4Element::GetElementTable();
   G4double AtomicNumber;
   size_t J;

   for ( J=0 ; J < G4Element::GetNumberOfElements(); J++ )
       {
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestEnergyLimit,HighestEnergyLimit,
                                           NumbBinTable );
        AtomicNumber = (*theElementTable)[J]->GetZ();

        for ( G4int i = 0 ; i < NumbBinTable ; i++ )
           {
             LowEdgeEnergy = ptrVector->GetLowEdgeEnergy(i);
             Value = ComputeCrossSectionPerAtom(LowEdgeEnergy, AtomicNumber);
             ptrVector->PutValue(i,Value);
           }

        theCrossSectionTable->insertAt( J , ptrVector ) ;

      }

// Build mean free path table for the Compton Scattering process

   if (theMeanFreePathTable) {
       theMeanFreePathTable->clearAndDestroy(); delete theMeanFreePathTable;}

   theMeanFreePathTable= new G4PhysicsTable(G4Material::GetNumberOfMaterials());
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   G4Material* material;

   for ( J=0 ; J < G4Material::GetNumberOfMaterials(); J++ )
       {
        //create physics vector then fill it ....
        ptrVector = new G4PhysicsLogVector(LowestEnergyLimit,HighestEnergyLimit,
                                           NumbBinTable ) ;
        material = (*theMaterialTable)[J];

        for ( G4int i = 0 ; i < NumbBinTable ; i++ )
           {
             LowEdgeEnergy = ptrVector->GetLowEdgeEnergy( i ) ;
             Value = ComputeMeanFreePath( LowEdgeEnergy, material);
             ptrVector->PutValue( i , Value ) ;
           }

        theMeanFreePathTable->insertAt( J , ptrVector ) ;
      }
*/
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* test31IonSpecialProcess::PostStepDoIt(const G4Track& aTrack,
                                                         const G4Step&  aStep)

{
   aParticleChange.Initialize(aTrack);

   return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31IonSpecialProcess::PrintInfoDefinition()
{

  G4cout << G4endl << GetProcessName() << " is built  " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
