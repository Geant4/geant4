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
// 
// G4PAIxSection.cc -- class implementation file
//
// GEANT 4 class implementation file
//
// For information related to this code, please, contact
// the Geant4 Collaboration.
//
// R&D: Vladimir.Grichine@cern.ch
//
// History:
//
// 13.05.03 V. Grichine, bug fixed for maxEnergyTransfer > max interval energy
// 28.05.01 V.Ivanchenko minor changes to provide ANSI -wall compilation 
// 17.05.01 V. Grichine, low energy extension down to 10*keV of proton
// 20.11.98 adapted to a new Material/SandiaTable interface, mma 
// 11.06.97 V. Grichine, 1st version 
//

#include "G4PAIxSection.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4Poisson.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4SandiaTable.hh"

using namespace std;

/* ******************************************************************

// Init  array of Lorentz factors

const G4double G4PAIxSection::fLorentzFactor[22] =
{
          0.0 ,     1.1 ,   1.2 ,   1.3 ,    1.5 ,    1.8 ,  2.0 ,
          2.5 ,     3.0 ,   4.0 ,   7.0 ,   10.0 ,   20.0 , 40.0 ,
         70.0 ,   100.0 , 300.0 , 600.0 , 1000.0 , 3000.0 ,
      10000.0 , 50000.0
};

const G4int G4PAIxSection::
fRefGammaNumber = 29;         // The number of gamma for creation of 
                               // spline (9)

***************************************************************** */ 

// Local class constants

const G4double G4PAIxSection::fDelta = 0.005; // 0.005 energy shift from interval border
const G4double G4PAIxSection::fError = 0.005; // 0.005 error in lin-log approximation

const G4int G4PAIxSection::fMaxSplineSize = 1000;  // Max size of output spline
                                                    // arrays
//////////////////////////////////////////////////////////////////
//
// Constructor
//

G4PAIxSection::G4PAIxSection()
{
  fSandia = nullptr;
  fMatSandiaMatrix = nullptr;
  fDensity = fElectronDensity = fNormalizationCof = fLowEnergyCof = 0.0;
  fIntervalNumber = fSplineNumber = 0;
  fVerbose = 0;
    
  fSplineEnergy          = G4DataVector(fMaxSplineSize,0.0);
  fRePartDielectricConst = G4DataVector(fMaxSplineSize,0.0);
  fImPartDielectricConst = G4DataVector(fMaxSplineSize,0.0);
  fIntegralTerm          = G4DataVector(fMaxSplineSize,0.0);
  fDifPAIxSection        = G4DataVector(fMaxSplineSize,0.0);
  fdNdxCerenkov          = G4DataVector(fMaxSplineSize,0.0);
  fdNdxPlasmon           = G4DataVector(fMaxSplineSize,0.0);
  fdNdxMM                = G4DataVector(fMaxSplineSize,0.0);
  fdNdxResonance         = G4DataVector(fMaxSplineSize,0.0);
  fIntegralPAIxSection   = G4DataVector(fMaxSplineSize,0.0);
  fIntegralPAIdEdx       = G4DataVector(fMaxSplineSize,0.0);
  fIntegralCerenkov      = G4DataVector(fMaxSplineSize,0.0);
  fIntegralPlasmon       = G4DataVector(fMaxSplineSize,0.0);
  fIntegralMM            = G4DataVector(fMaxSplineSize,0.0);
  fIntegralResonance     = G4DataVector(fMaxSplineSize,0.0);

  fMaterialIndex = 0;   

  for( G4int i = 0; i < 500; ++i ) 
  {
    for( G4int j = 0; j < 112; ++j )  fPAItable[i][j] = 0.0; 
  }
}

//////////////////////////////////////////////////////////////////
//
// Constructor
//

G4PAIxSection::G4PAIxSection(G4MaterialCutsCouple* matCC)
{
  fDensity       = matCC->GetMaterial()->GetDensity();
  G4int matIndex = (G4int)matCC->GetMaterial()->GetIndex();
  fMaterialIndex = matIndex;   

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  fSandia = (*theMaterialTable)[matIndex]->GetSandiaTable();

  fVerbose = 0;

  G4int i, j; 
  fMatSandiaMatrix = new G4OrderedTable();
 
  for (i = 0; i < fSandia->GetMaxInterval()-1; ++i)
  {
     fMatSandiaMatrix->push_back(new G4DataVector(5,0.));
  }                         
  for (i = 0; i < fSandia->GetMaxInterval()-1; ++i)
  {
    (*(*fMatSandiaMatrix)[i])[0] = fSandia->GetSandiaMatTable(i,0);

    for(j = 1; j < 5; ++j)
    {
      (*(*fMatSandiaMatrix)[i])[j] = fSandia->GetSandiaMatTable(i,j)*fDensity;
    }     
  }
  ComputeLowEnergyCof();                               
  //  fEnergyInterval = fA1 = fA2 = fA3 = fA4 = 0;
}

////////////////////////////////////////////////////////////////

G4PAIxSection::G4PAIxSection(G4int materialIndex,
                             G4double maxEnergyTransfer)
{
  fSandia = nullptr;
  fMatSandiaMatrix = nullptr;
  fVerbose = 0;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int i, j;   

  fMaterialIndex   = materialIndex;   
  fDensity                = (*theMaterialTable)[materialIndex]->GetDensity();
  fElectronDensity        = (*theMaterialTable)[materialIndex]->
                             GetElectronDensity();
  fIntervalNumber         = (*theMaterialTable)[materialIndex]->
                             GetSandiaTable()->GetMatNbOfIntervals();
  fIntervalNumber--;      
  // G4cout<<fDensity<<"\t"<<fElectronDensity<<"\t"<<fIntervalNumber<<G4endl;

  fEnergyInterval = G4DataVector(fIntervalNumber+2,0.0);
  fA1             = G4DataVector(fIntervalNumber+2,0.0);
  fA2             = G4DataVector(fIntervalNumber+2,0.0);
  fA3             = G4DataVector(fIntervalNumber+2,0.0);
  fA4             = G4DataVector(fIntervalNumber+2,0.0);

  for(i = 1; i <= fIntervalNumber; i++ )
    {
      if(((*theMaterialTable)[materialIndex]->
    GetSandiaTable()->GetSandiaCofForMaterial(i-1,0) >= maxEnergyTransfer) ||
              i > fIntervalNumber               )
        {
          fEnergyInterval[i] = maxEnergyTransfer;
          fIntervalNumber = i;
          break;
        }
         fEnergyInterval[i] = (*theMaterialTable)[materialIndex]->
                              GetSandiaTable()->GetSandiaCofForMaterial(i-1,0);
         fA1[i]             = (*theMaterialTable)[materialIndex]->
                              GetSandiaTable()->GetSandiaCofForMaterial(i-1,1);
         fA2[i]             = (*theMaterialTable)[materialIndex]->
                              GetSandiaTable()->GetSandiaCofForMaterial(i-1,2);
         fA3[i]             = (*theMaterialTable)[materialIndex]->
                              GetSandiaTable()->GetSandiaCofForMaterial(i-1,3);
         fA4[i]             = (*theMaterialTable)[materialIndex]->
                              GetSandiaTable()->GetSandiaCofForMaterial(i-1,4);
         // G4cout<<i<<"\t"<<fEnergyInterval[i]<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
         //                               <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }   
  if(fEnergyInterval[fIntervalNumber] != maxEnergyTransfer)
    {
         fIntervalNumber++;
         fEnergyInterval[fIntervalNumber] = maxEnergyTransfer;
    }

  // Now checking, if two borders are too close together

  for(i=1;i<fIntervalNumber;i++)
    {
        if(fEnergyInterval[i+1]-fEnergyInterval[i] >
           1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]))
        {
          continue;
        }
        else
        {
          for(j=i;j<fIntervalNumber;j++)
          {
            fEnergyInterval[j] = fEnergyInterval[j+1];
                        fA1[j] = fA1[j+1];
                        fA2[j] = fA2[j+1];
                        fA3[j] = fA3[j+1];
                        fA4[j] = fA4[j+1];
          }
          fIntervalNumber--;
          i--;
        }
    }


      /* *********************************

      fSplineEnergy          = new G4double[fMaxSplineSize];   
      fRePartDielectricConst = new G4double[fMaxSplineSize];   
      fImPartDielectricConst = new G4double[fMaxSplineSize];   
      fIntegralTerm          = new G4double[fMaxSplineSize];   
      fDifPAIxSection        = new G4double[fMaxSplineSize];   
      fIntegralPAIxSection   = new G4double[fMaxSplineSize];   
      
      for(i=0;i<fMaxSplineSize;i++)
      {
         fSplineEnergy[i]          = 0.0;   
         fRePartDielectricConst[i] = 0.0;   
         fImPartDielectricConst[i] = 0.0;   
         fIntegralTerm[i]          = 0.0;   
         fDifPAIxSection[i]        = 0.0;   
         fIntegralPAIxSection[i]   = 0.0;   
      }
      **************************************************  */   
      ComputeLowEnergyCof();      
      InitPAI();  // create arrays allocated above
      /*     
      delete[] fEnergyInterval;
      delete[] fA1;
      delete[] fA2;
      delete[] fA3;
      delete[] fA4; 
      */   
}

////////////////////////////////////////////////////////////////////////
//
// Constructor called from G4PAIPhotonModel !!!

G4PAIxSection::G4PAIxSection( G4int materialIndex,
                              G4double maxEnergyTransfer,
                              G4double betaGammaSq,
                              G4double** photoAbsCof, 
                              G4int intNumber                   )
{
  fSandia = nullptr;
  fDensity = fElectronDensity = fNormalizationCof = fLowEnergyCof = 0.0;
  fIntervalNumber = fSplineNumber = 0;
  fVerbose = 0;
    
  fSplineEnergy          = G4DataVector(500,0.0);
  fRePartDielectricConst = G4DataVector(500,0.0);
  fImPartDielectricConst = G4DataVector(500,0.0);
  fIntegralTerm          = G4DataVector(500,0.0);
  fDifPAIxSection        = G4DataVector(500,0.0);
  fdNdxCerenkov          = G4DataVector(500,0.0);
  fdNdxPlasmon           = G4DataVector(500,0.0);
  fdNdxMM                = G4DataVector(500,0.0);
  fdNdxResonance         = G4DataVector(500,0.0);
  fIntegralPAIxSection   = G4DataVector(500,0.0);
  fIntegralPAIdEdx       = G4DataVector(500,0.0);
  fIntegralCerenkov      = G4DataVector(500,0.0);
  fIntegralPlasmon       = G4DataVector(500,0.0);
  fIntegralMM            = G4DataVector(500,0.0);
  fIntegralResonance     = G4DataVector(500,0.0);

  for( G4int i = 0; i < 500; ++i ) 
  {
    for( G4int j = 0; j < 112; ++j )  fPAItable[i][j] = 0.0; 
  }

  fSandia = nullptr;
  fMatSandiaMatrix = nullptr;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int i, j; 
  
  fMaterialIndex   = materialIndex;      
  fDensity         = (*theMaterialTable)[materialIndex]->GetDensity();
  fElectronDensity = (*theMaterialTable)[materialIndex]->GetElectronDensity();

  fIntervalNumber         = intNumber;
  fIntervalNumber--;
  //   G4cout<<fDensity<<"\t"<<fElectronDensity<<"\t"<<fIntervalNumber<<G4endl;
  
  fEnergyInterval = G4DataVector(fIntervalNumber+2,0.0);
  fA1             = G4DataVector(fIntervalNumber+2,0.0);
  fA2             = G4DataVector(fIntervalNumber+2,0.0);
  fA3             = G4DataVector(fIntervalNumber+2,0.0);
  fA4             = G4DataVector(fIntervalNumber+2,0.0);


  /*
      fEnergyInterval = new G4double[fIntervalNumber+2];
      fA1             = new G4double[fIntervalNumber+2];
      fA2             = new G4double[fIntervalNumber+2];
      fA3             = new G4double[fIntervalNumber+2];
      fA4             = new G4double[fIntervalNumber+2];
  */
  for( i = 1; i <= fIntervalNumber; i++ )
    {
         if( ( photoAbsCof[i-1][0] >= maxEnergyTransfer ) ||
             i > fIntervalNumber )
         {
            fEnergyInterval[i] = maxEnergyTransfer;
            fIntervalNumber = i;
            break;
         }
         fEnergyInterval[i] = photoAbsCof[i-1][0];
         fA1[i]             = photoAbsCof[i-1][1];
         fA2[i]             = photoAbsCof[i-1][2];
         fA3[i]             = photoAbsCof[i-1][3];
         fA4[i]             = photoAbsCof[i-1][4];
         // G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
         //      <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }
      // G4cout<<"i last = "<<i<<"; "<<"fIntervalNumber = "<<fIntervalNumber<<G4endl; 
  
  if(fEnergyInterval[fIntervalNumber] != maxEnergyTransfer)
    {
         fIntervalNumber++;
         fEnergyInterval[fIntervalNumber] = maxEnergyTransfer;
    }
      // G4cout<<"after check of max energy transfer"<<G4endl;

  for( i = 1; i <= fIntervalNumber; i++ )
    {
        // G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
        //   <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }
      // Now checking, if two borders are too close together

  for( i = 1; i < fIntervalNumber; i++ )
    {
        if(fEnergyInterval[i+1]-fEnergyInterval[i] >
           1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]))
        {
          continue;
        }
        else
        {
          for(j=i;j<fIntervalNumber;j++)
          {
            fEnergyInterval[j] = fEnergyInterval[j+1];
                        fA1[j] = fA1[j+1];
                        fA2[j] = fA2[j+1];
                        fA3[j] = fA3[j+1];
                        fA4[j] = fA4[j+1];
          }
          fIntervalNumber--;
          i--;
        }
    }
  // G4cout<<"after check of close borders"<<G4endl;

  for( i = 1; i <= fIntervalNumber; i++ )
    {
        // G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
        //  <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }

  // Preparation of fSplineEnergy array corresponding to min ionisation, G~4

  ComputeLowEnergyCof();            
  G4double   betaGammaSqRef = 
    fLorentzFactor[fRefGammaNumber]*fLorentzFactor[fRefGammaNumber] - 1;

  NormShift(betaGammaSqRef);             
  SplainPAI(betaGammaSqRef);
      
  // Preparation of integral PAI cross section for input betaGammaSq
   
  for(i = 1; i <= fSplineNumber; i++)
    {
         fdNdxCerenkov[i]   = PAIdNdxCerenkov(i,betaGammaSq);
         fdNdxMM[i]   = PAIdNdxMM(i,betaGammaSq);
         fdNdxPlasmon[i]    = PAIdNdxPlasmon(i,betaGammaSq);
         fdNdxResonance[i]  = PAIdNdxResonance(i,betaGammaSq);
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);

         // G4cout<<i<<"; dNdxC = "<<fdNdxCerenkov[i]<<"; dNdxP = "<<fdNdxPlasmon[i]
         //    <<"; dNdxPAI = "<<fDifPAIxSection[i]<<G4endl;
    }
  IntegralCerenkov();
  IntegralMM();
  IntegralPlasmon();
  IntegralResonance();
  IntegralPAIxSection();
  /*      
      delete[] fEnergyInterval;
      delete[] fA1;
      delete[] fA2;
      delete[] fA3;
      delete[] fA4;
  */    
}

////////////////////////////////////////////////////////////////////////
//
// Test Constructor with beta*gamma square value

G4PAIxSection::G4PAIxSection( G4int materialIndex,
                              G4double maxEnergyTransfer,
                              G4double betaGammaSq          )
{
  fSandia = nullptr;
  fMatSandiaMatrix = nullptr;
  fVerbose = 0;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4int i, j, numberOfElements;   

  fMaterialIndex   = materialIndex;   
  fDensity         = (*theMaterialTable)[materialIndex]->GetDensity();
  fElectronDensity = (*theMaterialTable)[materialIndex]->GetElectronDensity();
  numberOfElements = (G4int)(*theMaterialTable)[materialIndex]->GetNumberOfElements();

  G4int* thisMaterialZ = new G4int[numberOfElements];
   
  for( i = 0; i < numberOfElements; ++i )
   {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[materialIndex]->
                                      GetElement(i)->GetZ();
   }
  // fSandia = new G4SandiaTable(materialIndex);
  fSandia = (*theMaterialTable)[materialIndex]->GetSandiaTable();
  G4SandiaTable     thisMaterialSandiaTable(materialIndex);
  fIntervalNumber = thisMaterialSandiaTable.SandiaIntervals(thisMaterialZ,
                                                            numberOfElements);
  fIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                      (*theMaterialTable)[materialIndex]->GetFractionVector() ,
                             numberOfElements,fIntervalNumber);

  fIntervalNumber--;

  fEnergyInterval = G4DataVector(fIntervalNumber+2,0.0);
  fA1             = G4DataVector(fIntervalNumber+2,0.0);
  fA2             = G4DataVector(fIntervalNumber+2,0.0);
  fA3             = G4DataVector(fIntervalNumber+2,0.0);
  fA4             = G4DataVector(fIntervalNumber+2,0.0);

  /*
      fEnergyInterval = new G4double[fIntervalNumber+2];
      fA1             = new G4double[fIntervalNumber+2];
      fA2             = new G4double[fIntervalNumber+2];
      fA3             = new G4double[fIntervalNumber+2];
      fA4             = new G4double[fIntervalNumber+2];
  */
  for( i = 1; i <= fIntervalNumber; i++ )
    {
  if((thisMaterialSandiaTable.GetPhotoAbsorpCof(i,0) >= maxEnergyTransfer) ||
          i > fIntervalNumber)
         {
            fEnergyInterval[i] = maxEnergyTransfer;
            fIntervalNumber = i;
            break;
         }
   fEnergyInterval[i] = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,0);
   fA1[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,1)*fDensity;
   fA2[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,2)*fDensity;
   fA3[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,3)*fDensity;
   fA4[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,4)*fDensity;

    }   
  if(fEnergyInterval[fIntervalNumber] != maxEnergyTransfer)
    {
         fIntervalNumber++;
         fEnergyInterval[fIntervalNumber] = maxEnergyTransfer;
         fA1[fIntervalNumber] = fA1[fIntervalNumber-1];
         fA2[fIntervalNumber] = fA2[fIntervalNumber-1];
         fA3[fIntervalNumber] = fA3[fIntervalNumber-1];
         fA4[fIntervalNumber] = fA4[fIntervalNumber-1];
    }
  for(i=1;i<=fIntervalNumber;i++)
    {
        // G4cout<<fEnergyInterval[i]<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
        //   <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }
  // Now checking, if two borders are too close together

  for( i = 1; i < fIntervalNumber; i++ )
    {
        if(fEnergyInterval[i+1]-fEnergyInterval[i] >
           1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]))
        {
          continue;
        }
        else
        {
          for( j = i; j < fIntervalNumber; j++ )
          {
            fEnergyInterval[j] = fEnergyInterval[j+1];
                        fA1[j] = fA1[j+1];
                        fA2[j] = fA2[j+1];
                        fA3[j] = fA3[j+1];
                        fA4[j] = fA4[j+1];
          }
          fIntervalNumber--;
          i--;
        }
    }

      /* *********************************
      fSplineEnergy          = new G4double[fMaxSplineSize];   
      fRePartDielectricConst = new G4double[fMaxSplineSize];   
      fImPartDielectricConst = new G4double[fMaxSplineSize];   
      fIntegralTerm          = new G4double[fMaxSplineSize];   
      fDifPAIxSection        = new G4double[fMaxSplineSize];   
      fIntegralPAIxSection   = new G4double[fMaxSplineSize];   
      
      for(i=0;i<fMaxSplineSize;i++)
      {
         fSplineEnergy[i]          = 0.0;   
         fRePartDielectricConst[i] = 0.0;   
         fImPartDielectricConst[i] = 0.0;   
         fIntegralTerm[i]          = 0.0;   
         fDifPAIxSection[i]        = 0.0;   
         fIntegralPAIxSection[i]   = 0.0;   
      }
      */ ////////////////////////

      // Preparation of fSplineEnergy array corresponding to min ionisation, G~4

  ComputeLowEnergyCof();      
  G4double   betaGammaSqRef = 
    fLorentzFactor[fRefGammaNumber]*fLorentzFactor[fRefGammaNumber] - 1;

  NormShift(betaGammaSqRef);             
  SplainPAI(betaGammaSqRef);
      
  // Preparation of integral PAI cross section for input betaGammaSq
   
  for(i = 1; i <= fSplineNumber; i++)
    {
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
         fdNdxCerenkov[i]   = PAIdNdxCerenkov(i,betaGammaSq);
         fdNdxMM[i]   = PAIdNdxMM(i,betaGammaSq);
         fdNdxPlasmon[i]    = PAIdNdxPlasmon(i,betaGammaSq);
         fdNdxResonance[i]  = PAIdNdxResonance(i,betaGammaSq);
    }
  IntegralPAIxSection();
  IntegralCerenkov();
  IntegralMM();
  IntegralPlasmon();
  IntegralResonance();    
}

////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4PAIxSection::~G4PAIxSection()
{
   /* ************************
   delete[] fSplineEnergy         ;   
   delete[] fRePartDielectricConst;   
   delete[] fImPartDielectricConst;   
   delete[] fIntegralTerm         ;   
   delete[] fDifPAIxSection       ;   
   delete[] fIntegralPAIxSection  ;
   */ ////////////////////////
  delete fMatSandiaMatrix;
}

G4double G4PAIxSection::GetLorentzFactor(G4int j) const
{
   return fLorentzFactor[j];
}

////////////////////////////////////////////////////////////////////////
//
// Constructor with beta*gamma square value called from G4PAIPhotModel/Data

void G4PAIxSection::Initialize( const G4Material* material,
                                G4double maxEnergyTransfer,
                                G4double betaGammaSq, 
                                G4SandiaTable* sandia)
{
  if(fVerbose > 0)
  {
    G4cout<<G4endl;
    G4cout<<"G4PAIxSection::Initialize(...,G4SandiaTable* sandia)"<<G4endl;
    G4cout<<G4endl;
  }
  G4int i, j;

  fSandia          = sandia;
  fIntervalNumber  = sandia->GetMaxInterval();
  fDensity         = material->GetDensity();
  fElectronDensity = material->GetElectronDensity();

  // fIntervalNumber--;

  if( fVerbose > 0 )
  {
    G4cout<<"fDensity = "<<fDensity<<"\t"<<fElectronDensity<<"\t fIntervalNumber = "<<fIntervalNumber<<G4endl;
  }  
  fEnergyInterval = G4DataVector(fIntervalNumber+2,0.0);
  fA1             = G4DataVector(fIntervalNumber+2,0.0);
  fA2             = G4DataVector(fIntervalNumber+2,0.0);
  fA3             = G4DataVector(fIntervalNumber+2,0.0);
  fA4             = G4DataVector(fIntervalNumber+2,0.0);

  for( i = 1; i <= fIntervalNumber; i++ ) 
  {
    if ( sandia->GetSandiaMatTablePAI(i-1,0) < 1.*eV && sandia->GetLowerI1() == false) 
    { 
      fIntervalNumber--;
      continue;
    }
    if( ( sandia->GetSandiaMatTablePAI(i-1,0) >= maxEnergyTransfer ) || i >= fIntervalNumber ) 
    {
      fEnergyInterval[i] = maxEnergyTransfer;
      fIntervalNumber = i;
      break;
    }
    fEnergyInterval[i] = sandia->GetSandiaMatTablePAI(i-1,0);
    fA1[i]             = sandia->GetSandiaMatTablePAI(i-1,1);
    fA2[i]             = sandia->GetSandiaMatTablePAI(i-1,2);
    fA3[i]             = sandia->GetSandiaMatTablePAI(i-1,3);
    fA4[i]             = sandia->GetSandiaMatTablePAI(i-1,4);

      if( fVerbose > 0 ) 
      {
        G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
             <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
      }
  }
  if( fVerbose > 0 ) G4cout<<"last i = "<<i<<"; "<<"fIntervalNumber = "<<fIntervalNumber<<G4endl;   

  if( fEnergyInterval[fIntervalNumber] != maxEnergyTransfer )
  {
      fIntervalNumber++;
      fEnergyInterval[fIntervalNumber] = maxEnergyTransfer;
  }
  if( fVerbose > 0 )
  {  
    for( i = 1; i <= fIntervalNumber; i++ )
    {
      G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
        <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }
  }  
  if( fVerbose > 0 )    G4cout<<"Now checking, if two borders are too close together"<<G4endl;

  for( i = 1; i < fIntervalNumber; i++ )
  {
    if( fEnergyInterval[i+1]-fEnergyInterval[i] >
         1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]) && fEnergyInterval[i] > 0.) continue;
    else
    {
      if( fVerbose > 0 )  G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fEnergyInterval[i+1]/keV;

      for( j = i; j < fIntervalNumber; j++ )
      {
              fEnergyInterval[j] = fEnergyInterval[j+1];
              fA1[j]             = fA1[j+1];
              fA2[j]             = fA2[j+1];
              fA3[j]             = fA3[j+1];
              fA4[j]             = fA4[j+1];
      }
      fIntervalNumber--;
      i--;
    }
  }
  if( fVerbose > 0 )
  {
    for( i = 1; i <= fIntervalNumber; i++ )
    {
      G4cout<<i<<"\t"<<fEnergyInterval[i]/keV<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
        <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl;
    }
  }
  // Preparation of fSplineEnergy array corresponding to min ionisation, G~4

  ComputeLowEnergyCof(material);
            
  G4double   betaGammaSqRef = 
    fLorentzFactor[fRefGammaNumber]*fLorentzFactor[fRefGammaNumber] - 1;

  NormShift(betaGammaSqRef);             
  SplainPAI(betaGammaSqRef);
      
  // Preparation of integral PAI cross section for input betaGammaSq
   
  for( i = 1; i <= fSplineNumber; i++ )
  {
     fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);


     fdNdxCerenkov[i]   = PAIdNdxCerenkov(i,betaGammaSq);
     fdNdxMM[i]   = PAIdNdxMM(i,betaGammaSq);
     fdNdxPlasmon[i]    = PAIdNdxPlasmon(i,betaGammaSq);
     fdNdxResonance[i]  = PAIdNdxResonance(i,betaGammaSq);
  }
  IntegralPAIxSection();   
  IntegralCerenkov();
  IntegralMM();
  IntegralPlasmon();
  IntegralResonance();
   
  for( i = 1; i <= fSplineNumber; i++ )
  {
    if(fVerbose>0) G4cout<<i<<"; w = "<<fSplineEnergy[i]/keV<<" keV; dN/dx_>w = "<<fIntegralPAIxSection[i]<<" 1/mm"<<G4endl;
  }
}


/////////////////////////////////////////////////////////////////////////
//
// Compute low energy cof. It reduces PAI xsc for Lorentz factors less than 4.
//

void G4PAIxSection::ComputeLowEnergyCof(const G4Material* material)
{    
  G4int i, numberOfElements = (G4int)material->GetNumberOfElements();
  G4double sumZ = 0., sumCof = 0.; 

  static const G4double p0 =  1.20923e+00; 
  static const G4double p1 =  3.53256e-01; 
  static const G4double p2 = -1.45052e-03; 
  
  G4double* thisMaterialZ   = new G4double[numberOfElements];
  G4double* thisMaterialCof = new G4double[numberOfElements];
   
  for( i = 0; i < numberOfElements; ++i )
  {
    thisMaterialZ[i] = material->GetElement(i)->GetZ();
    sumZ += thisMaterialZ[i];
    thisMaterialCof[i] = p0+p1*thisMaterialZ[i]+p2*thisMaterialZ[i]*thisMaterialZ[i];   
  }
  for( i = 0; i < numberOfElements; ++i )
  {
    sumCof += thisMaterialCof[i]*thisMaterialZ[i]/sumZ;
  }
  fLowEnergyCof = sumCof;
  delete [] thisMaterialZ;
  delete [] thisMaterialCof;
  // G4cout<<"fLowEnergyCof = "<<fLowEnergyCof<<G4endl;
}

/////////////////////////////////////////////////////////////////////////
//
// Compute low energy cof. It reduces PAI xsc for Lorentz factors less than 4.
//

void G4PAIxSection::ComputeLowEnergyCof()
{    
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  G4int i, numberOfElements = (G4int)(*theMaterialTable)[fMaterialIndex]->GetNumberOfElements();
  G4double sumZ = 0., sumCof = 0.; 

  const G4double p0 =  1.20923e+00; 
  const G4double p1 =  3.53256e-01; 
  const G4double p2 = -1.45052e-03; 
  
  G4double* thisMaterialZ   = new G4double[numberOfElements];
  G4double* thisMaterialCof = new G4double[numberOfElements];
   
  for( i = 0; i < numberOfElements; ++i )
  {
    thisMaterialZ[i] = (*theMaterialTable)[fMaterialIndex]->GetElement(i)->GetZ();
    sumZ += thisMaterialZ[i];
    thisMaterialCof[i] = p0+p1*thisMaterialZ[i]+p2*thisMaterialZ[i]*thisMaterialZ[i];   
  }
  for( i = 0; i < numberOfElements; ++i )
  {
    sumCof += thisMaterialCof[i]*thisMaterialZ[i]/sumZ;
  }
  fLowEnergyCof = sumCof;
  // G4cout<<"fLowEnergyCof = "<<fLowEnergyCof<<G4endl;
  delete [] thisMaterialZ;
  delete [] thisMaterialCof;
}

/////////////////////////////////////////////////////////////////////////
//
// General control function for class G4PAIxSection
//

void G4PAIxSection::InitPAI()
{    
   G4int i;
   G4double betaGammaSq = fLorentzFactor[fRefGammaNumber]*
                          fLorentzFactor[fRefGammaNumber] - 1;

   // Preparation of integral PAI cross section for reference gamma
   
   NormShift(betaGammaSq);             
   SplainPAI(betaGammaSq);

   IntegralPAIxSection();
   IntegralCerenkov();
   IntegralMM();
   IntegralPlasmon();
   IntegralResonance();

   for(i = 0; i<= fSplineNumber; i++)
   {
      fPAItable[i][fRefGammaNumber] = fIntegralPAIxSection[i];
      if(i != 0) 
      {
         fPAItable[i][0] = fSplineEnergy[i];
      }
   }
   fPAItable[0][0] = fSplineNumber;
   
   for(G4int j = 1; j < 112; j++)       // for other gammas
   {
      if( j == fRefGammaNumber ) continue;
      
      betaGammaSq = fLorentzFactor[j]*fLorentzFactor[j] - 1;
      
      for(i = 1; i <= fSplineNumber; i++)
      {
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
         fdNdxCerenkov[i]   = PAIdNdxCerenkov(i,betaGammaSq);
         fdNdxMM[i]   = PAIdNdxMM(i,betaGammaSq);
         fdNdxPlasmon[i]    = PAIdNdxPlasmon(i,betaGammaSq);
         fdNdxResonance[i]  = PAIdNdxResonance(i,betaGammaSq);
      }
      IntegralPAIxSection();
      IntegralCerenkov();
      IntegralMM();
      IntegralPlasmon();
      IntegralResonance();
      
      for(i = 0; i <= fSplineNumber; i++)
      {
         fPAItable[i][j] = fIntegralPAIxSection[i];
      }
   } 

}  

///////////////////////////////////////////////////////////////////////
//
// Shifting from borders to intervals Creation of first energy points
//

void G4PAIxSection::NormShift(G4double betaGammaSq)
{
  G4int i, j;

  if(fVerbose>0) G4cout<<"      G4PAIxSection::NormShift call "<<G4endl;


  for( i = 1; i <= fIntervalNumber-1; i++ )
  {
    for( j = 1; j <= 2; j++ )
    {
      fSplineNumber = (i-1)*2 + j;

      if( j == 1 ) fSplineEnergy[fSplineNumber] = fEnergyInterval[i  ]*(1+fDelta);
      else         fSplineEnergy[fSplineNumber] = fEnergyInterval[i+1]*(1-fDelta); 
      if(fVerbose>0) G4cout<<"cn = "<<fSplineNumber<<"; "<<"w = "<<fSplineEnergy[fSplineNumber]/keV<<" keV"<<G4endl;
    }
  }
  fIntegralTerm[1]=RutherfordIntegral(1,fEnergyInterval[1],fSplineEnergy[1]);

  j = 1;

  for( i = 2; i <= fSplineNumber; i++ )
  {
    if( fSplineEnergy[i]<fEnergyInterval[j+1] )
    {
         fIntegralTerm[i] = fIntegralTerm[i-1] + 
                            RutherfordIntegral(j,fSplineEnergy[i-1],
                                                 fSplineEnergy[i]   );
    }
    else
    {
       G4double x = RutherfordIntegral(j,fSplineEnergy[i-1],
                                           fEnergyInterval[j+1]   );
         j++;
         fIntegralTerm[i] = fIntegralTerm[i-1] + x + 
                            RutherfordIntegral(j,fEnergyInterval[j],
                                                 fSplineEnergy[i]    );
    }
   if(fVerbose>0)  G4cout<<i<<"  Shift: w = "<<fSplineEnergy[i]/keV<<" keV \t"<<fIntegralTerm[i]<<"\n"<<G4endl;
  } 
  fNormalizationCof = 2*pi*pi*hbarc*hbarc*fine_structure_const/electron_mass_c2;
  fNormalizationCof *= fElectronDensity/fIntegralTerm[fSplineNumber];

  // G4cout<<"fNormalizationCof = "<<fNormalizationCof<<G4endl;

          // Calculation of PAI differrential cross-section (1/(keV*cm))
          // in the energy points near borders of energy intervals

   for(G4int k = 1; k <= fIntervalNumber-1; k++ )
   {
      for( j = 1; j <= 2; j++ )
      {
         i = (k-1)*2 + j;
         fImPartDielectricConst[i] = fNormalizationCof*
                                     ImPartDielectricConst(k,fSplineEnergy[i]);
         fRePartDielectricConst[i] = fNormalizationCof*
                                     RePartDielectricConst(fSplineEnergy[i]);
         fIntegralTerm[i] *= fNormalizationCof;

         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
         fdNdxCerenkov[i]   = PAIdNdxCerenkov(i,betaGammaSq);
         fdNdxMM[i]   = PAIdNdxMM(i,betaGammaSq);
         fdNdxPlasmon[i]    = PAIdNdxPlasmon(i,betaGammaSq);
         fdNdxResonance[i]    = PAIdNdxResonance(i,betaGammaSq);
   if(fVerbose>0)  G4cout<<i<<"  Shift: w = "<<fSplineEnergy[i]/keV<<" keV, xsc = "<<fDifPAIxSection[i]<<"\n"<<G4endl;
      }
   }

}  // end of NormShift 

/////////////////////////////////////////////////////////////////////////
//
// Creation of new energy points as geometrical mean of existing
// one, calculation PAI_cs for them, while the error of logarithmic
// linear approximation would be smaller than 'fError'

void G4PAIxSection::SplainPAI(G4double betaGammaSq)
{
  G4int j, k = 1, i = 1;

  if(fVerbose>0) G4cout<<"                   G4PAIxSection::SplainPAI call "<<G4endl;

  while ( (i < fSplineNumber) && (fSplineNumber < fMaxSplineSize-1) )
  {
     // if( std::abs(fSplineEnergy[i+1]-fEnergyInterval[k+1]) > (fSplineEnergy[i+1]+fEnergyInterval[k+1])*5.e-7 )
     if( fSplineEnergy[i+1] > fEnergyInterval[k+1] )
     {
          k++;   // Here next energy point is in next energy interval
          i++;
          if(fVerbose>0) G4cout<<"                     in if: i = "<<i<<"; k = "<<k<<G4endl;
          continue;
     }
     if(fVerbose>0) G4cout<<"       out if: i = "<<i<<"; k = "<<k<<G4endl;

                        // Shifting of arrayes for inserting the geometrical 
                       // average of 'i' and 'i+1' energy points to 'i+1' place
     fSplineNumber++;

     for( j = fSplineNumber; j >= i+2; j-- )
     {
         fSplineEnergy[j]          = fSplineEnergy[j-1];
         fImPartDielectricConst[j] = fImPartDielectricConst[j-1];
         fRePartDielectricConst[j] = fRePartDielectricConst[j-1];
         fIntegralTerm[j]          = fIntegralTerm[j-1];

         fDifPAIxSection[j] = fDifPAIxSection[j-1];
         fdNdxCerenkov[j]   = fdNdxCerenkov[j-1];
         fdNdxMM[j]         = fdNdxMM[j-1];
         fdNdxPlasmon[j]    = fdNdxPlasmon[j-1];
         fdNdxResonance[j]  = fdNdxResonance[j-1];
     }
      G4double x1  = fSplineEnergy[i];
      G4double x2  = fSplineEnergy[i+1];
      G4double yy1 = fDifPAIxSection[i];
      G4double y2  = fDifPAIxSection[i+1];

      if(fVerbose>0) G4cout<<"Spline: x1 = "<<x1<<"; x2 = "<<x2<<", yy1 = "<<yy1<<"; y2 = "<<y2<<G4endl;


      G4double en1 = sqrt(x1*x2);
      // G4double    en1 = 0.5*(x1 + x2);


      fSplineEnergy[i+1] = en1;

                 // Calculation of logarithmic linear approximation
                 // in this (enr) energy point, which number is 'i+1' now

      G4double a = log10(y2/yy1)/log10(x2/x1);
      G4double b = log10(yy1) - a*log10(x1);
      G4double y = a*log10(en1) + b;

      y = pow(10.,y);

                 // Calculation of the PAI dif. cross-section at this point

      fImPartDielectricConst[i+1] = fNormalizationCof*
                                    ImPartDielectricConst(k,fSplineEnergy[i+1]);
      fRePartDielectricConst[i+1] = fNormalizationCof*
                                    RePartDielectricConst(fSplineEnergy[i+1]);
      fIntegralTerm[i+1] = fIntegralTerm[i] + fNormalizationCof*
                           RutherfordIntegral(k,fSplineEnergy[i],
                                                fSplineEnergy[i+1]);

      fDifPAIxSection[i+1] = DifPAIxSection(i+1,betaGammaSq);
      fdNdxCerenkov[i+1]   = PAIdNdxCerenkov(i+1,betaGammaSq);
      fdNdxMM[i+1]         = PAIdNdxMM(i+1,betaGammaSq);
      fdNdxPlasmon[i+1]    = PAIdNdxPlasmon(i+1,betaGammaSq);
      fdNdxResonance[i+1]  = PAIdNdxResonance(i+1,betaGammaSq);

                  // Condition for next division of this segment or to pass

    if(fVerbose>0) G4cout<<"Spline, a = "<<a<<"; b = "<<b<<"; new xsc = "<<y<<"; compxsc = "<<fDifPAIxSection[i+1]<<G4endl;

                  // to higher energies

      G4double x = 2*(fDifPAIxSection[i+1] - y)/(fDifPAIxSection[i+1] + y);

      G4double delta = 2.*(fSplineEnergy[i+1]-fSplineEnergy[i])/(fSplineEnergy[i+1]+fSplineEnergy[i]);

      if( x < 0 ) 
      {
         x = -x;
      }
      if( x > fError && fSplineNumber < fMaxSplineSize-1 && delta > 2.*fDelta )
      {
         continue;  // next division
      }
      i += 2;  // pass to next segment

      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  }   // close 'while'

}  // end of SplainPAI 


////////////////////////////////////////////////////////////////////
//
// Integration over electrons that could be considered
// quasi-free at energy transfer of interest

G4double G4PAIxSection::RutherfordIntegral( G4int k,
                                            G4double x1,
                                              G4double x2   )
{
   G4double  c1, c2, c3;
   // G4cout<<"RI: x1 = "<<x1<<"; "<<"x2 = "<<x2<<G4endl;   
   c1 = (x2 - x1)/x1/x2;
   c2 = (x2 - x1)*(x2 + x1)/x1/x1/x2/x2;
   c3 = (x2 - x1)*(x1*x1 + x1*x2 + x2*x2)/x1/x1/x1/x2/x2/x2;
   // G4cout<<" RI: c1 = "<<c1<<"; "<<"c2 = "<<c2<<"; "<<"c3 = "<<c3<<G4endl;   
   
   return  fA1[k]*log(x2/x1) + fA2[k]*c1 + fA3[k]*c2/2 + fA4[k]*c3/3;

}   // end of RutherfordIntegral 


/////////////////////////////////////////////////////////////////
//
// Imaginary part of dielectric constant
// (G4int k - interval number, G4double en1 - energy point)

G4double G4PAIxSection::ImPartDielectricConst( G4int    k ,
                                               G4double energy1 )
{
   G4double energy2,energy3,energy4,result;

   energy2 = energy1*energy1;
   energy3 = energy2*energy1;
   energy4 = energy3*energy1;
   
   result = fA1[k]/energy1+fA2[k]/energy2+fA3[k]/energy3+fA4[k]/energy4;  
   result *=hbarc/energy1;
   
   return result;

}  // end of ImPartDielectricConst 

/////////////////////////////////////////////////////////////////
//
// Returns lambda of photon with energy1 in current material 

G4double G4PAIxSection::GetPhotonRange( G4double energy1 )
{
  G4int i;
  G4double energy2, energy3, energy4, result, lambda;

  energy2 = energy1*energy1;
  energy3 = energy2*energy1;
  energy4 = energy3*energy1;

  // G4double* SandiaCof = fSandia->GetSandiaCofForMaterialPAI(energy1);
  // result = SandiaCof[0]/energy1+SandiaCof[1]/energy2+SandiaCof[2]/energy3+SandiaCof[3]/energy4;
  // result *= fDensity;

  for( i = 1; i <= fIntervalNumber; i++ )
  {
     if( energy1 < fEnergyInterval[i]) break;
  }
  i--;
  if(i == 0) i = 1;

  result = fA1[i]/energy1+fA2[i]/energy2+fA3[i]/energy3+fA4[i]/energy4;  

  if( result > DBL_MIN ) lambda = 1./result;
  else                   lambda = DBL_MAX;
   
  return lambda;
}  

/////////////////////////////////////////////////////////////////
//
// Return lambda of electron with energy1 in current material
// parametrisation from NIM A554(2005)474-493 

G4double G4PAIxSection::GetElectronRange( G4double energy )
{
  G4double range;
  /*
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  G4double Z = (*theMaterialTable)[fMaterialIndex]->GetIonisation()->GetZeffective();
  G4double A = (*theMaterialTable)[fMaterialIndex]->GetA();

  energy /= keV; // energy in keV in parametrised formula

  if (energy < 10.)
  {
    range = 3.872e-3*A/Z;
    range *= pow(energy,1.492);
  }
  else
  {
    range = 6.97e-3*pow(energy,1.6);
  }
  */
  // Blum&Rolandi Particle Detection with Drift Chambers, p. 7

  G4double cofA = 5.37e-4*g/cm2/keV;
  G4double cofB = 0.9815;
  G4double cofC = 3.123e-3/keV;
  // energy /= keV;

  range = cofA*energy*( 1 - cofB/(1 + cofC*energy) ); 

  // range *= g/cm2;
  range /= fDensity;

  return range;
}

//////////////////////////////////////////////////////////////////////////////
//
// Real part of dielectric constant minus unit: epsilon_1 - 1
// (G4double enb - energy point)
//

G4double G4PAIxSection::RePartDielectricConst(G4double enb)
{       
   G4double x0, x02, x03, x04, x05, x1, x2, xx1 ,xx2 , xx12,
            c1, c2, c3, cof1, cof2, xln1, xln2, xln3, result;

   x0 = enb;
   result = 0;
   
   for(G4int i=1;i<=fIntervalNumber-1;i++)
   {
      x1 = fEnergyInterval[i];
      x2 = fEnergyInterval[i+1];
      xx1 = x1 - x0;
      xx2 = x2 - x0;
      xx12 = xx2/xx1;
      
      if(xx12<0)
      {
         xx12 = -xx12;
      }
      xln1 = log(x2/x1);
      xln2 = log(xx12);
      xln3 = log((x2 + x0)/(x1 + x0));
      x02 = x0*x0;
      x03 = x02*x0;
      x04 = x03*x0;
      x05 = x04*x0;
      c1  = (x2 - x1)/x1/x2;
      c2  = (x2 - x1)*(x2 +x1)/x1/x1/x2/x2;
      c3  = (x2 -x1)*(x1*x1 + x1*x2 + x2*x2)/x1/x1/x1/x2/x2/x2;

      result -= (fA1[i]/x02 + fA3[i]/x04)*xln1;
      result -= (fA2[i]/x02 + fA4[i]/x04)*c1;
      result -= fA3[i]*c2/2/x02;
      result -= fA4[i]*c3/3/x02;

      cof1 = fA1[i]/x02 + fA3[i]/x04;
      cof2 = fA2[i]/x03 + fA4[i]/x05;

      result += 0.5*(cof1 +cof2)*xln2;
      result += 0.5*(cof1 - cof2)*xln3;
   } 
   result *= 2*hbarc/pi;
   
   return result;

}   // end of RePartDielectricConst 

//////////////////////////////////////////////////////////////////////
//
// PAI differential cross-section in terms of
// simplified Allison's equation
//

G4double G4PAIxSection::DifPAIxSection( G4int              i ,
                                        G4double betaGammaSq  )
{        
   G4double cof,x1,x2,x3,x4,x5,x6,x7,x8,result;

   G4double betaBohr  = fine_structure_const;
   // G4double betaBohr2 = fine_structure_const*fine_structure_const;
   // G4double betaBohr3 = betaBohr*betaBohr2; // *4.0;

   G4double be2  = betaGammaSq/(1 + betaGammaSq);
   G4double beta = sqrt(be2);
   // G4double be3 = beta*be2;

   cof = 1.;
   x1  = log(2*electron_mass_c2/fSplineEnergy[i]);

   if( betaGammaSq < 0.01 ) x2 = log(be2);
   else
   {
     x2 = -log( (1/betaGammaSq - fRePartDielectricConst[i])*
                (1/betaGammaSq - fRePartDielectricConst[i]) + 
                fImPartDielectricConst[i]*fImPartDielectricConst[i] )/2;
   }
   if( fImPartDielectricConst[i] == 0.0 ||betaGammaSq < 0.01 )
   {
     x6 = 0.;
   }
   else
   {
     x3 = -fRePartDielectricConst[i] + 1/betaGammaSq;
     x5 = -1 - fRePartDielectricConst[i] +
          be2*((1 +fRePartDielectricConst[i])*(1 + fRePartDielectricConst[i]) +
          fImPartDielectricConst[i]*fImPartDielectricConst[i]);

     x7 = atan2(fImPartDielectricConst[i],x3);
     x6 = x5 * x7;
   }
    // if(fImPartDielectricConst[i] == 0) x6 = 0.;
   
   x4 = ((x1 + x2)*fImPartDielectricConst[i] + x6)/hbarc;

   //   if( x4 < 0.0 ) x4 = 0.0;

   x8 = (1 + fRePartDielectricConst[i])*(1 + fRePartDielectricConst[i]) + 
        fImPartDielectricConst[i]*fImPartDielectricConst[i];

   result = (x4 + cof*fIntegralTerm[i]/fSplineEnergy[i]/fSplineEnergy[i]);

   if( result < 1.0e-8 ) result = 1.0e-8;

   result *= fine_structure_const/be2/pi;

   // low energy correction

   G4double lowCof = fLowEnergyCof; // 6.0 ; // Ar ~ 4.; -> fLowCof as f(Z1,Z2)? 

   result *= (1 - exp(-beta/betaBohr/lowCof));


   // result *= (1 - exp(-be2/betaBohr2/lowCof));

   // result *= (1 - exp(-be3/betaBohr3/lowCof)); // ~ be for be<<betaBohr

   // result *= (1 - exp(-be4/betaBohr4/lowCof));

   if(fDensity >= 0.1)
   { 
      result /= x8;
   }
   return result;

} // end of DifPAIxSection 

//////////////////////////////////////////////////////////////////////////
//
// Calculation od dN/dx of collisions with creation of Cerenkov pseudo-photons

G4double G4PAIxSection::PAIdNdxCerenkov( G4int    i ,
                                         G4double betaGammaSq  )
{        
   G4double logarithm, x3, x5, argument, modul2, dNdxC; 
   G4double be2, betaBohr2, cofBetaBohr;

   cofBetaBohr = 4.0;
   betaBohr2   = fine_structure_const*fine_structure_const;
   G4double betaBohr4   = betaBohr2*betaBohr2*cofBetaBohr;

   be2 = betaGammaSq/(1 + betaGammaSq);
   G4double be4 = be2*be2;

   if( betaGammaSq < 0.01 ) logarithm = log(1.0+betaGammaSq); // 0.0;
   else
   {
     logarithm  = -log( (1/betaGammaSq - fRePartDielectricConst[i])*
                        (1/betaGammaSq - fRePartDielectricConst[i]) + 
                        fImPartDielectricConst[i]*fImPartDielectricConst[i] )*0.5;
     logarithm += log(1+1.0/betaGammaSq);
   }

   if( fImPartDielectricConst[i] == 0.0 || betaGammaSq < 0.01 )
   {
     argument = 0.0;
   }
   else
   {
     x3 = -fRePartDielectricConst[i] + 1.0/betaGammaSq;
     x5 = -1.0 - fRePartDielectricConst[i] +
          be2*((1.0 +fRePartDielectricConst[i])*(1.0 + fRePartDielectricConst[i]) +
          fImPartDielectricConst[i]*fImPartDielectricConst[i]);
     if( x3 == 0.0 ) argument = 0.5*pi;
     else            argument = atan2(fImPartDielectricConst[i],x3);
     argument *= x5 ;
   }   
   dNdxC = ( logarithm*fImPartDielectricConst[i] + argument )/hbarc;
  
   if(dNdxC < 1.0e-8) dNdxC = 1.0e-8;

   dNdxC *= fine_structure_const/be2/pi;

   dNdxC *= (1-exp(-be4/betaBohr4));

   if(fDensity >= 0.1)
   { 
      modul2 = (1.0 + fRePartDielectricConst[i])*(1.0 + fRePartDielectricConst[i]) + 
                    fImPartDielectricConst[i]*fImPartDielectricConst[i];
      dNdxC /= modul2;
   }
   return dNdxC;

} // end of PAIdNdxCerenkov 

//////////////////////////////////////////////////////////////////////////
//
// Calculation od dN/dx of collisions of MM with creation of Cerenkov pseudo-photons

G4double G4PAIxSection::PAIdNdxMM( G4int    i ,
                                         G4double betaGammaSq  )
{        
   G4double logarithm, x3, x5, argument, dNdxC; 
   G4double be2, be4, betaBohr2,betaBohr4,cofBetaBohr;

   cofBetaBohr = 4.0;
   betaBohr2   = fine_structure_const*fine_structure_const;
   betaBohr4   = betaBohr2*betaBohr2*cofBetaBohr;

   be2 = betaGammaSq/(1 + betaGammaSq);
   be4 = be2*be2;

   if( betaGammaSq < 0.01 ) logarithm = log(1.0+betaGammaSq); // 0.0;
   else
   {
     logarithm  = -log( (1/betaGammaSq - fRePartDielectricConst[i])*
                        (1/betaGammaSq - fRePartDielectricConst[i]) + 
                        fImPartDielectricConst[i]*fImPartDielectricConst[i] )*0.5;
     logarithm += log(1+1.0/betaGammaSq);
   }

   if( fImPartDielectricConst[i] == 0.0 || betaGammaSq < 0.01 )
   {
     argument = 0.0;
   }
   else
   {
     x3 = -fRePartDielectricConst[i] + 1.0/betaGammaSq;
     x5 = be2*( 1.0 + fRePartDielectricConst[i] ) - 1.0;
     if( x3 == 0.0 ) argument = 0.5*pi;
     else            argument = atan2(fImPartDielectricConst[i],x3);
     argument *= x5 ;
   }   
   dNdxC = ( logarithm*fImPartDielectricConst[i]*be2 + argument )/hbarc;
  
   if(dNdxC < 1.0e-8) dNdxC = 1.0e-8;

   dNdxC *= fine_structure_const/be2/pi;

   dNdxC *= (1-exp(-be4/betaBohr4));
   return dNdxC;

} // end of PAIdNdxMM 

//////////////////////////////////////////////////////////////////////////
//
// Calculation od dN/dx of collisions with creation of longitudinal EM
// excitations (plasmons, delta-electrons)

G4double G4PAIxSection::PAIdNdxPlasmon( G4int    i ,
                                        G4double betaGammaSq  )
{        
   G4double resonance, modul2, dNdxP, cof = 1.;
   G4double be2, betaBohr;
  
   betaBohr   = fine_structure_const;
   be2 = betaGammaSq/(1 + betaGammaSq);

   G4double beta = sqrt(be2);
 
   resonance = log(2*electron_mass_c2*be2/fSplineEnergy[i]);  
   resonance *= fImPartDielectricConst[i]/hbarc;


   dNdxP = ( resonance + cof*fIntegralTerm[i]/fSplineEnergy[i]/fSplineEnergy[i] );

   if( dNdxP < 1.0e-8 ) dNdxP = 1.0e-8;

   dNdxP *= fine_structure_const/be2/pi;

   dNdxP  *= (1 - exp(-beta/betaBohr/fLowEnergyCof));

   // dNdxP *= (1-exp(-be4/betaBohr4));

   if( fDensity >= 0.1 )
   { 
     modul2 = (1 + fRePartDielectricConst[i])*(1 + fRePartDielectricConst[i]) + 
        fImPartDielectricConst[i]*fImPartDielectricConst[i];
     dNdxP /= modul2;
   }
   return dNdxP;

} // end of PAIdNdxPlasmon 

//////////////////////////////////////////////////////////////////////////
//
// Calculation od dN/dx of collisions with creation of longitudinal EM
// resonance excitations (plasmons, delta-electrons)

G4double G4PAIxSection::PAIdNdxResonance( G4int    i ,
                                        G4double betaGammaSq  )
{        
   G4double resonance, modul2, dNdxP;
   G4double be2, be4, betaBohr2, betaBohr4, cofBetaBohr;

   cofBetaBohr = 4.0;
   betaBohr2   = fine_structure_const*fine_structure_const;
   betaBohr4   = betaBohr2*betaBohr2*cofBetaBohr;

   be2 = betaGammaSq/(1 + betaGammaSq);
   be4 = be2*be2;
 
   resonance = log(2*electron_mass_c2*be2/fSplineEnergy[i]);  
   resonance *= fImPartDielectricConst[i]/hbarc;


   dNdxP = resonance;

   if( dNdxP < 1.0e-8 ) dNdxP = 1.0e-8;

   dNdxP *= fine_structure_const/be2/pi;
   dNdxP *= (1-exp(-be4/betaBohr4));

   if( fDensity >= 0.1 )
   { 
     modul2 = (1 + fRePartDielectricConst[i])*(1 + fRePartDielectricConst[i]) + 
        fImPartDielectricConst[i]*fImPartDielectricConst[i];
     dNdxP /= modul2;
   }
   return dNdxP;

} // end of PAIdNdxResonance 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI integral cross-section
// fIntegralPAIxSection[1] = specific primary ionisation, 1/cm
// and fIntegralPAIxSection[0] = mean energy loss per cm  in keV/cm

void G4PAIxSection::IntegralPAIxSection()
{
  fIntegralPAIxSection[fSplineNumber] = 0;
  fIntegralPAIdEdx[fSplineNumber]     = 0;
  fIntegralPAIxSection[0]             = 0;
  G4int i, k = fIntervalNumber -1;

  for( i = fSplineNumber-1; i >= 1; i--)
  {
    if(fSplineEnergy[i] >= fEnergyInterval[k])
    {
      fIntegralPAIxSection[i] = fIntegralPAIxSection[i+1] + SumOverInterval(i);
      fIntegralPAIdEdx[i] = fIntegralPAIdEdx[i+1] + SumOverIntervaldEdx(i);
    }
    else
    {
      fIntegralPAIxSection[i] = fIntegralPAIxSection[i+1] + 
                                   SumOverBorder(i+1,fEnergyInterval[k]);
      fIntegralPAIdEdx[i] = fIntegralPAIdEdx[i+1] + 
                                   SumOverBorderdEdx(i+1,fEnergyInterval[k]);
      k--;
    }
    if(fVerbose>0) G4cout<<"i = "<<i<<"; k = "<<k<<"; intPAIxsc[i] = "<<fIntegralPAIxSection[i]<<G4endl;
  }
}   // end of IntegralPAIxSection 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI Cerenkov integral cross-section
// fIntegralCrenkov[1] = specific Crenkov ionisation, 1/cm
// and fIntegralCerenkov[0] = mean Cerenkov loss per cm  in keV/cm

void G4PAIxSection::IntegralCerenkov()
{
  G4int i, k;
   fIntegralCerenkov[fSplineNumber] = 0;
   fIntegralCerenkov[0] = 0;
   k = fIntervalNumber -1;

   for( i = fSplineNumber-1; i >= 1; i-- )
   {
      if(fSplineEnergy[i] >= fEnergyInterval[k])
      {
        fIntegralCerenkov[i] = fIntegralCerenkov[i+1] + SumOverInterCerenkov(i);
        // G4cout<<"int: i = "<<i<<"; sumC = "<<fIntegralCerenkov[i]<<G4endl;
      }
      else
      {
        fIntegralCerenkov[i] = fIntegralCerenkov[i+1] + 
                                   SumOverBordCerenkov(i+1,fEnergyInterval[k]);
        k--;
        // G4cout<<"bord: i = "<<i<<"; sumC = "<<fIntegralCerenkov[i]<<G4endl;
      }
   }

}   // end of IntegralCerenkov 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI MM-Cerenkov integral cross-section
// fIntegralMM[1] = specific MM-Cerenkov ionisation, 1/cm
// and fIntegralMM[0] = mean MM-Cerenkov loss per cm  in keV/cm

void G4PAIxSection::IntegralMM()
{
  G4int i, k;
   fIntegralMM[fSplineNumber] = 0;
   fIntegralMM[0] = 0;
   k = fIntervalNumber -1;

   for( i = fSplineNumber-1; i >= 1; i-- )
   {
      if(fSplineEnergy[i] >= fEnergyInterval[k])
      {
        fIntegralMM[i] = fIntegralMM[i+1] + SumOverInterMM(i);
        // G4cout<<"int: i = "<<i<<"; sumC = "<<fIntegralMM[i]<<G4endl;
      }
      else
      {
        fIntegralMM[i] = fIntegralMM[i+1] + 
                                   SumOverBordMM(i+1,fEnergyInterval[k]);
        k--;
        // G4cout<<"bord: i = "<<i<<"; sumC = "<<fIntegralMM[i]<<G4endl;
      }
   }

}   // end of IntegralMM 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI Plasmon integral cross-section
// fIntegralPlasmon[1] = splasmon primary ionisation, 1/cm
// and fIntegralPlasmon[0] = mean plasmon loss per cm  in keV/cm

void G4PAIxSection::IntegralPlasmon()
{
   fIntegralPlasmon[fSplineNumber] = 0;
   fIntegralPlasmon[0] = 0;
   G4int k = fIntervalNumber -1;
   for(G4int i=fSplineNumber-1;i>=1;i--)
   {
      if(fSplineEnergy[i] >= fEnergyInterval[k])
      {
        fIntegralPlasmon[i] = fIntegralPlasmon[i+1] + SumOverInterPlasmon(i);
      }
      else
      {
        fIntegralPlasmon[i] = fIntegralPlasmon[i+1] + 
                                   SumOverBordPlasmon(i+1,fEnergyInterval[k]);
        k--;
      }
   }

}   // end of IntegralPlasmon

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI resonance integral cross-section
// fIntegralResonance[1] = resonance primary ionisation, 1/cm
// and fIntegralResonance[0] = mean resonance loss per cm  in keV/cm

void G4PAIxSection::IntegralResonance()
{
   fIntegralResonance[fSplineNumber] = 0;
   fIntegralResonance[0] = 0;
   G4int k = fIntervalNumber -1;
   for(G4int i=fSplineNumber-1;i>=1;i--)
   {
      if(fSplineEnergy[i] >= fEnergyInterval[k])
      {
        fIntegralResonance[i] = fIntegralResonance[i+1] + SumOverInterResonance(i);
      }
      else
      {
        fIntegralResonance[i] = fIntegralResonance[i+1] + 
                                   SumOverBordResonance(i+1,fEnergyInterval[k]);
        k--;
      }
   }

}   // end of IntegralResonance

//////////////////////////////////////////////////////////////////////
//
// Calculation the PAI integral cross-section inside
// of interval of continuous values of photo-ionisation
// cross-section. Parameter  'i' is the number of interval.

G4double G4PAIxSection::SumOverInterval( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,c,result;

   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   if(fVerbose>0) G4cout<<"SumOverInterval i= " << i << " x0 = "<<x0<<"; x1 = "<<x1<<G4endl;

   if( x1+x0 <= 0.0 || std::abs( 2.*(x1-x0)/(x1+x0) ) < 1.e-6) return 0.;

   y0 = fDifPAIxSection[i];
   yy1 = fDifPAIxSection[i+1];

   if(fVerbose>0) G4cout<<"x0 = "<<x0<<"; x1 = "<<x1<<", y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;

   c = x1/x0;
   a = log10(yy1/y0)/log10(c);

   if(fVerbose>0) G4cout<<"SumOverInterval, a = "<<a<<"; c = "<<c<<G4endl;

   // b = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);
   a += 1.;
   if( std::abs(a) < 1.e-6 ) 
   {
      result = b*log(x1/x0);
   }
   else
   {
      result = y0*(x1*pow(c,a-1) - x0)/a;
   }
   a += 1.;
   if( std::abs(a) < 1.e-6 ) 
   {
      fIntegralPAIxSection[0] += b*log(x1/x0);
   }
   else
   {
      fIntegralPAIxSection[0] += y0*(x1*x1*pow(c,a-2) - x0*x0)/a;
   }
   if(fVerbose>0) G4cout<<"SumOverInterval, result = "<<result<<G4endl;
   return result;

} //  end of SumOverInterval

/////////////////////////////////

G4double G4PAIxSection::SumOverIntervaldEdx( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,c,result;

   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];

   if(x1+x0 <= 0.0 || std::abs( 2.*(x1-x0)/(x1+x0) ) < 1.e-6) return 0.;

   y0 = fDifPAIxSection[i];
   yy1 = fDifPAIxSection[i+1];
   c = x1/x0;
   a = log10(yy1/y0)/log10(c);
   // b = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);
   a += 2;
   if(a == 0) 
   {
     result = b*log(x1/x0);
   }
   else
   {
     result = y0*(x1*x1*pow(c,a-2) - x0*x0)/a;
   }
   return result;

} //  end of SumOverInterval

//////////////////////////////////////////////////////////////////////
//
// Calculation the PAI Cerenkov integral cross-section inside
// of interval of continuous values of photo-ionisation Cerenkov
// cross-section. Parameter  'i' is the number of interval.

G4double G4PAIxSection::SumOverInterCerenkov( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,c,result;

   x0  = fSplineEnergy[i];
   x1  = fSplineEnergy[i+1];

   if(x1+x0 <= 0.0 || std::abs( 2.*(x1-x0)/(x1+x0) ) < 1.e-6) return 0.;

   y0  = fdNdxCerenkov[i];
   yy1 = fdNdxCerenkov[i+1];
   // G4cout<<"SumC, i = "<<i<<"; x0 ="<<x0<<"; x1 = "<<x1
   //   <<"; y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;

   c = x1/x0;
   a = log10(yy1/y0)/log10(c);
   b = y0/pow(x0,a);

   a += 1.0;
   if(a == 0) result = b*log(c);
   else       result = y0*(x1*pow(c,a-1) - x0)/a;   
   a += 1.0;

   if( a == 0 ) fIntegralCerenkov[0] += b*log(x1/x0);
   else         fIntegralCerenkov[0] += y0*(x1*x1*pow(c,a-2) - x0*x0)/a;
   //  G4cout<<"a = "<<a<<"; b = "<<b<<"; result = "<<result<<G4endl;   
   return result;

} //  end of SumOverInterCerenkov

//////////////////////////////////////////////////////////////////////
//
// Calculation the PAI MM-Cerenkov integral cross-section inside
// of interval of continuous values of photo-ionisation Cerenkov
// cross-section. Parameter  'i' is the number of interval.

G4double G4PAIxSection::SumOverInterMM( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,c,result;

   x0  = fSplineEnergy[i];
   x1  = fSplineEnergy[i+1];

   if(x1+x0 <= 0.0 || std::abs( 2.*(x1-x0)/(x1+x0) ) < 1.e-6) return 0.;

   y0  = fdNdxMM[i];
   yy1 = fdNdxMM[i+1];
   //G4cout<<"SumC, i = "<<i<<"; x0 ="<<x0<<"; x1 = "<<x1
   //   <<"; y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;

   c = x1/x0;
   //G4cout<<" c = "<<c<< " yy1/y0= " << yy1/y0 <<G4endl;   
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   b = y0/pow(x0,a);

   a += 1.0;
   if(a == 0) result = b*log(c);
   else       result = y0*(x1*pow(c,a-1) - x0)/a;   
   a += 1.0;

   if( a == 0 ) fIntegralMM[0] += b*log(c);
   else         fIntegralMM[0] += y0*(x1*x1*pow(c,a-2) - x0*x0)/a;
   //G4cout<<"a = "<<a<<"; b = "<<b<<"; result = "<<result<<G4endl;   
   return result;

} //  end of SumOverInterMM

//////////////////////////////////////////////////////////////////////
//
// Calculation the PAI Plasmon integral cross-section inside
// of interval of continuous values of photo-ionisation Plasmon
// cross-section. Parameter  'i' is the number of interval.

G4double G4PAIxSection::SumOverInterPlasmon( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,c,result;

   x0  = fSplineEnergy[i];
   x1  = fSplineEnergy[i+1];

   if(x1+x0 <= 0.0 || std::abs( 2.*(x1-x0)/(x1+x0) ) < 1.e-6) return 0.;

   y0  = fdNdxPlasmon[i];
   yy1 = fdNdxPlasmon[i+1];
   c =x1/x0;
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   // b = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);

   a += 1.0;
   if(a == 0) result = b*log(x1/x0);
   else       result = y0*(x1*pow(c,a-1) - x0)/a;   
   a += 1.0;

   if( a == 0 ) fIntegralPlasmon[0] += b*log(x1/x0);
   else         fIntegralPlasmon[0] += y0*(x1*x1*pow(c,a-2) - x0*x0)/a;
   
   return result;

} //  end of SumOverInterPlasmon

//////////////////////////////////////////////////////////////////////
//
// Calculation the PAI resonance integral cross-section inside
// of interval of continuous values of photo-ionisation resonance
// cross-section. Parameter  'i' is the number of interval.

G4double G4PAIxSection::SumOverInterResonance( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,c,result;

   x0  = fSplineEnergy[i];
   x1  = fSplineEnergy[i+1];

   if(x1+x0 <= 0.0 || std::abs( 2.*(x1-x0)/(x1+x0) ) < 1.e-6) return 0.;

   y0  = fdNdxResonance[i];
   yy1 = fdNdxResonance[i+1];
   c =x1/x0;
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   // b = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);

   a += 1.0;
   if(a == 0) result = b*log(x1/x0);
   else       result = y0*(x1*pow(c,a-1) - x0)/a;   
   a += 1.0;

   if( a == 0 ) fIntegralResonance[0] += b*log(x1/x0);
   else         fIntegralResonance[0] += y0*(x1*x1*pow(c,a-2) - x0*x0)/a;
   
   return result;

} //  end of SumOverInterResonance

///////////////////////////////////////////////////////////////////////////////
//
// Integration of PAI cross-section for the case of
// passing across border between intervals

G4double G4PAIxSection::SumOverBorder( G4int      i , 
                                       G4double en0    )
{               
  G4double x0,x1,y0,yy1,a,b,/*c,*/d,e0,result;

   e0 = en0;
   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   y0 = fDifPAIxSection[i];
   yy1 = fDifPAIxSection[i+1];

   //c = x1/x0;
   d = e0/x0;   
   a = log10(yy1/y0)/log10(x1/x0);
   if(a > 10.0) return 0.;  

   if(fVerbose>0) G4cout<<"SumOverBorder, a = "<<a<<G4endl;

   // b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);  // pow(10.,b);
   
   a += 1.;
   if( std::abs(a) < 1.e-6 )
   {
      result = b*log(x0/e0);
   }
   else
   {
      result = y0*(x0 - e0*pow(d,a-1))/a;
   }
   a += 1.;
   if( std::abs(a) < 1.e-6 )
   {
      fIntegralPAIxSection[0] += b*log(x0/e0);
   }
   else 
   {
      fIntegralPAIxSection[0] += y0*(x0*x0 - e0*e0*pow(d,a-2))/a;
   }
   x0 = fSplineEnergy[i - 1];
   x1 = fSplineEnergy[i - 2];
   y0 = fDifPAIxSection[i - 1];
   yy1 = fDifPAIxSection[i - 2];

   //c = x1/x0;
   d = e0/x0;   
   a = log10(yy1/y0)/log10(x1/x0);
   //  b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);
   a += 1.;
   if( std::abs(a) < 1.e-6 )
   {
      result += b*log(e0/x0);
   }
   else
   {
      result += y0*(e0*pow(d,a-1) - x0)/a;
   }
   a += 1.;
   if( std::abs(a) < 1.e-6 ) 
   {
      fIntegralPAIxSection[0] += b*log(e0/x0);
   }
   else
   {
      fIntegralPAIxSection[0] += y0*(e0*e0*pow(d,a-2) - x0*x0)/a;
   }
   return result;

} 

///////////////////////////////////////////////////////////////////////

G4double G4PAIxSection::SumOverBorderdEdx( G4int      i , 
                                       G4double en0    )
{               
  G4double x0,x1,y0,yy1,a,b,/*c,*/d,e0,result;

   e0 = en0;
   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   y0 = fDifPAIxSection[i];
   yy1 = fDifPAIxSection[i+1];

   //c = x1/x0;
   d = e0/x0;   
   a = log10(yy1/y0)/log10(x1/x0);
   if(a > 10.0) return 0.;  
   // b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);  // pow(10.,b);
   
   a += 2;
   if(a == 0)
   {
      result = b*log(x0/e0);
   }
   else 
   {
      result = y0*(x0*x0 - e0*e0*pow(d,a-2))/a;
   }
   x0 = fSplineEnergy[i - 1];
   x1 = fSplineEnergy[i - 2];
   y0 = fDifPAIxSection[i - 1];
   yy1 = fDifPAIxSection[i - 2];

   // c = x1/x0;
   d = e0/x0;   
   a = log10(yy1/y0)/log10(x1/x0);
   //  b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);
   a += 2;
   if(a == 0) 
   {
      result += b*log(e0/x0);
   }
   else
   {
      result += y0*(e0*e0*pow(d,a-2) - x0*x0)/a;
   }
   return result;

} 

///////////////////////////////////////////////////////////////////////////////
//
// Integration of Cerenkov cross-section for the case of
// passing across border between intervals

G4double G4PAIxSection::SumOverBordCerenkov( G4int      i , 
                                             G4double en0    )
{               
   G4double x0,x1,y0,yy1,a,b,e0,c,d,result;

   e0 = en0;
   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   y0 = fdNdxCerenkov[i];
   yy1 = fdNdxCerenkov[i+1];

   //  G4cout<<G4endl;
   //  G4cout<<"SumBordC, i = "<<i<<"; en0 = "<<en0<<"; x0 ="<<x0<<"; x1 = "<<x1
   //     <<"; y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;
   c = x1/x0;
   d = e0/x0;
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   // b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a); // pow(10.,b0);   
   
   a += 1.0;
   if( a == 0 ) result = b*log(x0/e0);
   else         result = y0*(x0 - e0*pow(d,a-1))/a;   
   a += 1.0;

   if( a == 0 ) fIntegralCerenkov[0] += b*log(x0/e0);
   else         fIntegralCerenkov[0] += y0*(x0*x0 - e0*e0*pow(d,a-2))/a;

// G4cout<<"a = "<<a<<"; b0 = "<<b0<<"; b = "<<b<<"; result = "<<result<<G4endl;
   
   x0  = fSplineEnergy[i - 1];
   x1  = fSplineEnergy[i - 2];
   y0  = fdNdxCerenkov[i - 1];
   yy1 = fdNdxCerenkov[i - 2];

   // G4cout<<"x0 ="<<x0<<"; x1 = "<<x1
   //    <<"; y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;

   c = x1/x0;
   d = e0/x0;
   a  = log10(yy1/y0)/log10(x1/x0);
   // b0 = log10(y0) - a*log10(x0);
   b  =  y0/pow(x0,a);  // pow(10.,b0);

   a += 1.0;
   if( a == 0 ) result += b*log(e0/x0);
   else         result += y0*(e0*pow(d,a-1) - x0 )/a;
   a += 1.0;

   if( a == 0 )   fIntegralCerenkov[0] += b*log(e0/x0);
   else           fIntegralCerenkov[0] += y0*(e0*e0*pow(d,a-2) - x0*x0)/a;

   // G4cout<<"a = "<<a<<"; b0 = "<<b0<<"; b = "
   // <<b<<"; result = "<<result<<G4endl;    

   return result;

} 

///////////////////////////////////////////////////////////////////////////////
//
// Integration of MM-Cerenkov cross-section for the case of
// passing across border between intervals

G4double G4PAIxSection::SumOverBordMM( G4int      i , 
                                             G4double en0    )
{               
   G4double x0,x1,y0,yy1,a,b,e0,c,d,result;

   e0 = en0;
   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   y0 = fdNdxMM[i];
   yy1 = fdNdxMM[i+1];

   //  G4cout<<G4endl;
   //  G4cout<<"SumBordC, i = "<<i<<"; en0 = "<<en0<<"; x0 ="<<x0<<"; x1 = "<<x1
   //     <<"; y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;
   c = x1/x0;
   d = e0/x0;
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   // b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a); // pow(10.,b0);   
   
   a += 1.0;
   if( a == 0 ) result = b*log(x0/e0);
   else         result = y0*(x0 - e0*pow(d,a-1))/a;   
   a += 1.0;

   if( a == 0 ) fIntegralMM[0] += b*log(x0/e0);
   else         fIntegralMM[0] += y0*(x0*x0 - e0*e0*pow(d,a-2))/a;

// G4cout<<"a = "<<a<<"; b0 = "<<b0<<"; b = "<<b<<"; result = "<<result<<G4endl;
   
   x0  = fSplineEnergy[i - 1];
   x1  = fSplineEnergy[i - 2];
   y0  = fdNdxMM[i - 1];
   yy1 = fdNdxMM[i - 2];

   // G4cout<<"x0 ="<<x0<<"; x1 = "<<x1
   //    <<"; y0 = "<<y0<<"; yy1 = "<<yy1<<G4endl;

   c = x1/x0;
   d = e0/x0;
   a  = log10(yy1/y0)/log10(x1/x0);
   // b0 = log10(y0) - a*log10(x0);
   b  =  y0/pow(x0,a);  // pow(10.,b0);

   a += 1.0;
   if( a == 0 ) result += b*log(e0/x0);
   else         result += y0*(e0*pow(d,a-1) - x0 )/a;
   a += 1.0;

   if( a == 0 )   fIntegralMM[0] += b*log(e0/x0);
   else           fIntegralMM[0] += y0*(e0*e0*pow(d,a-2) - x0*x0)/a;

   // G4cout<<"a = "<<a<<"; b0 = "<<b0<<"; b = "
   // <<b<<"; result = "<<result<<G4endl;    

   return result;

} 

///////////////////////////////////////////////////////////////////////////////
//
// Integration of Plasmon cross-section for the case of
// passing across border between intervals

G4double G4PAIxSection::SumOverBordPlasmon( G4int      i , 
                                             G4double en0    )
{               
   G4double x0,x1,y0,yy1,a,b,c,d,e0,result;

   e0 = en0;
   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   y0 = fdNdxPlasmon[i];
   yy1 = fdNdxPlasmon[i+1];

   c = x1/x0;
   d = e0/x0;   
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   //  b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a); //pow(10.,b);
   
   a += 1.0;
   if( a == 0 ) result = b*log(x0/e0);
   else         result = y0*(x0 - e0*pow(d,a-1))/a;   
   a += 1.0;

   if( a == 0 ) fIntegralPlasmon[0] += b*log(x0/e0);
   else         fIntegralPlasmon[0] += y0*(x0*x0 - e0*e0*pow(d,a-2))/a;
   
   x0 = fSplineEnergy[i - 1];
   x1 = fSplineEnergy[i - 2];
   y0 = fdNdxPlasmon[i - 1];
   yy1 = fdNdxPlasmon[i - 2];

   c = x1/x0;
   d = e0/x0;
   a = log10(yy1/y0)/log10(c);
   // b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);// pow(10.,b0);

   a += 1.0;
   if( a == 0 ) result += b*log(e0/x0);
   else         result += y0*(e0*pow(d,a-1) - x0)/a;
   a += 1.0;

   if( a == 0 )   fIntegralPlasmon[0] += b*log(e0/x0);
   else           fIntegralPlasmon[0] += y0*(e0*e0*pow(d,a-2) - x0*x0)/a;
   
   return result;

} 

///////////////////////////////////////////////////////////////////////////////
//
// Integration of resonance cross-section for the case of
// passing across border between intervals

G4double G4PAIxSection::SumOverBordResonance( G4int      i , 
                                             G4double en0    )
{               
   G4double x0,x1,y0,yy1,a,b,c,d,e0,result;

   e0 = en0;
   x0 = fSplineEnergy[i];
   x1 = fSplineEnergy[i+1];
   y0 = fdNdxResonance[i];
   yy1 = fdNdxResonance[i+1];

   c = x1/x0;
   d = e0/x0;   
   a = log10(yy1/y0)/log10(c);
   if(a > 10.0) return 0.;  
   //  b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a); //pow(10.,b);
   
   a += 1.0;
   if( a == 0 ) result = b*log(x0/e0);
   else         result = y0*(x0 - e0*pow(d,a-1))/a;   
   a += 1.0;

   if( a == 0 ) fIntegralResonance[0] += b*log(x0/e0);
   else         fIntegralResonance[0] += y0*(x0*x0 - e0*e0*pow(d,a-2))/a;
   
   x0 = fSplineEnergy[i - 1];
   x1 = fSplineEnergy[i - 2];
   y0 = fdNdxResonance[i - 1];
   yy1 = fdNdxResonance[i - 2];

   c = x1/x0;
   d = e0/x0;
   a = log10(yy1/y0)/log10(c);
   // b0 = log10(y0) - a*log10(x0);
   b = y0/pow(x0,a);// pow(10.,b0);

   a += 1.0;
   if( a == 0 ) result += b*log(e0/x0);
   else         result += y0*(e0*pow(d,a-1) - x0)/a;
   a += 1.0;

   if( a == 0 )   fIntegralResonance[0] += b*log(e0/x0);
   else           fIntegralResonance[0] += y0*(e0*e0*pow(d,a-2) - x0*x0)/a;
   
   return result;

} 

/////////////////////////////////////////////////////////////////////////
//
// Returns random PAI-total energy loss over step

G4double G4PAIxSection::GetStepEnergyLoss( G4double step )
{  
  G4long numOfCollisions;
  G4double meanNumber, loss = 0.0;

  // G4cout<<" G4PAIxSection::GetStepEnergyLoss "<<G4endl;

  meanNumber = fIntegralPAIxSection[1]*step;
  numOfCollisions = G4Poisson(meanNumber);

  //   G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl;

  while(numOfCollisions)
  {
    loss += GetEnergyTransfer();
    numOfCollisions--;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  }
  // G4cout<<"PAI energy loss = "<<loss/keV<<" keV"<<G4endl; 

  return loss;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns random PAI-total energy transfer in one collision

G4double G4PAIxSection::GetEnergyTransfer()
{  
  G4int iTransfer ;

  G4double energyTransfer, position;

  position = fIntegralPAIxSection[1]*G4UniformRand();

  for( iTransfer = 1; iTransfer <= fSplineNumber; iTransfer++ )
  {
        if( position >= fIntegralPAIxSection[iTransfer] ) break;
  }
  if(iTransfer > fSplineNumber) iTransfer--;
 
  energyTransfer = fSplineEnergy[iTransfer];

  if(iTransfer > 1)
  {
    energyTransfer -= (fSplineEnergy[iTransfer]-fSplineEnergy[iTransfer-1])*G4UniformRand();
  }
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns random Cerenkov energy loss over step

G4double G4PAIxSection::GetStepCerenkovLoss( G4double step )
{  
  G4long numOfCollisions;
  G4double meanNumber, loss = 0.0;

  // G4cout<<" G4PAIxSection::GetStepCerenkovLoss "<<G4endl;

  meanNumber = fIntegralCerenkov[1]*step;
  numOfCollisions = G4Poisson(meanNumber);

  //   G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl;

  while(numOfCollisions)
  {
    loss += GetCerenkovEnergyTransfer();
    numOfCollisions--;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  }
  // G4cout<<"PAI Cerenkov loss = "<<loss/keV<<" keV"<<G4endl; 

  return loss;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns random MM-Cerenkov energy loss over step

G4double G4PAIxSection::GetStepMMLoss( G4double step )
{  
  G4long numOfCollisions;
  G4double meanNumber, loss = 0.0;

  // G4cout<<" G4PAIxSection::GetStepMMLoss "<<G4endl;

  meanNumber = fIntegralMM[1]*step;
  numOfCollisions = G4Poisson(meanNumber);

  //   G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl;

  while(numOfCollisions)
  {
    loss += GetMMEnergyTransfer();
    numOfCollisions--;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  }
  // G4cout<<"PAI MM-Cerenkov loss = "<<loss/keV<<" keV"<<G4endl; 

  return loss;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns Cerenkov energy transfer in one collision

G4double G4PAIxSection::GetCerenkovEnergyTransfer()
{  
  G4int iTransfer ;

  G4double energyTransfer, position;

  position = fIntegralCerenkov[1]*G4UniformRand();

  for( iTransfer = 1; iTransfer <= fSplineNumber; iTransfer++ )
  {
        if( position >= fIntegralCerenkov[iTransfer] ) break;
  }
  if(iTransfer > fSplineNumber) iTransfer--;
 
  energyTransfer = fSplineEnergy[iTransfer];

  if(iTransfer > 1)
  {
    energyTransfer -= (fSplineEnergy[iTransfer]-fSplineEnergy[iTransfer-1])*G4UniformRand();
  }
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns MM-Cerenkov energy transfer in one collision

G4double G4PAIxSection::GetMMEnergyTransfer()
{  
  G4int iTransfer ;

  G4double energyTransfer, position;

  position = fIntegralMM[1]*G4UniformRand();

  for( iTransfer = 1; iTransfer <= fSplineNumber; iTransfer++ )
  {
        if( position >= fIntegralMM[iTransfer] ) break;
  }
  if(iTransfer > fSplineNumber) iTransfer--;
 
  energyTransfer = fSplineEnergy[iTransfer];

  if(iTransfer > 1)
  {
    energyTransfer -= (fSplineEnergy[iTransfer]-fSplineEnergy[iTransfer-1])*G4UniformRand();
  }
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns random plasmon energy loss over step

G4double G4PAIxSection::GetStepPlasmonLoss( G4double step )
{  
  G4long numOfCollisions;
  G4double  meanNumber, loss = 0.0;

  // G4cout<<" G4PAIxSection::GetStepPlasmonLoss "<<G4endl;

  meanNumber = fIntegralPlasmon[1]*step;
  numOfCollisions = G4Poisson(meanNumber);

  //   G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl;

  while(numOfCollisions)
  {
    loss += GetPlasmonEnergyTransfer();
    numOfCollisions--;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  }
  // G4cout<<"PAI Plasmon loss = "<<loss/keV<<" keV"<<G4endl; 

  return loss;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns plasmon energy transfer in one collision

G4double G4PAIxSection::GetPlasmonEnergyTransfer()
{  
  G4int iTransfer ;

  G4double energyTransfer, position;

  position = fIntegralPlasmon[1]*G4UniformRand();

  for( iTransfer = 1; iTransfer <= fSplineNumber; iTransfer++ )
  {
        if( position >= fIntegralPlasmon[iTransfer] ) break;
  }
  if(iTransfer > fSplineNumber) iTransfer--;
 
  energyTransfer = fSplineEnergy[iTransfer];

  if(iTransfer > 1)
  {
    energyTransfer -= (fSplineEnergy[iTransfer]-fSplineEnergy[iTransfer-1])*G4UniformRand();
  }
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////////////
//
// Returns random resonance energy loss over step

G4double G4PAIxSection::GetStepResonanceLoss( G4double step )
{  
  G4long numOfCollisions;
  G4double meanNumber, loss = 0.0;

  // G4cout<<" G4PAIxSection::GetStepCreLosnkovs "<<G4endl;

  meanNumber = fIntegralResonance[1]*step;
  numOfCollisions = G4Poisson(meanNumber);

  //   G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl;

  while(numOfCollisions)
  {
    loss += GetResonanceEnergyTransfer();
    numOfCollisions--;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  }
  // G4cout<<"PAI resonance loss = "<<loss/keV<<" keV"<<G4endl; 

  return loss;
}


/////////////////////////////////////////////////////////////////////////
//
// Returns resonance energy transfer in one collision

G4double G4PAIxSection::GetResonanceEnergyTransfer()
{  
  G4int iTransfer ;

  G4double energyTransfer, position;

  position = fIntegralResonance[1]*G4UniformRand();

  for( iTransfer = 1; iTransfer <= fSplineNumber; iTransfer++ )
  {
        if( position >= fIntegralResonance[iTransfer] ) break;
  }
  if(iTransfer > fSplineNumber) iTransfer--;
 
  energyTransfer = fSplineEnergy[iTransfer];

  if(iTransfer > 1)
  {
    energyTransfer -= (fSplineEnergy[iTransfer]-fSplineEnergy[iTransfer-1])*G4UniformRand();
  }
  return energyTransfer;
}


/////////////////////////////////////////////////////////////////////////
//
// Returns Rutherford energy transfer in one collision

G4double G4PAIxSection::GetRutherfordEnergyTransfer()
{  
  G4int iTransfer ;

  G4double energyTransfer, position;

  position = (fIntegralPlasmon[1]-fIntegralResonance[1])*G4UniformRand();

  for( iTransfer = 1; iTransfer <= fSplineNumber; iTransfer++ )
  {
        if( position >= (fIntegralPlasmon[iTransfer]-fIntegralResonance[iTransfer]) ) break;
  }
  if(iTransfer > fSplineNumber) iTransfer--;
 
  energyTransfer = fSplineEnergy[iTransfer];

  if(iTransfer > 1)
  {
    energyTransfer -= (fSplineEnergy[iTransfer]-fSplineEnergy[iTransfer-1])*G4UniformRand();
  }
  return energyTransfer;
}

/////////////////////////////////////////////////////////////////////////////
//

void G4PAIxSection::CallError(G4int i, const G4String& methodName) const
{
  G4String head = "G4PAIxSection::" + methodName + "()";
  G4ExceptionDescription ed;
  ed << "Wrong index " << i << " fSplineNumber= " << fSplineNumber;
  G4Exception(head,"pai001",FatalException,ed);
}

/////////////////////////////////////////////////////////////////////////////
//
// Init  array of Lorentz factors
//

G4int G4PAIxSection::fNumberOfGammas = 111;

const G4double G4PAIxSection::fLorentzFactor[112] =     // fNumberOfGammas+1
{
0.0,
1.094989e+00, 1.107813e+00, 1.122369e+00, 1.138890e+00, 1.157642e+00,
1.178925e+00, 1.203082e+00, 1.230500e+00, 1.261620e+00, 1.296942e+00, // 10
1.337032e+00, 1.382535e+00, 1.434181e+00, 1.492800e+00, 1.559334e+00,
1.634850e+00, 1.720562e+00, 1.817845e+00, 1.928263e+00, 2.053589e+00, // 20
2.195835e+00, 2.357285e+00, 2.540533e+00, 2.748522e+00, 2.984591e+00,
3.252533e+00, 3.556649e+00, 3.901824e+00, 4.293602e+00, 4.738274e+00, // 30
5.242981e+00, 5.815829e+00, 6.466019e+00, 7.203990e+00, 8.041596e+00,
8.992288e+00, 1.007133e+01, 1.129606e+01, 1.268614e+01, 1.426390e+01, // 40
1.605467e+01, 1.808721e+01, 2.039417e+01, 2.301259e+01, 2.598453e+01,
2.935771e+01, 3.318630e+01, 3.753180e+01, 4.246399e+01, 4.806208e+01, // 50
5.441597e+01, 6.162770e+01, 6.981310e+01, 7.910361e+01, 8.964844e+01,
1.016169e+02, 1.152013e+02, 1.306197e+02, 1.481198e+02, 1.679826e+02, // 60
1.905270e+02, 2.161152e+02, 2.451581e+02, 2.781221e+02, 3.155365e+02,
3.580024e+02, 4.062016e+02, 4.609081e+02, 5.230007e+02, 5.934765e+02, // 70
6.734672e+02, 7.642575e+02, 8.673056e+02, 9.842662e+02, 1.117018e+03,
1.267692e+03, 1.438709e+03, 1.632816e+03, 1.853128e+03, 2.103186e+03, // 80
2.387004e+03, 2.709140e+03, 3.074768e+03, 3.489760e+03, 3.960780e+03,
4.495394e+03, 5.102185e+03, 5.790900e+03, 6.572600e+03, 7.459837e+03, // 90
8.466860e+03, 9.609843e+03, 1.090714e+04, 1.237959e+04, 1.405083e+04,
1.594771e+04, 1.810069e+04, 2.054434e+04, 2.331792e+04, 2.646595e+04, // 100
3.003901e+04, 3.409446e+04, 3.869745e+04, 4.392189e+04, 4.985168e+04,
5.658206e+04, 6.422112e+04, 7.289153e+04, 8.273254e+04, 9.390219e+04, // 110
1.065799e+05
};

///////////////////////////////////////////////////////////////////////
//
// The number of gamma for creation of  spline (near ion-min , G ~ 4 )
//

const
G4int G4PAIxSection::fRefGammaNumber = 29; 

   
//   
// end of G4PAIxSection implementation file 
//
////////////////////////////////////////////////////////////////////////////

