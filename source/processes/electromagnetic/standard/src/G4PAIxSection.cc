// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIxSection.cc,v 1.4 1999-12-15 14:51:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4PAIxSection.cc -- class implementation file
//
// GEANT 4 class implementation file --- Copyright CERN 1995
// CERN Geneva Switzerland
//
// For information related to this code, please, contact
// CERN, CN Division, ASD Group
//
// History:
// 1st version 11.06.97 V. Grichine

// 20.11.98 adapted to a new Material/SandiaTable interface, mma 

#include "G4ios.hh"
#include <math.h>
#include "G4PAIxSection.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"

/* ******************************************************************

// Init  array of Lorentz factors

const G4double G4PAIxSection::fLorentzFactor[22] =
{
          0.0 ,     1.1 ,   1.2 ,   1.3 ,    1.5 ,    1.8 ,  2.0 ,
          2.5 ,     3.0 ,   4.0 ,   7.0 ,   10.0 ,   20.0 , 40.0 ,
         70.0 ,   100.0 , 300.0 , 600.0 , 1000.0 , 3000.0 ,
      10000.0 , 50000.0
} ;

const G4int G4PAIxSection::
fRefGammaNumber = 29 ;         // The number of gamma for creation of 
                               // spline (9)

***************************************************************** */ 

// Local class constants

const G4double G4PAIxSection::fDelta = 0.005 ; // energy shift from interval border
const G4double G4PAIxSection::fError = 0.005 ; // error in lin-log approximation

const G4int G4PAIxSection::fMaxSplineSize = 500 ;  // Max size of output spline
                                                    // arrays

//////////////////////////////////////////////////////////////////
//
// Constructor
//

G4PAIxSection::G4PAIxSection(G4int materialIndex,
			     G4double maxEnergyTransfer)
{
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
   G4int i, j ;   

      fDensity                = (*theMaterialTable)[materialIndex]->GetDensity() ;
      fElectronDensity        = (*theMaterialTable)[materialIndex]->
                             GetElectronDensity() ;
      fIntervalNumber         = (*theMaterialTable)[materialIndex]->
                             GetSandiaTable()->GetMatNbOfIntervals() ;
      G4cout<<fDensity<<"\t"<<fElectronDensity<<"\t"<<fIntervalNumber<<G4endl ;
      // G4double maxEnergyTransfer  = 100*keV ;

      fEnergyInterval = new G4double[fIntervalNumber+2] ;
      fA1             = new G4double[fIntervalNumber+2] ;
      fA2             = new G4double[fIntervalNumber+2] ;
      fA3             = new G4double[fIntervalNumber+2] ;
      fA4             = new G4double[fIntervalNumber+2] ;
      for(i=1;i<=fIntervalNumber;i++)
      {
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
	 // G4cout<<fEnergyInterval[i]<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
	 //                               <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl ;
         if(fEnergyInterval[i] >= maxEnergyTransfer)
         {
            fEnergyInterval[i] = maxEnergyTransfer ;
	    fIntervalNumber = i ;
	    break;
         }
      }   
      if(fEnergyInterval[fIntervalNumber] != maxEnergyTransfer)
      {
         fIntervalNumber++;
         fEnergyInterval[fIntervalNumber] = maxEnergyTransfer ;
      }

      // Now checking, if two borders are too close together

      for(i=1;i<fIntervalNumber;i++)
      {
        if(fEnergyInterval[i+1]-fEnergyInterval[i] >
           1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]))
	{
          continue ;
	}
        else
	{
          for(j=i;j<fIntervalNumber;j++)
	  {
            fEnergyInterval[j] = fEnergyInterval[j+1] ;
                        fA1[j] = fA1[j+1] ;
                        fA2[j] = fA2[j+1] ;
                        fA3[j] = fA3[j+1] ;
                        fA4[j] = fA4[j+1] ;
	  }
          fIntervalNumber-- ;
          i-- ;
	}
      }


      /* *********************************

      fSplineEnergy          = new G4double[fMaxSplineSize] ;   
      fRePartDielectricConst = new G4double[fMaxSplineSize] ;   
      fImPartDielectricConst = new G4double[fMaxSplineSize] ;   
      fIntegralTerm          = new G4double[fMaxSplineSize] ;   
      fDifPAIxSection        = new G4double[fMaxSplineSize] ;   
      fIntegralPAIxSection   = new G4double[fMaxSplineSize] ;   
      
      for(i=0;i<fMaxSplineSize;i++)
      {
         fSplineEnergy[i]          = 0.0 ;   
         fRePartDielectricConst[i] = 0.0 ;   
         fImPartDielectricConst[i] = 0.0 ;   
         fIntegralTerm[i]          = 0.0 ;   
         fDifPAIxSection[i]        = 0.0 ;   
         fIntegralPAIxSection[i]   = 0.0 ;   
      }
      **************************************************  */   

      InitPAI() ;  // create arrays allocated above
      
      delete[] fEnergyInterval ;
      delete[] fA1 ;
      delete[] fA2 ;
      delete[] fA3 ;
      delete[] fA4 ;    
}

////////////////////////////////////////////////////////////////////////
//
// Constructor with beta*gamma square value

G4PAIxSection::G4PAIxSection( G4int materialIndex,
			      G4double maxEnergyTransfer,
			      G4double betaGammaSq,
                              G4double** photoAbsCof, 
                              G4int intNumber                   )
{
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   G4int i, j ;   
   
      fDensity                = (*theMaterialTable)[materialIndex]->GetDensity();
      fElectronDensity        = (*theMaterialTable)[materialIndex]->
                             GetElectronDensity() ;

      fIntervalNumber         = intNumber ;

// (*theMaterialTable)[materialIndex]->GetSandiaTable()->GetMatNbOfIntervals() ;

      // G4cout<<fDensity<<"\t"<<fElectronDensity<<"\t"<<fIntervalNumber<<G4endl ;
      // G4double maxEnergyTransfer  = 100*keV ;

      fEnergyInterval = new G4double[fIntervalNumber+2] ;
      fA1             = new G4double[fIntervalNumber+2] ;
      fA2             = new G4double[fIntervalNumber+2] ;
      fA3             = new G4double[fIntervalNumber+2] ;
      fA4             = new G4double[fIntervalNumber+2] ;
      for(i=1;i<=fIntervalNumber;i++)
      {
         fEnergyInterval[i] = photoAbsCof[i-1][0] ;
	   // (*theMaterialTable)[materialIndex]->
	   //                   GetSandiaTable()->GetSandiaCofForMaterial(i-1,0);
         fA1[i]             = photoAbsCof[i-1][1] ;
	   //(*theMaterialTable)[materialIndex]->
	   //     GetSandiaTable()->GetSandiaCofForMaterial(i-1,1);
         fA2[i]             = photoAbsCof[i-1][2] ;
	   //(*theMaterialTable)[materialIndex]->
	   //    GetSandiaTable()->GetSandiaCofForMaterial(i-1,2);
         fA3[i]             = photoAbsCof[i-1][3] ;
	   //(*theMaterialTable)[materialIndex]->
	   //      GetSandiaTable()->GetSandiaCofForMaterial(i-1,3);
         fA4[i]             = photoAbsCof[i-1][4] ;
	   //(*theMaterialTable)[materialIndex]->
	   //       GetSandiaTable()->GetSandiaCofForMaterial(i-1,4);
         if( i == 1 || i == fIntervalNumber)
	 {
	   // G4cout<<fEnergyInterval[i]<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
	   //         <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl ;
	 }
         if(fEnergyInterval[i] >= maxEnergyTransfer)
         {
            fEnergyInterval[i] = maxEnergyTransfer ;
	    fIntervalNumber = i ;
	    break;
         }
      }   
      if(fEnergyInterval[fIntervalNumber] != maxEnergyTransfer)
      {
         fIntervalNumber++;
         fEnergyInterval[fIntervalNumber] = maxEnergyTransfer ;
      }

      // Now checking, if two borders are too close together

      for(i=1;i<fIntervalNumber;i++)
      {
        if(fEnergyInterval[i+1]-fEnergyInterval[i] >
           1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]))
	{
          continue ;
	}
        else
	{
          for(j=i;j<fIntervalNumber;j++)
	  {
            fEnergyInterval[j] = fEnergyInterval[j+1] ;
                        fA1[j] = fA1[j+1] ;
                        fA2[j] = fA2[j+1] ;
                        fA3[j] = fA3[j+1] ;
                        fA4[j] = fA4[j+1] ;
	  }
          fIntervalNumber-- ;
          i-- ;
	}
      }

      /* *********************************
      fSplineEnergy          = new G4double[fMaxSplineSize] ;   
      fRePartDielectricConst = new G4double[fMaxSplineSize] ;   
      fImPartDielectricConst = new G4double[fMaxSplineSize] ;   
      fIntegralTerm          = new G4double[fMaxSplineSize] ;   
      fDifPAIxSection        = new G4double[fMaxSplineSize] ;   
      fIntegralPAIxSection   = new G4double[fMaxSplineSize] ;   
      
      for(i=0;i<fMaxSplineSize;i++)
      {
         fSplineEnergy[i]          = 0.0 ;   
         fRePartDielectricConst[i] = 0.0 ;   
         fImPartDielectricConst[i] = 0.0 ;   
         fIntegralTerm[i]          = 0.0 ;   
         fDifPAIxSection[i]        = 0.0 ;   
         fIntegralPAIxSection[i]   = 0.0 ;   
      }
      */ ////////////////////////

      // Preparation of fSplineEnergy array corresponding to min ionisation, G~4
      
      G4double   betaGammaSqRef = 
      fLorentzFactor[fRefGammaNumber]*fLorentzFactor[fRefGammaNumber] - 1;

      NormShift(betaGammaSqRef) ;             
      SplainPAI(betaGammaSqRef) ;
      
      // Preparation of integral PAI cross section for input betaGammaSq
   
      for(i = 1 ; i <= fSplineNumber ; i++)
      {
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
      }
      IntegralPAIxSection() ;
      
      delete[] fEnergyInterval ;
      delete[] fA1 ;
      delete[] fA2 ;
      delete[] fA3 ;
      delete[] fA4 ;    
}

////////////////////////////////////////////////////////////////////////
//
// Test Constructor with beta*gamma square value

G4PAIxSection::G4PAIxSection( G4int materialIndex,
			      G4double maxEnergyTransfer,
			      G4double betaGammaSq          )
{
   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   G4int i, j, numberOfElements ;   
   
   fDensity         = (*theMaterialTable)[materialIndex]->GetDensity();
   fElectronDensity = (*theMaterialTable)[materialIndex]->GetElectronDensity() ;

   G4SandiaTable thisMaterialSandiaTable(materialIndex) ;
   numberOfElements = (*theMaterialTable)[materialIndex]->GetNumberOfElements() ;

   G4int* thisMaterialZ = new G4int[numberOfElements] ;
   for(i=0;i<numberOfElements;i++)
   {
         thisMaterialZ[i] = (G4int)(*theMaterialTable)[materialIndex]->
                                      GetElement(i)->GetZ() ;
   }
   fIntervalNumber = thisMaterialSandiaTable.SandiaIntervals
                           (thisMaterialZ,numberOfElements) ;
   
   fIntervalNumber = thisMaterialSandiaTable.SandiaMixing
                           ( thisMaterialZ ,
                      (*theMaterialTable)[materialIndex]->GetFractionVector() ,
        		     numberOfElements,fIntervalNumber) ;


      fEnergyInterval = new G4double[fIntervalNumber+2] ;
      fA1             = new G4double[fIntervalNumber+2] ;
      fA2             = new G4double[fIntervalNumber+2] ;
      fA3             = new G4double[fIntervalNumber+2] ;
      fA4             = new G4double[fIntervalNumber+2] ;
      for(i=1;i<=fIntervalNumber;i++)
      {
         fEnergyInterval[i] = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,0) ;

   fA1[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,1)*fDensity ;
   fA2[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,2)*fDensity ;
   fA3[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,3)*fDensity ;
   fA4[i]             = thisMaterialSandiaTable.GetPhotoAbsorpCof(i,4)*fDensity ;

         if( i == 1 || i == fIntervalNumber)
	 {
	   // G4cout<<fEnergyInterval[i]<<"\t"<<fA1[i]<<"\t"<<fA2[i]<<"\t"
	   //         <<fA3[i]<<"\t"<<fA4[i]<<"\t"<<G4endl ;
	 }
         if(fEnergyInterval[i] >= maxEnergyTransfer)
         {
            fEnergyInterval[i] = maxEnergyTransfer ;
	    fIntervalNumber = i ;
	    break;
         }
      }   
      if(fEnergyInterval[fIntervalNumber] != maxEnergyTransfer)
      {
         fIntervalNumber++;
         fEnergyInterval[fIntervalNumber] = maxEnergyTransfer ;
      }

      // Now checking, if two borders are too close together

      for(i=1;i<fIntervalNumber;i++)
      {
        if(fEnergyInterval[i+1]-fEnergyInterval[i] >
           1.5*fDelta*(fEnergyInterval[i+1]+fEnergyInterval[i]))
	{
          continue ;
	}
        else
	{
          for(j=i;j<fIntervalNumber;j++)
	  {
            fEnergyInterval[j] = fEnergyInterval[j+1] ;
                        fA1[j] = fA1[j+1] ;
                        fA2[j] = fA2[j+1] ;
                        fA3[j] = fA3[j+1] ;
                        fA4[j] = fA4[j+1] ;
	  }
          fIntervalNumber-- ;
          i-- ;
	}
      }

      /* *********************************
      fSplineEnergy          = new G4double[fMaxSplineSize] ;   
      fRePartDielectricConst = new G4double[fMaxSplineSize] ;   
      fImPartDielectricConst = new G4double[fMaxSplineSize] ;   
      fIntegralTerm          = new G4double[fMaxSplineSize] ;   
      fDifPAIxSection        = new G4double[fMaxSplineSize] ;   
      fIntegralPAIxSection   = new G4double[fMaxSplineSize] ;   
      
      for(i=0;i<fMaxSplineSize;i++)
      {
         fSplineEnergy[i]          = 0.0 ;   
         fRePartDielectricConst[i] = 0.0 ;   
         fImPartDielectricConst[i] = 0.0 ;   
         fIntegralTerm[i]          = 0.0 ;   
         fDifPAIxSection[i]        = 0.0 ;   
         fIntegralPAIxSection[i]   = 0.0 ;   
      }
      */ ////////////////////////

      // Preparation of fSplineEnergy array corresponding to min ionisation, G~4
      
      G4double   betaGammaSqRef = 
      fLorentzFactor[fRefGammaNumber]*fLorentzFactor[fRefGammaNumber] - 1;

      NormShift(betaGammaSqRef) ;             
      SplainPAI(betaGammaSqRef) ;
      
      // Preparation of integral PAI cross section for input betaGammaSq
   
      for(i = 1 ; i <= fSplineNumber ; i++)
      {
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
      }
      IntegralPAIxSection() ;
      
      //   delete[] fEnergyInterval ;
      delete[] fA1 ;
      delete[] fA2 ;
      delete[] fA3 ;
      delete[] fA4 ;    
}


////////////////////////////////////////////////////////////////////////////
//
// Destructor

G4PAIxSection::~G4PAIxSection()
{
   /* ************************
   delete[] fSplineEnergy          ;   
   delete[] fRePartDielectricConst ;   
   delete[] fImPartDielectricConst ;   
   delete[] fIntegralTerm          ;   
   delete[] fDifPAIxSection        ;   
   delete[] fIntegralPAIxSection   ;
   */ ////////////////////////
}

/////////////////////////////////////////////////////////////////////////
//
// General control function for class G4PAIxSection
//

void G4PAIxSection::InitPAI()
{    
   G4int i ;
   G4double betaGammaSq = fLorentzFactor[fRefGammaNumber]*
                          fLorentzFactor[fRefGammaNumber] - 1;

   // Preparation of integral PAI cross section for reference gamma
   
   NormShift(betaGammaSq) ;             
   SplainPAI(betaGammaSq) ;
   IntegralPAIxSection() ;

   for(i = 0 ; i<=fSplineNumber ; i++)
   {
      fPAItable[i][fRefGammaNumber] = fIntegralPAIxSection[i] ;
      if(i != 0) 
      {
	 fPAItable[i][0] = fSplineEnergy[i] ;
      }
   }
   fPAItable[0][0] = fSplineNumber ;
   
   for(G4int j = 1 ; j < 112 ; j++)       // for other gammas
   {
      if(j == fRefGammaNumber)
      {
	 continue ;
      }
      betaGammaSq = fLorentzFactor[j]*fLorentzFactor[j] - 1 ;
      
      for(i = 1 ; i <= fSplineNumber ; i++)
      {
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
      }
      IntegralPAIxSection() ;
      
      for(i = 0 ; i <= fSplineNumber ; i++)
      {
         fPAItable[i][j] = fIntegralPAIxSection[i] ;
      }
   } 

}  

///////////////////////////////////////////////////////////////////////
//
// Shifting from borders to intervals Creation of first energy points
//

void G4PAIxSection::NormShift(G4double betaGammaSq)
{
   G4int i,j;
   for(i=1;i<=fIntervalNumber-1;i++)
   {
      for(j=1;j<=2;j++)
      {
         fSplineNumber = (i-1)*2 + j ;

         if(j==1)
	 {
	    fSplineEnergy[fSplineNumber]=fEnergyInterval[i]*(1+fDelta);
	 }
         else 
	 {
	    fSplineEnergy[fSplineNumber]=fEnergyInterval[i+1]*(1-fDelta);
	 }
      }
   }
   fIntegralTerm[1]=RutherfordIntegral(1,fEnergyInterval[1],fSplineEnergy[1]);
   j=1;
   for(i=2;i<=fSplineNumber;i++)
   {
      if(fSplineEnergy[i]<fEnergyInterval[j+1])
      {
         fIntegralTerm[i] = fIntegralTerm[i-1] + 
	                    RutherfordIntegral(j,fSplineEnergy[i-1],
                                                 fSplineEnergy[i]   ) ;
      }
      else
      {
         G4double x = RutherfordIntegral(j,fSplineEnergy[i-1],
                                           fEnergyInterval[j+1]   ) ;
         j++;
         fIntegralTerm[i] = fIntegralTerm[i-1] + x + 
	                    RutherfordIntegral(j,fEnergyInterval[j],
                                                 fSplineEnergy[i]    ) ;
      }
      // G4cout<<i<<"\t"<<fSplineEnergy[i]<<"\t"<<fIntegralTerm[i]<<"\n"<<G4endl;
   } 
   fNormalizationCof = 2*pi*pi*hbarc*hbarc*fine_structure_const/electron_mass_c2 ;
   fNormalizationCof *= fElectronDensity/fIntegralTerm[fSplineNumber] ;

   // G4cout<<"fNormalizationCof = "<<fNormalizationCof<<G4endl ;

	  // Calculation of PAI differrential cross-section (1/(keV*cm))
	  // in the energy points near borders of energy intervals

   for(G4int k=1;k<=fIntervalNumber-1;k++)
   {
      for(j=1;j<=2;j++)
      {
         i = (k-1)*2 + j ;
         fImPartDielectricConst[i] = fNormalizationCof*
	                             ImPartDielectricConst(k,fSplineEnergy[i]);
         fRePartDielectricConst[i] = fNormalizationCof*
	                             RePartDielectricConst(fSplineEnergy[i]);
         fIntegralTerm[i] *= fNormalizationCof;
         fDifPAIxSection[i] = DifPAIxSection(i,betaGammaSq);
      }
   }

}  // end of NormShift 

/////////////////////////////////////////////////////////////////////////
//
// Creation of new energy points as geometrical mean of existing
// one, calculation PAI_cs for them, while the error of logarithmic
// linear approximation would be smaller than 'fError'

void
   G4PAIxSection::SplainPAI(G4double betaGammaSq)
{
   G4int k = 1 ;
   G4int i = 1 ;

   while ( (i < fSplineNumber) && (fSplineNumber < fMaxSplineSize-1) )
   {
      if(fSplineEnergy[i+1] > fEnergyInterval[k+1])
      {
          k++ ;   // Here next energy point is in next energy interval
	  i++;
          continue;
      }
 	               // Shifting of arrayes for inserting the geometrical 
		       // average of 'i' and 'i+1' energy points to 'i+1' place
      fSplineNumber++;

      for(G4int j=fSplineNumber;j>=i+2;j--)
      {
         fSplineEnergy[j] = fSplineEnergy[j-1];
         fImPartDielectricConst[j] = fImPartDielectricConst[j-1];
	 fRePartDielectricConst[j] = fRePartDielectricConst[j-1];
	 fIntegralTerm[j] = fIntegralTerm[j-1];
	 fDifPAIxSection[j] = fDifPAIxSection[j-1];
      }
      G4double x1  = fSplineEnergy[i];
      G4double x2  = fSplineEnergy[i+1];
      G4double yy1 = fDifPAIxSection[i];
      G4double y2  = fDifPAIxSection[i+1];

      G4double en1 = sqrt(x1*x2);
      fSplineEnergy[i+1] = en1;

		 // Calculation of logarithmic linear approximation
		 // in this (enr) energy point, which number is 'i+1' now

      G4double a = log10(y2/yy1)/log10(x2/x1);
      G4double b = log10(yy1) - a*log10(x1);
      G4double y = a*log10(en1) + b ;
      y = pow(10,y);

		 // Calculation of the PAI dif. cross-section at this point

      fImPartDielectricConst[i+1] = fNormalizationCof*
	                            ImPartDielectricConst(k,fSplineEnergy[i+1]);
      fRePartDielectricConst[i+1] = fNormalizationCof*
	                            RePartDielectricConst(fSplineEnergy[i+1]);
      fIntegralTerm[i+1] = fIntegralTerm[i] + fNormalizationCof*
	                   RutherfordIntegral(k,fSplineEnergy[i],
                                                fSplineEnergy[i+1]);
      fDifPAIxSection[i+1] = DifPAIxSection(i+1,betaGammaSq);

		  // Condition for next division of this segment or to pass
		  // to higher energies

      G4double x = 2*(fDifPAIxSection[i+1] - y)/(fDifPAIxSection[i+1] + y);

      if( x < 0 ) 
      {
	 x = -x ;
      }
      if( x > fError && fSplineNumber < fMaxSplineSize-1 )
      {
	 continue;  // next division
      }
      i += 2;  // pass to next segment

   }   // close 'while'

}  // end of SplainPAI 

////////////////////////////////////////////////////////////////////////
//
// Calculation of the PAI integral cross-section
// fIntegralPAIxSection[1] = specific primary ionisation, 1/cm
// and fIntegralPAIxSection[0] = mean energy loss per cm  in keV/cm

void G4PAIxSection::IntegralPAIxSection()
{
   fIntegralPAIxSection[fSplineNumber] = 0 ;
   fIntegralPAIxSection[0] = 0 ;
   G4int k = fIntervalNumber -1 ;
   for(G4int i=fSplineNumber-1;i>=1;i--)
   {
      if(fSplineEnergy[i] >= fEnergyInterval[k])
      {
        fIntegralPAIxSection[i] = fIntegralPAIxSection[i+1] + SumOverInterval(i) ;
      }
      else
      {
        fIntegralPAIxSection[i] = fIntegralPAIxSection[i+1] + 
	                           SumOverBorder(i+1,fEnergyInterval[k]) ;
	k-- ;
      }
   }

}   // end of IntegralPAIxSection 

////////////////////////////////////////////////////////////////////
//
// Integration over electrons that could be considered
// quasi-free at energy transfer of interest

G4double G4PAIxSection::RutherfordIntegral( G4int k,
				            G4double x1,
			  	            G4double x2   )
{
   G4double  c1, c2, c3 ;
   
   c1 = (x2 - x1)/x1/x2 ;
   c2 = (x2 - x1)*(x2 + x1)/x1/x1/x2/x2 ;
   c3 = (x2 - x1)*(x1*x1 + x1*x2 + x2*x2)/x1/x1/x1/x2/x2/x2 ;
   
   return  fA1[k]*log(x2/x1) + fA2[k]*c1 + fA3[k]*c2/2 + fA4[k]*c3/3 ;

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
   
   result = fA1[k]/energy1+fA2[k]/energy2+fA3[k]/energy3+fA4[k]/energy4 ;  
   result *=hbarc/energy1 ;
   
   return result ;

}  // end of ImPartDielectricConst 


//////////////////////////////////////////////////////////////////////////////
//
// Real part of dielectric constant minus unit
// (G4double enb - energy point)
//

G4double G4PAIxSection::RePartDielectricConst(G4double enb)
{       
   G4double x0, x02, x03, x04, x05, x1, x2, xx1 ,xx2 , xx12,
            c1, c2, c3, cof1, cof2, xln1, xln2, xln3, result ;

   x0 = enb ;
   result = 0 ;
   
   for(G4int i=1;i<=fIntervalNumber-1;i++)
   {
      x1 = fEnergyInterval[i] ;
      x2 = fEnergyInterval[i+1] ;
      xx1 = x1 - x0 ;
      xx2 = x2 - x0 ;
      xx12 = xx2/xx1 ;
      
      if(xx12<0)
      {
	 xx12 = -xx12;
      }
      xln1 = log(x2/x1) ;
      xln2 = log(xx12) ;
      xln3 = log((x2 + x0)/(x1 + x0)) ;
      x02 = x0*x0 ;
      x03 = x02*x0 ;
      x04 = x03*x0 ;
      x05 = x04*x0;
      c1  = (x2 - x1)/x1/x2 ;
      c2  = (x2 - x1)*(x2 +x1)/x1/x1/x2/x2 ;
      c3  = (x2 -x1)*(x1*x1 + x1*x2 + x2*x2)/x1/x1/x1/x2/x2/x2 ;

      result -= (fA1[i]/x02 + fA3[i]/x04)*xln1 ;
      result -= (fA2[i]/x02 + fA4[i]/x04)*c1 ;
      result -= fA3[i]*c2/2/x02 ;
      result -= fA4[i]*c3/3/x02 ;

      cof1 = fA1[i]/x02 + fA3[i]/x04 ;
      cof2 = fA2[i]/x03 + fA4[i]/x05 ;

      result += 0.5*(cof1 +cof2)*xln2 ;
      result += 0.5*(cof1 - cof2)*xln3 ;
   } 
   result *= 2*hbarc/pi ;
   
   return result ;

}   // end of RePartDielectricConst 

//////////////////////////////////////////////////////////////////////
//
// PAI differential cross-section in terms of
// simplified Allison's equation
//

G4double G4PAIxSection::DifPAIxSection( G4int              i ,
                                        G4double betaGammaSq  )
{        
   G4double be2,cof,x1,x2,x3,x4,x5,x6,x7,x8,result ;

   be2 = betaGammaSq/(1 + betaGammaSq) ;
   cof = 1 ;
   x1 = log(2*electron_mass_c2/fSplineEnergy[i]) ;
   x2 = -log((1/betaGammaSq - fRePartDielectricConst[i])*
	     (1/betaGammaSq - fRePartDielectricConst[i]) + 
	     fImPartDielectricConst[i]*fImPartDielectricConst[i])/2 ;
   
   x3 = -fRePartDielectricConst[i] + 1/betaGammaSq ;
   x5 = -1 - fRePartDielectricConst[i] +
         be2*((1 +fRePartDielectricConst[i])*(1 + fRePartDielectricConst[i]) +
	 fImPartDielectricConst[i]*fImPartDielectricConst[i]) ;

   if(fImPartDielectricConst[i]==0)
   {
      x6=0 ;
   }
   else
   {
      x7 = atan2(fImPartDielectricConst[i],x3) ;
      x6 = x5 * x7 ;
   }
    // if(fImPartDielectricConst[i] == 0) x6 = 0 ;
   
   x4 = ((x1 + x2)*fImPartDielectricConst[i] + x6)/hbarc ;
   x8 = (1 + fRePartDielectricConst[i])*(1 + fRePartDielectricConst[i]) + 
        fImPartDielectricConst[i]*fImPartDielectricConst[i] ;

   result = (x4 + cof*fIntegralTerm[i]/fSplineEnergy[i]/fSplineEnergy[i])*
            fine_structure_const/be2/pi ;

   if(fDensity >= 0.1)
   { 
      result /= x8 ;
   }
   return result ;

} // end of DifPAIxSection 

//////////////////////////////////////////////////////////////////////
//
// Calculation the PAI integral cross-section inside
// of interval of continuous values of photo-ionisation
// cross-section. Parameter  'i' is the number of interval.

G4double G4PAIxSection::SumOverInterval( G4int i )
{         
   G4double x0,x1,y0,yy1,a,b,result ;

   x0 = fSplineEnergy[i] ;
   x1 = fSplineEnergy[i+1] ;
   y0 = fDifPAIxSection[i] ;
   yy1 = fDifPAIxSection[i+1];

   a = log10(yy1/y0)/log10(x1/x0) ;
   b = log10(y0) - a*log10(x0) ;
   b = pow(10.0,b) ;
   a += 1 ;
   if(a == 0) 
   {
      result = b*log(x1/x0) ;
   }
   else
   {
      result = b*(pow(x1,a) - pow(x0,a))/a ;
   }
   a++;
   if(a == 0) 
   {
      fIntegralPAIxSection[0] += b*log(x1/x0) ;
   }
   else
   {
      fIntegralPAIxSection[0] += b*(pow(x1,a) - pow(x0,a))/a ;
   }
   return result ;

} //  end of SumOverInterval

///////////////////////////////////////////////////////////////////////////////
//
// Integration of PAI cross-section for the case of
// passing across border between intervals

G4double G4PAIxSection::SumOverBorder( G4int      i , 
                                       G4double en0    )
{               
   G4double x0,x1,y0,yy1,a,b,e0,result ;

   e0 = en0 ;
   x0 = fSplineEnergy[i] ;
   x1 = fSplineEnergy[i+1] ;
   y0 = fDifPAIxSection[i] ;
   yy1 = fDifPAIxSection[i+1] ;
   
   a = log10(yy1/y0)/log10(x1/x0) ;
   b = log10(y0) - a*log10(x0) ;
   b = pow(10,b) ;
   
   a += 1 ;
   if(a == 0)
   {
      result = b*log(x0/e0) ;
   }
   else
   {
      result = b*(pow(x0,a) - pow(e0,a))/a ;
   }
   a++ ;
   if(a == 0)
   {
      fIntegralPAIxSection[0] += b*log(x0/e0) ;
   }
   else 
   {
      fIntegralPAIxSection[0] += b*(pow(x0,a) - pow(e0,a))/a ;
   }
   x0 = fSplineEnergy[i - 1] ;
   x1 = fSplineEnergy[i - 2] ;
   y0 = fDifPAIxSection[i - 1] ;
   yy1 = fDifPAIxSection[i - 2] ;

   a = log10(yy1/y0)/log10(x1/x0) ;
   b = log10(y0) - a*log10(x0) ;
   b = pow(10,b) ;
   a += 1 ;
   if(a == 0)
   {
      result += b*log(e0/x0) ;
   }
   else
   {
      result += b*(pow(e0,a) - pow(x0,a))/a ;
   }
   a++ ;
   if(a == 0) 
   {
      fIntegralPAIxSection[0] += b*log(e0/x0) ;
   }
   else
   {
      fIntegralPAIxSection[0] += b*(pow(e0,a) - pow(x0,a))/a ;
   }
   return result ;

} 

/////////////////////////////////////////////////////////////////////////
//
//

G4double G4PAIxSection::GetStepEnergyLoss( G4double step )
{  
  G4int iTransfer  ;
  G4long numOfCollisions ;
  G4double loss = 0.0 ;
  G4double meanNumber, position ;

  // G4cout<<" G4PAIxSection::GetStepEnergyLoss "<<G4endl ;



  meanNumber = fIntegralPAIxSection[1]*step ;
  numOfCollisions = RandPoisson::shoot(meanNumber) ;

  //   G4cout<<"numOfCollisions = "<<numOfCollisions<<G4endl ;

  while(numOfCollisions)
  {
    position = fIntegralPAIxSection[1]*G4UniformRand() ;

    for( iTransfer=1 ; iTransfer<=fSplineNumber ; iTransfer++ )
    {
        if( position >= fIntegralPAIxSection[iTransfer] ) break ;
    }
    loss += fSplineEnergy[iTransfer]  ;
    numOfCollisions-- ;
  }
  // G4cout<<"PAI energy loss = "<<loss/keV<<" keV"<<G4endl ; 

  return loss ;
}



/////////////////////////////////////////////////////////////////////////////
//
// Init  array of Lorentz factors
//

G4int G4PAIxSection::fNumberOfGammas = 111 ;

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
} ;

///////////////////////////////////////////////////////////////////////
//
// The number of gamma for creation of  spline (near ion-min , G ~ 4 )
//

const
G4int G4PAIxSection::fRefGammaNumber = 29 ; 

   
//   
// end of G4PAIxSection implementation file 
//
////////////////////////////////////////////////////////////////////////////

