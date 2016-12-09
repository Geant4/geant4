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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 09-May-06 fix in Sample by T. Koi
// 080318 Fix Compilation warnings - gcc-4.3.0 by T. Koi
//        (This fix has a real effect to the code.) 
// 080409 Fix div0 error with G4FPE by T. Koi
// 080612 Fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #1
// 080714 Limiting the sum of energy of secondary particles by T. Koi
// 080801 Fix div0 error wiht G4FPE and memory leak by T. Koi
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPContAngularPar.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleHPLegendreStore.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4ParticleHPVector.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleHPKallbachMannSyst.hh"
#include "G4IonTable.hh"
#include <set>
 
G4ParticleHPContAngularPar::G4ParticleHPContAngularPar( G4ParticleDefinition* projectile )
{  
  theAngular = 0;
  if ( fCache.Get() == NULL ) cacheInit();
  fCache.Get()->currentMeanEnergy = -2;
  fCache.Get()->fresh = true;
  adjustResult = true;
  if ( getenv( "G4PHP_DO_NOT_ADJUST_FINAL_STATE" ) ) adjustResult = false;

  theMinEner = DBL_MAX;
  theMaxEner = -DBL_MAX;
  theProjectile = projectile; 

  theEnergy = 0.0;
  nEnergies = 0;
  nDiscreteEnergies = 0;
  nAngularParameters = 0;
}

  void G4ParticleHPContAngularPar::Init(std::istream & aDataFile, G4ParticleDefinition* projectile)
  { 
    adjustResult = true;
    if ( getenv( "G4PHP_DO_NOT_ADJUST_FINAL_STATE" ) ) adjustResult = false;

    theProjectile = projectile;

    aDataFile >> theEnergy >> nEnergies >> nDiscreteEnergies >> nAngularParameters;
    /*if( getenv("G4PHPTEST") )*/
    theEnergy *= eV;
    theAngular = new G4ParticleHPList [nEnergies];
    for(G4int i=0; i<nEnergies; i++)
    {
      G4double sEnergy;
      aDataFile >> sEnergy;
      sEnergy*=eV;
      theAngular[i].SetLabel(sEnergy);
      theAngular[i].Init(aDataFile, nAngularParameters, 1.);
      theMinEner = std::min(theMinEner,sEnergy);
      theMaxEner = std::max(theMaxEner,sEnergy);
    }
  }

  G4ReactionProduct * 
  G4ParticleHPContAngularPar::Sample(G4double anEnergy, G4double massCode, G4double /*targetMass*/, 
                                    G4int angularRep, G4int /*interpolE*/ )
  {
    if( getenv("G4PHPTEST") ) G4cout << "  G4ParticleHPContAngularPar::Sample " << anEnergy << " " << massCode << " " << angularRep << G4endl; //GDEB
    if ( fCache.Get() == NULL ) cacheInit();
    G4ReactionProduct * result = new G4ReactionProduct;
    G4int Z = static_cast<G4int>(massCode/1000);
    G4int A = static_cast<G4int>(massCode-1000*Z);
    if(massCode==0)
    {
      result->SetDefinition(G4Gamma::Gamma());
    }
    else if(A==0)
    {
      result->SetDefinition(G4Electron::Electron());     
      if(Z==1) result->SetDefinition(G4Positron::Positron());
    }
    else if(A==1)
    {
      result->SetDefinition(G4Neutron::Neutron());
      if(Z==1) result->SetDefinition(G4Proton::Proton());
    }
    else if(A==2)
    {
      result->SetDefinition(G4Deuteron::Deuteron());      
    }
    else if(A==3)
    {
      result->SetDefinition(G4Triton::Triton());  
      if(Z==2) result->SetDefinition(G4He3::He3());
    }
    else if(A==4)
    {
      result->SetDefinition(G4Alpha::Alpha());
      if(Z!=2) throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPContAngularPar: Unknown ion case 1");    
    }
    else
    {
       result->SetDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0));
    }
    G4int i(0);
    G4int it(0);
    G4double fsEnergy(0);
    G4double cosTh(0);

   if( angularRep == 1 )
   {
// 080612 Fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #1
       //if (interpolE == 2)
//110609 above was wrong interupition, pointed out by E.Mendoza and D.Cano (CIMAT)
//Following are reviesd version written by T.Koi (SLAC)
      if ( nDiscreteEnergies != 0 )
      {

//1st check remaining_energy 
//	if this is the first set it. (How?)
         if ( fCache.Get()->fresh == true ) 
         { 
            //Discrete Lines, larger energies come first 
            //Continues Emssions, low to high                                      LAST  
            fCache.Get()->remaining_energy = std::max ( theAngular[0].GetLabel() , theAngular[nEnergies-1].GetLabel() );
            fCache.Get()->fresh = false; 
         }

         //Cheating for small remaining_energy 
         //TEMPORAL SOLUTION
         if ( nDiscreteEnergies == nEnergies )
         {
            fCache.Get()->remaining_energy = std::max ( fCache.Get()->remaining_energy , theAngular[nDiscreteEnergies-1].GetLabel() ); //Minimum Line
         }
         else
         {
            //G4double cont_min = theAngular[nDiscreteEnergies].GetLabel();   
            //if ( theAngular[nDiscreteEnergies].GetLabel() == 0.0 ) cont_min = theAngular[nDiscreteEnergies+1].GetLabel();   
            G4double cont_min=0.0; 
            for ( G4int j = nDiscreteEnergies ; j < nEnergies ; j++ )
            {
               cont_min = theAngular[j].GetLabel();   
               if ( theAngular[j].GetValue(0) != 0.0 ) break;  
            }
            fCache.Get()->remaining_energy = std::max ( fCache.Get()->remaining_energy , std::min ( theAngular[nDiscreteEnergies-1].GetLabel() , cont_min ) );   //Minimum Line or grid 
         }
//
	 G4double random = G4UniformRand();

	 G4double * running = new G4double[nEnergies+1];
	 running[0] = 0.0;

         for ( G4int j = 0 ; j < nDiscreteEnergies ; j++ ) 
         {
            G4double delta = 0.0;
            if ( theAngular[j].GetLabel() <= fCache.Get()->remaining_energy ) delta = theAngular[i].GetValue(0);
            running[j+1] = running[j] + delta;
         }
         G4double tot_prob_DIS = running[ nDiscreteEnergies ];
 
         for ( G4int j = nDiscreteEnergies ; j < nEnergies ; j++ ) 
         {
            G4double delta = 0.0;
            G4double e_low = 0.0;
            G4double e_high = 0.0;
            if ( theAngular[j].GetLabel() <= fCache.Get()->remaining_energy ) delta = theAngular[j].GetValue(0);

            //To calculate Prob. e_low and e_high should be in eV 
            //There are two case
            //1:theAngular[nDiscreteEnergies].GetLabel() != 0.0
            //   delta should be used between j-1 and j 
            //   At j = nDiscreteEnergies (the first) e_low should be set explicitly  
            if ( theAngular[j].GetLabel() != 0 )
            {
               if ( j == nDiscreteEnergies ) {
                  e_low = 0.0/eV;
               } else {
                  e_low = theAngular[j-1].GetLabel()/eV;
               }
               e_high = theAngular[j].GetLabel()/eV;
            }
            //2:theAngular[nDiscreteEnergies].GetLabel() == 0.0
            //   delta should be used between j and j+1 
            if ( theAngular[j].GetLabel() == 0.0 ) {
               e_low = theAngular[j].GetLabel()/eV;
               if ( j != nEnergies-1 ) {
                  e_high = theAngular[j+1].GetLabel()/eV;
               } else {
                  e_high = theAngular[j].GetLabel()/eV;
                  if ( theAngular[j].GetValue(0) != 0.0 ) {
                     throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPContAngularPar: Unexpected non zero value of theAngular[nEnergies-1].GetValue(0)");    
                  }
               }
            }

            running[j+1] = running[j] + ( ( e_high - e_low ) * delta );
         }
         G4double tot_prob_CON = running[ nEnergies ] - running[ nDiscreteEnergies ];

/*
         For FPE debugging 
         if (tot_prob_DIS + tot_prob_CON == 0 ) { 
            G4cout << "TKDB tot_prob_DIS + tot_prob_CON " << tot_prob_DIS + tot_prob_CON << G4endl;
            G4cout << "massCode " << massCode << G4endl;
            G4cout << "nDiscreteEnergies " << nDiscreteEnergies << " nEnergies " << nEnergies << G4endl;
            for ( int j = nDiscreteEnergies ; j < nEnergies ; j++ ) {
               G4cout << j << " " << theAngular[j].GetLabel() << " " << theAngular[j].GetValue(0) << G4endl;
            }
          }
*/
         // Normalize random 
         random *= (tot_prob_DIS + tot_prob_CON);
//2nd Judge Discrete or not             This shoudl be relatively close to 1  For safty 
         if ( random <= ( tot_prob_DIS / ( tot_prob_DIS + tot_prob_CON ) ) || nDiscreteEnergies == nEnergies )      
         {
//          Discrete Emission 
            for ( G4int j = 0 ; j < nDiscreteEnergies ; j++ )
	    {
               //Here we should use i+1
	       if ( random < running[ j+1 ] ) 
               {
                  it = j; 
                  break;
               }
            }
            fsEnergy = theAngular[ it ].GetLabel();

 	    G4ParticleHPLegendreStore theStore(1);
	    theStore.Init(0,fsEnergy,nAngularParameters);
	    for (G4int j=0;j<nAngularParameters;j++)
	    {
	       theStore.SetCoeff(0,j,theAngular[it].GetValue(j));
	    }
	    // use it to sample.
	    cosTh = theStore.SampleMax(fsEnergy);
         //Done 
         }
         else
         {
//          Continuous Emission
            for ( G4int j = nDiscreteEnergies ; j < nEnergies ; j++ )
	    {
               //Here we should use i
	       if ( random < running[ j ] ) 
               {
                  it = j; 
                  break;
               }
            }

            G4double x1 = running[it-1];
            G4double x2 = running[it];

            G4double y1 = 0.0;
            if ( it != nDiscreteEnergies ) 
                y1 = theAngular[it-1].GetLabel();
            G4double y2 = theAngular[it].GetLabel();

            fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                         random,x1,x2,y1,y2);

            G4ParticleHPLegendreStore theStore(2);
            theStore.Init(0,y1,nAngularParameters);
            theStore.Init(1,y2,nAngularParameters);
            theStore.SetManager(theManager);
            for (G4int j=0;j<nAngularParameters;j++)
            {
               G4int itt = it;
               if ( it == nDiscreteEnergies ) itt = it+1; //"This case "it-1" has data for Discrete, so we will use an extrpolate values it and it+1
               if ( it == 0 ) 
               {
                  //Safty for unexpected it = 0;
                  //G4cout << "110611 G4ParticleHPContAngularPar::Sample it = 0; invetigation required " << G4endl;
                  itt = it+1; 
               }
               theStore.SetCoeff(0,j,theAngular[itt-1].GetValue(j));
               theStore.SetCoeff(1,j,theAngular[itt].GetValue(j));
            }
            // use it to sample.
            cosTh = theStore.SampleMax(fsEnergy);

        //Done 
        }

         //TK080711
	 if( adjustResult )  fCache.Get()->remaining_energy -= fsEnergy;
         //TK080711

         //080801b
	 delete[] running;
         //080801b
      } 
      else 
      {
         // Only continue, TK will clean up 

         //080714 
         if ( fCache.Get()->fresh == true )
         {
            fCache.Get()->remaining_energy = theAngular[ nEnergies-1 ].GetLabel();
            fCache.Get()->fresh = false;
         }
         //080714 
         G4double random = G4UniformRand();
         G4double * running = new G4double[nEnergies];
         running[0]=0;
         G4double weighted = 0;
         for(i=1; i<nEnergies; i++)
         {
/*
           if(i!=0) 
           {
             running[i]=running[i-1];
           }
           running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                                theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                                theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
           weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                                theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                                theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
*/

             running[i]=running[i-1];
             if ( fCache.Get()->remaining_energy >= theAngular[i].GetLabel() )
             {
                running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                                 theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                                 theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
                weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                                 theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                                 theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
             }
         }
         // cash the mean energy in this distribution
         //080409 TKDB
         if ( nEnergies == 1 || running[nEnergies-1] == 0 )  
            fCache.Get()->currentMeanEnergy = 0.0;
         else
         { 
            fCache.Get()->currentMeanEnergy = weighted/running[nEnergies-1];
         }
         
         //080409 TKDB
         if ( nEnergies == 1 ) it = 0; 

         //080729
         if ( running[nEnergies-1] != 0 )  
         {
            for ( i = 1 ; i < nEnergies ; i++ )
            {
               it = i;
               if ( random < running [ i ] / running [ nEnergies-1 ] ) break;
            } 
         }

         //080714
         if ( running [ nEnergies-1 ] == 0 ) it = 0;
         //080714

         if (it<nDiscreteEnergies||it==0) 
         {
           if(it == 0)
           {
             fsEnergy = theAngular[it].GetLabel();
             G4ParticleHPLegendreStore theStore(1);
             theStore.Init(0,fsEnergy,nAngularParameters);
             for(i=0;i<nAngularParameters;i++)
             {
               theStore.SetCoeff(0,i,theAngular[it].GetValue(i));
             }
             // use it to sample.
             cosTh = theStore.SampleMax(fsEnergy);
           }
           else
           {
             G4double e1, e2;
             e1 = theAngular[it-1].GetLabel();
             e2 = theAngular[it].GetLabel();
             fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                           random,
                                           running[it-1]/running[nEnergies-1], 
                                           running[it]/running[nEnergies-1],
                                           e1, e2);
             // fill a Legendrestore
             G4ParticleHPLegendreStore theStore(2);
             theStore.Init(0,e1,nAngularParameters);
             theStore.Init(1,e2,nAngularParameters);
             for(i=0;i<nAngularParameters;i++)
             {
               theStore.SetCoeff(0,i,theAngular[it-1].GetValue(i));
               theStore.SetCoeff(1,i,theAngular[it].GetValue(i));
             }
             // use it to sample.
             theStore.SetManager(theManager);
             cosTh = theStore.SampleMax(fsEnergy);
           }
         }
         else // continuum contribution
         {
           G4double x1 = running[it-1]/running[nEnergies-1];
           G4double x2 = running[it]/running[nEnergies-1];
           G4double y1 = theAngular[it-1].GetLabel();
           G4double y2 = theAngular[it].GetLabel();
           fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                         random,x1,x2,y1,y2);
           G4ParticleHPLegendreStore theStore(2);
           theStore.Init(0,y1,nAngularParameters);
           theStore.Init(1,y2,nAngularParameters);
           theStore.SetManager(theManager);
           for(i=0;i<nAngularParameters;i++)
           {
             theStore.SetCoeff(0,i,theAngular[it-1].GetValue(i));
             theStore.SetCoeff(1,i,theAngular[it].GetValue(i));
           }
           // use it to sample.
           cosTh = theStore.SampleMax(fsEnergy);
         }
         delete [] running;

         //080714
	 if( adjustResult )  fCache.Get()->remaining_energy -= fsEnergy;
         //080714
      }
   }
    else if(angularRep==2)
    {
      // first get the energy (already the right for this incoming energy)
      G4int j;
      G4double * running = new G4double[nEnergies];
      running[0]=0;
      G4double weighted = 0;
      if( getenv("G4PHPTEST") ) G4cout << "  G4ParticleHPContAngularPar::Sample nEnergies " << nEnergies << G4endl;
      for(j=1; j<nEnergies; j++)
      {
        if(j!=0) running[j]=running[j-1];
        running[j] += theInt.GetBinIntegral(theManager.GetScheme(j-1),
                             theAngular[j-1].GetLabel(), theAngular[j].GetLabel(),
                             theAngular[j-1].GetValue(0), theAngular[j].GetValue(0));
        weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(j-1),
                             theAngular[j-1].GetLabel(), theAngular[j].GetLabel(),
                             theAngular[j-1].GetValue(0), theAngular[j].GetValue(0));
	if( getenv("G4PHPTEST") ) G4cout << "  G4ParticleHPContAngularPar::Sample " << j << " running " << running[j] 
	       << " " << theManager.GetScheme(j-1) << " " << theAngular[j-1].GetLabel() << " " <<  theAngular[j].GetLabel() << " " << theAngular[j-1].GetValue(0) << " " <<  theAngular[j].GetValue(0) << G4endl; //GDEB
      }
      // cash the mean energy in this distribution
      //080409 TKDB
      //currentMeanEnergy = weighted/running[nEnergies-1];
      if ( nEnergies == 1 )
         fCache.Get()->currentMeanEnergy = 0.0;
      else
        fCache.Get()->currentMeanEnergy = weighted/running[nEnergies-1];
      
      G4int itt(0);
      G4double randkal = G4UniformRand();
      //080409 TKDB
      //for(i=0; i<nEnergies; i++)
      for(j=1; j<nEnergies; j++)
      {
        itt = j;
        if(randkal<running[j]/running[nEnergies-1]) break;
      }
      
      // interpolate the secondary energy.
      G4double x, x1,x2,y1,y2;
      if(itt==0) itt=1;
      x = randkal*running[nEnergies-1];
      x1 = running[itt-1];
      x2 = running[itt];
      G4double compoundFraction;
      // interpolate energy
      y1 = theAngular[itt-1].GetLabel();
      y2 = theAngular[itt].GetLabel();
      fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(itt-1), 
                                    x, x1,x2,y1,y2);
      if( getenv("G4PHPTEST") ) G4cout << itt << " G4particleHPContAngularPar fsEnergy " << fsEnergy << " " << theManager.GetInverseScheme(itt-1) << " x " << x << " " << x1 << " " << x2 << " y " << y1 << " " << y2 << G4endl; //GDEB
      // for theta interpolate the compoundFractions
      G4double cLow = theAngular[itt-1].GetValue(1);
      G4double cHigh = theAngular[itt].GetValue(1);
      compoundFraction = theInt.Interpolate(theManager.GetScheme(itt),
                                            fsEnergy, y1, y2, cLow,cHigh);
      if( getenv("G4PHPTEST") )  G4cout << itt << " G4particleHPContAngularPar compoundFraction " << compoundFraction << " E " << fsEnergy << " " << theManager.GetScheme(itt) << " ener " << fsEnergy << " y " << y1 << " " << y2 << " cLH " << cLow << " " << cHigh << G4endl; //GDEB
      delete [] running;
      
      // get cosTh
      G4double incidentEnergy = anEnergy;
      G4double incidentMass = theProjectile->GetPDGMass();
      G4double productEnergy = fsEnergy;
      G4double productMass = result->GetMass();
      G4int targetZ = G4int(fCache.Get()->theTargetCode/1000);
      G4int targetA = G4int(fCache.Get()->theTargetCode-1000*targetZ);
      // To correspond to natural composition (-nat-) data files. 
      if ( targetA == 0 ) 
         targetA = G4int ( fCache.Get()->theTarget->GetMass()/amu_c2 + 0.5 );
      G4double targetMass = fCache.Get()->theTarget->GetMass();
      G4int residualA = targetA+1-A;
      G4int residualZ = targetZ-Z;
      G4double residualMass =  residualZ*G4Proton::Proton()->GetPDGMass();
               residualMass +=(residualA-residualZ)*theProjectile->GetPDGMass();
               residualMass -= G4NucleiProperties::GetBindingEnergy( residualA , residualZ );
      G4ParticleHPKallbachMannSyst theKallbach(compoundFraction,
                                              incidentEnergy, incidentMass,
                                              productEnergy, productMass,
                                              residualMass, residualA, residualZ,
                                              targetMass, targetA, targetZ);
      cosTh = theKallbach.Sample(anEnergy);
      if( getenv("G4PHPTEST") ) G4cout << " G4ParticleHPKallbachMannSyst::Sample resulttest " << cosTh << G4endl; //GDEB
    }
    else if(angularRep>10&&angularRep<16)
    {
      G4double random = G4UniformRand();
      G4double * running = new G4double[nEnergies];
      running[0]=0;      
      G4double weighted = 0;
      for(i=1; i<nEnergies; i++)
      {
        if(i!=0) running[i]=running[i-1];
        running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
        weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                             theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                             theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
      }
       // cash the mean energy in this distribution
      //currentMeanEnergy = weighted/running[nEnergies-1];
      if ( nEnergies == 1 )  
         fCache.Get()->currentMeanEnergy = 0.0;
      else
         fCache.Get()->currentMeanEnergy = weighted/running[nEnergies-1];
      
      //080409 TKDB
      if ( nEnergies == 1 ) it = 0; 
      //for(i=0; i<nEnergies; i++)
      for(i=1; i<nEnergies; i++)
      {
        it = i;
        if(random<running[i]/running[nEnergies-1]) break;
      }
      if(it<nDiscreteEnergies||it==0) 
      {
        if(it==0)
        {
          fsEnergy = theAngular[0].GetLabel();          
          G4ParticleHPVector theStore; 
	  G4int aCounter = 0;
          for(G4int j=1; j<nAngularParameters; j+=2) 
          {
            theStore.SetX(aCounter, theAngular[0].GetValue(j));
            theStore.SetY(aCounter, theAngular[0].GetValue(j+1));
	    aCounter++;	    
          }
          G4InterpolationManager aMan;
          aMan.Init(angularRep-10, nAngularParameters-1);
          theStore.SetInterpolationManager(aMan);
          cosTh = theStore.Sample();
        }
        else 
        {
          fsEnergy = theAngular[it].GetLabel();
          G4ParticleHPVector theStore; 
          G4InterpolationManager aMan;
          aMan.Init(angularRep-10, nAngularParameters-1);
          theStore.SetInterpolationManager(aMan); // Store interpolates f(costh)
          G4InterpolationScheme currentScheme = theManager.GetInverseScheme(it);
	  G4int aCounter = 0;
          for(G4int j=1; j<nAngularParameters; j+=2) 
          {
            theStore.SetX(aCounter, theAngular[it].GetValue(j));
            theStore.SetY(aCounter, theInt.Interpolate(currentScheme, 
                                       random,
                                       running[it-1]/running[nEnergies-1],
                                       running[it]/running[nEnergies-1],
                                       theAngular[it-1].GetValue(j+1),
                                       theAngular[it].GetValue(j+1)));
	    aCounter++;	    
          }
          cosTh = theStore.Sample();
        }
      }
      else
      {
        G4double x1 = running[it-1]/running[nEnergies-1];
        G4double x2 = running[it]/running[nEnergies-1];
        G4double y1 = theAngular[it-1].GetLabel();
        G4double y2 = theAngular[it].GetLabel();
        fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                      random,x1,x2,y1,y2);
        G4ParticleHPVector theBuff1;
        G4ParticleHPVector theBuff2;
        G4InterpolationManager aMan;
        aMan.Init(angularRep-10, nAngularParameters-1);
//        theBuff1.SetInterpolationManager(aMan); // Store interpolates f(costh)
//        theBuff2.SetInterpolationManager(aMan); // Store interpolates f(costh)
//      Bug Report #1366 from L. Russell 
        //for(i=0; i<nAngularParameters; i++) // i=1 ist wichtig!
        //{
        //  theBuff1.SetX(i, theAngular[it-1].GetValue(i));
        //  theBuff1.SetY(i, theAngular[it-1].GetValue(i+1));
        //  theBuff2.SetX(i, theAngular[it].GetValue(i));
        //  theBuff2.SetY(i, theAngular[it].GetValue(i+1));
        //  i++;
        //}
        {
        G4int j;
        for(i=0,j=1; i<nAngularParameters; i++,j+=2) 
        {
          theBuff1.SetX(i, theAngular[it-1].GetValue(j));
          theBuff1.SetY(i, theAngular[it-1].GetValue(j+1));
          theBuff2.SetX(i, theAngular[it].GetValue(j));
          theBuff2.SetY(i, theAngular[it].GetValue(j+1));
        }
        }
        G4ParticleHPVector theStore;
        theStore.SetInterpolationManager(aMan); // Store interpolates f(costh)        
        x1 = y1;
        x2 = y2;
        G4double x, y;
        //for(i=0;i<theBuff1.GetVectorLength(); i++);
        for(i=0;i<theBuff1.GetVectorLength(); i++)
        {
          x = theBuff1.GetX(i); // costh binning identical
          y1 = theBuff1.GetY(i);
          y2 = theBuff2.GetY(i);
          y = theInt.Interpolate(theManager.GetScheme(it),
                                 fsEnergy, theAngular[it-1].GetLabel(), 
                                 theAngular[it].GetLabel(), y1, y2);
          theStore.SetX(i, x);
          theStore.SetY(i, y);
        }
        cosTh = theStore.Sample();
      }
      delete [] running;
    }
    else
    {
      throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPContAngularPar::Sample: Unknown angular representation");
    }
    result->SetKineticEnergy(fsEnergy);
    G4double phi = twopi*G4UniformRand();
    G4double theta = std::acos(cosTh);
    G4double sinth = std::sin(theta);
    G4double mtot = result->GetTotalMomentum();
    G4ThreeVector tempVector(mtot*sinth*std::cos(phi), mtot*sinth*std::sin(phi), mtot*std::cos(theta) );
    result->SetMomentum(tempVector);
//  return the result.    
    return result;
  }


#define MERGE_NEW

void G4ParticleHPContAngularPar::PrepareTableInterpolation(const G4ParticleHPContAngularPar* angParPrev)
{

  //----- Discrete energies: store own energies in a map for faster searching
  G4int ie;
  for(ie=0; ie<nDiscreteEnergies; ie++) {
    theDiscreteEnergiesOwn[theAngular[ie].GetLabel()] = ie;
  }
  if( !angParPrev ) return;

  //----- Discrete energies: use energies that appear in one or another
  for(ie=0; ie<nDiscreteEnergies; ie++) {
    theDiscreteEnergies.insert(theAngular[ie].GetLabel());
  }
  G4int nDiscreteEnergiesPrev = angParPrev->GetNDiscreteEnergies();
  for(ie=0; ie<nDiscreteEnergiesPrev; ie++) {
    theDiscreteEnergies.insert(angParPrev->theAngular[ie].GetLabel());
  }
  
  //--- Get the values for which interpolation will be done : all energies of this and previous ContAngularPar
  for(ie=nDiscreteEnergies; ie<nEnergies; ie++) {
    G4double ener = theAngular[ie].GetLabel();
    G4double enerT = (ener-theMinEner)/(theMaxEner-theMinEner);
    theEnergiesTransformed.insert(enerT);
    //-    if( getenv("G4PHPTEST2") ) G4cout <<this << " G4ParticleHPContAngularPar::PrepareTableInterpolation  theEnergiesTransformed1 " << enerT << G4endl; //GDEB
  } 
  G4int nEnergiesPrev = angParPrev->GetNEnergies();
  G4double minEnerPrev = angParPrev->GetMinEner();
  G4double maxEnerPrev = angParPrev->GetMaxEner();
  for(ie=nDiscreteEnergiesPrev; ie<nEnergiesPrev; ie++) {
    G4double ener = angParPrev->theAngular[ie].GetLabel();
    G4double enerT = (ener-minEnerPrev)/(maxEnerPrev-minEnerPrev);
    theEnergiesTransformed.insert(enerT);
    //-    if( getenv("G4PHPTEST2") ) G4cout << this << " G4ParticleHPContAngularPar::PrepareTableInterpolation  theEnergiesTransformed2 " << enerT << G4endl; //GDEB
  }
  // add the maximum energy
  theEnergiesTransformed.insert(1.);

}

void G4ParticleHPContAngularPar::BuildByInterpolation(G4double anEnergy, G4InterpolationScheme aScheme, 
							G4ParticleHPContAngularPar & angpar1, 
							G4ParticleHPContAngularPar & angpar2) 
{
  G4int ie,ie1,ie2, ie1Prev, ie2Prev;
  nAngularParameters = angpar1.nAngularParameters;
  theManager = angpar1.theManager;
  theEnergy = anEnergy;

  nDiscreteEnergies = theDiscreteEnergies.size();
  std::set<G4double>::const_iterator itede;
  std::map<G4double,G4int> discEnerOwn1 = angpar1.GetDiscreteEnergiesOwn();
  std::map<G4double,G4int> discEnerOwn2 = angpar2.GetDiscreteEnergiesOwn();
  std::map<G4double,G4int>::const_iterator itedeo;
  ie = 0;
  for( itede = theDiscreteEnergies.begin(); itede != theDiscreteEnergies.end(); itede++, ie++ ) {
    G4double discEner = *itede;
    itedeo = discEnerOwn1.find(discEner);
    if( itedeo == discEnerOwn1.end() ) {
      ie1 = 0;
    } else {
      ie1 = -1;
    }
    itedeo = discEnerOwn2.find(discEner);
    if( itedeo == discEnerOwn2.end() ) {
      ie2 = 0;
    } else {
      ie2 = -1;
    }

    theAngular[ie].SetLabel(discEner);
    G4double val1, val2;
    for(G4int ip=0; ip<nAngularParameters; ip++) {
      if( ie1 != -1 ) {
	val1 = angpar1.theAngular[ie1].GetValue(ip);
      } else {
	val1 = 0.;
      }
      if( ie2 != -1 ) {
	val2 = angpar2.theAngular[ie2].GetValue(ip);
      } else {
	val2 = 0.;
      }
      
      G4double value = theInt.Interpolate(aScheme, anEnergy, 
					  angpar1.theEnergy, angpar2.theEnergy,
					  val1,
					  val2);
      if( getenv("G4PHPTEST2") ) G4cout << ie << " " << ip << " G4ParticleHPContAngularPar::Merge DiscreteEnergies  val1 " << val1 << " val2 " << val2 << " value " << value << G4endl; //GDEB
      
      theAngular[ie].SetValue(ip, value);
    }
  }
  
  if(theAngular != 0) delete [] theAngular;
  nEnergies = nDiscreteEnergies + angpar2.GetNEnergiesTransformed();
  theAngular = new G4ParticleHPList [nEnergies];

  //---- Get minimum and maximum energy interpolating
  theMinEner = angpar1.GetMinEner() + (theEnergy-angpar1.GetEnergy()) * (angpar2.GetMinEner()-angpar1.GetMinEner())/(angpar2.GetEnergy()-angpar1.GetEnergy());
  theMaxEner = angpar1.GetMaxEner() + (theEnergy-angpar1.GetEnergy()) * (angpar2.GetMaxEner()-angpar1.GetMaxEner())/(angpar2.GetEnergy()-angpar1.GetEnergy());

  if( getenv("G4PHPTEST2") )  G4cout << " G4ParticleHPContAngularPar::Merge E " << anEnergy << " minmax " << theMinEner << " " << theMaxEner << G4endl; //GDEB

  //--- Loop to energies of new set
  std::set<G4double> energiesTransformed = angpar2.GetEnergiesTransformed();
  std::set<G4double>::const_iterator iteet = energiesTransformed.begin();
  G4int nEnergies1 = angpar1.GetNEnergies();
  G4int nDiscreteEnergies1 = angpar1.GetNDiscreteEnergies();
  G4double minEner1 = angpar1.GetMinEner();
  G4double maxEner1 = angpar1.GetMaxEner();
  G4int nEnergies2 = angpar2.GetNEnergies();
  G4int nDiscreteEnergies2 = angpar2.GetNDiscreteEnergies();
  G4double minEner2 = angpar2.GetMinEner();
  G4double maxEner2 = angpar2.GetMaxEner();
  for(ie=nDiscreteEnergies; ie<nEnergies; ie++,iteet++) { 
    G4double eT = (*iteet);

    //--- Use eT1 = eT: Get energy and parameters of angpar1 for this eT
    G4double e1 = (maxEner1-minEner1) * eT + minEner1;
    //----- Get parameter value corresponding to this e1
    for(ie1=nDiscreteEnergies1; ie1<nEnergies1; ie1++) {
      if( (angpar1.theAngular[ie1].GetLabel() - e1) > 1.E-10*e1 ) break;
    }
    ie1Prev = ie1 - 1;
    if( ie1 == 0 ) ie1Prev++; 
    if( ie1 == nEnergies1 ) {
      ie1--;
      ie1Prev = ie1;
    }
    //--- Use eT2 = eT: Get energy and parameters of angpar2 for this eT
    G4double e2 = (maxEner2-minEner2) * eT + minEner2;
    //----- Get parameter value corresponding to this e2
    for(ie2=nDiscreteEnergies2; ie2<nEnergies2; ie2++) {
      //      G4cout << " GET IE2 " << ie2 << " - " << angpar2.theAngular[ie2].GetLabel() - e2 << " " << angpar2.theAngular[ie2].GetLabel() << " " << e2 <<G4endl;
      if( (angpar2.theAngular[ie2].GetLabel() - e2) > 1.E-10*e2 ) break;
    }
    ie2Prev = ie2 - 1;
    if( ie2 == 0 ) ie2Prev++; 
    if( ie2 == nEnergies2 ) {
      ie2--;
      ie2Prev = ie2;
    }

    //---- Energy corresponding to energy transformed    
    G4double eN = (theMaxEner-theMinEner) * eT + theMinEner;
    if( getenv("G4PHPTEST2") )  G4cout << ie << " " << ie1 << " " << ie2 << " G4ParticleHPContAngularPar::loop eT " << eT << " -> eN " << eN << " e1 " << e1 << " e2 " << e2 << G4endl; //GDEB
    
    theAngular[ie].SetLabel(eN);
    
    for(G4int ip=0; ip<nAngularParameters; ip++) {
      G4double val1 = theInt.Interpolate2(theManager.GetScheme(ie),
					  e1,
					  angpar1.theAngular[ie1Prev].GetLabel(),
					  angpar1.theAngular[ie1].GetLabel(),
					  angpar1.theAngular[ie1Prev].GetValue(ip),
					  angpar1.theAngular[ie1].GetValue(ip)) * (maxEner1-minEner1);  
      G4double val2 = theInt.Interpolate2(theManager.GetScheme(ie),
					  e2,
					  angpar2.theAngular[ie2Prev].GetLabel(),
					  angpar2.theAngular[ie2].GetLabel(),
					  angpar2.theAngular[ie2Prev].GetValue(ip),
					  angpar2.theAngular[ie2].GetValue(ip)) * (maxEner2-minEner2);  
      
      G4double value = theInt.Interpolate(aScheme, anEnergy, 
					  angpar1.theEnergy, angpar2.theEnergy,
					  val1,
					  val2);
      //value /= (theMaxEner-theMinEner); 
      if ( theMaxEner != theMinEner ) {
         value /= (theMaxEner-theMinEner); 
      } else if ( value != 0 ) {
         throw G4HadronicException(__FILE__, __LINE__, "G4ParticleHPContAngularPar::PrepareTableInterpolation theMaxEner == theMinEner and  value != 0.");
      }
      if( getenv("G4PHPTEST2") ) G4cout << ie << " " << ip << " G4ParticleHPContAngularPar::Merge val1 " << val1 << " val2 " << val2 << " value " << value << G4endl; //GDEB
      //-      val1 = angpar1.theAngular[ie1-1].GetValue(ip) * (maxEner1-minEner1); 
      //-      val2 = angpar2.theAngular[ie2-1].GetValue(ip) * (maxEner2-minEner2); 
      //-      if( getenv("G4PHPTEST2") ) G4cout << ie << " " << ip << " G4ParticleHPContAngularPar::MergeOLD val1 " << val1 << " val2 " << val2 << " value " << value << G4endl; //GDEB
      
      theAngular[ie].SetValue(ip, value);
    }
  }

  if( getenv("G4PHPTEST2") ) {
    G4cout << " G4ParticleHPContAngularPar::Merge ANGPAR1 " << G4endl; //GDEB
    angpar1.Dump();
    G4cout << " G4ParticleHPContAngularPar::Merge ANGPAR2 " << G4endl;
    angpar2.Dump();
    G4cout << " G4ParticleHPContAngularPar::Merge ANGPARNEW " << G4endl;
    Dump();
  }   
}

void G4ParticleHPContAngularPar::Dump()
{
  G4cout << theEnergy << " " << nEnergies << " " << nDiscreteEnergies << " " << nAngularParameters << G4endl;

  for(G4int ii=0; ii<nEnergies; ii++) {
    theAngular[ii].Dump();
  }

}
