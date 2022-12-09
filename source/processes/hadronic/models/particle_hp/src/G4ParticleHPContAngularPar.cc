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
// June-2019 - E. Mendoza --> redefinition of the residual mass to consider incident particles different than neutrons.

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
#include "G4ParticleHPManager.hh"
#include <set>
#include <vector>
 
G4ParticleHPContAngularPar::G4ParticleHPContAngularPar(G4ParticleDefinition* projectile)
{  
  theAngular = 0;
  if (fCache.Get() == 0) cacheInit();
  fCache.Get()->currentMeanEnergy = -2;
  fCache.Get()->fresh = true;
  adjustResult = true;
  if (G4ParticleHPManager::GetInstance()->GetDoNotAdjustFinalState() )
    adjustResult = false;
  
  theMinEner = DBL_MAX;
  theMaxEner = -DBL_MAX;
  theProjectile = projectile; 

  theEnergy = 0.0;
  nEnergies = 0;
  nDiscreteEnergies = 0;
  nAngularParameters = 0;
}


void
G4ParticleHPContAngularPar::Init(std::istream& aDataFile, G4ParticleDefinition* projectile)
{ 
  adjustResult = true;
  if (G4ParticleHPManager::GetInstance()->GetDoNotAdjustFinalState() )
    adjustResult = false;

  theProjectile = projectile;

  aDataFile >> theEnergy >> nEnergies >> nDiscreteEnergies >> nAngularParameters;
  theEnergy *= eV;
  theAngular = new G4ParticleHPList [nEnergies];
  G4double sEnergy;
  for (G4int i = 0; i < nEnergies; i++) {
    aDataFile >> sEnergy;
    sEnergy *= eV;
    theAngular[i].SetLabel(sEnergy);
    theAngular[i].Init(aDataFile, nAngularParameters, 1.);
    theMinEner = std::min(theMinEner,sEnergy);
    theMaxEner = std::max(theMaxEner,sEnergy);
  }
}


G4ReactionProduct* 
G4ParticleHPContAngularPar::Sample(G4double anEnergy, G4double massCode,
                                   G4double /*targetMass*/, 
                                   G4int angularRep, G4int /*interpolE*/ )
{
  // The following line is needed because it may change between runs by UI command
  if (G4ParticleHPManager::GetInstance()->GetDoNotAdjustFinalState() )
    adjustResult = false;
 
  if (fCache.Get() == 0 ) cacheInit();
  G4ReactionProduct * result = new G4ReactionProduct;
  G4int Z = static_cast<G4int>(massCode/1000);
  G4int A = static_cast<G4int>(massCode-1000*Z);
  if (massCode == 0) {
    result->SetDefinition(G4Gamma::Gamma());
  } else if (A == 0) {
    result->SetDefinition(G4Electron::Electron());     
    if (Z == 1) result->SetDefinition(G4Positron::Positron());
  } else if (A == 1) {
    result->SetDefinition(G4Neutron::Neutron());
    if (Z == 1) result->SetDefinition(G4Proton::Proton());
  } else if (A == 2) {
    result->SetDefinition(G4Deuteron::Deuteron());      
  } else if (A == 3) {
    result->SetDefinition(G4Triton::Triton());  
    if(Z == 2) result->SetDefinition(G4He3::He3());
  } else if (A == 4) {
    result->SetDefinition(G4Alpha::Alpha());
    if (Z != 2) throw G4HadronicException(__FILE__, __LINE__,
                                          "G4ParticleHPContAngularPar: Unknown ion case 1");    
  } else {
    result->SetDefinition(G4IonTable::GetIonTable()->GetIon(Z,A,0) );
  }

  G4int i(0);
  G4int it(0);
  G4double fsEnergy(0);
  G4double cosTh(0);

  if (angularRep == 1) {
    if (nDiscreteEnergies != 0) {
      // 1st check remaining_energy 
      // if this is the first set it. (How?)
      if (fCache.Get()->fresh == true) { 
        // Discrete Lines, larger energies come first 
        // Continues Emssions, low to high                                 LAST  
        fCache.Get()->remaining_energy = std::max (theAngular[0].GetLabel(),
                                                   theAngular[nEnergies-1].GetLabel() );
        fCache.Get()->fresh = false; 
      }

      // Cheating for small remaining_energy 
      // Temporary solution
      if (nDiscreteEnergies == nEnergies) {
        fCache.Get()->remaining_energy = std::max(fCache.Get()->remaining_energy,
                                                  theAngular[nDiscreteEnergies-1].GetLabel() ); //Minimum Line
      } else {
        G4double cont_min = 0.0; 
        for (G4int j = nDiscreteEnergies; j < nEnergies; j++) {
          cont_min = theAngular[j].GetLabel();   
          if (theAngular[j].GetValue(0) != 0.0 ) break;  
        }
        fCache.Get()->remaining_energy =
          std::max (fCache.Get()->remaining_energy,
                    std::min(theAngular[nDiscreteEnergies-1].GetLabel(), cont_min) );   //Minimum Line or grid 
      }

      G4double random = G4UniformRand();
      G4double* running = new G4double[nEnergies+1];
      running[0] = 0.0;

      G4double delta;
      for (G4int j = 0; j < nDiscreteEnergies; j++) {
        delta = 0.0;
        if (theAngular[j].GetLabel() <= fCache.Get()->remaining_energy)
          delta = theAngular[j].GetValue(0);
        running[j+1] = running[j] + delta;
      }

      G4double tot_prob_DIS = std::max( running[nDiscreteEnergies], 0.0 );

      G4double delta1;
      for (G4int j = nDiscreteEnergies; j < nEnergies; j++) {
        delta1 = 0.0;
        G4double e_low = 0.0;
        G4double e_high = 0.0;
        if (theAngular[j].GetLabel() <= fCache.Get()->remaining_energy)
          delta1 = theAngular[j].GetValue(0);

        // To calculate Prob. e_low and e_high should be in eV 
        // There are two cases:
        // 1: theAngular[nDiscreteEnergies].GetLabel() != 0.0
        //    delta1 should be used between j-1 and j 
        //    At j = nDiscreteEnergies (the first) e_low should be set explicitly  
        if (theAngular[j].GetLabel() != 0) {
          if (j == nDiscreteEnergies) {
            e_low = 0.0/eV;
          } else {
	    if ( j < 1 ) j = 1;  // Protection against evaluation of arrays at index j-1
            e_low = theAngular[j-1].GetLabel()/eV;
          }
          e_high = theAngular[j].GetLabel()/eV;
        }

        // 2: theAngular[nDiscreteEnergies].GetLabel() == 0.0
        //    delta1 should be used between j and j+1 
        if (theAngular[j].GetLabel() == 0.0) {
          e_low = theAngular[j].GetLabel()/eV;
          if (j != nEnergies-1) {
            e_high = theAngular[j+1].GetLabel()/eV;
          } else {
            e_high = theAngular[j].GetLabel()/eV;
            if (theAngular[j].GetValue(0) != 0.0) {
              throw G4HadronicException(__FILE__, __LINE__,
                      "G4ParticleHPContAngularPar: Unexpected non zero value of theAngular[nEnergies-1].GetValue(0)");    
            }
          }
        }

        running[j+1] = running[j] + ( ( e_high - e_low ) * delta1);
      }
      G4double tot_prob_CON = std::max( running[ nEnergies ] - running[ nDiscreteEnergies ], 0.0 );

      // Give up in the pathological case of null probabilities 
      if ( tot_prob_DIS == 0.0 && tot_prob_CON == 0.0 ) return result;

      // Normalize random 
      random *= (tot_prob_DIS + tot_prob_CON);
      // 2nd Judge Discrete or not

      // This should be relatively close to 1  For safty 
      if (random <= (tot_prob_DIS/(tot_prob_DIS + tot_prob_CON) ) || nDiscreteEnergies == nEnergies) {
        // Discrete Emission 
        for (G4int j = 0; j < nDiscreteEnergies; j++) {
          // Here we should use i+1
	  if (random < running[ j+1 ] ) {
            it = j; 
            break;
          }
        }
        fsEnergy = theAngular[ it ].GetLabel();

        G4ParticleHPLegendreStore theStore(1);
        theStore.Init(0,fsEnergy,nAngularParameters);
        for (G4int j = 0; j < nAngularParameters; j++) {
          theStore.SetCoeff(0,j,theAngular[it].GetValue(j));
        }
        // use it to sample.
        cosTh = theStore.SampleMax(fsEnergy);
        // Done

      } else {
        // Continuous emission
        for (G4int j = nDiscreteEnergies; j < nEnergies; j++) {
          // Here we should use i
          if (random < running[ j ] ) {
            it = j; 
            break;
          }
        }

	if ( it < 1 ) it = 1;  // Protection against evaluation of arrays at index it-1

        G4double x1 = running[it-1];
        G4double x2 = running[it];

        G4double y1 = 0.0;
        if (it != nDiscreteEnergies) y1 = theAngular[it-1].GetLabel();
        G4double y2 = theAngular[it].GetLabel();

        fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(it),
                                      random,x1,x2,y1,y2);

        G4ParticleHPLegendreStore theStore(2);
        theStore.Init(0, y1, nAngularParameters);
        theStore.Init(1, y2, nAngularParameters);
        theStore.SetManager(theManager);
        G4int itt;
        for (G4int j = 0; j < nAngularParameters; j++) {
          itt = it;
          if (it == nDiscreteEnergies) itt = it+1;
          // "This case "it-1" has data for Discrete, so we will use an extrpolated values it and it+1
          theStore.SetCoeff(0, j, theAngular[itt-1].GetValue(j) );
          theStore.SetCoeff(1, j, theAngular[itt].GetValue(j) );
        }
        // use it to sample.
        cosTh = theStore.SampleMax(fsEnergy);

        //Done 
      }

      // The remaining energy needs to be lowered by the photon energy in *any* case.
      // Otherwise additional photons with too high energy will be produced - therefore the
	 // adjustResult condition has been removed
      fCache.Get()->remaining_energy -= fsEnergy;
      delete[] running;

      // end (nDiscreteEnergies != 0) branch 

    } else {
      // Only continue, TK will clean up 
      if (fCache.Get()->fresh == true) {
        fCache.Get()->remaining_energy = theAngular[ nEnergies-1 ].GetLabel();
        fCache.Get()->fresh = false;
      }

      G4double random = G4UniformRand();
      G4double* running = new G4double[nEnergies];
      running[0] = 0;
      G4double weighted = 0;
      for (i = 1; i < nEnergies; i++) {
        running[i]=running[i-1];
        if (fCache.Get()->remaining_energy >= theAngular[i].GetLabel() ) {
          running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                                              theAngular[i-1].GetLabel(),
                                              theAngular[i].GetLabel(),
                                              theAngular[i-1].GetValue(0),
                                              theAngular[i].GetValue(0) );
          weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                                                    theAngular[i-1].GetLabel(),
                                                    theAngular[i].GetLabel(),
                                                    theAngular[i-1].GetValue(0),
                                                    theAngular[i].GetValue(0));
        }
      }

      // Cache the mean energy in this distribution
      if (nEnergies == 1 || running[nEnergies-1] == 0) { 
        fCache.Get()->currentMeanEnergy = 0.0;
      } else { 
        fCache.Get()->currentMeanEnergy = weighted/running[nEnergies-1];
      }
         
      if (nEnergies == 1) it = 0; 
      if (running[nEnergies-1] != 0) {
        for (i = 1; i < nEnergies; i++) {
          it = i;
          if (random < running [i]/running[nEnergies-1] ) break;
        }
      }

      if (running[nEnergies-1] == 0) it = 0;
      if (it < nDiscreteEnergies || it == 0) {
        if (it == 0) {
          fsEnergy = theAngular[it].GetLabel();
          G4ParticleHPLegendreStore theStore(1);
          theStore.Init(0,fsEnergy,nAngularParameters);
          for (i = 0; i < nAngularParameters; i++) {
            theStore.SetCoeff(0, i, theAngular[it].GetValue(i) );
          }
          // use it to sample.
          cosTh = theStore.SampleMax(fsEnergy);

        } else {
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
          for (i = 0; i < nAngularParameters; i++) {
            theStore.SetCoeff(0, i, theAngular[it-1].GetValue(i) );
            theStore.SetCoeff(1, i, theAngular[it].GetValue(i) );
          }
          // use it to sample.
          theStore.SetManager(theManager);
          cosTh = theStore.SampleMax(fsEnergy);
        }

      } else { // continuum contribution
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
        for (i = 0; i < nAngularParameters; i++) {
          theStore.SetCoeff(0, i, theAngular[it-1].GetValue(i));
          theStore.SetCoeff(1, i, theAngular[it].GetValue(i));
        }
        // use it to sample.
        cosTh = theStore.SampleMax(fsEnergy);
      }
      delete [] running;

      // The remaining energy needs to be lowered by the photon energy in
      // *any* case.  Otherwise additional photons with too much energy will be
      // produced - therefore the  adjustResult condition has been removed

      fCache.Get()->remaining_energy -= fsEnergy;
      // end if (nDiscreteEnergies != 0)
    }
    // end of (angularRep == 1) branch

  } else if (angularRep == 2) {
    // first get the energy (already the right for this incoming energy)
    G4int j;
    G4double* running = new G4double[nEnergies];
    running[0] = 0;
    G4double weighted = 0;
    for (j = 1; j < nEnergies; j++) {
      if (j != 0) running[j] = running[j-1];
      running[j] += theInt.GetBinIntegral(theManager.GetScheme(j-1),
                           theAngular[j-1].GetLabel(), theAngular[j].GetLabel(),
                           theAngular[j-1].GetValue(0), theAngular[j].GetValue(0));
      weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(j-1),
                           theAngular[j-1].GetLabel(), theAngular[j].GetLabel(),
                           theAngular[j-1].GetValue(0), theAngular[j].GetValue(0));
    }

    // Cache the mean energy in this distribution
    if (nEnergies == 1)
      fCache.Get()->currentMeanEnergy = 0.0;
    else
      fCache.Get()->currentMeanEnergy = weighted/running[nEnergies-1];
      
    G4int itt(0);
    G4double randkal = G4UniformRand();
    for (j = 1; j < nEnergies; j++) {
      itt = j;
      if (randkal < running[j]/running[nEnergies-1] ) break;
    }
      
    // Interpolate the secondary energy
    G4double x, x1, x2, y1, y2;
    if (itt == 0) itt = 1;
    x = randkal*running[nEnergies-1];
    x1 = running[itt-1];
    x2 = running[itt];
    G4double compoundFraction;
    // interpolate energy
    y1 = theAngular[itt-1].GetLabel();
    y2 = theAngular[itt].GetLabel();
    fsEnergy = theInt.Interpolate(theManager.GetInverseScheme(itt-1), 
                                  x, x1, x2, y1, y2);

    // For theta, interpolate the compoundFractions
    G4double cLow = theAngular[itt-1].GetValue(1);
    G4double cHigh = theAngular[itt].GetValue(1);
    compoundFraction = theInt.Interpolate(theManager.GetScheme(itt),
                                          fsEnergy, y1, y2, cLow, cHigh);

    if (compoundFraction > 1.0) compoundFraction = 1.0;  // Protection against unphysical interpolation

    delete [] running;
      
    // get cosTh
    G4double incidentEnergy = anEnergy;
    G4double incidentMass = theProjectile->GetPDGMass();
    G4double productEnergy = fsEnergy;
    G4double productMass = result->GetMass();
    G4int targetZ = G4int(fCache.Get()->theTargetCode/1000);
    G4int targetA = G4int(fCache.Get()->theTargetCode-1000*targetZ);

    // To correspond to natural composition (-nat-) data files. 
    if (targetA == 0) 
      targetA = G4int ( fCache.Get()->theTarget->GetMass()/amu_c2 + 0.5 );
    G4double targetMass = fCache.Get()->theTarget->GetMass();
    G4int incidentA=G4int(incidentMass/amu_c2 + 0.5 );
    G4int incidentZ=G4int(theProjectile->GetPDGCharge()+ 0.5 );
    G4int residualA = targetA+incidentA-A;
    G4int residualZ = targetZ+incidentZ-Z;
    G4double residualMass =G4NucleiProperties::GetNuclearMass( residualA , residualZ );
    G4ParticleHPKallbachMannSyst theKallbach(compoundFraction,
                                             incidentEnergy, incidentMass,
                                             productEnergy, productMass,
                                             residualMass, residualA, residualZ,
                                             targetMass, targetA, targetZ,
                                             incidentA,incidentZ,A,Z);
    cosTh = theKallbach.Sample(anEnergy);
    // end (angularRep == 2) branch

  } else if (angularRep > 10 && angularRep < 16) {
    G4double random = G4UniformRand();
    G4double* running = new G4double[nEnergies];
    running[0]=0;      
    G4double weighted = 0;
    for (i = 1; i < nEnergies; i++) {
      if (i != 0) running[i] = running[i-1];
      running[i] += theInt.GetBinIntegral(theManager.GetScheme(i-1),
                           theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                           theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
      weighted += theInt.GetWeightedBinIntegral(theManager.GetScheme(i-1),
                           theAngular[i-1].GetLabel(), theAngular[i].GetLabel(),
                           theAngular[i-1].GetValue(0), theAngular[i].GetValue(0));
    }

    // Cache the mean energy in this distribution
    if (nEnergies == 1)  
      fCache.Get()->currentMeanEnergy = 0.0;
    else
      fCache.Get()->currentMeanEnergy = weighted/running[nEnergies-1];
      
    if (nEnergies == 1) it = 0; 
    for (i = 1; i < nEnergies; i++) {
      it = i;
      if(random<running[i]/running[nEnergies-1]) break;
    }

    if (it < nDiscreteEnergies || it == 0) {
      if (it == 0) {
        fsEnergy = theAngular[0].GetLabel();          
        G4ParticleHPVector theStore; 
        G4int aCounter = 0;
        for (G4int j=1; j<nAngularParameters; j+=2) {
          theStore.SetX(aCounter, theAngular[0].GetValue(j));
          theStore.SetY(aCounter, theAngular[0].GetValue(j+1));
          aCounter++;	    
        }
        G4InterpolationManager aMan;
        aMan.Init(angularRep-10, nAngularParameters-1);
        theStore.SetInterpolationManager(aMan);
        cosTh = theStore.Sample();
      } else {
        fsEnergy = theAngular[it].GetLabel();
        G4ParticleHPVector theStore; 
        G4InterpolationManager aMan;
        aMan.Init(angularRep-10, nAngularParameters-1);
        theStore.SetInterpolationManager(aMan); // Store interpolates f(costh)
        G4InterpolationScheme currentScheme = theManager.GetInverseScheme(it);
        G4int aCounter = 0;
        for (G4int j=1; j<nAngularParameters; j+=2) {
          theStore.SetX(aCounter, theAngular[it].GetValue(j));
          theStore.SetY(aCounter,
                        theInt.Interpolate(currentScheme, random,
                                           running[it-1]/running[nEnergies-1],
                                           running[it]/running[nEnergies-1],
                                           theAngular[it-1].GetValue(j+1),
                                           theAngular[it].GetValue(j+1) ) );
          aCounter++;	    
        }
        cosTh = theStore.Sample();
      }

    } else {
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

      G4int j;
      for (i = 0, j = 1; i < nAngularParameters; i++, j+=2) {
        theBuff1.SetX(i, theAngular[it-1].GetValue(j));
        theBuff1.SetY(i, theAngular[it-1].GetValue(j+1));
        theBuff2.SetX(i, theAngular[it].GetValue(j));
        theBuff2.SetY(i, theAngular[it].GetValue(j+1));
      }

      G4ParticleHPVector theStore;
      theStore.SetInterpolationManager(aMan); // Store interpolates f(costh)        
      x1 = y1;
      x2 = y2;
      G4double x, y;
      for (i = 0; i < theBuff1.GetVectorLength(); i++) {
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

  } else {
    throw G4HadronicException(__FILE__, __LINE__,
                   "G4ParticleHPContAngularPar::Sample: Unknown angular representation");
  }

  result->SetKineticEnergy(fsEnergy);
  G4double phi = twopi*G4UniformRand();
  G4double theta = std::acos(cosTh);
  G4double sinth = std::sin(theta);
  G4double mtot = result->GetTotalMomentum();
  G4ThreeVector tempVector(mtot*sinth*std::cos(phi), mtot*sinth*std::sin(phi), mtot*std::cos(theta) );
  result->SetMomentum(tempVector);
  return result;
}


#define MERGE_NEW

void G4ParticleHPContAngularPar::PrepareTableInterpolation()
{
  // Discrete energies: store own energies in a map for faster searching
  //
  // The data files sometimes have identical discrete energies (likely typos)
  // which would lead to overwriting the already existing index and hence
  // creating a hole in the lookup table.
  // No attempt is made here to correct for the energies - rather an epsilon
  // is subtracted from the energy in order to uniquely identify the line

  for (G4int ie = 0; ie < nDiscreteEnergies; ie++) {
    // check if energy is already present and subtract epsilon if that's the case
    G4double myE = theAngular[ie].GetLabel();
    while (theDiscreteEnergiesOwn.find(myE) != theDiscreteEnergiesOwn.end() ) {
      myE -= 1e-6;
    }
    theDiscreteEnergiesOwn[myE] = ie;
  }

  /*
   * the approach here makes no sense. It would work only for two sets that
   * have identical min and max energy. If the 2 sets differ in min, max or
   * both, the energy inserted would be normalized to its original set but
   * interpreted with the new - which is not correct.
   *
   * Disable the code for now and simply return ...
   */

  return;

  /*
   *
   
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
  }
 
  G4int nEnergiesPrev = angParPrev->GetNEnergies();
  G4double minEnerPrev = angParPrev->GetMinEner();
  G4double maxEnerPrev = angParPrev->GetMaxEner();
  for(ie=nDiscreteEnergiesPrev; ie<nEnergiesPrev; ie++) {
    G4double ener = angParPrev->theAngular[ie].GetLabel();
    G4double enerT = (ener-minEnerPrev)/(maxEnerPrev-minEnerPrev);
    theEnergiesTransformed.insert(enerT);
  }
  // add the maximum energy
  //theEnergiesTransformed.insert(1.);

  *
  */
}


void
G4ParticleHPContAngularPar::BuildByInterpolation(G4double anEnergy,
                                                 G4InterpolationScheme aScheme, 
                                                 G4ParticleHPContAngularPar& angpar1, 
                                                 G4ParticleHPContAngularPar& angpar2) 
{
  G4int ie, ie1, ie2, ie1Prev, ie2Prev;
  // Only rebuild the interpolation table if there is a new interaction.
  // For several subsequent samplings of final state particles in the same
  // interaction the existing table should be used
  if (fCache.Get()->fresh != true) return;

  // Make copies of angpar1 and angpar2. Since these are given by reference
  // it can not be excluded that one of them is "this". Hence this code uses
  // potentially the old "this" for creating the new this - which leads to
  // memory corruption if the old is not stored as separarte object for lookup
  const G4ParticleHPContAngularPar copyAngpar1(angpar1),copyAngpar2(angpar2);
  
  nAngularParameters = copyAngpar1.nAngularParameters;
  theManager = copyAngpar1.theManager;
  theEnergy = anEnergy;
  theMinEner = DBL_MAX; // min and max will be re-calculated after interpolation
  theMaxEner = -DBL_MAX;

  // The two discrete sets must be merged. A vector holds the temporary data to
  // be copied to the array in the end.  Since the G4ParticleHPList class
  // contains pointers, can't simply assign elements of this type. Each member
  // needs to call the explicit Set() method instead.

  // First, average probabilities for those lines that are in both sets
  const std::map<G4double,G4int> discEnerOwn1 = copyAngpar1.GetDiscreteEnergiesOwn();
  const std::map<G4double,G4int> discEnerOwn2 = copyAngpar2.GetDiscreteEnergiesOwn();
  std::map<G4double,G4int>::const_iterator itedeo1;
  std::map<G4double,G4int>::const_iterator itedeo2;
  std::vector<G4ParticleHPList*> vAngular(discEnerOwn1.size() );
  G4double discEner1;
  for (itedeo1 = discEnerOwn1.cbegin(); itedeo1 != discEnerOwn1.cend(); ++itedeo1) {
    discEner1 = itedeo1->first;
    if (discEner1 < theMinEner) {
      theMinEner = discEner1;
    }
    if (discEner1 > theMaxEner) {
      theMaxEner = discEner1;
    }
    ie1 = itedeo1->second;
    itedeo2 = discEnerOwn2.find(discEner1);
    if( itedeo2 == discEnerOwn2.cend() ) {
      ie2 = -1;
    } else {
      ie2 = itedeo2->second;
    }
    vAngular[ie1] = new G4ParticleHPList();
    vAngular[ie1]->SetLabel(copyAngpar1.theAngular[ie1].GetLabel());
    G4double val1, val2;
    for (G4int ip = 0; ip < nAngularParameters; ++ip) {
      val1 = copyAngpar1.theAngular[ie1].GetValue(ip);
      if (ie2 != -1) {
        val2 = copyAngpar2.theAngular[ie2].GetValue(ip);
      } else {
        val2 = 0.;
      }
      G4double value = theInt.Interpolate(aScheme, anEnergy, 
					  copyAngpar1.theEnergy, copyAngpar2.theEnergy,
                                          val1, val2);
      vAngular[ie1]->SetValue(ip, value);
    }
  } // itedeo1 loop

  // Add the ones in set2 but not in set1
  std::vector<G4ParticleHPList*>::const_iterator itv;
  G4double discEner2;
  for (itedeo2 = discEnerOwn2.cbegin(); itedeo2 != discEnerOwn2.cend(); ++itedeo2) {
    discEner2 = itedeo2->first;
    ie2 = itedeo2->second;
    G4bool notFound = true;
    itedeo1 = discEnerOwn1.find(discEner2);
    if (itedeo1 != discEnerOwn1.cend() ) {
      notFound = false;
    }
    if (notFound) {
      // not yet in list
      if (discEner2 < theMinEner) {
        theMinEner = discEner2;
      }
      if (discEner2 > theMaxEner) {
        theMaxEner = discEner2;
      }
      // find position to insert
      G4bool isInserted = false;
      ie = 0;
      for (itv = vAngular.cbegin(); itv != vAngular.cend(); ++itv,++ie) {
        if (discEner2 > (*itv)->GetLabel() ) {
          itv = vAngular.insert(itv, new G4ParticleHPList);
          (*itv)->SetLabel(copyAngpar2.theAngular[ie2].GetLabel());
          isInserted = true;
          break;
        }
      }
      if (!isInserted) {
        ie=(G4int)vAngular.size();
        vAngular.push_back(new G4ParticleHPList);
        vAngular[ie]->SetLabel(copyAngpar2.theAngular[ie2].GetLabel());
        isInserted = true;
      }

      G4double val1, val2;
      for (G4int ip = 0; ip < nAngularParameters; ++ip) {
        val1 = 0;
        val2 = copyAngpar2.theAngular[ie2].GetValue(ip);
        G4double value = theInt.Interpolate(aScheme, anEnergy, 
                                            copyAngpar1.theEnergy,
                                            copyAngpar2.theEnergy,
                                            val1, val2);
        vAngular[ie]->SetValue(ip, value);
      }
    }  // end if(notFound) 
  } // end loop on itedeo2
  
  // Store new discrete list 
  nDiscreteEnergies = (G4int)vAngular.size();
  if (theAngular != 0) delete [] theAngular;
  theAngular = 0;
  if (nDiscreteEnergies > 0) {
    theAngular = new G4ParticleHPList [nDiscreteEnergies];
  }
  theDiscreteEnergiesOwn.clear();
  theDiscreteEnergies.clear();
  for (ie = 0; ie < nDiscreteEnergies; ++ie) {
    theAngular[ie].SetLabel(vAngular[ie]->GetLabel() );
    for (G4int ip = 0; ip < nAngularParameters; ++ip) {
      theAngular[ie].SetValue(ip, vAngular[ie]->GetValue(ip));
    }
    theDiscreteEnergiesOwn[theAngular[ie].GetLabel()] = ie;
    theDiscreteEnergies.insert(theAngular[ie].GetLabel());
  }

  // The continuous energies need to be made from scratch like the discrete
  // ones. Therefore the re-assignemnt of theAngular needs to be done
  // after the continuous energy set is also finalized. Only then the
  // total number of nEnergies is known and the array can be allocated.

  // Get minimum and maximum energy interpolating
  // Don't use theMinEner or theMaxEner here, since the transformed energies
  // need the interpolated range from the original Angpar 
  G4double interMinEner = copyAngpar1.GetMinEner() + (theEnergy-copyAngpar1.GetEnergy() )
                                                   * (copyAngpar2.GetMinEner() - copyAngpar1.GetMinEner() )
                                                   / (copyAngpar2.GetEnergy()-copyAngpar1.GetEnergy() );
  G4double interMaxEner = copyAngpar1.GetMaxEner() + (theEnergy-copyAngpar1.GetEnergy() )
                                                   * (copyAngpar2.GetMaxEner()-copyAngpar1.GetMaxEner() )
                                                   / (copyAngpar2.GetEnergy()-copyAngpar1.GetEnergy() );

  // Loop to energies of new set
  theEnergiesTransformed.clear();

  G4int nEnergies1 = copyAngpar1.GetNEnergies();
  G4int nDiscreteEnergies1 = copyAngpar1.GetNDiscreteEnergies();
  G4double minEner1 = copyAngpar1.GetMinEner();
  G4double maxEner1 = copyAngpar1.GetMaxEner();
  G4int nEnergies2 = copyAngpar2.GetNEnergies();
  G4int nDiscreteEnergies2 = copyAngpar2.GetNDiscreteEnergies();
  G4double minEner2 = copyAngpar2.GetMinEner();
  G4double maxEner2 = copyAngpar2.GetMaxEner();

  // First build the list of transformed energies normalized
  // to the new min max by assuming that the min-max range of
  // each set would be scalable to the new, interpolated min
  // max range

  G4double e1(0.);
  G4double eTNorm1(0.);
  for (ie1 = nDiscreteEnergies1; ie1 < nEnergies1; ++ie1) {
    e1 = copyAngpar1.theAngular[ie1].GetLabel();
    eTNorm1 = (e1 - minEner1);
    if (maxEner1 != minEner1) eTNorm1 /= (maxEner1-minEner1);
    if (eTNorm1 >= 0 && eTNorm1 <= 1) theEnergiesTransformed.insert(eTNorm1);
  }

  G4double e2(0.);
  G4double eTNorm2(0.);
  for (ie2 = nDiscreteEnergies2; ie2 < nEnergies2; ++ie2) {
    e2 = copyAngpar2.theAngular[ie2].GetLabel();
    eTNorm2 = (e2 - minEner2);
    if (maxEner2 != minEner2) eTNorm2 /= (maxEner2-minEner2);
    if (eTNorm2 >= 0 && eTNorm2 <= 1) theEnergiesTransformed.insert(eTNorm2);
  }

  // Now the list of energies is complete
  nEnergies = nDiscreteEnergies+(G4int)theEnergiesTransformed.size();
  
  // Create final array of angular parameters
  G4ParticleHPList* theNewAngular = new G4ParticleHPList [nEnergies];

  // Copy discrete energies and interpolated parameters to new array
  
  if (theAngular != 0) {
    for (ie = 0; ie < nDiscreteEnergies; ++ie) { 
      theNewAngular[ie].SetLabel(theAngular[ie].GetLabel());
      for (G4int ip = 0; ip < nAngularParameters; ++ip) {
        theNewAngular[ie].SetValue(ip,theAngular[ie].GetValue(ip));
      }
    }
    delete [] theAngular;
  }
  theAngular = theNewAngular;
  
  // Interpolate the continuous energies for new array
  std::set<G4double>::const_iterator iteet = theEnergiesTransformed.begin();

  G4double e1Interp(0.);
  G4double e2Interp(0.);
  for (ie = nDiscreteEnergies; ie < nEnergies; ++ie, ++iteet) { 
    G4double eT = (*iteet);

    //--- Use eT1 = eT: Get energy and parameters of copyAngpar1 for this eT
    e1Interp = (maxEner1 - minEner1) * eT + minEner1;
    //----- Get parameter value corresponding to this e1Interp
    for (ie1 = nDiscreteEnergies1; ie1 < nEnergies1; ++ie1) {
      if ((copyAngpar1.theAngular[ie1].GetLabel() - e1Interp) > 1.E-10*e1Interp) break;
    }
    ie1Prev = ie1 - 1;
    if (ie1 == 0) ++ie1Prev; 
    if (ie1 == nEnergies1) {
      ie1--;
      ie1Prev = ie1;
    }

    //--- Use eT2 = eT: Get energy and parameters of copyAngpar2 for this eT
    e2Interp = (maxEner2-minEner2) * eT + minEner2;
    //----- Get parameter value corresponding to this e2Interp
    for (ie2 = nDiscreteEnergies2; ie2 < nEnergies2; ++ie2) {
      if ((copyAngpar2.theAngular[ie2].GetLabel() - e2Interp) > 1.E-10*e2Interp) break;
    }
    ie2Prev = ie2 - 1;
    if (ie2 == 0) ++ie2Prev; 
    if (ie2 == nEnergies2) {
      ie2--;
      ie2Prev = ie2;
    }

    //---- Energy corresponding to energy transformed    
    G4double eN = (interMaxEner-interMinEner) * eT + interMinEner;
    
    theAngular[ie].SetLabel(eN);
    if (eN < theMinEner) {
      theMinEner = eN;
    }
    if (eN > theMaxEner) {
      theMaxEner = eN;
    }
   
    G4double val1(0.);
    G4double val2(0.);
    G4double value(0.);
    for (G4int ip = 0; ip < nAngularParameters; ++ip) {
      val1 = theInt.Interpolate2(theManager.GetScheme(ie),
                                 e1Interp,
                                 copyAngpar1.theAngular[ie1Prev].GetLabel(),
                                 copyAngpar1.theAngular[ie1].GetLabel(),
                                 copyAngpar1.theAngular[ie1Prev].GetValue(ip),
                                 copyAngpar1.theAngular[ie1].GetValue(ip)) * (maxEner1-minEner1);  
      val2 = theInt.Interpolate2(theManager.GetScheme(ie),
                                 e2Interp,
                                 copyAngpar2.theAngular[ie2Prev].GetLabel(),
                                 copyAngpar2.theAngular[ie2].GetLabel(),
                                 copyAngpar2.theAngular[ie2Prev].GetValue(ip),
                                 copyAngpar2.theAngular[ie2].GetValue(ip)) * (maxEner2-minEner2);

      value = theInt.Interpolate(aScheme, anEnergy, 
                                 copyAngpar1.theEnergy, copyAngpar2.theEnergy,
                                 val1, val2);
      if (interMaxEner != interMinEner) {
        value /= (interMaxEner-interMinEner); 
      } else if (value != 0) {
         throw G4HadronicException(__FILE__, __LINE__,
                                   "G4ParticleHPContAngularPar::PrepareTableInterpolation interMaxEner == interMinEner and  value != 0.");
      }
      theAngular[ie].SetValue(ip, value);
    }
  }  // end loop on nDiscreteEnergies

  for (itv = vAngular.cbegin(); itv != vAngular.cend(); ++itv) delete (*itv);
  
}


void G4ParticleHPContAngularPar::Dump() const
{
  G4cout << theEnergy << " " << nEnergies << " " << nDiscreteEnergies
         << " " << nAngularParameters << G4endl;

  for (G4int ii = 0; ii < nEnergies; ii++) theAngular[ii].Dump();
}

