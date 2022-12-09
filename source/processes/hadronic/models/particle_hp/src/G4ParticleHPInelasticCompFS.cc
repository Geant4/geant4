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
// particle_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 070606 bug fix and migrate to enable to Partial cases by T. Koi 
// 080603 bug fix for Hadron Hyper News #932 by T. Koi 
// 080612 bug fix contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #4,6
// 080717 bug fix of calculation of residual momentum by T. Koi
// 080801 protect negative available energy by T. Koi
//        introduce theNDLDataA,Z which has A and Z of NDL data by T. Koi
// 081024 G4NucleiPropertiesTable:: to G4NucleiProperties::
// 090514 Fix bug in IC electron emission case 
//        Contribution from Chao Zhang (Chao.Zhang@usd.edu) and Dongming Mei(Dongming.Mei@usd.edu)
// 100406 "nothingWasKnownOnHadron=1" then sample mu isotropic in CM 
//        add two_body_reaction
// 100909 add safty 
// 101111 add safty for _nat_ data case in Binary reaction, but break conservation  
// 110430 add Reaction Q value and break up flag (MF3::QI and LR)
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPInelasticCompFS.hh"
#include "G4ParticleHPManager.hh"
#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4ParticleHPDataUsed.hh"
#include "G4IonTable.hh"
#include "G4Pow.hh"
#include "G4SystemOfUnits.hh"

#include "G4NRESP71M03.hh" // nresp71_m03.hh and nresp71_m02.hh are alike. The only difference between m02 and m03 is in the total carbon cross section that is properly included in the latter. These data are not used in nresp71_m0*.hh.

// June-2019 - E. Mendoza - re-build "two_body_reaction", to be used by incident charged particles (now isotropic emission in the CMS). Also restrict nresp use below 20 MeV (for future developments). Add photon emission when no data available.

void G4ParticleHPInelasticCompFS::InitGammas(G4double AR, G4double ZR)
{
   std::ostringstream ost;
   ost <<gammaPath<<"z"<<ZR<<".a"<<AR;
   G4String aName = ost.str();
   std::ifstream from(aName, std::ios::in);

   if(!from) return; // no data found for this isotope
   std::ifstream theGammaData(aName, std::ios::in);
    
   theGammas.Init(theGammaData);
}

void G4ParticleHPInelasticCompFS::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & aFSType, G4ParticleDefinition*)
{
  gammaPath = "/Inelastic/Gammas/";  //only in neutron data base 
    if(!G4FindDataDir("G4NEUTRONHPDATA"))
       throw G4HadronicException(__FILE__, __LINE__, "Please setenv G4NEUTRONHPDATA to point to the neutron cross-section files where Inelastic/Gammas data is found.");
  G4String tBase = G4FindDataDir("G4NEUTRONHPDATA");
  gammaPath = tBase+gammaPath;
  G4String tString = dirName;
  G4bool dbool;
  G4ParticleHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, tString, aFSType, dbool);
  G4String filename = aFile.GetName();
#ifdef G4PHPDEBUG
  if( std::getenv("G4ParticleHPDebug") ) G4cout << " G4ParticleHPInelasticCompFS::Init FILE " << filename << G4endl;
#endif

  SetAZMs( A, Z, M, aFile ); 

  if ( !dbool || ( Z<2.5 && ( std::abs(theNDLDataZ - Z)>0.0001 || std::abs(theNDLDataA - A)>0.0001)) )
  {
#ifdef G4PHPDEBUG
    if(std::getenv("G4ParticleHPDebug_NamesLogging")) G4cout << "Skipped = "<< filename <<" "<<A<<" "<<Z<<G4endl;
#endif
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    return;
  }
  std::istringstream theData(std::ios::in);
  G4ParticleHPManager::GetInstance()->GetDataStream(filename,theData);
  if(!theData) //"!" is a operator of ios
  {
    hasAnyData = false;
    hasFSData = false; 
    hasXsec = false;
    //    theData.close();
    return;
  }
  // here we go
  G4int infoType, dataType, dummy;
  G4int sfType, it;
  hasFSData = false; 
  while (theData >> infoType) // Loop checking, 11.05.2015, T. Koi
  {
    hasFSData = true; 
    theData >> dataType;
    theData >> sfType >> dummy;
    it = 50;
    if(sfType>=600||(sfType<100&&sfType>=50)) it = sfType%50;
    if(dataType==3) 
    {
      //theData >> dummy >> dummy;
      //TK110430
      // QI and LR introudced since G4NDL3.15
      G4double dqi;
      G4int ilr;
      theData >> dqi >> ilr;

      QI[ it ] = dqi*CLHEP::eV;
      LR[ it ] = ilr;
      theXsection[it] = new G4ParticleHPVector;
      G4int total;
      theData >> total;
      theXsection[it]->Init(theData, total, CLHEP::eV);
      //std::cout << theXsection[it]->GetXsec(1*MeV) << std::endl;
    }
    else if(dataType==4)
    {
      theAngularDistribution[it] = new G4ParticleHPAngular;
      theAngularDistribution[it]->Init(theData);
    }
    else if(dataType==5)
    {
      theEnergyDistribution[it] = new G4ParticleHPEnergyDistribution;
      theEnergyDistribution[it]->Init(theData); 
    }
    else if(dataType==6)
    {
      theEnergyAngData[it] = new G4ParticleHPEnAngCorrelation(theProjectile);
      //      G4cout << this << " CompFS theEnergyAngData " << it << theEnergyAngData[it] << G4endl; //GDEB
      theEnergyAngData[it]->Init(theData);
    }
    else if(dataType==12)
    {
      theFinalStatePhotons[it] = new G4ParticleHPPhotonDist;
      theFinalStatePhotons[it]->InitMean(theData);
    }
    else if(dataType==13)
    {
      theFinalStatePhotons[it] = new G4ParticleHPPhotonDist;
      theFinalStatePhotons[it]->InitPartials(theData, theXsection[50]);
    }
    else if(dataType==14)
    {
      theFinalStatePhotons[it]->InitAngular(theData);
    }
    else if(dataType==15)
    {
      theFinalStatePhotons[it]->InitEnergies(theData);
    }
    else
    {
      throw G4HadronicException(__FILE__, __LINE__, "Data-type unknown to G4ParticleHPInelasticCompFS");
    }
  }
  //  theData.close();
}

G4int G4ParticleHPInelasticCompFS::SelectExitChannel(G4double eKinetic)
{
  G4double running[50];
  running[0] = 0;
  G4int i;
  for(i=0; i<50; ++i)
  {
    if(i!=0) running[i]=running[i-1];
    if(theXsection[i] != 0) 
    {
      running[i] += std::max(0., theXsection[i]->GetXsec(eKinetic));
    }
  }
  G4double random = G4UniformRand();
  G4double sum = running[49];
  G4int it = 50;
  if(0!=sum)
  {
    G4int i0;
    for(i0=0; i0<50; ++i0)
    {
      it = i0;
      if(random < running[i0]/sum) break;
    }
  }
  return it;
}

// n,p,d,t,he3,a
void G4ParticleHPInelasticCompFS::CompositeApply(const G4HadProjectile& theTrack,
                                                 G4ParticleDefinition* aDefinition)
{

    // prepare neutron
    if ( theResult.Get() == nullptr ) theResult.Put( new G4HadFinalState );
    theResult.Get()->Clear();
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4HadProjectile *hadProjectile = &theTrack;
    G4ReactionProduct incidReactionProduct( const_cast<G4ParticleDefinition *>(hadProjectile->GetDefinition()) ); // incidReactionProduct
    incidReactionProduct.SetMomentum( hadProjectile->Get4Momentum().vect() );
    incidReactionProduct.SetKineticEnergy( eKinetic );

    // prepare target
    G4int i;
    for(i=0; i<50; ++i)
    { if(theXsection[i] != 0) { break; } } 

    G4double targetMass=0;
    G4double eps = 0.0001;
    targetMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(theBaseA+eps), static_cast<G4int>(theBaseZ+eps));
#ifdef G4PHPDEBUG
    if( std::getenv("G4ParticleHPDebug"))
      G4cout <<this <<" G4ParticleHPInelasticCompFS::CompositeApply A "
             <<theBaseA <<" Z " <<theBaseZ <<" incident "
             <<hadProjectile->GetDefinition()->GetParticleName() <<G4endl;
#endif
    G4ReactionProduct theTarget; 
    G4Nucleus aNucleus;
    //G4ThreeVector neuVelo = (1./hadProjectile->GetDefinition()->GetPDGMass())*incidReactionProduct.GetMomentum();
    //theTarget = aNucleus.GetBiasedThermalNucleus( targetMass/hadProjectile->GetDefinition()->GetPDGMass() , neuVelo, theTrack.GetMaterial()->GetTemperature());
    //G4Nucleus::GetBiasedThermalNucleus requests normalization of mass and velocity in neutron mass
    G4ThreeVector neuVelo = ( 1./G4Neutron::Neutron()->GetPDGMass() )*incidReactionProduct.GetMomentum();
    theTarget = aNucleus.GetBiasedThermalNucleus( targetMass/G4Neutron::Neutron()->GetPDGMass()
                                                , neuVelo, theTrack.GetMaterial()->GetTemperature() );
    
    theTarget.SetDefinition( G4IonTable::GetIonTable()->GetIon( G4int(theBaseZ), G4int(theBaseA) , 0.0 ) );  //XX

    // prepare the residual mass
    G4double residualMass=0;
    G4double residualZ = theBaseZ + theProjectile->GetPDGCharge() - aDefinition->GetPDGCharge();
    G4double residualA = theBaseA + theProjectile->GetBaryonNumber() - aDefinition->GetBaryonNumber();
    residualMass = G4NucleiProperties::GetNuclearMass(static_cast<G4int>(residualA+eps), static_cast<G4int>(residualZ+eps));

    // prepare energy in target rest frame
    G4ReactionProduct boosted;
    boosted.Lorentz(incidReactionProduct, theTarget);
    eKinetic = boosted.GetKineticEnergy();
  
    // select exit channel for composite FS class.
    G4int it = SelectExitChannel( eKinetic );
  
    //E. Mendoza (2018) -- to use JENDL/AN-2005
    if(theEnergyDistribution[it]==0 && theAngularDistribution[it]==0 && theEnergyAngData[it]==0){
      if(theEnergyDistribution[50]!=0 || theAngularDistribution[50]!=0 || theEnergyAngData[50]!=0){
        it=50;
      }
    }
 
    // set target and neutron in the relevant exit channel
    InitDistributionInitialState(incidReactionProduct, theTarget, it);    

   //---------------------------------------------------------------------//
   // Hook for NRESP71MODEL
   if ( G4ParticleHPManager::GetInstance()->GetUseNRESP71Model() && eKinetic<20*MeV) {
      if ( (G4int)(theBaseZ+0.1) == 6 ) // If the reaction is with Carbon...
      {
         if ( theProjectile == G4Neutron::Definition() ) {
            if ( use_nresp71_model( aDefinition , it , theTarget , boosted ) ) return;
         }
      }
   }
   //---------------------------------------------------------------------//

    G4ReactionProductVector * thePhotons = 0;
    G4ReactionProductVector * theParticles = 0;
    G4ReactionProduct aHadron;
    aHadron.SetDefinition(aDefinition); // what if only cross-sections exist ==> Na 23 11 @@@@    
    G4double availableEnergy = incidReactionProduct.GetKineticEnergy() + incidReactionProduct.GetMass() - aHadron.GetMass() +
                             (targetMass - residualMass);

    if ( availableEnergy < 0 )
    {
       availableEnergy = 0; 
    }
    G4int nothingWasKnownOnHadron = 0;
    G4int dummy;
    G4double eGamm = 0;
    G4int iLevel=it-1;

    // without photon has it = 0
    if( 50 == it ) 
    {
      // Excitation level is not determined
      iLevel=-1;
      aHadron.SetKineticEnergy(availableEnergy*residualMass/
                               (aHadron.GetMass()+residualMass));

      //aHadron.SetMomentum(incidReactionProduct.GetMomentum()*(1./incidReactionProduct.GetTotalMomentum())*
      //                  std::sqrt(aHadron.GetTotalEnergy()*aHadron.GetTotalEnergy()-
      //                            aHadron.GetMass()*aHadron.GetMass()));

      //TK add safty 100909
      G4double p2 = ( aHadron.GetTotalEnergy()*aHadron.GetTotalEnergy() - aHadron.GetMass()*aHadron.GetMass() );
      G4double p = 0.0;
      if ( p2 > 0.0 )
        p = std::sqrt( p );
      aHadron.SetMomentum(incidReactionProduct.GetMomentum()*(1./incidReactionProduct.GetTotalMomentum())*p );
    }
    else
    {
      while ( iLevel!=-1 && theGammas.GetLevel(iLevel) == 0 ) { iLevel--; } // Loop checking, 11.05.2015, T. Koi
    }

    if ( theAngularDistribution[it] != 0 ) // MF4
    {
      if(theEnergyDistribution[it]!=0) // MF5
      {
	//************************************************************
	/*
        aHadron.SetKineticEnergy(theEnergyDistribution[it]->Sample(eKinetic, dummy));
        G4double eSecN = aHadron.GetKineticEnergy();
	*/
	//************************************************************
	//EMendoza --> maximum allowable energy should be taken into account.
        G4double dqi = 0.0;
        if ( QI[it] < 0 || 849 < QI[it] ) dqi = QI[it]; //For backword compatibility QI introduced since G4NDL3.15
	G4double MaxEne=eKinetic+dqi;
	G4double eSecN=0.;

        G4int icounter=0;
        G4int icounter_max=1024;
        do {
           ++icounter;
           if ( icounter > icounter_max ) {
	      G4cout << "Loop-counter exceeded the threshold value at " << __LINE__ << "th line of " << __FILE__ << "." << G4endl;
              break;
           }
	   eSecN=theEnergyDistribution[it]->Sample(eKinetic, dummy);
	} while(eSecN>MaxEne); // Loop checking, 11.05.2015, T. Koi
	aHadron.SetKineticEnergy(eSecN);
	//************************************************************
        eGamm = eKinetic-eSecN;
        for(iLevel=theGammas.GetNumberOfLevels()-1; iLevel>=0; --iLevel)
        {
          if(theGammas.GetLevelEnergy(iLevel)<eGamm) break;
        }
        G4double random = 2*G4UniformRand();
        iLevel+=G4int(random);
        if(iLevel>theGammas.GetNumberOfLevels()-1)iLevel = theGammas.GetNumberOfLevels()-1;
      }
      else
      {
        G4double eExcitation = 0;
        if(iLevel>=0) eExcitation = theGammas.GetLevel(iLevel)->GetLevelEnergy();    
        while (eKinetic-eExcitation < 0 && iLevel>0) // Loop checking, 11.05.2015, T. Koi
	{
	  --iLevel;
	  eExcitation = theGammas.GetLevel(iLevel)->GetLevelEnergy();    
	}

        // Use QI value for calculating excitation energy of residual.
        G4bool useQI=false;
        G4double dqi = QI[it]; 
        if (dqi < 0 || 849 < dqi) useQI = true; // Former libraries do not have values in this range
 
        if (useQI) {
          // QI introudced since G4NDL3.15
          // G4double QM=(incidReactionProduct.GetMass()+targetMass)-(aHadron.GetMass()+residualMass);
          // eExcitation = QM-QI[it];
          // eExcitation = QI[0] - QI[it];   // Bug fix #1838
          // if(eExcitation < 20*CLHEP::keV) eExcitation = 0;

          eExcitation = std::max(0.,QI[0] - QI[it]);  // Bug fix 2333

          // Re-evaluate iLevel based on this eExcitation
          iLevel = 0;
          G4bool find = false;
          G4int imaxEx = 0;
          G4double level_tolerance = 1.0*CLHEP::keV;

          while( theGammas.GetLevel(iLevel+1) != 0 ) { // Loop checking, 11.05.2015, T. Koi
            G4double maxEx = 0.0;
            if (maxEx < theGammas.GetLevel(iLevel)->GetLevelEnergy() ) {
              maxEx = theGammas.GetLevel(iLevel)->GetLevelEnergy();
              imaxEx = iLevel;
            }

            // Fix bug 1789 DHW - first if-branch added because gamma data come from ENSDF
            //                    and do not necessarily match the excitations used in ENDF-B.VII
            //                    Compromise solution: use 1 keV tolerance suggested by T. Koi
            if (std::abs(eExcitation - theGammas.GetLevel(iLevel)->GetLevelEnergy() ) < level_tolerance) {
              find = true;
              break;

            } else if (eExcitation < theGammas.GetLevel(iLevel)->GetLevelEnergy() ) {
              find = true;
              --iLevel;
              // very small eExcitation, iLevel becomes -1, this is protected below
              if (theTrack.GetDefinition() == aDefinition) { // this line added as part of fix #1838
                if (iLevel == -1) iLevel = 0;
              }
              break;
            }
            ++iLevel;
          }

          // If proper level cannot be found, use the maximum level
          if (!find) iLevel = imaxEx;
        }
	
	if(std::getenv("G4ParticleHPDebug") && eKinetic-eExcitation < 0) 
	{
	  throw G4HadronicException(__FILE__, __LINE__, "SEVERE: InelasticCompFS: Consistency of data not good enough, please file report");
	}
	if(eKinetic-eExcitation < 0) eExcitation = 0;
	if(iLevel!= -1) aHadron.SetKineticEnergy(eKinetic - eExcitation);
	
      }
      theAngularDistribution[it]->SampleAndUpdate(aHadron);

      if( theFinalStatePhotons[it] == 0 )
      {
	thePhotons = theGammas.GetDecayGammas(iLevel);
	eGamm -= theGammas.GetLevelEnergy(iLevel);
	if(eGamm>0) // @ ok for now, but really needs an efficient way of correllated sampling @
	{
          G4ReactionProduct * theRestEnergy = new G4ReactionProduct;
          theRestEnergy->SetDefinition(G4Gamma::Gamma());
          theRestEnergy->SetKineticEnergy(eGamm);
          G4double costh = 2.*G4UniformRand()-1.;
          G4double phi = CLHEP::twopi*G4UniformRand();
          theRestEnergy->SetMomentum(eGamm*std::sin(std::acos(costh))*std::cos(phi), 
                                     eGamm*std::sin(std::acos(costh))*std::sin(phi),
                                     eGamm*costh);
          if(thePhotons == 0) { thePhotons = new G4ReactionProductVector; }
          thePhotons->push_back(theRestEnergy);
	}
      }
    }
    else if(theEnergyAngData[it] != 0) // MF6  
    {

      theParticles = theEnergyAngData[it]->Sample(eKinetic);

      // Adjust A and Z in the case of miss much between selected data and target nucleus 
      if ( theParticles != NULL ) {
         G4int sumA = 0;
         G4int sumZ = 0;
         G4int maxA = 0;
         G4int jAtMaxA = 0;
         for ( G4int j = 0 ; j != (G4int)theParticles->size() ; ++j ) {
            if ( theParticles->at(j)->GetDefinition()->GetBaryonNumber() > maxA ) {
               maxA = theParticles->at(j)->GetDefinition()->GetBaryonNumber(); 
               jAtMaxA = j; 
            }
            sumA += theParticles->at(j)->GetDefinition()->GetBaryonNumber();
            sumZ += G4int( theParticles->at(j)->GetDefinition()->GetPDGCharge() + eps );
         }
         G4int dA = (G4int)theBaseA + hadProjectile->GetDefinition()->GetBaryonNumber() - sumA;
         G4int dZ = (G4int)theBaseZ + G4int( hadProjectile->GetDefinition()->GetPDGCharge() + eps ) - sumZ;
         if ( dA < 0 || dZ < 0 ) {
            G4int newA = theParticles->at(jAtMaxA)->GetDefinition()->GetBaryonNumber() + dA ;
            G4int newZ = G4int( theParticles->at(jAtMaxA)->GetDefinition()->GetPDGCharge() + eps ) + dZ;
            G4ParticleDefinition* pd = G4IonTable::GetIonTable()->GetIon ( newZ , newA );
            theParticles->at( jAtMaxA )->SetDefinition( pd );
         }
      }
    }
    else
    {
      // @@@ what to do, if we have photon data, but no info on the hadron itself
      nothingWasKnownOnHadron = 1;
    }

    if ( theFinalStatePhotons[it] != 0 ) 
    {
       // the photon distributions are in the Nucleus rest frame.
       // TK residual rest frame
      G4ReactionProduct boosted_tmp;
      boosted_tmp.Lorentz(incidReactionProduct, theTarget);
      G4double anEnergy = boosted_tmp.GetKineticEnergy();
      thePhotons = theFinalStatePhotons[it]->GetPhotons(anEnergy);
      G4double aBaseEnergy = theFinalStatePhotons[it]->GetLevelEnergy();
      G4double testEnergy = 0;
      if(thePhotons!=0 && thePhotons->size()!=0)
      { aBaseEnergy-=thePhotons->operator[](0)->GetTotalEnergy(); }
      if(theFinalStatePhotons[it]->NeedsCascade())
      {
	while(aBaseEnergy>0.01*CLHEP::keV) // Loop checking, 11.05.2015, T. Koi
        {
          // cascade down the levels
	  G4bool foundMatchingLevel = false;
          G4int closest = 2;
	  G4double deltaEold = -1;
	  for(G4int j=1; j<it; ++j)
          {
            if(theFinalStatePhotons[j]!=0) 
            {
              testEnergy = theFinalStatePhotons[j]->GetLevelEnergy();
            }
            else
            {
              testEnergy = 0;
            }
	    G4double deltaE = std::abs(testEnergy-aBaseEnergy);
            if(deltaE<0.1*CLHEP::keV)
            {
              G4ReactionProductVector * theNext = 
        	theFinalStatePhotons[j]->GetPhotons(anEnergy);
              if ( thePhotons != nullptr )
                thePhotons->push_back(theNext->operator[](0));
              aBaseEnergy = testEnergy-theNext->operator[](0)->GetTotalEnergy();
              delete theNext;
	      foundMatchingLevel = true;
              break; // ===>
            }
	    if(theFinalStatePhotons[j]!=0 && ( deltaE<deltaEold||deltaEold<0.) )
	    {
	      closest = j;
	      deltaEold = deltaE;     
	    }
          } // <=== the break goes here.
	  if(!foundMatchingLevel)
	  {
            G4ReactionProductVector * theNext = 
               theFinalStatePhotons[closest]->GetPhotons(anEnergy);
            if ( thePhotons != nullptr )
              thePhotons->push_back(theNext->operator[](0));
            aBaseEnergy = aBaseEnergy-theNext->operator[](0)->GetTotalEnergy();
            delete theNext;
	  }
        } 
      }
    }
    unsigned int i0;
    if(thePhotons!=0)
    {
      for(i0=0; i0<thePhotons->size(); ++i0)
      {
	// back to lab
	thePhotons->operator[](i0)->Lorentz(*(thePhotons->operator[](i0)), -1.*theTarget);
      }
    }
    if (nothingWasKnownOnHadron)
    {
      // In this case, hadron should be isotropic in CM
      // Next 12 lines are Emilio's replacement 
      // G4double QM=(incidReactionProduct.GetMass()+targetMass)-(aHadron.GetMass()+residualMass);
      // G4double eExcitation = QM-QI[it];
      // G4double eExcitation = QI[0] - QI[it];  // Fix of bug #1838
      // if(eExcitation<20*CLHEP::keV){eExcitation=0;}

      G4double eExcitation = std::max(0.,QI[0] - QI[it]);  // Fix of bug #2333

      two_body_reaction(&incidReactionProduct,&theTarget,&aHadron,eExcitation);
      if(thePhotons==0 && eExcitation>0){
        for(iLevel=theGammas.GetNumberOfLevels()-1; iLevel>=0; --iLevel)
        {
          if(theGammas.GetLevelEnergy(iLevel)<eExcitation+5*keV) break; // 5 keV tolerance
        }
        thePhotons = theGammas.GetDecayGammas(iLevel);
      }
    }

    // fill the result
    // Beware - the recoil is not necessarily in the particles...
    // Can be calculated from momentum conservation?
    // The idea is that the particles ar emitted forst, and the gammas only once the
    // recoil is on the residual; assumption is that gammas do not contribute to 
    // the recoil.
    // This needs more design @@@

    G4bool needsSeparateRecoil = false;
    G4int totalBaryonNumber = 0;
    G4int totalCharge = 0;
    G4ThreeVector totalMomentum(0);
    if(theParticles != 0) 
    {
      const G4ParticleDefinition * aDef;
      std::size_t ii0;
      for(ii0=0; ii0<theParticles->size(); ++ii0)
      {
        aDef = theParticles->operator[](ii0)->GetDefinition();
	totalBaryonNumber+=aDef->GetBaryonNumber();
	totalCharge+=G4int(aDef->GetPDGCharge()+eps);
        totalMomentum += theParticles->operator[](ii0)->GetMomentum();
      } 
      if(totalBaryonNumber!=G4int(theBaseA+eps+hadProjectile->GetDefinition()->GetBaryonNumber())) 
      {
        needsSeparateRecoil = true;
	residualA = G4int(theBaseA+eps+hadProjectile->GetDefinition()->GetBaryonNumber()
	                  -totalBaryonNumber);
	residualZ = G4int(theBaseZ+eps+hadProjectile->GetDefinition()->GetPDGCharge()
	                  -totalCharge);
      }
    }
    
    std::size_t nPhotons = 0;
    if(thePhotons!=0) { nPhotons = thePhotons->size(); }
        
    G4DynamicParticle * theSec;
    
    if( theParticles==0 )
    {
      theSec = new G4DynamicParticle;   
      theSec->SetDefinition(aHadron.GetDefinition());
      theSec->SetMomentum(aHadron.GetMomentum());
      theResult.Get()->AddSecondary(theSec, secID);    
#ifdef G4PHPDEBUG
      if( std::getenv("G4ParticleHPDebug"))
        G4cout << this << " G4ParticleHPInelasticCompFS::BaseApply  add secondary1 "
               << theSec->GetParticleDefinition()->GetParticleName() << " E= "
               << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
      
      aHadron.Lorentz(aHadron, theTarget);
      G4ReactionProduct theResidual;   
      theResidual.SetDefinition(G4IonTable::GetIonTable()
				->GetIon(static_cast<G4int>(residualZ), static_cast<G4int>(residualA), 0));  
      theResidual.SetKineticEnergy(aHadron.GetKineticEnergy()*aHadron.GetMass()/theResidual.GetMass());
      
      //080612TK contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #6
      //theResidual.SetMomentum(-1.*aHadron.GetMomentum());
      G4ThreeVector incidentNeutronMomentum = incidReactionProduct.GetMomentum();
      theResidual.SetMomentum(incidentNeutronMomentum - aHadron.GetMomentum());
      
      theResidual.Lorentz(theResidual, -1.*theTarget);
      G4ThreeVector totalPhotonMomentum(0,0,0);
      if(thePhotons!=0)
      {
        for(i=0; i<(G4int)nPhotons; ++i)
        {
          totalPhotonMomentum += thePhotons->operator[](i)->GetMomentum();
        }
      }
      theSec = new G4DynamicParticle;   
      theSec->SetDefinition(theResidual.GetDefinition());
      theSec->SetMomentum(theResidual.GetMomentum()-totalPhotonMomentum);
      theResult.Get()->AddSecondary(theSec, secID);    
#ifdef G4PHPDEBUG
      if( std::getenv("G4ParticleHPDebug"))
        G4cout << this
               << " G4ParticleHPInelasticCompFS::BaseApply add secondary2 "
               << theSec->GetParticleDefinition()->GetParticleName() << " E= "
               << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
    }
    else
    {
      for(i0=0; i0<theParticles->size(); ++i0)
      {
        theSec = new G4DynamicParticle; 
        theSec->SetDefinition(theParticles->operator[](i0)->GetDefinition());
        theSec->SetMomentum(theParticles->operator[](i0)->GetMomentum());
        theResult.Get()->AddSecondary(theSec, secID); 
#ifdef G4PHPDEBUG
      if( std::getenv("G4ParticleHPDebug"))
        G4cout << this
               << " G4ParticleHPInelasticCompFS::BaseApply add secondary3 "
               << theSec->GetParticleDefinition()->GetParticleName() << " E= "
               << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
        delete theParticles->operator[](i0); 
      } 
      delete theParticles;
      if(needsSeparateRecoil && residualZ!=0)
      {
        G4ReactionProduct theResidual;   
        theResidual.SetDefinition(G4IonTable::GetIonTable()
	                          ->GetIon(static_cast<G4int>(residualZ), static_cast<G4int>(residualA), 0));  
        G4double resiualKineticEnergy  = theResidual.GetMass()*theResidual.GetMass();
                 resiualKineticEnergy += totalMomentum*totalMomentum;
  	         resiualKineticEnergy  = std::sqrt(resiualKineticEnergy) - theResidual.GetMass();
	theResidual.SetKineticEnergy(resiualKineticEnergy);

        //080612TK contribution from Benoit Pirard and Laurent Desorgher (Univ. Bern) #4
        //theResidual.SetMomentum(-1.*totalMomentum);
	//G4ThreeVector incidentNeutronMomentum = incidReactionProduct.GetMomentum();
        //theResidual.SetMomentum(incidentNeutronMomentum - aHadron.GetMomentum());
        //080717 TK Comment still do NOT include photon's mometum which produce by thePhotons
        theResidual.SetMomentum( incidReactionProduct.GetMomentum() + theTarget.GetMomentum() - totalMomentum );

        theSec = new G4DynamicParticle;   
        theSec->SetDefinition(theResidual.GetDefinition());
        theSec->SetMomentum(theResidual.GetMomentum());
        theResult.Get()->AddSecondary(theSec, secID);  
#ifdef G4PHPDEBUG
      if( std::getenv("G4ParticleHPDebug"))
        G4cout << this
               << " G4ParticleHPInelasticCompFS::BaseApply add secondary4 "
               << theSec->GetParticleDefinition()->GetParticleName() << " E= "
               << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif
      }  
    }
    if(thePhotons!=0)
    {
      for(i=0; i<(G4int)nPhotons; ++i)
      {
        theSec = new G4DynamicParticle;    
        //Bug reported Chao Zhang (Chao.Zhang@usd.edu), Dongming Mei(Dongming.Mei@usd.edu) Feb. 25, 2009 
        //theSec->SetDefinition(G4Gamma::Gamma());
        theSec->SetDefinition( thePhotons->operator[](i)->GetDefinition() );
        //But never cause real effect at least with G4NDL3.13 TK
        theSec->SetMomentum(thePhotons->operator[](i)->GetMomentum());
        theResult.Get()->AddSecondary(theSec, secID); 
#ifdef G4PHPDEBUG
      if( std::getenv("G4ParticleHPDebug"))
        G4cout << this
               << " G4ParticleHPInelasticCompFS::BaseApply add secondary5 "
               << theSec->GetParticleDefinition()->GetParticleName() << " E= "
               << theSec->GetKineticEnergy() << " NSECO "
               << theResult.Get()->GetNumberOfSecondaries() << G4endl;
#endif

        delete thePhotons->operator[](i);
      }
      // some garbage collection
      delete thePhotons;
    }

    G4ParticleDefinition* targ_pd = G4IonTable::GetIonTable()->GetIon ( (G4int)theBaseZ , (G4int)theBaseA , 0.0 );
    G4LorentzVector targ_4p_lab ( theTarget.GetMomentum() , std::sqrt( targ_pd->GetPDGMass()*targ_pd->GetPDGMass() + theTarget.GetMomentum().mag2() ) );
    G4LorentzVector proj_4p_lab = theTrack.Get4Momentum();
    G4LorentzVector init_4p_lab = proj_4p_lab + targ_4p_lab;
    adjust_final_state ( init_4p_lab ); 

    // clean up the primary neutron
    theResult.Get()->SetStatusChange( stopAndKill );
}

//Re-implemented by E. Mendoza (2019). Isotropic emission in the CMS:
// proj: projectile in target-rest-frame (input)
// targ: target in target-rest-frame (input)
// product: secondary particle in target-rest-frame (output)
// resExcitationEnergy: excitation energy of the residual nucleus
//
void G4ParticleHPInelasticCompFS::two_body_reaction(G4ReactionProduct* proj,
                                                    G4ReactionProduct* targ,
                                                    G4ReactionProduct* product, 
                                                    G4double resExcitationEnergy)
{
  //CMS system:
  G4ReactionProduct theCMS= *proj+ *targ;

  //Residual definition:
  G4int resZ=(G4int)(proj->GetDefinition()->GetPDGCharge()+targ->GetDefinition()->GetPDGCharge()-product->GetDefinition()->GetPDGCharge()+0.1);
  G4int resA=proj->GetDefinition()->GetBaryonNumber()+targ->GetDefinition()->GetBaryonNumber()-product->GetDefinition()->GetBaryonNumber();
  G4ReactionProduct theResidual;
  theResidual.SetDefinition(G4IonTable::GetIonTable()->GetIon(resZ,resA,0.0));

  //CMS system:
  G4ReactionProduct theCMSproj;
  G4ReactionProduct theCMStarg;
  theCMSproj.Lorentz(*proj,theCMS);
  theCMStarg.Lorentz(*targ,theCMS);
  //final Momentum in the CMS:
  G4double totE=std::sqrt(theCMSproj.GetMass()*theCMSproj.GetMass()+theCMSproj.GetTotalMomentum()*theCMSproj.GetTotalMomentum())+std::sqrt(theCMStarg.GetMass()*theCMStarg.GetMass()+theCMStarg.GetTotalMomentum()*theCMStarg.GetTotalMomentum());
  G4double prodmass=product->GetMass();
  G4double resmass=theResidual.GetMass()+resExcitationEnergy;
  G4double fmomsquared=1./4./totE/totE*(totE*totE-(prodmass-resmass)*(prodmass-resmass))*(totE*totE-(prodmass+resmass)*(prodmass+resmass));
  G4double fmom=0;
  if(fmomsquared>0){
    fmom=std::sqrt(fmomsquared);
  }

  //random (isotropic direction):
  G4double cosTh = 2.*G4UniformRand()-1.;
  G4double phi = CLHEP::twopi*G4UniformRand();
  G4double theta = std::acos(cosTh);
  G4double sinth = std::sin(theta);
  product->SetMomentum(fmom*sinth*std::cos(phi),fmom*sinth*std::sin(phi),fmom*cosTh); //CMS
  product->SetTotalEnergy(std::sqrt(prodmass*prodmass+fmom*fmom)); //CMS
  //Back to the LAB system:
  product->Lorentz(*product,-1.*theCMS);
}

G4bool G4ParticleHPInelasticCompFS::use_nresp71_model( const G4ParticleDefinition* aDefinition , const G4int it , const G4ReactionProduct& theTarget , G4ReactionProduct& boosted )
{
  if ( aDefinition == G4Neutron::Definition() ) // If the outgoing particle is a neutron...
  {
    // LR: flag LR in ENDF. It indicates whether there is breakup of the residual nucleus or not.
    // it: exit channel (index of the carbon excited state)

    // Added by A. R. Garcia (CIEMAT) to include the physics of C(N,N'3A) reactions from NRESP71.

    if ( LR[it] > 0 ) // If there is breakup of the residual nucleus LR(flag LR in ENDF)>0 (i.e. Z=6 MT=52-91 (it=MT-50)).
    {
      // Defining carbon as the target in the reference frame at rest.
      G4ReactionProduct theCarbon(theTarget);

      theCarbon.SetMomentum(G4ThreeVector());
      theCarbon.SetKineticEnergy(0.);

      // Creating four reaction products.
      G4ReactionProduct theProds[4];

      // Applying C(N,N'3A) reaction mechanisms in the target rest frame.
      if ( it == 41 )
      {
        // QI=QM=-7.275 MeV for C-0(N,N')C-C(3A) in ENDF/B-VII.1.
        // This is not the value of the QI of the first step according
        // to the model. So we don't take it. Instead, we set the one
        // we have calculated: QI=(mn+m12C)-(ma+m9Be+Ex9Be)=-8.130 MeV.
        nresp71_model.ApplyMechanismI_NBeA2A(boosted, theCarbon, theProds, -8.130/*QI[it]*/);
        // N+C --> A[0]+9BE* | 9BE* --> N[1]+8BE | 8BE --> 2*A[2,3].
      }
      else
      {
        nresp71_model.ApplyMechanismII_ACN2A(boosted, theCarbon, theProds, QI[it]);
        // N+C --> N'[0]+C* | C* --> A[1]+8BE | 8BE --> 2*A[2,3].
      }

      // Returning to the reference frame where the target was in motion.
      for ( G4int j=0; j<4; ++j )
      {
        theProds[j].Lorentz(theProds[j], -1.*theTarget);
        theResult.Get()->AddSecondary(new G4DynamicParticle(theProds[j].GetDefinition(), theProds[j].GetMomentum()), secID);
      }

      // Killing the primary neutron.
      theResult.Get()->SetStatusChange(stopAndKill);
 
      return true;
    }
  }
  else if ( aDefinition == G4Alpha::Definition() ) // If the outgoing particle is an alpha, ...
  {
    // Added by A. R. Garcia (CIEMAT) to include the physics of C(N,A)9BE reactions from NRESP71.

    if ( LR[it] == 0 ) // If Z=6, an alpha particle is emitted and there is no breakup of the residual nucleus LR(flag LR in ENDF)==0.
    {
      // Defining carbon as the target in the reference frame at rest.
      G4ReactionProduct theCarbon(theTarget);
      theCarbon.SetMomentum(G4ThreeVector());
      theCarbon.SetKineticEnergy(0.);

      // Creating four reaction products.
      G4ReactionProduct theProds[2];

      // Applying C(N,A)9BE reaction mechanism.
      nresp71_model.ApplyMechanismABE(boosted, theCarbon, theProds);
      // N+C --> A[0]+9BE[1].

      for ( G4int j=0; j<2; ++j )
      {
        // Returning to the system of reference where the target was in motion.
        theProds[j].Lorentz(theProds[j], -1.*theTarget);
        theResult.Get()->AddSecondary(new G4DynamicParticle(theProds[j].GetDefinition(), theProds[j].GetMomentum()), secID);
      }

      // Killing the primary neutron.
      theResult.Get()->SetStatusChange(stopAndKill);;

      return true;
    }
    else
    {
      G4Exception("G4ParticleHPInelasticCompFS::CompositeApply()",
                  "G4ParticleInelasticCompFS.cc", FatalException,
                  "Alpha production with LR!=0.");
    }
  }
  return false;
}
