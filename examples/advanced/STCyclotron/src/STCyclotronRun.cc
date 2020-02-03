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
// Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronRun.cc

#include "STCyclotronRun.hh"
#include "STCyclotronAnalysis.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

STCyclotronRun::STCyclotronRun()
  : G4Run(),fTotalEnergyDepositTarget(0.),fTotalEnergyDepositFoil(0.),fParticleTarget(0),fTargetThickness(0.),fTargetDiameter(0.),fFoilThickness(0.),fTargetVolume(0.),fFoilVolume(0.),fPrimariesPerEvent(0),fTimePerEvent(0),fBeamName(""),fBeamCurrent(0.),fBeamEnergy(0.)
{ }

STCyclotronRun::~STCyclotronRun()
{ }

void STCyclotronRun::Merge(const G4Run* aRun)
{

  const STCyclotronRun* localRun = static_cast<const STCyclotronRun*>(aRun);

  //Merging cumulable variables
  fTotalEnergyDepositTarget += localRun->fTotalEnergyDepositTarget;
  fTotalEnergyDepositFoil += localRun->fTotalEnergyDepositFoil;
  fParticleTarget += localRun->fParticleTarget;

  //Constant over the different runs
  if(localRun->fTargetVolume!=0)fTargetVolume = localRun->fTargetVolume;
  if(localRun->fFoilVolume!=0)fFoilVolume = localRun->fFoilVolume;
  if(localRun->fPrimariesPerEvent!=0)fPrimariesPerEvent = localRun->fPrimariesPerEvent;
  if(localRun->fTimePerEvent!=0)fTimePerEvent = localRun->fTimePerEvent;
  if(localRun->fTargetThickness!=0)fTargetThickness = localRun->fTargetThickness;
  if(localRun->fTargetDiameter!=0)fTargetDiameter = localRun->fTargetDiameter;
  if(localRun->fFoilThickness!=0)fFoilThickness = localRun->fFoilThickness;
  fBeamName = localRun->fBeamName;
  if(localRun->fBeamCurrent!=0.)fBeamCurrent = localRun->fBeamCurrent;
  if(localRun->fBeamEnergy!=0.)fBeamEnergy = localRun->fBeamEnergy;

 
  
  //<<<----toMerge
  std::map<G4String,G4int>::iterator itSI;
  std::map<G4String,G4double>::iterator itSD;
  std::map<G4String,G4String>::iterator itSS;
  std::map<G4int,G4String>::iterator itIS;

  //----Merging results for primary isotopes
  std::map<G4String,G4int> locPrimaryIsotopeCountTarget = localRun->fPrimaryIsotopeCountTarget;
  for (itSI = locPrimaryIsotopeCountTarget.begin(); itSI != locPrimaryIsotopeCountTarget.end(); itSI++)
    { 
      G4String name = itSI->first;
      G4int count   = itSI->second;
      fPrimaryIsotopeCountTarget[name] += count;
    }
  
  std::map<G4String,G4double> locPrimaryIsotopeTimeTarget = localRun->fPrimaryIsotopeTimeTarget;
  for (itSD = locPrimaryIsotopeTimeTarget.begin(); itSD != locPrimaryIsotopeTimeTarget.end(); itSD++)
    { 
      G4String name = itSD->first;
      G4double time   = itSD->second;
      fPrimaryIsotopeTimeTarget[name] = time;
    }

  //----Merging results for decay isotopes
  // std::map<G4int,G4String>    fIsotopeIDTarget;
  std::map<G4String,G4String> locDecayIsotopeCountTarget = localRun->fDecayIsotopeCountTarget;
  for (itSS = locDecayIsotopeCountTarget.begin(); itSS != locDecayIsotopeCountTarget.end(); itSS++)
    { 
      G4String nameDaughter = itSS->first;
      G4String mum   = itSS->second;
      fDecayIsotopeCountTarget[nameDaughter] = mum;
    }
  std::map<G4String,G4double> locDecayIsotopeTimeTarget = localRun->fDecayIsotopeTimeTarget;
   for (itSD = locDecayIsotopeTimeTarget.begin(); itSD != locDecayIsotopeTimeTarget.end(); itSD++)
    { 
      G4String nameDaughter = itSD->first;
      G4double time   = itSD->second;
      fDecayIsotopeTimeTarget[nameDaughter] = time;
    }
   std::map<G4String,G4String> locParticleParent = localRun->fParticleParent;
   for (itSS = locParticleParent.begin(); itSS != locParticleParent.end(); itSS++)
     { 
       G4String nameDaughter = itSS->first;
       G4String parent   = itSS->second;
       fParticleParent[nameDaughter] = parent;
     }
   std::map<G4int,G4String> locIsotopeIDTarget = localRun->fIsotopeIDTarget;
   for (itIS = locIsotopeIDTarget.begin(); itIS != locIsotopeIDTarget.end(); itIS++)
     { 
       G4int ID = itIS->first;
       G4String name   = itIS->second;
       fIsotopeIDTarget[ID] = name;
     }
 
  
  //----Merging results for stable isotopes
  std::map<G4String,G4int> locStableIsotopeCountTarget = localRun->fStableIsotopeCountTarget;
  for (itSI = locStableIsotopeCountTarget.begin(); itSI != locStableIsotopeCountTarget.end(); itSI++)
    { 
      G4String name = itSI->first;
      G4int count   = itSI->second;
      fStableIsotopeCountTarget[name] += count;
    }
 
  //----Merging results for particles
  std::map<G4String,G4int> locParticleCountTarget = localRun->fParticleCountTarget;
  for (itSI = locParticleCountTarget.begin(); itSI != locParticleCountTarget.end(); itSI++)
    { 
      G4String name = itSI->first;
      G4int count   = itSI->second;
      fParticleCountTarget[name] += count;
    }
  
  
  G4Run::Merge(aRun); 

} 

void STCyclotronRun::EndOfRun(G4double irradiationTime) 
{

  G4int nbEvents = GetNumberOfEvent();  

  if (nbEvents == 0) return;

  //------------------------------------------------------
  //  Opening the ASCII file
  //------------------------------------------------------
  fOutPut.open("Output_General.txt",std::ofstream::out);
  fOutPut1.open("Output_ParentIsotopes.txt",std::ofstream::out);
  fOutPut2.open("Output_DaughterIsotopes.txt",std::ofstream::out);
  fOutPut3.open("Output_OtherParticles.txt",std::ofstream::out);
  fOutPut4.open("Output_StableIsotopes.txt",std::ofstream::out); 
  
  //------------------------------------------------------
  //  Calculates the equivalent time for a given run
  //------------------------------------------------------
  G4double timePerEvent = fTimePerEvent; //in seconds
  G4double timeForARun = nbEvents*timePerEvent; //in seconds
  G4double minDecay = 0.0001; //in seconds
  G4double maxDecay = 1000000.; //in seconds
  
  //------------------------------------------------------
  //    Rescale the value of the beam current to account for
  //    the loss of primary particles due to the foil.
  //------------------------------------------------------
  
  G4int totalPrimaries = fPrimariesPerEvent*nbEvents;
  G4double currentFactor; 
  if(fParticleTarget>0.) currentFactor =(fParticleTarget*1.)/(totalPrimaries*1.);
  else currentFactor = 0.;


  fOutPut << "//-----------------------------------//" << G4endl;
  fOutPut << "//    Parameters of the simulation:  //" << G4endl;
  fOutPut << "//-----------------------------------//" << G4endl;
  fOutPut << "Beam parameters: " << G4endl;
  fOutPut << fBeamName << " - Name of beam primary particles." << G4endl;
  fOutPut << fBeamEnergy << " - Energy of beam primary particles (MeV)." << G4endl;
  fOutPut << fBeamCurrent << " - Beam current (Ampere)." << G4endl;
  fOutPut << irradiationTime << " - Irradiation time in hour(s)." << G4endl;
  fOutPut << currentFactor << " - Current factor." << G4endl;
  fOutPut << "//-----------------------------------//" << G4endl;
  fOutPut << "Simulation parameters: " << G4endl;
  fOutPut << timePerEvent << " - Equivalent time per event (s)." << G4endl;
  fOutPut << nbEvents << " - Number of events" << G4endl;
  fOutPut << fPrimariesPerEvent << " - Primaries per event" << G4endl;
  fOutPut << fPrimariesPerEvent*nbEvents << " - Total number of particles sent." << G4endl;
  fOutPut << "//-----------------------------------//" << G4endl;
  fOutPut << "Geometry parameters: " << G4endl;
  fOutPut << fTargetThickness << " - target thickness (mm)." << G4endl;
  fOutPut << fTargetDiameter << " - target diameter (mm)." << G4endl;
  fOutPut << fFoilThickness << " - foil thickness (mm)." << G4endl;
  //Add particle type, particle energy, beam diameter, beam current

  //target material???
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////Calculation of the number of isotopes at the end of the irradiation and the activity generated/////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  
  //Maps to fill
  std::map<G4String,G4double> fPrimaryIsotopeEOBTarget;
  std::map<G4String,G4double> fPrimaryActivityTarget;
  std::map<G4String,G4double> fDecayIsotopeEOBTarget;
  std::map<G4String,G4double> fDecayActivityTarget;

  
  //----------------------------------------------
  //        CASE 1 : Parent isotopes
  //----------------------------------------------
 
  
  std::map<G4String,G4int>::iterator it;  
  G4double primaryActivityTotal = 0.;
  G4double decayActivityTotal=0.;

  fOutPut1 << "//-----------------------------------//\n"
    << "//       Data for parent isotopes    //\n"
	  << "//-----------------------------------//\n" << G4endl;
 
  for (it = fPrimaryIsotopeCountTarget.begin(); it != fPrimaryIsotopeCountTarget.end(); it++)
    {

     
      G4String name = it->first;
      G4double count   = (it->second)*currentFactor;
      G4double halfLifeTime = fPrimaryIsotopeTimeTarget[name]*10E-10/3600.*std::log(2.);
      G4String process = fParticleParent[name];

      //Only store isotopes with a life time between minDecay and maxDecay.
      G4bool store;
      if(halfLifeTime > minDecay && halfLifeTime < maxDecay) store = true;
      else store = false;

   
      //Calculation of the yield (s-1)
      G4double decayConstant = 1/(fPrimaryIsotopeTimeTarget[name]*10E-10);
 

      //----------------------------------------------
      //       Number of particles per second
      //----------------------------------------------

      G4double particlesPerSecond = fPrimaryIsotopeCountTarget[name]*currentFactor/timeForARun;

      //----------------------------------------------
      //           Calculation yield EOB
      //----------------------------------------------

      fPrimaryIsotopeEOBTarget[name] = particlesPerSecond/decayConstant * (1. - std::exp(-irradiationTime*3600*decayConstant));

      //----------------------------------------------
      //       Calculation of the activity
      //       conversion factor Bq to mCi 
      //----------------------------------------------

      G4double conv = 2.7E-8;
      fPrimaryActivityTarget[name]= fPrimaryIsotopeEOBTarget[name]*decayConstant*conv;
  
      if(store)
	{

	  //----------------------------------------------
	  //    Incrementation for total primary activity
	  //----------------------------------------------

	  primaryActivityTotal = primaryActivityTotal + fPrimaryActivityTarget[name];
	}

       
      //---------------------------//
      //   Printing out results    //
      //---------------------------//
      
      if(store)
	{
	  fOutPut1 << name << " - name of parent isotope." << G4endl;
	  fOutPut1 << count/currentFactor << " - number of isotopes created during the simulation." << G4endl;
	  fOutPut1 << decayConstant << " - decay constant in s-1." << G4endl; 
	  fOutPut1 << halfLifeTime << " - half life time in hour(s)." << G4endl; 
	  fOutPut1 << process << " - creation process." << G4endl; 
	  fOutPut1 << particlesPerSecond << " - isotope per sec." << G4endl; 
	  fOutPut1 << fPrimaryIsotopeEOBTarget[name] << " - yield EOB." << G4endl;
	  fOutPut1 << fPrimaryActivityTarget[name] << " - activity (mCi) at the EOB." << G4endl;
	  fOutPut1 << "------------------------" << G4endl; 
	}        
 
    }


  //----------------------------------------------
  //  CASE 2 : isotopes from primary isotopes decay
  //----------------------------------------------

  fOutPut2 << "//-----------------------------------//\n"
	  << "//     Data for daughter isotopes    //\n"
	  << "//-----------------------------------//\n" << G4endl;

  std::map<G4String,G4String>::iterator it1;  

  for (it1 = fDecayIsotopeCountTarget.begin(); it1 != fDecayIsotopeCountTarget.end(); it1++)
    {

      
      G4String nameDaughter = it1->first;
      G4String nameMum   = it1->second;

      G4double halfLifeTimeMum = fDecayIsotopeTimeTarget[nameDaughter]*10E10/3600;
      G4double halfLifeTimeDaughter = fPrimaryIsotopeTimeTarget[nameMum]*10E10/3600;

      G4bool store;
      if(halfLifeTimeMum > minDecay && halfLifeTimeMum < maxDecay &&
	 halfLifeTimeDaughter > minDecay && halfLifeTimeDaughter < maxDecay){store=true;}
      else{store=false;}

  
   
      //----------------------------------------------
      //       Calculation of the yield
      //       fParticleTime[name] is the time
      //       life of the particle, divided by ln(2), in nS
      //----------------------------------------------

      G4double decayConstantMum = 1/(fPrimaryIsotopeTimeTarget[nameMum]*10.E-10);
      G4double decayConstantDaughter = 1/(fDecayIsotopeTimeTarget[nameDaughter]*10.E-10);

  
      //----------------------------------------------
      //         Number of particles per second
      //----------------------------------------------

      G4double particlesPerSecond = fPrimaryIsotopeCountTarget[nameMum]*currentFactor/timeForARun;

      //----------------------------------------------
      //        Number of particles at the EOB
      //----------------------------------------------  

      fDecayIsotopeEOBTarget[nameDaughter] = particlesPerSecond*((1 - std::exp(-irradiationTime*3600*decayConstantDaughter))/decayConstantDaughter + (std::exp(-irradiationTime*3600*decayConstantDaughter) - std::exp(-irradiationTime*3600*decayConstantMum))/(decayConstantDaughter-decayConstantMum));

 
      //----------------------------------------------
      //       Calculation of activity
      //       conversion factor Bq to mCu 
      //----------------------------------------------

      G4double conv = 2.7E-8;
      fDecayActivityTarget[nameDaughter]= fDecayIsotopeEOBTarget[nameDaughter]*decayConstantDaughter*conv;

      if(store)
	{
	  decayActivityTotal = decayActivityTotal + fDecayActivityTarget[nameDaughter];
	}

      if(store)
	{
	  fOutPut2 << nameDaughter << " - name of daughter isotope." << G4endl;
	  fOutPut2 << nameMum << " - name of parent isotope." << G4endl;
	  fOutPut2 << decayConstantDaughter << " - decay constant of daughter in s-1." << G4endl; 
	  fOutPut2 << decayConstantMum << " - decay constant of mum in s-1." << G4endl; 
	  fOutPut2 << halfLifeTimeDaughter << " - half life time of daughter in hour(s)." << G4endl; 
	  fOutPut2 << halfLifeTimeMum << " - half life time of mum in hour(s)." << G4endl; 
	  fOutPut2 << particlesPerSecond << "  - isotope per sec." << G4endl;
	  fOutPut2 << fDecayIsotopeEOBTarget[nameDaughter] << " - yield at the EOB." << G4endl;
	  fOutPut2 << fDecayActivityTarget[nameDaughter] << " - activity (mCi) at the EOB." << G4endl;
	  fOutPut2 << "------------------------" << G4endl; 
	}

    }

  //----------------------------------------------
  //    Particles created, other than nuclei
  //----------------------------------------------

  fOutPut3 << "//-----------------------------------//\n"
	  << "//      Data for other particles     //\n"
	  << "//-----------------------------------//" << G4endl;
  
  std::map<G4String, G4int>::iterator it3;
  for(it3=fParticleCountTarget.begin(); it3!= fParticleCountTarget.end(); it3++)
    {
      G4String name = it3->first;
      G4double number = it3->second;
      fOutPut3 << name << " - name of the particle" << G4endl;
      fOutPut3 << number << " - number of particles" << G4endl;
      fOutPut3 << "------------------------" << G4endl; 
    }

  


  fOutPut4 << "//-----------------------------------//\n"
	  << "//      Data for stable isotopes     //\n"
	  << "//-----------------------------------//\n" << G4endl;
    
  std::map<G4String, G4int>::iterator it6;
  for(it6=fStableIsotopeCountTarget.begin();it6!=fStableIsotopeCountTarget.end();it6++)
    {
      G4String isotope = it6 ->first;
      G4int number = it6 -> second;
      fOutPut4 << isotope << " - name of the isotope" << G4endl;
      fOutPut4 << number << " - number of isotopes" << G4endl;
      fOutPut4 << "------------------------" << G4endl; 
    }
 


  //Clear the maps
  fPrimaryIsotopeEOBTarget.clear();
  fPrimaryActivityTarget.clear();
  fDecayIsotopeEOBTarget.clear(); 
  fDecayActivityTarget.clear(); 

  
  //Clear the maps
  fPrimaryIsotopeCountTarget.clear();
  fPrimaryIsotopeTimeTarget.clear();
  fDecayIsotopeCountTarget.clear();
  fDecayIsotopeTimeTarget.clear();
  fParticleParent.clear();
  fParticleCountTarget.clear();
  fStableIsotopeCountTarget.clear();
  fIsotopeIDTarget.clear(); 

 
  //-----------------------------
  //     Calculation of heat
  //-----------------------------

  G4double totalEnergyDepositTargetEOB = fTotalEnergyDepositTarget/timeForARun * irradiationTime * 3600.;
  G4double totalEnergyDepositTargetPerSecond = fTotalEnergyDepositTarget/timeForARun;
  //Heat calculation in W/mm3
  G4double heatTarget = totalEnergyDepositTargetPerSecond/fTargetVolume * 1.60E-13;
  G4double heatFoil = fTotalEnergyDepositFoil / fFoilVolume * 1.60E-13; 

  
  //Output data in a .txt file
  fOutPut << "//-------------------------------------------------//\n"
	 << "//     Heating, total activity and process data    //\n"
	 << "//-------------------------------------------------//" << G4endl;

  fOutPut << "Total heating in the target : "
	 << heatTarget << " W/mm3" << G4endl;
  fOutPut << "The total heating during the irradiation is " << totalEnergyDepositTargetEOB << "J/mm3" << G4endl;
  fOutPut << "Total heating in the foil : " << heatFoil << " W/mm3" << G4endl;
  

  fOutPut.close();
  fOutPut1.close();
  fOutPut2.close();
  fOutPut3.close();
  fOutPut4.close();
  
}




// Accumulation functions for maps used at the end of run action

void STCyclotronRun::PrimaryIsotopeCountTarget(G4String name,G4double time)
{
  fPrimaryIsotopeCountTarget[name]++;
  fPrimaryIsotopeTimeTarget[name]=time;
}
//------------------------------------
void STCyclotronRun::CountStableIsotopes(G4String name)
{
  fStableIsotopeCountTarget[name]++;
}
//------------------------------------
void STCyclotronRun::DecayIsotopeCountTarget(G4String nameDaughter,G4String mum, G4double time)
{
  fDecayIsotopeCountTarget[nameDaughter]=mum;
  fDecayIsotopeTimeTarget[nameDaughter]=time;  
}
//------------------------------------
void STCyclotronRun::ParticleParent(G4String isotope, G4String parent)
{
  fParticleParent[isotope]=parent;
}
//
//-----> Count other particles
//------------------------------------
void STCyclotronRun::ParticleCountTarget(G4String name)
{
  fParticleCountTarget[name]++;
}



//-------------------------------------------------------------------------------------------------------------
// Accumulation functions for maps used only during the run
//
//-----> Isotope ID to obtain the "mother isotope" in SensitiveTarget()
//------------------------------------
void STCyclotronRun::StoreIsotopeID(G4int ID, G4String name)
{
  fIsotopeIDTarget[ID]=name;   
}
//
std::map<G4int,G4String> STCyclotronRun::GetIsotopeID()
{
  return fIsotopeIDTarget;
}


void STCyclotronRun::EnergyDepositionTarget(G4double edep)
{
  fTotalEnergyDepositTarget += edep;
}

void STCyclotronRun::EnergyDepositionFoil(G4double edep)
{
  fTotalEnergyDepositFoil += edep;
}

void STCyclotronRun::CountParticlesTarget()
{
  fParticleTarget++;
}

void STCyclotronRun::SetFoilVolume(G4double foilVolume)
{
  fFoilVolume = foilVolume;
}

void STCyclotronRun::SetFoilThickness(G4double foilThickness)
{
  fFoilThickness = foilThickness;
}

void STCyclotronRun::SetTargetVolume(G4double targetVolume)
{
  fTargetVolume = targetVolume;
}

void STCyclotronRun::SetTargetThickness(G4double targetThickness)
{
  fTargetThickness = targetThickness;
}

void STCyclotronRun::SetTargetDiameter(G4double targetDiameter)
{
  fTargetDiameter = targetDiameter;
}

void STCyclotronRun::SetPrimariesPerEvent(G4int primaries)
{
  fPrimariesPerEvent = primaries;
}

void STCyclotronRun::SetTimePerEvent(G4double timePerEvent)
{
  fTimePerEvent = timePerEvent;
}

void STCyclotronRun::SetBeamName(G4String beamName)
{
  fBeamName = beamName;
}

void STCyclotronRun::SetBeamCurrent(G4double beamCurrent)
{
  fBeamCurrent = beamCurrent;
}

void STCyclotronRun::SetBeamEnergy(G4double beamEnergy)
{
  fBeamEnergy = beamEnergy;
}
