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
// $Id: G4VelocityTable.cc 107603 2017-11-27 07:21:59Z gcosmo $
//
//
//---------------------------------------------------------------
//
//  G4VelocityTable.cc
//
//  class description:
//    This class keeps a table of velocity as a function of
//    the ratio kinetic erngy and mass
//
//---------------------------------------------------------------
//   created                     17.Aug.  2011 H.Kurashige
//

#include "G4VelocityTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#include "G4ios.hh" 

G4ThreadLocal G4VelocityTable* G4VelocityTable::theInstance = nullptr;

////////////////
G4VelocityTable::G4VelocityTable()
////////////////
  : edgeMin(0.), edgeMax(0.), numberOfNodes(0),
    dBin(0.), baseBin(0.),
    lastEnergy(-DBL_MAX), lastValue(0.), lastBin(0),
    maxT( 1000.0 ), minT( 0.0001 ), NbinT( 500 )
{
  PrepareVelocityTable();
} 

/////////////////
G4VelocityTable::~G4VelocityTable()
/////////////////
{
  dataVector.clear();
  binVector.clear();
}


///////////////////
void G4VelocityTable::PrepareVelocityTable()
///////////////////
{
  //const G4double g4log10 = std::log(10.); 
  
  dataVector.clear();
  binVector.clear();
  dBin     =  G4Log(maxT/minT)/NbinT;
  baseBin  =  G4Log(minT)/dBin;

  numberOfNodes = NbinT + 1;
  dataVector.reserve(numberOfNodes);
  binVector.reserve(numberOfNodes);

  binVector.push_back(minT);
  dataVector.push_back(0.0);

  for (size_t i=1; i<numberOfNodes-1; i++){
    binVector.push_back(G4Exp((baseBin+i)*dBin));
    dataVector.push_back(0.0);
  }
  binVector.push_back(maxT);
  dataVector.push_back(0.0);
  
  edgeMin = binVector[0];
  edgeMax = binVector[numberOfNodes-1];

  for (G4int i=0; i<=NbinT; i++){
    G4double T = binVector[i];
    dataVector[i]= c_light*std::sqrt(T*(T+2.))/(T+1.0);    
  }

  return;
} 

size_t G4VelocityTable::FindBinLocation(G4double theEnergy) const
{
  // For G4PhysicsLogVector, FindBinLocation is implemented using
  // a simple arithmetic calculation.
  //
  // Because this is a virtual function, it is accessed through a
  // pointer to the G4PhyiscsVector object for most usages. In this
  // case, 'inline' will not be invoked. However, there is a possibility 
  // that the user access to the G4PhysicsLogVector object directly and 
  // not through pointers or references. In this case, the 'inline' will
  // be invoked. (See R.B.Murray, "C++ Strategies and Tactics", Chap.6.6)

  //const G4double g4log10 = G4Log(10.); 
  return size_t( G4Log(theEnergy)/dBin - baseBin );
}

G4double G4VelocityTable::Value(G4double theEnergy) 
{
  // Use cache for speed up - check if the value 'theEnergy' is same as the 
  // last call. If it is same, then use the last bin location. Also the
  // value 'theEnergy' lies between the last energy and low edge of of the 
  // bin of last call, then the last bin location is used.

  if( theEnergy == lastEnergy ) {

  } else if( theEnergy < lastEnergy
        &&   theEnergy >= binVector[lastBin]) {
     lastEnergy = theEnergy;
     lastValue = Interpolation();

  } else if( theEnergy <= edgeMin ) {
     lastBin = 0;
     lastEnergy = edgeMin;
     lastValue  = dataVector[0];

  } else if( theEnergy >= edgeMax ){
     lastBin = numberOfNodes-1;
     lastEnergy = edgeMax;
     lastValue  = dataVector[lastBin];

  } else {
     lastBin = (size_t)( G4Log(theEnergy)/dBin - baseBin ); 
     if(lastBin == numberOfNodes) { --lastBin; } // VI: fix possible precision lost
     lastEnergy = theEnergy;
     lastValue = Interpolation();
     
  }
  return lastValue;        
}

///////////////////////////////////////////////////////
// Static methods

///////////////////
G4VelocityTable* G4VelocityTable::GetVelocityTable()
///////////////////
{
  if (!theInstance)  { 
    static G4ThreadLocalSingleton<G4VelocityTable> inst;
    theInstance = inst.Instance();
  }
  return theInstance;
}

///////////////////
void G4VelocityTable::SetVelocityTableProperties(G4double t_max, G4double t_min, G4int nbin)
///////////////////
{
  if (theInstance == nullptr)  { GetVelocityTable(); }

  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  
  // check if state is outside event loop
  if(!(currentState==G4State_Idle||currentState==G4State_PreInit)){ 
    G4Exception("G4VelocityTable::SetVelocityTableProperties",
                "Track101", JustWarning,
                "Can modify only in PreInit or Idle state : Method ignored.");
    return;
  }

  if (nbin > 100 )  theInstance->NbinT = nbin;
  if ((t_min < t_max)&&(t_min>0.))  {
    theInstance->minT = t_min; 
    theInstance->maxT = t_max;
  } 
  theInstance->PrepareVelocityTable();
}

///////////////////
G4double G4VelocityTable::GetMaxTOfVelocityTable()
///////////////////
{
  return GetVelocityTable()->maxT; 
}

///////////////////
G4double G4VelocityTable::GetMinTOfVelocityTable() 
///////////////////
{
  return GetVelocityTable()->minT; 
}

///////////////////
G4int    G4VelocityTable::GetNbinOfVelocityTable() 
///////////////////
{
  return GetVelocityTable()->NbinT; 
}
