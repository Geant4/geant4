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
// $Id: G4MolecularDissociationTable.cc 84858 2014-10-21 16:08:22Z gcosmo $
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// ----------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation by Alfonso Mantero 4 Mar 2009
//
// ----------------------------------------------------------------

#include "G4MolecularDissociationTable.hh"
#include "G4MolecularDissociationChannel.hh"

using namespace std;

G4MolecularDissociationTable::G4MolecularDissociationTable()
{
  ;
}

G4MolecularDissociationTable::~G4MolecularDissociationTable()
{
  channelsMap::iterator it_map = fDecayChannelsMap.begin();

  for (; it_map != fDecayChannelsMap.end(); it_map++)
  {
    vector<const G4MolecularDissociationChannel*>& decayChannels = it_map
        ->second;
    if (!decayChannels.empty())
    {
      for (int i = 0; i < (int) decayChannels.size(); i++)
      {
        if (decayChannels[i])
        {
          delete decayChannels[i];
          decayChannels[i] = 0;
        }
      }
      decayChannels.clear();
    }
  }
  fDecayChannelsMap.clear();
}

G4MolecularDissociationTable::G4MolecularDissociationTable(const G4MolecularDissociationTable& right)
{
  *this = right;
}

G4MolecularDissociationTable& G4MolecularDissociationTable::operator=(const G4MolecularDissociationTable& aMolecularDecayTable)
{
  fExcitedStatesMap = aMolecularDecayTable.fExcitedStatesMap;
  fDecayChannelsMap = channelsMap(aMolecularDecayTable.GetDecayChannelsMap());
  return *this;
}

const vector<const G4MolecularDissociationChannel*>* G4MolecularDissociationTable::GetDecayChannels(const G4ElectronOccupancy* conf) const
{
  statesMap::const_iterator it_exstates = fExcitedStatesMap.find(*conf);
  if (it_exstates == fExcitedStatesMap.end()) return 0;
  channelsMap::const_iterator it_decchannel = fDecayChannelsMap.find(
      it_exstates->second);
  if (it_decchannel == fDecayChannelsMap.end()) return 0;
  return &(it_decchannel->second);
}

const vector<const G4MolecularDissociationChannel*>* G4MolecularDissociationTable::GetDecayChannels(const G4String& exState) const
{
  channelsMap::const_iterator it_decchannel = fDecayChannelsMap.find(exState);
  if (it_decchannel == fDecayChannelsMap.end()) return 0;
  return &(it_decchannel->second);
}

const G4String& G4MolecularDissociationTable::GetExcitedState(const G4ElectronOccupancy* conf) const
{
  statesMap::const_iterator it_exstates = fExcitedStatesMap.find(*conf);

  if (it_exstates == fExcitedStatesMap.end())
  {
    G4String errMsg = "Excited state not found";
    G4Exception(
        "G4MolecularDecayTable::GetExcitedState(const G4ElectronOccupancy*)",
        "G4MolecularDecayTable001", FatalErrorInArgument, errMsg);
//        return *(new G4String("IM FAKE"));  // fake return statement
  }

  return it_exstates->second;
}

const G4ElectronOccupancy& G4MolecularDissociationTable::GetElectronOccupancy(const G4String& exState) const
{
  statesMap::const_iterator statesIter;
  const G4ElectronOccupancy* conf(0);
  for (statesIter = fExcitedStatesMap.begin();
      statesIter != fExcitedStatesMap.end(); statesIter++)
  {
    if (exState == statesIter->second) conf = &(statesIter->first);
  }

  if (statesIter == fExcitedStatesMap.end())
  {
    G4String errMsg = "Excited state" + exState + " not found";
    G4Exception("G4MolecularDecayTable::GetElectronOccupancy(const G4String&)",
                "G4MolecularDecayTable002", FatalErrorInArgument, errMsg);
  }
  return *conf;
}

void G4MolecularDissociationTable::AddExcitedState(const G4String& label)
{
  channelsMap::iterator channelsIter = fDecayChannelsMap.find(label);
  if (channelsIter != fDecayChannelsMap.end())
  {
    G4String errMsg = "Excited state" + label
                      + " already registered in the decay table.";
    G4Exception("G4MolecularDecayTable::AddExcitedState",
                "G4MolecularDecayTable003", FatalErrorInArgument, errMsg);
    return;
  }
  fDecayChannelsMap[label];
}

void G4MolecularDissociationTable::AddeConfToExcitedState(const G4String& label,
                                                          const G4ElectronOccupancy& conf)
{
  statesMap::iterator statesIter = fExcitedStatesMap.find(conf);

  if (statesIter == fExcitedStatesMap.end())
  {
    fExcitedStatesMap[conf] = label;
  }
  else
  {
    G4Exception(
        "G4MolecularDecayTable::AddExcitedState", "G4MolecularDecayTable004",
        FatalErrorInArgument,
        "Electronic configuration already registered in the decay table");
  }
}

void G4MolecularDissociationTable::AddDecayChannel(const G4String& label,
                                                   const G4MolecularDissociationChannel* channel)
{
  fDecayChannelsMap[label].push_back(channel);
}

void G4MolecularDissociationTable::CheckDataConsistency()
{
  channelsMap::iterator channelsIter;

  //Let's check probabilities

  for (channelsIter = fDecayChannelsMap.begin();
      channelsIter != fDecayChannelsMap.end(); channelsIter++)
  {

    vector<const G4MolecularDissociationChannel*>& decayVect = channelsIter
        ->second;
    G4double sum = 0;

    G4double max = decayVect.size();

    for (size_t i = 0; i < max; i++)
    {
      const G4MolecularDissociationChannel* decay = decayVect[i];
      const G4double prob = decay->GetProbability();
      sum += prob;
    }

    if (sum != 1)
    {
      G4String errMsg = "Deexcitation Channels probabilities in "
          + channelsIter->first + "excited state don't sum up to 1";
      G4Exception("G4MolecularDecayTable::CheckDataConsistency",
                  "G4MolecularDecayTable005", FatalErrorInArgument, errMsg);
    }
  }

}

