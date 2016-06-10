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
// $Id: G4MolecularDissociationTable.cc 93883 2015-11-03 08:25:04Z gcosmo $
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
#include "G4MolecularConfiguration.hh"

using namespace std;
using namespace G4DNA;

//______________________________________________________________________________

G4MolecularDissociationTable::G4MolecularDissociationTable()
{
  ;
}

//______________________________________________________________________________

G4MolecularDissociationTable::~G4MolecularDissociationTable()
{
  ChannelMap::iterator it_map = fDissociationChannels.begin();

  for (; it_map != fDissociationChannels.end(); it_map++)
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
  fDissociationChannels.clear();
}

//______________________________________________________________________________

G4MolecularDissociationTable::
  G4MolecularDissociationTable(const G4MolecularDissociationTable& right)
{
  *this = right;
}

//______________________________________________________________________________

G4MolecularDissociationTable&
G4MolecularDissociationTable::operator=(const G4MolecularDissociationTable& right)
{
  if(this == &right) return *this;
  fDissociationChannels = right.fDissociationChannels;
  return *this;
}

//______________________________________________________________________________

const vector<const G4MolecularDissociationChannel*>*
G4MolecularDissociationTable::
  GetDecayChannels(const G4MolecularConfiguration* conf) const
{
  ChannelMap::const_iterator it_exstates = fDissociationChannels.find(conf);
  if (it_exstates == fDissociationChannels.end()) return 0;
  return &(it_exstates->second);
}

//______________________________________________________________________________

const vector<const G4MolecularDissociationChannel*>*
G4MolecularDissociationTable::GetDecayChannels(const G4String& exState) const
{
  for(ChannelMap::const_iterator it = fDissociationChannels.begin() ;
      it!=fDissociationChannels.end() ; ++it
      )
  {
    if(it->first->GetLabel() == exState) return &(it->second);
  }
  return 0;
}

//______________________________________________________________________________

//void G4MolecularDissociationTable::
//  AddExcitedState(const G4String& label,
//                  const G4MolecularConfiguration* molConf);
//{
//
//}

//______________________________________________________________________________

void G4MolecularDissociationTable::
  AddChannel(const G4MolecularConfiguration* molConf,
             const G4MolecularDissociationChannel* channel)
{
  fDissociationChannels[molConf].push_back(channel);
}

//______________________________________________________________________________

void G4MolecularDissociationTable::CheckDataConsistency() const
{
  ChannelMap::const_iterator channelsIter;

  for(channelsIter = fDissociationChannels.begin();
      channelsIter != fDissociationChannels.end(); ++channelsIter)
  {

    const vector<const G4MolecularDissociationChannel*>& decayVect =
        channelsIter->second;
    G4double sum = 0;

    G4double max = decayVect.size();

    for(size_t i = 0; i < max; i++)
    {
      const G4MolecularDissociationChannel* decay = decayVect[i];
      const G4double prob = decay->GetProbability();
      sum += prob;
    }

    if(sum != 1)
    {
      G4ExceptionDescription errMsg;
      errMsg << "The probabilities for deecitation of molecular configuration "
             << channelsIter->first->GetName()
             << " with label :" <<  channelsIter->first->GetLabel()
             << " don't sum up to 1";
      G4Exception("G4MolecularDissociationTable::CheckDataConsistency",
                  "BRANCHING_RATIOS_CONSISTENCY",
                  FatalErrorInArgument,
                  errMsg);
    }
  }
}

void G4MolecularDissociationTable::Serialize(std::ostream& /*char_traits*/)
{
  // TODO
}
