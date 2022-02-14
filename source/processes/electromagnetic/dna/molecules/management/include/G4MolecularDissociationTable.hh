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
// Contact: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 
//
// ----------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Alfonso Mantero 4 Mar 2009
// ----------------------------------------------------------------
//

#ifndef G4MolecularDecayTable_h
#define G4MolecularDecayTable_h 1

#include <map>
#include <vector>
#include "globals.hh"

class G4MolecularDissociationChannel;
class G4MolecularConfiguration;

//______________________________________________________________________________

namespace G4DNA
{
typedef std::map<const G4MolecularConfiguration*,
    std::vector<const G4MolecularDissociationChannel*> > ChannelMap;
}
//______________________________________________________________________________

/** Class Description
 *  G4MolecularDecayTable operates as a container of deexcitation modes
 *  for excited or ionized molecules
 */

class G4MolecularDissociationTable
{

public:

  G4MolecularDissociationTable();
  ~G4MolecularDissociationTable();
  G4MolecularDissociationTable(const G4MolecularDissociationTable&);
  G4MolecularDissociationTable & operator=(const G4MolecularDissociationTable &right);

public:

  //____________________________________________________________________________
  // Create the table
  // methods to construct the table "interactively"

  void AddChannel(const G4MolecularConfiguration* molConf,
                  const G4MolecularDissociationChannel* channel);

  void CheckDataConsistency() const;
  // Checks that probabilities sum up to 100% for each excited state

  //____________________________________________________________________________
  // Get methods to retrieve data

  const std::vector<const G4MolecularDissociationChannel*>*
    GetDecayChannels(const G4MolecularConfiguration*) const;
  const std::vector<const G4MolecularDissociationChannel*>*
    GetDecayChannels(const G4String& excitedStateLabel) const;

  void Serialize(std::ostream&);

private:
  G4DNA::ChannelMap fDissociationChannels;
};

#endif

