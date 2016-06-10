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
#include "G4ElectronOccupancy.hh"

class G4MolecularDecayChannel;

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
struct eOccComp
{
    bool operator() (const G4ElectronOccupancy& occ1, const G4ElectronOccupancy& occ2) const
    {

        if (occ1.GetTotalOccupancy() != occ2.GetTotalOccupancy())
        {
            return occ1.GetTotalOccupancy()<occ2.GetTotalOccupancy();
        }
        else
        {

            for (G4int i=0; i<occ1.GetSizeOfOrbit();)
            {

                if (occ1.GetOccupancy(i) != occ2.GetOccupancy(i))
                {
                    return occ1.GetOccupancy(i) < occ2.GetOccupancy(i);
                }
                else
                {

                    i++;
                    if (i >= occ1.GetSizeOfOrbit()) return false;
                }

            }
        }
        return false;
    }
};

//°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
typedef std::map<G4ElectronOccupancy, G4String, eOccComp>  statesMap;
typedef std::map<G4String, std::vector<const G4MolecularDecayChannel*>, std::less <G4String> > channelsMap;

/** Class Description
 *  G4MolecularDecayTable operates as a container of deexcitation modes
 *  for excited or ionized molecules
 */

class G4MolecularDecayTable
{

public:

    G4MolecularDecayTable();
    ~G4MolecularDecayTable();
    G4MolecularDecayTable(const G4MolecularDecayTable&);
    G4MolecularDecayTable & operator=(const G4MolecularDecayTable &right);

public:

    //°°°°°°°°°°°°°°°°°°
    // Create the table
    ////
    // methods to construct the table "interactively"
    void  AddExcitedState(const G4String&); // creates an empty excited state

    // creates or adds to an existing excited state an electronic configuration
    void  AddeConfToExcitedState(const G4String&,const G4ElectronOccupancy&);
    void  AddDecayChannel(const G4String&,const G4MolecularDecayChannel*);
    void  CheckDataConsistency() ;
    // Checks that probabilities sum up to 100% for each excited state

//------------Inline fuctions------------
// Get methods to retrieve data

    const std::vector<const G4MolecularDecayChannel*>* GetDecayChannels(const G4ElectronOccupancy*) const ;
    const std::vector<const G4MolecularDecayChannel*>* GetDecayChannels(const G4String&) const ;

    const G4String& GetExcitedState(const G4ElectronOccupancy*) const ;
    const G4ElectronOccupancy& GetElectronOccupancy (const G4String&) const ;

    inline const statesMap& GetExcitedStateMaps() const;
    inline const channelsMap& GetDecayChannelsMap() const;

private:
    statesMap fExcitedStatesMap ;
    channelsMap fDecayChannelsMap ;
};

inline const statesMap&  G4MolecularDecayTable::GetExcitedStateMaps() const
{
    return fExcitedStatesMap;
}

inline const channelsMap&  G4MolecularDecayTable::GetDecayChannelsMap() const
{
    return fDecayChannelsMap;
}

#endif








