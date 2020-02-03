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
/// \file medical/DICOM/include/DicomRun.hh
/// \brief Definition of the DicomRun class
//

#ifndef DicomRun_h
#define DicomRun_h 1

#include "G4Run.hh"
#include "G4Event.hh"

#include "G4THitsMap.hh"
#include <vector>

//---------------------------------------------------------------------
/// DicomRun class
///
/// Example implementation for multi-functional-detector and 
/// primitive scorer.
/// This DicomRun class has collections which accumulate
/// a event information into a run information.
//---------------------------------------------------------------------

class DicomRun : public G4Run
{

  public:

    // constructor and destructor.
    //  vector of multifunctionaldetector name has to given to constructor.
    DicomRun();
    DicomRun(const std::vector<G4String> mfdName);
    virtual ~DicomRun();

    // virtual method from G4Run. 
    // The method is overriden in this class for scoring.
    virtual void RecordEvent(const G4Event*);

    // Access methods for scoring information.
    // - Number of HitsMap for this RUN. 
    //   This is equal to number of collections.
    size_t GetNumberOfHitsMap() const {return fRunMap.size();}
    // - Get HitsMap of this RUN.
    //   by sequential number, by multifucntional name and collection name,
    //   and by collection name with full path.
    G4THitsMap<G4double>* GetHitsMap(G4int i) const {return fRunMap[i];}
    G4THitsMap<G4double>* GetHitsMap(const G4String& detName, 
                                    const G4String& colName) const;
    G4THitsMap<G4double>* GetHitsMap(const G4String& fullName) const;

    void ConstructMFD(const std::vector<G4String>&);

    virtual void Merge(const G4Run*);

  private:

    std::vector<G4String> fCollName;
    std::vector<G4int> fCollID;
    std::vector<G4THitsMap<G4double>*> fRunMap;
};

//==========================================================================
//          Generic Functions to help with merge
//==========================================================================
template <typename T>
inline void Copy(std::vector<T>& main, const std::vector<T>& data)
{
    for(size_t i = main.size(); i < data.size(); ++i) {
        main.push_back(data.at(i));
    }
}
//==========================================================================
template <typename T>
inline size_t Copy(std::vector<T*>& main, const std::vector<T*>& data)
{
    size_t size_diff = data.size() - main.size();
    for(size_t i = main.size(); i < data.size(); ++i) {
        main.push_back(new T(*data.at(i)));
    }
    return size_diff;
}
//==========================================================================
template <typename T>
inline void Print(const std::vector<T>& data)
{
    G4cout << G4endl;
    for(size_t i = 0; i < data.size(); ++i) {
        G4cout << "\t\t" << i << " \t\t " << data.at(i) << G4endl;
    }
    G4cout << G4endl;
}
//==========================================================================

#endif
