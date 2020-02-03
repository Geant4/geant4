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
/// \file Dicom2Run.cc
/// \brief Implementation of the Dicom2Run class

//=====================================================================
///
///  (Description)
///    Dicom2Run Class is for accumulating scored quantities which is
///  scored using G4MutiFunctionalDetector and G4VPrimitiveScorer.
///  Accumulation is done using G4THitsVector object.
///
///    The constructor Dicom2Run(const std::vector<G4String> mfdName)
///  needs a vector filled with MultiFunctionalDetector names which
///  was assigned at instantiation of MultiFunctionalDetector(MFD).
///  Then Dicom2Run constructor automatically scans primitive scorers
///  in the MFD, and obtains collectionIDs of all collections associated
///  to those primitive scorers. Futhermore, the G4THitsVector objects
///  for accumulating during a RUN are automatically created too.
///  (*) Collection Name is same as primitive scorer name.
///
///    The resultant information is kept inside Dicom2Run objects as
///  data members.
///  std::vector<G4String> fCollName;            // Collection Name,
///  std::vector<G4int> fCollID;                 // Collection ID,
///  std::vector<Dicom2RunVector*> fRunMap; // HitsVector for RUN.
///
///  The resualtant HitsVector objects are obtain using access method,
///  GetHitsVector(..).
///
//=====================================================================

#include "Dicom2Run.hh"
#include "DicomDetectorConstruction.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include <cstdint>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Constructor.
Dicom2Run::Dicom2Run()
: DicomRun()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Constructor.
//   (The vector of MultiFunctionalDetector name has to given.)
Dicom2Run::Dicom2Run(const std::vector<G4String> mfdName)
: DicomRun()
{
    ConstructMFD(mfdName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Destructor
//    clear all data members.
Dicom2Run::~Dicom2Run()
{
    //--- Clear HitsVector for RUN
    for(std::size_t i = 0; i < fRunMap.size(); ++i)
    {
         if(fRunMap[i])
            fRunMap[i]->clear();
    }
    fCollName.clear();
    fCollID.clear();
    fRunMap.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Destructor
//    clear all data members.
void Dicom2Run::ConstructMFD(const std::vector<G4String>& mfdName)
{
    DicomRun::ConstructMFD(mfdName);

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //=================================================
    //  Initalize RunMaps for accumulation.
    //  Get CollectionIDs for HitCollections.
    //=================================================
    for(std::size_t idet = 0; idet < mfdName.size(); ++idet)
    {
        // Loop for all MFD.
        G4String detName = mfdName[idet];
        //--- Seek and Obtain MFD objects from SDmanager.
        G4MultiFunctionalDetector* mfd =
                (G4MultiFunctionalDetector*)
                (SDman->FindSensitiveDetector(detName));
        //
        if(mfd)
        {
            //--- Loop over the registered primitive scorers.
            for (G4int icol = 0; icol < mfd->GetNumberOfPrimitives(); ++icol)
            {
                // Get Primitive Scorer object.
                G4VPrimitiveScorer* scorer = mfd->GetPrimitive(icol);
                // collection name and collectionID for HitsCollection,
                // where type of HitsCollection is G4THitsVector in case
                // of primitive scorer.
                // The collection name is given by :
                //  <MFD name>/<Primitive Scorer name>.
                G4String collectionName = scorer->GetName();
                G4String fullCollectionName = detName+"/"+collectionName;
                G4int    collectionID = SDman->GetCollectionID(fullCollectionName);
                //
                if(collectionID >= 0)
                {
                    G4cout << "++ "<<fullCollectionName<< " id " << collectionID
                           << G4endl;
                    // Store obtained HitsCollection information into data
                    // members. qnd creates new G4THitsVector for accumulating
                    // quantities during RUN.
                    fCollName.push_back(fullCollectionName);
                    fCollID.push_back(collectionID);
                    fRunMap.push_back(new Dicom2RunVector(detName,
                                                          collectionName));
                }
                else
                {
                    G4cout << "** collection " << fullCollectionName << " not found. "
                           <<G4endl;
                }
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void Dicom2Run::RecordEvent(const G4Event* aEvent)
{
    DicomRun::RecordEvent(aEvent);

    //G4cout << "Dicom Run :: Recording event " << aEvent->GetEventID()
    //<< "..." << G4endl;
    //=============================
    // HitsCollection of This Event
    //============================
    G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
    if (!HCE)
        return;

    //=======================================================
    // Sum up HitsVector of this Event  into HitsVector of this RUN
    //=======================================================
    for(std::size_t i = 0; i < fCollID.size(); ++i)
    {
        // Loop over HitsCollection
        G4THitsMap<G4double>* EvtMap = nullptr;
        if(fCollID[i] >= 0)
        {
            // Collection is attached to HCE
            EvtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID[i]));
        }
        else
        {
            G4cout <<" Error EvtMap Not Found "<< i << G4endl;
        }

        // if valid pointer, add the pointer
        if(EvtMap)
        {
            //for(auto itr = EvtMap->begin(); itr != EvtMap->end(); ++itr)
            //    G4cout << "adding " << *EvtMap->GetObject(itr) << G4endl;
            //=== Sum up HitsVector of this event to HitsVector of RUN.===
            *fRunMap[i] += *EvtMap;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Merge hits map from threads
void Dicom2Run::Merge(const G4Run* aRun)
{
    DicomRun::Merge(aRun);

    const Dicom2Run* localRun = static_cast<const Dicom2Run*>(aRun);

    Copy(fCollName, localRun->fCollName);
    Copy(fCollID, localRun->fCollID);
    G4int ncopies = G4int(Copy(fRunMap, localRun->fRunMap));
    // copy function returns the fRunMap size if all data is copied
    // so this loop isn't executed the first time around
    for(G4int i = ncopies; i < G4int(fRunMap.size()); ++i)
        *fRunMap[i] += *localRun->fRunMap[i];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//=================================================================
//  Access method for HitsVector of the RUN
//
//-----
// Access HitsVector.
//  By  MultiFunctionalDetector name and Collection Name.
Dicom2Run::Dicom2RunVector*
Dicom2Run::GetHitsVector(const G4String& detName,
                         const G4String& colName) const
{
    G4String fullName = detName+"/"+colName;
    return GetHitsVector(fullName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Access HitsVector.
//  By full description of collection name, that is
//    <MultiFunctional Detector Name>/<Primitive Scorer Name>
Dicom2Run::Dicom2RunVector*
Dicom2Run::GetHitsVector(const G4String& fullName) const
{

    std::size_t Ncol = fCollName.size();
    for(std::size_t i = 0; i < Ncol; ++i)
    {
        if(fCollName[i] == fullName)
        {
            return fRunMap[i];
        }
    }

    G4Exception("Dicom2Run", fullName.c_str(), JustWarning,
                "GetHitsVector failed to locate the requested HitsVector");
    return nullptr;
}
