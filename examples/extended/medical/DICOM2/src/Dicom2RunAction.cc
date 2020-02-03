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
/// \file Dicom2RunAction.cc
/// \brief Implementation of the Dicom2RunAction class
//

#include "Dicom2RunAction.hh"
#include "Dicom2Run.hh"

//-- In order to obtain detector information.
#include <fstream>
#include <iomanip>
#include "G4THitsMap.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4StatAnalysis.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/// Constructor
Dicom2RunAction::Dicom2RunAction()
: DicomRunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/// Destructor.
Dicom2RunAction::~Dicom2RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Run* Dicom2RunAction::GenerateRun()
{
    // Generate new RUN object, which is specially
    // dedicated for MultiFunctionalDetector scheme.
    //  Detail description can be found in Dicom2Run.hh/cc.
    return fDcmrun = new Dicom2Run(fSDName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dicom2RunAction::EndOfRunAction(const G4Run* aRun)
{
    // Lock the output because of the external calls to DicomRunAction
    // otherwise, the output gets too confusing
    G4AutoLock l(G4TypeMutex<Dicom2RunAction>());

    G4cout << G4endl;
    G4cout << "[==========================================================="
           << " DICOM "
           << "===========================================================]"
           << G4endl;
    G4cout << G4endl;

    DicomRunAction::EndOfRunAction(aRun);

    G4cout << G4endl;
    G4cout << "[==========================================================="
           << " DICOM2 "
           << "==========================================================]"
           << G4endl;
    G4cout << G4endl;

    G4int nofEvents = aRun->GetNumberOfEvent();

    G4StatAnalysis local_total_dose;

    const Dicom2Run* dcm2Run = static_cast<const Dicom2Run*>(aRun);
    //--- Dump all scored quantities involved in Dicom2Run.
    for(uintmax_t i = 0; i < fSDName.size(); i++)
    {
        //
        //---------------------------------------------
        // Dump accumulated quantities for this RUN.
        //  (Display only central region of x-y plane)
        //      0       ConcreteSD/DoseDeposit
        //---------------------------------------------
        Dicom2RunVector* DoseDeposit =
                dcm2Run->GetHitsVector(fSDName[i]+"/DoseDeposit");

        if(DoseDeposit && DoseDeposit->size() != 0 )
        {
            for(auto itr = DoseDeposit->begin(); itr != DoseDeposit->end(); ++itr)
            {
                // this will sometimes return null pointers
                if(!DoseDeposit->GetObject(itr))
                    continue;
                local_total_dose += (*DoseDeposit->GetObject(itr));
            }
        }
    }

    if(IsMaster())
    {
        G4cout << " ###### EndOfRunAction ###### " << G4endl;
        //- Dicom2Run object.
        const Dicom2Run* re02Run = static_cast<const Dicom2Run*>(aRun);
        //--- Dump all scored quantities involved in Dicom2Run.

        for(uintmax_t i = 0; i < fSDName.size(); i++)
        {
            //
            //---------------------------------------------
            // Dump accumulated quantities for this RUN.
            //  (Display only central region of x-y plane)
            //      0       ConcreteSD/DoseDeposit
            //---------------------------------------------
            Dicom2RunVector* DoseDeposit =
                    re02Run->GetHitsVector(fSDName[i]+"/DoseDeposit");

            G4cout << "============================================================="
                   <<G4endl;
            G4cout << " Number of event processed : "
                   << aRun->GetNumberOfEvent() << G4endl;
            G4cout << "============================================================="
                   <<G4endl;

            std::ofstream fileout;
            G4String fname = "dicom2-vector.out";
            fileout.open(fname);
            G4cout << " opened file " << fname << " for dose output" << G4endl;

            if(DoseDeposit && DoseDeposit->size() != 0)
            {
                std::ostream *myout = &G4cout;
                PrintHeader(myout);
                for(auto itr = DoseDeposit->begin(); itr != DoseDeposit->end();
                    ++itr)
                {
                    auto _idx = DoseDeposit->GetIndex(itr);
                    G4StatAnalysis* _stat = DoseDeposit->GetObject(itr);
                    if(_stat && _stat->GetHits() > 0)
                    {
                        G4StatAnalysis _tmp_stat = *_stat;
                        _tmp_stat /= CLHEP::gray;
                        fileout << _idx << "     "  << (*_stat) << G4endl;
                    }
                }
                G4cout << "============================================="<<G4endl;
            }
            else
            {
                G4Exception("Dicom2RunAction", "000", JustWarning,
                            "DoseDeposit HitsMap is either a null pointer "
                            "of the HitsMap was empty");
            }
            fileout.close();
            G4cout << " closed file " << fname << " for dose output" << G4endl;

        }
    }


    if (IsMaster())
    {
        // convert to units of Gy
        local_total_dose /= gray;
        G4cout << "--------------------End of Global Run-----------------------"
               << G4endl;
        G4cout << " The run was " << nofEvents << " events " << G4endl;
        G4cout << "      TOTAL DOSE : \t" << local_total_dose << " Gy" << G4endl;
        if(nofEvents > 0)
        {
            local_total_dose /= nofEvents;
            G4cout << " TOTAL DOSE/Bq-s : \t" << local_total_dose << " Gy/Bq-s"
                   << G4endl;
        }
    }
    else
    {
        // convert to units of Gy
        local_total_dose /= gray;
        G4cout << "--------------------End of Local Run------------------------"
               << G4endl;
        G4cout << " The run was " << nofEvents << " events" << G4endl;
        G4cout << "LOCAL TOTAL DOSE : \t" << local_total_dose << " Gy" << G4endl;
        if(nofEvents > 0)
        {
            local_total_dose /= nofEvents;
            G4cout << " LOCAL DOSE/Bq-s : \t" << local_total_dose << " Gy/Bq-s"
                   << G4endl;
        }
    }

    G4cout << G4endl;
    G4cout << "Finished : End of Run Action " << aRun->GetRunID() << "\n" << G4endl;

}
