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
/// \file parallel/ThreadsafeScorers/src/TSRunAction.cc
/// \brief Implementation of the TSRunAction class
//
//
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "TSRunAction.hh"
#include "TSActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4Timer.hh"

#include "TSRun.hh"
#include "TSDetectorConstruction.hh"
#include "G4StatAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSRunAction::TSRunAction()
: fDetector(TSDetectorConstruction::Instance()),
  fName(fDetector->GetMFDName())
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSRunAction::~TSRunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* TSRunAction::GenerateRun()
{
    return new TSRun(fName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSRunAction::BeginOfRunAction(const G4Run* aRun)
{
    G4int evts_to_process = aRun->GetNumberOfEventToBeProcessed();
    G4RunManager::GetRunManager()->SetPrintProgress((evts_to_process > 100)
                                                    ? evts_to_process/100
                                                    : 1);
    if(IsMaster())
        G4PrintEnv();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSRunAction::EndOfRunAction(const G4Run* aRun)
{

    if(IsMaster())
    {
<<<<<<< HEAD
        G4cout << " ###### EndOfTSRunAction ###### " << G4endl;

        aRun->GetNumberOfEvent();
        std::ofstream fileout;
        G4String fname = "";
        std::stringstream separator;

        separator
        << "============================================================";

        typedef std::set<G4int> IDSet_t;
        IDSet_t IDs;

        //- TSRun object.
        const TSRun* tsRun = static_cast<const TSRun*>(aRun);
        //--- Dump all scored quantities involved in TSRun.

        //---------------------------------------------
        // Dump accumulated quantities for this RUN.
        //---------------------------------------------
        std::vector<G4String> primScorerNames { "EnergyDeposit",
                                                "NumberOfSteps" };
        std::vector<G4String> fnames { "mfd_tl", "mfd_tg" };
        std::vector<G4double> units { CLHEP::eV, CLHEP::keV, 1, 1 };
        std::vector<G4String> unitstr { "keV", "steps" };

        for(unsigned i = 0; i < primScorerNames.size(); ++i)
        {
          for(unsigned j = 0; j < fnames.size(); ++j)
          {
            //----------------------------------------------------------------//
            auto print = [] (std::ostream& fout,
                G4int first, G4double second,
                G4double unit1, G4double unit2, G4String unit2str)
            {
              if(fout)
                  fout <<  first
                  << "     "  << second/unit1
                  << G4endl;

              G4cout
              << "    " << std::setw(10) << first
              << "    " << std::setw(15) << std::setprecision(6)
              << std::fixed << second/unit2 << " " << unit2str
              << G4endl;
              G4cout.unsetf(std::ios::fixed);
            };
            //----------------------------------------------------------------//

=======
      G4cout << " ###### EndOfTSRunAction ###### " << G4endl;

      aRun->GetNumberOfEvent();
      std::ofstream fileout;
      G4String fname = "";
      std::stringstream separator;

      separator
          << "============================================================";

      typedef std::set<G4int> IDSet_t;
      IDSet_t IDs;

      //- TSRun object.
      const TSRun* tsRun = static_cast<const TSRun*>(aRun);
      //--- Dump all scored quantities involved in TSRun.

      //---------------------------------------------
      // Dump accumulated quantities for this RUN.
      //---------------------------------------------
      std::vector<G4String> primScorerNames { "EnergyDeposit",
                                              "NumberOfSteps" };
      std::vector<G4String> fnames { "mfd_tl", "mfd_tg" };
      std::vector<G4double> units { CLHEP::eV, CLHEP::keV, 1, 1 };
      std::vector<G4String> unitstr { "keV", "steps" };

      //----------------------------------------------------------------------//
      // lambda to print double value
      auto print = [] (std::ostream& fout,
              G4int first, G4double second,
              G4double unit1, G4double unit2, G4String unit2str)
      {
          if(fout)
              fout <<  first
                    << "     "  << second/unit1
                    << G4endl;

          G4cout
                  << "    " << std::setw(10) << first
                  << "    " << std::setw(15) << std::setprecision(6)
                  << std::fixed << second/unit2 << " " << unit2str
                  << G4endl;
          G4cout.unsetf(std::ios::fixed);
      };
      //----------------------------------------------------------------------//
      // lambda to print statistics
      auto stat_print = [] (std::ostream& fout,
              G4int first, G4StatAnalysis* stat, G4ConvergenceTester* conv,
              G4double unit1, G4double unit2, G4String unit2str)
      {
          if(!stat || !conv)
              return;
          auto fsecond = (*stat);
          auto psecond = (*stat);
          fsecond /= unit1;
          psecond /= unit2;
          if(fout)
          {
              fout << first << "     "  << fsecond << G4endl;
              conv->ShowResult(fout);
          }
          std::stringstream ss;
          ss << "    " << std::setw(10) << first
             << "    " << std::setw(15) << std::setprecision(6)
             << std::fixed << psecond << " " << unit2str;
          // skip print of ConvergenceTester to stdout
          G4cout << ss.str() << G4endl;
      };
      //----------------------------------------------------------------------//

      for(unsigned i = 0; i < primScorerNames.size(); ++i)
      {
        for(unsigned j = 0; j < fnames.size(); ++j)
        {
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
            fname = fnames.at(j) + "_" + primScorerNames.at(i) + ".out";
            fileout.open(fname);
            G4cout << separator.str() << G4endl;
            G4cout << " opened file " << fname << " for output" << G4endl;
            G4cout << separator.str() << G4endl;

            G4bool valid = true;
            if(j == 0)
            {
<<<<<<< HEAD
              G4THitsMap<G4double>* hitmap
                  = tsRun->GetHitsMap(fName + "/" + primScorerNames.at(i));
              if(hitmap && hitmap->GetMap()->size() != 0)
              {
                for(auto itr = hitmap->GetMap()->begin();
                    itr != hitmap->GetMap()->end(); itr++)
                {
                  IDs.insert(itr->first);
                  std::get<0>(fTypeCompare[primScorerNames.at(i)][itr->first])
                      = *itr->second/units.at(i);
                  print(fileout, itr->first, *itr->second,
                        units.at(i), units.at(i+1), unitstr.at(i));
                }
              } else {
                valid = false;
              }
            } else {
              G4TAtomicHitsMap<G4double>* hitmap
                  = tsRun->GetAtomicHitsMap(fName + "/" +
                                              primScorerNames.at(i));
              if(hitmap && hitmap->GetMap()->size() != 0)
              {
                for(auto itr = hitmap->GetMap()->begin();
                    itr != hitmap->GetMap()->end(); itr++)
                {
                  IDs.insert(itr->first);
                  std::get<1>(fTypeCompare[primScorerNames.at(i)][itr->first])
                      = *itr->second/units.at(i);
                  print(fileout, itr->first, *itr->second,
                        units.at(i), units.at(i+1), unitstr.at(i));
                }
              } else {
                valid = false;
              }
              if(!valid)
              {
                G4Exception("TSRunAction", "000", JustWarning,
                            G4String(primScorerNames.at(i) +
                            " HitsMap is either not "
                            "created or the HitsMap was empty").c_str());
              }
            }

            fileout.close();
            G4cout << separator.str() << G4endl;
            G4cout << " closed file " << fname << " for output" << G4endl;
          }
          // add the mutex data
          TSRun::MutexHitsMap_t* hitmap
              = tsRun->GetMutexHitsMap(fName + "/" +
                                         primScorerNames.at(i));
          if(hitmap && hitmap->size() != 0)
          {
            for(auto itr = hitmap->begin();
                itr != hitmap->end(); itr++)
            {
              IDs.insert(itr->first);
              std::get<2>(fTypeCompare[primScorerNames.at(i)][itr->first])
                  = itr->second/units.at(i);
=======
                G4THitsMap<G4double>* hitmap
                        = tsRun->GetHitsMap(fName + "/" + primScorerNames.at(i));
                G4StatContainer<G4StatAnalysis>* statmap
                        = tsRun->GetStatMap(fName + "/" + primScorerNames.at(i));
                G4StatContainer<G4ConvergenceTester>* convmap
                        = tsRun->GetConvMap(fName + "/" + primScorerNames.at(i));

                if(hitmap && hitmap->size() != 0)
                {
                    for(auto itr = hitmap->begin(); itr != hitmap->end(); itr++)
                    {
                        if(!hitmap->GetObject(itr))
                            continue;
                        IDs.insert(itr->first);
                        std::get<0>(fTypeCompare[primScorerNames.at(i)][itr->first])
                                = *itr->second/units.at(i);
                        print(fileout, itr->first, *itr->second,
                              units.at(i), units.at(i+1), unitstr.at(i));
                    }
                }
                else
                {
                    valid = false;
                }

                if(statmap && statmap->size() != 0 &&
                   convmap && convmap->size() != 0)
                {
                    auto stat_fname = "stat_" + fname;
                    std::ofstream statout;
                    statout.open(stat_fname);
                    for(auto itr = statmap->begin(); itr != statmap->end(); itr++)
                    {
                        G4int _f = statmap->GetIndex(itr);
                        G4StatAnalysis* _s = statmap->GetObject(itr);
                        G4ConvergenceTester* _c = convmap->GetObject(_f);
                        stat_print(statout, _f, _s, _c,
                                   units.at(i), units.at(i+1), unitstr.at(i));
                    }
                    statout.close();
                }
                else
                {
                    std::stringstream ss;
                    ss << " StatMap/ConvMap is either not "
                       << "created or the StatMap/ConvMap was empty";
                    if(statmap)
                        ss << " (StatMap size == " << statmap->size() << ")";
                    if(convmap)
                        ss << " (ConvMap size == " << convmap->size() << ")";

                    G4Exception("TSRunAction", "002", JustWarning,
                                G4String(primScorerNames.at(i) +
                                         ss.str()).c_str());
                }

                if(!valid)
                {
                    G4Exception("TSRunAction", "000", JustWarning,
                                G4String(primScorerNames.at(i) +
                                         " HitsMap is either not "
                                         "created or the HitsMap was empty").c_str());
                }
            }
            else
            {
                G4TAtomicHitsMap<G4double>* hitmap
                        = tsRun->GetAtomicHitsMap(fName + "/" +
                                                  primScorerNames.at(i));
                if(hitmap && hitmap->size() != 0)
                {
                    for(auto itr = hitmap->begin(); itr != hitmap->end(); itr++)
                    {
                        IDs.insert(itr->first);
                        std::get<1>(fTypeCompare[primScorerNames.at(i)][itr->first])
                                = *itr->second/units.at(i);
                        print(fileout, itr->first, *itr->second,
                              units.at(i), units.at(i+1), unitstr.at(i));
                    }
                }
                else
                {
                    valid = false;
                }

                if(!valid)
                {
                    G4Exception("TSRunAction", "001", JustWarning,
                                G4String(primScorerNames.at(i) +
                                         " HitsMap is either not "
                                         "created or the HitsMap was empty").c_str());
                }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
            }

<<<<<<< HEAD
        }

        //--------------------------------------------------------------------//
        // Check that the values are equivalent and there are no
        // IDs in one container that aren't in another
        //--------------------------------------------------------------------//

        fname = "mfd_diff.out";
        fileout.open(fname);

        G4cout << separator.str() << G4endl;
        G4cout << " opened file " << fname << " for difference output" << G4endl;
        G4cout << separator.str() << G4endl;

        fileout
        << "    " << std::setw(10) << "ID"
        << "    "
        << std::setw(30) << std::setprecision(12) << std::fixed
        << "MFD value"
        << "    "
        << std::setw(30) << std::setprecision(12) << std::fixed
        << "Atomic Hits Map value"
        << "    "
        << std::setw(30) << std::setprecision(12) << std::fixed
        << "Difference"
        << "    "
        << std::setw(30) << std::setprecision(12) << std::fixed
        << "Diff (MFD - MUTEXED)"
        << "    "
        << std::setw(30) << std::setprecision(12) << std::fixed
        << "Diff (ATOM_HIT_MAP - MUTEXED)"
        << G4endl << G4endl;

        for(auto itr1 = fTypeCompare.begin();
            itr1 != fTypeCompare.end(); ++itr1)
        {
=======
            fileout.close();
            G4cout << separator.str() << G4endl;
            G4cout << " closed file " << fname << " for output" << G4endl;
        }
        // add the mutex data
        TSRun::MutexHitsMap_t* hitmap
                = tsRun->GetMutexHitsMap(fName + "/" +
                                         primScorerNames.at(i));
        if(hitmap && hitmap->size() != 0)
        {
            for(auto itr = hitmap->begin();
                itr != hitmap->end(); itr++)
            {
                IDs.insert(itr->first);
                std::get<2>(fTypeCompare[primScorerNames.at(i)][itr->first])
                        = itr->second/units.at(i);
            }
        }

      }

      //--------------------------------------------------------------------//
      // Check that the values are equivalent and there are no
      // IDs in one container that aren't in another
      //--------------------------------------------------------------------//

      fname = "mfd_diff.out";
      fileout.open(fname);

      G4cout << separator.str() << G4endl;
      G4cout << " opened file " << fname << " for difference output" << G4endl;
      G4cout << separator.str() << G4endl;

      fileout << "    " << std::setw(10) << "ID"
              << "    "
              << std::setw(30) << std::setprecision(12) << std::fixed
              << "MFD value"
              << "    "
              << std::setw(30) << std::setprecision(12) << std::fixed
              << "Atomic Hits Map value"
              << "    "
              << std::setw(30) << std::setprecision(8) << std::scientific
              << "Difference"
              << "    "
              << std::setw(30) << std::setprecision(8) << std::scientific
              << "Diff (MFD - MUTEXED)"
              << "    "
              << std::setw(30) << std::setprecision(8) << std::scientific
              << "Diff (ATOM_HIT_MAP - MUTEXED)"
              << G4endl << G4endl;

      for(auto itr1 = fTypeCompare.begin();
          itr1 != fTypeCompare.end(); ++itr1)
      {
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
          fileout << "\n\nType = " << itr1->first << "\n" << G4endl;
          for(auto itr2 = itr1->second.begin();
              itr2 != itr1->second.end(); ++itr2)
          {
<<<<<<< HEAD
            fileout
            << "    " << std::setw(10) << itr2->first
            << "    "
            << std::setw(30) << std::setprecision(12) << std::fixed
            << std::get<0>(itr2->second)
            << "    "
            << std::setw(30) << std::setprecision(12) << std::fixed
            << std::get<1>(itr2->second)
            << "    "
            << std::setw(30) << std::setprecision(12) << std::fixed
            << (std::fabs(std::get<0>(itr2->second) - std::get<1>(itr2->second)))
            << "    "
            << std::setw(30) << std::setprecision(12) << std::fixed
            << (std::fabs(std::get<0>(itr2->second) - std::get<2>(itr2->second)))
            << "    "
            << std::setw(30) << std::setprecision(12) << std::fixed
            << (std::fabs(std::get<1>(itr2->second) - std::get<2>(itr2->second)))
            << G4endl;

          }
        }
=======
              G4double d01
                      = std::fabs(std::get<0>(itr2->second) -
                                  std::get<1>(itr2->second));
              G4double d02
                      = std::fabs(std::get<0>(itr2->second) -
                                  std::get<2>(itr2->second));
              G4double d03
                      = std::fabs(std::get<1>(itr2->second) -
                                  std::get<2>(itr2->second));


              auto _print_diff = [&] (const G4double& _dval)
              {
                  if(_dval > 0.0)
                      fileout << std::setprecision(8) << std::scientific
                              << std::setw(30) << _dval << "    ";
                  else
                      fileout << std::setprecision(1) << std::fixed
                              << std::setw(30) << _dval << "    ";
              };

              fileout
                      << "    " << std::setw(10) << itr2->first
                      << "    "
                      << std::setw(30) << std::setprecision(12) << std::fixed
                      << std::get<0>(itr2->second)
                      << "    "
                      << std::setw(30) << std::setprecision(12) << std::fixed
                      << std::get<1>(itr2->second)
                      << "    ";

              _print_diff(d01);
              _print_diff(d02);
              _print_diff(d03);

              fileout << G4endl;

          }
      }
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

        fileout.close();
        G4cout << " closed file " << fname << " for difference output" << G4endl;
        G4cout << separator.str() << G4endl;

    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




