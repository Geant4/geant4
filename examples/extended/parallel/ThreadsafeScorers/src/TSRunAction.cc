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
#include "G4TaskRunManager.hh"

#include "TSRun.hh"
#include "TSDetectorConstruction.hh"
#include "G4StatAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSRunAction::TSRunAction()
  : fDetector(TSDetectorConstruction::Instance())
  , fName(fDetector->GetMFDName())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSRunAction::~TSRunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* TSRunAction::GenerateRun() { return new TSRun(fName); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int evts_to_process = aRun->GetNumberOfEventToBeProcessed();
  G4RunManager::GetRunManager()->SetPrintProgress(
    (evts_to_process > 100) ? evts_to_process / 100 : 1);
  if(IsMaster())
    G4PrintEnv();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TSRunAction::EndOfRunAction(const G4Run* aRun)
{
  if(IsMaster())
  {
    G4cout << " ###### EndOfTSRunAction ###### " << G4endl;

    aRun->GetNumberOfEvent();
    std::ofstream fileout;
    G4String fname = "";
    std::stringstream separator;

    separator << "============================================================";

    typedef std::set<G4int> IDSet_t;
    IDSet_t IDs;

    //- TSRun object.
    const TSRun* tsRun = static_cast<const TSRun*>(aRun);
    //--- Dump all scored quantities involved in TSRun.

    //---------------------------------------------
    // Dump accumulated quantities for this RUN.
    //---------------------------------------------
    std::vector<G4String> primScorerNames{ "EnergyDeposit", "NumberOfSteps" };
    std::vector<G4String> fnames{ "mfd_tl", "mfd_tg" };
    std::vector<G4double> units{ CLHEP::eV, CLHEP::keV, 1, 1 };
    std::vector<G4String> unitstr{ "keV", "steps" };

    //----------------------------------------------------------------------//
    // lambda to print double value
    auto print = [](std::ostream& fout, G4int first, G4double second,
                    G4double unit1, G4double unit2, G4String unit2str) {
      if(fout)
        fout << first << "     " << second / unit1 << G4endl;

      G4cout << "    " << std::setw(10) << first << "    " << std::setw(15)
             << std::setprecision(6) << std::fixed << second / unit2 << " "
             << unit2str << G4endl;
      G4cout.unsetf(std::ios::fixed);
    };
    //----------------------------------------------------------------------//
    // lambda to print statistics
    auto stat_print = [](std::ostream& fout, G4int first, G4StatAnalysis* stat,
                         G4ConvergenceTester* conv, G4double unit1,
                         G4double unit2, G4String unit2str) {
      if(!stat || !conv)
        return;
      auto fsecond = (*stat);
      auto psecond = (*stat);
      fsecond /= unit1;
      psecond /= unit2;
      if(fout)
      {
        fout << first << "     " << fsecond << G4endl;
        conv->ShowResult(fout);
      }
      std::stringstream ss;
      ss << "    " << std::setw(10) << first << "    " << std::setw(15)
         << std::setprecision(6) << std::fixed << psecond << " " << unit2str;
      // skip print of ConvergenceTester to stdout
      G4cout << ss.str() << G4endl;
    };
    //----------------------------------------------------------------------//

    for(unsigned i = 0; i < primScorerNames.size(); ++i)
    {
      for(unsigned j = 0; j < fnames.size(); ++j)
      {
        fname = fnames.at(j) + "_" + primScorerNames.at(i) + ".out";
        fileout.open(fname);
        G4cout << separator.str() << G4endl;
        G4cout << " opened file " << fname << " for output" << G4endl;
        G4cout << separator.str() << G4endl;

        G4bool valid = true;
        if(j == 0)
        {
          G4THitsMap<G4double>* hitmap =
            tsRun->GetHitsMap(fName + "/" + primScorerNames.at(i));
          G4StatContainer<G4StatAnalysis>* statmap =
            tsRun->GetStatMap(fName + "/" + primScorerNames.at(i));
          G4StatContainer<G4ConvergenceTester>* convmap =
            tsRun->GetConvMap(fName + "/" + primScorerNames.at(i));

          if(hitmap && hitmap->size() != 0)
          {
            for(auto itr = hitmap->begin(); itr != hitmap->end(); itr++)
            {
              if(!hitmap->GetObject(itr))
                continue;
              IDs.insert(itr->first);
              std::get<0>(fTypeCompare[primScorerNames.at(i)][itr->first]) =
                *itr->second / units.at(i);
              print(fileout, itr->first, *itr->second, units.at(i),
                    units.at(i + 1), unitstr.at(i));
            }
          }
          else
          {
            valid = false;
          }

          if(statmap && statmap->size() != 0 && convmap && convmap->size() != 0)
          {
            auto stat_fname = "stat_" + fname;
            std::ofstream statout;
            statout.open(stat_fname);
            for(auto itr = statmap->begin(); itr != statmap->end(); itr++)
            {
              G4int _f                = statmap->GetIndex(itr);
              G4StatAnalysis* _s      = statmap->GetObject(itr);
              G4ConvergenceTester* _c = convmap->GetObject(_f);
              stat_print(statout, _f, _s, _c, units.at(i), units.at(i + 1),
                         unitstr.at(i));
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
                        G4String(primScorerNames.at(i) + ss.str()).c_str());
          }

          if(!valid)
          {
            G4Exception("TSRunAction", "000", JustWarning,
                        G4String(primScorerNames.at(i) +
                                 " HitsMap is either not "
                                 "created or the HitsMap was empty")
                          .c_str());
          }
        }
        else
        {
          G4TAtomicHitsMap<G4double>* hitmap =
            tsRun->GetAtomicHitsMap(fName + "/" + primScorerNames.at(i));
          if(hitmap && hitmap->size() != 0)
          {
            for(auto itr = hitmap->begin(); itr != hitmap->end(); itr++)
            {
              IDs.insert(itr->first);
              std::get<1>(fTypeCompare[primScorerNames.at(i)][itr->first]) =
                *itr->second / units.at(i);
              print(fileout, itr->first, *itr->second, units.at(i),
                    units.at(i + 1), unitstr.at(i));
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
                                 "created or the HitsMap was empty")
                          .c_str());
          }
        }

        fileout.close();
        G4cout << separator.str() << G4endl;
        G4cout << " closed file " << fname << " for output" << G4endl;
      }
      // add the mutex data
      TSRun::MutexHitsMap_t* hitmap =
        tsRun->GetMutexHitsMap(fName + "/" + primScorerNames.at(i));
      if(hitmap && hitmap->size() != 0)
      {
        for(auto itr = hitmap->begin(); itr != hitmap->end(); itr++)
        {
          IDs.insert(itr->first);
          std::get<2>(fTypeCompare[primScorerNames.at(i)][itr->first]) =
            itr->second / units.at(i);
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
            << "    " << std::setw(30) << std::setprecision(12) << std::fixed
            << "MFD value"
            << "    " << std::setw(30) << std::setprecision(12) << std::fixed
            << "Atomic Hits Map value"
            << "    " << std::setw(30) << std::setprecision(8)
            << std::scientific << "Difference"
            << "    " << std::setw(30) << std::setprecision(8)
            << std::scientific << "Diff (MFD - MUTEXED)"
            << "    " << std::setw(30) << std::setprecision(8)
            << std::scientific << "Diff (ATOM_HIT_MAP - MUTEXED)" << G4endl
            << G4endl;

    //----------------------------------------------------------------------//
    //
    //    Example of using tasking in the user-application. Although this
    //    is sort of trivial case and might not result in any speed-up
    //    it is a good validation test because the order that strings
    //    are joined matters. Although the tasks are executed asynchronously
    //    and may complete at different times, the return values from
    //    the tasks are stored in futures and are "joined" in the order
    //    that they were submitted to the task-group
    //
    //----------------------------------------------------------------------//
    // do not directly call G4TaskManager::GetInstance() as this will generate
    // an instance
    auto tm = dynamic_cast<G4TaskRunManager*>(G4RunManager::GetRunManager());
    // Get the thread-pool if available
    auto tp = (tm) ? tm->GetThreadPool() : nullptr;
    // write a join algorithm which combines the strings from the tasks
    auto join_output = [](std::string& lhs, std::string&& rhs) {
      return (lhs += rhs);
    };

    // this is the outer-loop of tasks
    auto report_type_comparison = [=](const G4String& id,
                                      const IDcompare_t& comp) {
      // the 'report_type_comparison' generates more tasks
      auto report_subtype_comparison = [](const G4int& idx,
                                          const Compare_t& value) {
        std::stringstream streamout;
        G4double d01 = std::fabs(std::get<0>(value) - std::get<1>(value));
        G4double d02 = std::fabs(std::get<0>(value) - std::get<2>(value));
        G4double d03 = std::fabs(std::get<1>(value) - std::get<2>(value));

        auto _print_diff = [&](const G4double& _dval) {
          if(_dval > 0.0)
            streamout << std::setprecision(8) << std::scientific
                      << std::setw(30) << _dval << "    ";
          else
            streamout << std::setprecision(1) << std::fixed << std::setw(30)
                      << _dval << "    ";
        };

        streamout << "    " << std::setw(10) << idx << "    " << std::setw(30)
                  << std::setprecision(12) << std::fixed << std::get<0>(value)
                  << "    " << std::setw(30) << std::setprecision(12)
                  << std::fixed << std::get<1>(value) << "    ";

        _print_diff(d01);
        _print_diff(d02);
        _print_diff(d03);

        streamout << G4endl;
        return streamout.str();
      };

      std::stringstream streamout;
      streamout << "\n\nType = " << id << "\n" << G4endl;
      if(tp)
      {
        // create a task group (nested inside the 'report_type_comparison' task)
        G4TaskGroup<std::string> tg(join_output, tp);
        // create the tasks in the task-group
        for(auto titr = comp.begin(); titr != comp.end(); ++titr)
          tg.exec(report_subtype_comparison, titr->first, titr->second);
        // wait on the tasks to finish and execute the join function
        // this will block the outer task from completing until all the inner
        // tasks have been completed
        streamout << tg.join();
      }
      else
      {
        // if there isn't a tasking thread-pool then we make traditional
        // function call on this thread
        for(auto titr = comp.begin(); titr != comp.end(); ++titr)
          streamout << report_subtype_comparison(titr->first, titr->second);
      }
      // this is the completion of the outer tasks
      return streamout.str();
    };

    G4String tasking_result = "";
    if(tp)
    {
      G4cout << "\n\nGenerating diff output via tasking... ";
      // create a task group to
      G4TaskGroup<std::string> tg(join_output, tp);
      for(auto itr = fTypeCompare.begin(); itr != fTypeCompare.end(); ++itr)
        tg.exec(report_type_comparison, itr->first, itr->second);
      // wait on the tasks to finish and execute the join function
      tasking_result = tg.join();
    }

    // if thread-pool was available, lets validate that tasking did what was
    // expected
    if(tp)
    {
      // generate the output serially
      G4String serial_result = "";
      for(auto itr = fTypeCompare.begin(); itr != fTypeCompare.end(); ++itr)
        serial_result += report_type_comparison(itr->first, itr->second);

      // write the tasking result even if it was bad so that it can viewed
      fileout << tasking_result;

      // compare the strings -- should be the same
      if(serial_result != tasking_result)
      {
        G4Exception("TSRunAction", "003", JustWarning,
                    "Output written via tasking did not match output written "
                    "serially. Appending serial result to output file");
        fileout
          << "\n\n#================CORRECT_SERIAL_OUTPUT================#\n\n";
        fileout << serial_result;
      }
    }
    else
    {
      // if thread-pool was not available, then just write serially
      for(auto itr = fTypeCompare.begin(); itr != fTypeCompare.end(); ++itr)
        fileout << report_type_comparison(itr->first, itr->second);
    }

    fileout.close();
    G4cout << " closed file " << fname << " for difference output" << G4endl;
    G4cout << separator.str() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
