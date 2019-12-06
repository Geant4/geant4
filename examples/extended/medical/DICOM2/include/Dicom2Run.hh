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
/// \file Dicom2Run.hh
/// \brief Definition of the Dicom2Run class
//

#ifndef Dicom2Run_h
#define Dicom2Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "DicomRun.hh"

#include "G4THitsMap.hh"
#include "G4THitsVector.hh"
#include "G4StatAnalysis.hh"

#include <vector>

//---------------------------------------------------------------------
/// Dicom2Run class
///
/// Demonstrates how to use G4StatAnalysis as a vector of objects
/// (instead of vector of pointers to objects) to reduce memory
/// consumption
//---------------------------------------------------------------------

class Dicom2Run : public DicomRun
{
  public:

    typedef G4VTHitsVector<G4StatAnalysis, std::vector<G4StatAnalysis>> Dicom2RunVector;

    // constructor and destructor.
    //  vector of multifunctionaldetector name has to given to constructor.
    Dicom2Run();
    Dicom2Run(const std::vector<G4String> mfdName);
    virtual ~Dicom2Run();

    // virtual method from G4Run.
    // The method is overriden in this class for scoring.
    virtual void RecordEvent(const G4Event*);

    // Access methods for scoring information.
    // - Number of HitsMap for this RUN.
    //   This is equal to number of collections.
    size_t GetNumberOfHitsMap() const { return fRunMap.size(); }
    // - Get HitsMap of this RUN.
    //   by sequential number, by multifucntional name and collection name,
    //   and by collection name with full path.
    Dicom2RunVector* GetHitsVector(G4int i) const {return fRunMap[i];}
    Dicom2RunVector* GetHitsVector(const G4String& detName,
                                   const G4String& colName) const;
    Dicom2RunVector* GetHitsVector(const G4String& fullName) const;

    void ConstructMFD(const std::vector<G4String>&);

    virtual void Merge(const G4Run*);

  private:

    std::vector<G4String> fCollName;
    std::vector<G4int> fCollID;
    std::vector<Dicom2RunVector*> fRunMap;
};

#endif
