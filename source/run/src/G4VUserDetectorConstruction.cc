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
// $Id$
//

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VUserParallelWorld.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VSensitiveDetector.hh"
#include "G4FieldManager.hh"
#include "G4SDManager.hh"
#include <assert.h>

G4VUserDetectorConstruction::G4VUserDetectorConstruction()
{;}

G4VUserDetectorConstruction::~G4VUserDetectorConstruction()
{;}

void G4VUserDetectorConstruction::RegisterParallelWorld(G4VUserParallelWorld* aPW)
{
  std::vector<G4VUserParallelWorld*>::iterator pwItr;
  for(pwItr=parallelWorld.begin();pwItr!=parallelWorld.end();pwItr++)
  {
    if((*pwItr)->GetName()==aPW->GetName())
    {
      G4String eM = "A parallel world <";
      eM += aPW->GetName();
      eM += "> is already registered to the user detector construction.";
      G4Exception("G4VUserDetectorConstruction::RegisterParallelWorld",
                  "Run0051",FatalErrorInArgument,eM);
    }
  }
  parallelWorld.push_back(aPW);
}

G4int G4VUserDetectorConstruction::ConstructParallelGeometries()
{
  G4int nP = 0;
  std::vector<G4VUserParallelWorld*>::iterator pwItr;
  for(pwItr=parallelWorld.begin();pwItr!=parallelWorld.end();pwItr++)
  {
    (*pwItr)->Construct();
    nP++;
  }
  return nP;
}

G4int G4VUserDetectorConstruction::GetNumberOfParallelWorld() const
{ return parallelWorld.size(); }

G4VUserParallelWorld* G4VUserDetectorConstruction::GetParallelWorld(G4int i) const
{
  if(i<0||i>=GetNumberOfParallelWorld()) return 0;
  return parallelWorld[i];
}

void G4VUserDetectorConstruction::ConstructSDandField()
{
    G4ExceptionDescription msg;
	msg << "User-derived class does not implement ConstructSDandFiled method:\n"
    << "workers will not have SD and fields!\n"
    << "If user defined SD and Field classes implement cloning, you can\n"
    << "re-implement this method as:\n"
    << "void UserDerivedClass::ConstructSDandField() { CloneSD(); CloneF(); }\n";
    G4Exception("G4VUserDetectorConstruction::ConstructSDandField", "Run0052", JustWarning,msg);
}

#include <map>
//TODO: Evaluate if we want to delegate this functionality of cloning to the correc class (the G4FieldManager and/or G4LogicalVolume
//      The issue is that it is G4LogicalVolume class that knows what are the data memebers that need to be cloned/prepared
void G4VUserDetectorConstruction::CloneF()
{
    //TODO: For moment G4FieldManager is per thread variable (i.e. belongs to the splitted part
    //of G4LogivalVolume. However we may want to review this and share G4FieldManager and have
    //Thread private only the Field itself and the mutable parts
    typedef std::map<G4FieldManager*,G4FieldManager*> FMtoFMmap;
    typedef std::pair<G4FieldManager*,G4FieldManager*> FMpair;
    FMtoFMmap masterToWorker;
    G4LogicalVolumeStore* const logVolStore = G4LogicalVolumeStore::GetInstance();
    assert( logVolStore != NULL );
    for ( G4LogicalVolumeStore::const_iterator it = logVolStore->begin() ; it != logVolStore->end() ; ++it )
    {
        G4LogicalVolume *g4LogicalVolume = *it;
        //Use shadow of master to get instance of FM
        //TODO: understand logic with Xin
        //TODO: use getter when available
        G4FieldManager* masterFM = 0;//g4LogicalVolume->fFieldManager;
        G4FieldManager* clonedFM = 0;
        if ( masterFM )
        {
            FMtoFMmap::iterator fmFound = masterToWorker.find(masterFM);
            if ( fmFound == masterToWorker.end() )
            {
                //First time we see this SD, let's clone and remember...
                try {
                    std::pair<FMtoFMmap::iterator,bool> insertedEl = masterToWorker.insert( FMpair(masterFM, masterFM->Clone()) );
                    clonedFM = (insertedEl.first)->second;
                }
                catch (...)
                {
                    G4ExceptionDescription msg;
                    msg << "Cloning of G4FieldManager failed."
                    << " But derived class does not implement cloning. Cannot continue.";
                    G4Exception("G4VUserDetectorConstruction::CloneSD", "Run0053", FatalException,msg);

                }
            }
            else
            {
                // We have already seen this SD attached to a fifferent LogicalVolume, let's re-use previous clone
                clonedFM = (*fmFound).second;
            }
        }// masterFM != 0
        //Note that we do not push FM to doughters (false argument), however, since we area looping on all
        //logical volumes and we implemented the "trick" of the map master<->cloned the final
        //effect is the same as using here the correct boolean flag: log-volumes that originally were sharing
        //the same FM they will have cloned ones
        g4LogicalVolume->SetFieldManager(clonedFM, false);
    }
}

void G4VUserDetectorConstruction::CloneSD()
{
    //Loop on ALL logial volumes to search for attached SD
    G4LogicalVolumeStore* const logVolStore = G4LogicalVolumeStore::GetInstance();
    assert( logVolStore != NULL );
    
    typedef std::map<G4VSensitiveDetector*,G4VSensitiveDetector*> SDtoSDmap;
    typedef std::pair<G4VSensitiveDetector*,G4VSensitiveDetector*> SDpair;
    SDtoSDmap masterToWorker;

    for ( G4LogicalVolumeStore::const_iterator it = logVolStore->begin() ; it != logVolStore->end() ; ++it )
    {
        G4LogicalVolume *g4LogicalVolume = *it;
        //Use shadow of master to get the instance of SD
        //TODO: understand logic with Xin
        //TODO: use getter when available
        G4VSensitiveDetector* masterSD = 0;//g4LogicalVolume->fSensitiveDetector;
        G4VSensitiveDetector* clonedSD = 0;
        if ( masterSD )
        {
            SDtoSDmap::iterator sdFound = masterToWorker.find(masterSD);
            if ( sdFound == masterToWorker.end() )
            {
                //First time we see this SD, let's clone and remember...
                try {
                    std::pair<SDtoSDmap::iterator,bool> insertedEl = masterToWorker.insert( SDpair(masterSD,masterSD->Clone()) );
                    clonedSD = (insertedEl.first)->second;
                    }
                    catch (...)
                    {
                        G4ExceptionDescription msg;
                        msg << "Cloning of G4VSensitiveDetector requested for:"<<masterSD->GetName()<<" (full path name: "<<masterSD->GetFullPathName()<<")."
                        << " But derived class does not implement cloning. Cannot continue.";
                        G4Exception("G4VUserDetectorConstruction::CloneSD", "Run0053", FatalException,msg);
                    }

            }
            else
            {
                // We have already seen this SD attached to a fifferent LogicalVolume, let's re-use previous clone
                clonedSD = (*sdFound).second;
                
            }
        }// masterSD!=0
        g4LogicalVolume->SetSensitiveDetector(clonedSD);

    }
}

void G4VUserDetectorConstruction::SetSensitiveDetector
(const G4String& logVolName, G4VSensitiveDetector* aSD, G4bool multi)
{ 
  G4bool found = false;
  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++)
  {
    if((*pos)->GetName()==logVolName)
    {
      if(found && !multi)
      {
        G4String eM = "More than one logical volumes of the name <";
        eM += (*pos)->GetName();
        eM += "> are found and thus the sensitive detector <";
        eM += aSD->GetName();
        eM += "> cannot be uniquely assigned.";
        G4Exception("G4VUserDetectorConstruction::SetSensitiveDetector",
                    "Run0052",FatalErrorInArgument,eM);
      }
      found = true;
      G4SDManager::GetSDMpointer()->AddNewDetector(aSD);
      (*pos)->SetSensitiveDetector(aSD);
    }
  } 
  if(!found)
  {
    G4String eM2 = "No logical volume of the name <";
    eM2 += logVolName;
    eM2 += "> is found. The specified sensitive detector <";
    eM2 += aSD->GetName();
    eM2 += "> couldn't be assigned to any volume.";
    G4Exception("G4VUserDetectorConstruction::SetSensitiveDetector",
                    "Run0053",FatalErrorInArgument,eM2);
  }
} 
