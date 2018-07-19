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
// $Id: G4VReadOutGeometry.cc 94771 2015-12-09 09:44:05Z gcosmo $
//

#include "G4VReadOutGeometry.hh"
#include "G4Navigator.hh"


G4VReadOutGeometry::G4VReadOutGeometry()
  :ROworld(nullptr),fincludeList(nullptr),
   fexcludeList(nullptr),touchableHistory(nullptr)
{
  name = "unknown";
  ROnavigator = new G4Navigator();
  G4ExceptionDescription ed;
  ed<<"The concept and the functionality of Readout Geometry has been merged\n"
    <<"into Parallel World. This G4VReadOutGeometry is kept for the sake of\n"
    <<"not breaking the commonly-used interface in the sensitive detector class.\n"
    <<"But this functionality of G4VReadOutGeometry class is no longer tested\n"
    <<"and thus may not be working well. We strongly recommend our customers to\n"
    <<"migrate to Parallel World scheme.";
  G4Exception("G4VReadOutGeometry","DIGIHIT1001",JustWarning,ed);
}

G4VReadOutGeometry::G4VReadOutGeometry(const G4VReadOutGeometry &right)
{
  fincludeList = nullptr;
  fexcludeList = nullptr;
  name = right.name;
  ROworld = right.ROworld;
  touchableHistory = nullptr;
  ROnavigator = new G4Navigator();
  // COPY CONSTRUCTOR NOT STRAIGHT FORWARD: need to copy the touchabelHistory
  // VALUE, same for navigator and same for the World+Geom hierachy
}

G4VReadOutGeometry::G4VReadOutGeometry(G4String n) 
  :ROworld(nullptr),fincludeList(nullptr),
   fexcludeList(nullptr),name(n),touchableHistory(nullptr)
{
  ROnavigator = new G4Navigator();
  G4ExceptionDescription ed;
  ed<<"The concept and the functionality of Readout Geometry has been merged\n"
    <<"into Parallel World. This G4VReadOutGeometry is kept for the sake of\n"
    <<"not breaking the commonly-used interface in the sensitive detector class.\n"
    <<"But this functionality of G4VReadOutGeometry class is no longer tested\n"
    <<"and thus may not be working well. We strongly recommend our customers to\n"
    <<"migrate to Parallel World scheme.";
  G4Exception("G4VReadOutGeometry","DIGIHIT1001",JustWarning,ed);
}

G4VReadOutGeometry::~G4VReadOutGeometry()
{ 
  //if(ROworld) delete ROworld; //should we do ? will it delete the goem tree also ?
  if(fincludeList)     delete fincludeList;
  if(fexcludeList)     delete fexcludeList;
  if(touchableHistory) delete touchableHistory;
  if(ROnavigator)      delete ROnavigator;
}

G4VReadOutGeometry & G4VReadOutGeometry::operator=(const G4VReadOutGeometry &right)
{
  if ( this == &right ) return *this;
  delete fincludeList; fincludeList     = nullptr;
  delete fexcludeList; fexcludeList     = nullptr;
  name             = right.name;
  ROworld          = right.ROworld;
  delete touchableHistory; touchableHistory = nullptr;
  delete ROnavigator; ROnavigator = new G4Navigator();
  return *this;
}

G4int G4VReadOutGeometry::operator==(const G4VReadOutGeometry &right) const
{ return (this == (G4VReadOutGeometry *) &right); }

G4int G4VReadOutGeometry::operator!=(const G4VReadOutGeometry &right) const
{ return (this != (G4VReadOutGeometry *) &right); }

void G4VReadOutGeometry::BuildROGeometry()
{
  ROworld = Build();
  ROnavigator->SetWorldVolume(ROworld);
}

G4bool G4VReadOutGeometry::CheckROVolume(G4Step*currentStep,G4TouchableHistory*& ROhist)
{
  ROhist = nullptr;
  G4bool incFlg = true;
  G4VPhysicalVolume* PV = currentStep->GetPreStepPoint()->GetPhysicalVolume();
  if((fexcludeList)&&(fexcludeList->CheckPV(PV)))
    { incFlg = false; }
  else if ((fincludeList)&&(fincludeList->CheckPV(PV)))
    { incFlg = true; }
  else if((fexcludeList)&&(fexcludeList->CheckLV(PV->GetLogicalVolume())))
    { incFlg = false; }
  else if((fincludeList)&&(fincludeList->CheckLV(PV->GetLogicalVolume())))
    { incFlg = true; }
  if(!incFlg) return false;
  
  if(ROworld)
    { incFlg = FindROTouchable(currentStep); }
  if(incFlg)
    { ROhist = touchableHistory; }
  return incFlg;
}

G4bool G4VReadOutGeometry::FindROTouchable(G4Step*currentStep)
{
  // Update G4TouchableHistory object (touchableHistory)
  // using the parallel readout world (ROworld)
  // Return false in case the current Step is outside of the
  // sensitive volume of the readout world.

  // At first invokation, creates the touchable history. Note
  // that default value (false) of Locate method is used.
  //  ---------> But the default Value is TRUE <-------------------- J.A. 
  if(!touchableHistory)
    {
      touchableHistory = new G4TouchableHistory();
      ROnavigator->LocateGlobalPointAndUpdateTouchable(
		         currentStep->GetPreStepPoint()->GetPosition(),
		         currentStep->GetPreStepPoint()->GetMomentumDirection(),
		         touchableHistory);
    }
  else
    {
      ROnavigator->LocateGlobalPointAndUpdateTouchable(
		         currentStep->GetPreStepPoint()->GetPosition(),
		         currentStep->GetPreStepPoint()->GetMomentumDirection(),
		         touchableHistory,
		         true);
    }
  // Can the above be improved by the use of an isotropic safety
  // in order to avoid LocateGlobalPointAndUpdateTouchable
  // at each Step ?
  // Should require that an RO geometry is notified at the
  // starting of a track to avoid possible confusion looking
  // at the safety value only.
  
  // checks if volume is sensitive:
  G4VPhysicalVolume* currentVolume = touchableHistory->GetVolume();
  // checks first if a physical volume exists here:
  if ( currentVolume )
    {
      return currentVolume->GetLogicalVolume()->
	GetSensitiveDetector() != 0;
    }
  // no sensitive volume found: returns false
  return false;
}

