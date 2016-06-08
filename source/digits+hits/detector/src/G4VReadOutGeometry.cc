// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VReadOutGeometry.cc,v 1.2.6.1 1999/12/07 20:47:45 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#include "G4VReadOutGeometry.hh"
#include "G4Navigator.hh"


G4VReadOutGeometry::G4VReadOutGeometry()
  :ROworld(NULL),touchableHistory(NULL),
   fincludeList(NULL),fexcludeList(NULL)
{
  name = "unknown";
  ROnavigator = new G4Navigator();
}

G4VReadOutGeometry::G4VReadOutGeometry(const G4VReadOutGeometry &right)
{
  fincludeList = right.fincludeList;
  fexcludeList = right.fexcludeList;
  name = right.name;
  ROworld = right.ROworld;
  touchableHistory = NULL;
  // COPY CONSTRUCTOR NOT STRAIGHT FORWARD: need to copy the touchabelHistory
  // VALUE, same foe navigator and same for the World+Geom hierachy ??
}

G4VReadOutGeometry::G4VReadOutGeometry(G4String n) 
  :name(n),ROworld(NULL),touchableHistory(NULL),
   fincludeList(NULL),fexcludeList(NULL)
{
  ROnavigator = new G4Navigator();
}

G4VReadOutGeometry::~G4VReadOutGeometry()
{ 
  //if(ROworld) delete ROworld; //should we do ? will it delete the goem tree also ?
  if(fincludeList)     delete fincludeList;
  if(fexcludeList)     delete fexcludeList;
  if(touchableHistory) delete touchableHistory;
  if(ROnavigator)      delete ROnavigator;
}

const G4VReadOutGeometry & G4VReadOutGeometry::operator=(const G4VReadOutGeometry &right)
{
  fincludeList     = right.fincludeList;
  fexcludeList     = right.fexcludeList;
  name             = right.name;
  ROworld          = right.ROworld;
  touchableHistory = NULL;
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
  ROhist = NULL;
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
  G4VPhysicalVolume* currentVolume;
  // checks first if a physical volume exists here:
  if ( currentVolume = touchableHistory->GetVolume() )
    {
      return currentVolume->GetLogicalVolume()->
	GetSensitiveDetector() != NULL;
    }
  // no sensitive volume found: returns false
  return false;
}

