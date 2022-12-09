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
// Jane Tinslay August 2006
//
#include "G4TrajectoryDrawByAttribute.hh"
#include "G4AttDef.hh"
#include "G4AttFilterUtils.hh"
#include "G4AttUtils.hh"
#include "G4AttValue.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VAttValueFilter.hh"
#include "G4VisTrajContext.hh"
#include "G4VTrajectory.hh"
#include <assert.h>

G4TrajectoryDrawByAttribute::G4TrajectoryDrawByAttribute(const G4String& name, G4VisTrajContext* context)
  :G4VTrajectoryModel(name, context)
  ,fAttName("")
  ,fFirst(true)
  ,fWarnedMissingAttribute(false)
  ,filter(0)
{}

G4TrajectoryDrawByAttribute::~G4TrajectoryDrawByAttribute() 
{
  ContextMap::iterator iter = fContextMap.begin();
  
  while (iter != fContextMap.end()) {
    delete iter->second;
    iter++;
  }
  
  delete filter;
}

void
G4TrajectoryDrawByAttribute::Draw(const G4VTrajectory& object, 
				  const G4bool& /*visible*/) const
{
  // Return if attribute name has not been set. Just print one warning
  if (fAttName.empty()) {

    if (!fWarnedMissingAttribute) {
      G4ExceptionDescription ed;
      ed<<"Null attribute name";
      G4Exception("G4TrajectoryDrawByAttribute::Draw",
		  "modeling0116",
		  JustWarning, ed);
      fWarnedMissingAttribute = true;
    }
    
    return;
  }
  
  // Basically cache data loaded filter for efficiency
  if (fFirst) {
    
    fFirst = false;
    
    // Get attribute definition
    G4AttDef attDef;
    
    // Expect definition to exist    
    if (!G4AttUtils::ExtractAttDef(object, fAttName, attDef)) {
      static G4bool warnedUnableToExtract = false;
      if (!warnedUnableToExtract) {
	G4ExceptionDescription ed;
	ed <<"Unable to extract attribute definition named "<<fAttName << '\n'
        << "Available attributes:\n"
        << *object.GetAttDefs();
	G4Exception
	  ("G4TrajectoryDrawByAttribute::Draw",
	   "modeling0117", JustWarning, ed, ". Invalid attribute name");
	warnedUnableToExtract = true;
      }
      return;
    }
    
    // Get new G4AttValue filter
    filter = G4AttFilterUtils::GetNewFilter(attDef);
    assert (0 != filter);
    
    // Load both interval and single valued data. Single valued data should
    // override interval data.
    ContextMap::const_iterator iter = fContextMap.begin();
    
    while (iter != fContextMap.end()) {
      if (iter->first.second == G4TrajectoryDrawByAttribute::Interval) {
	filter->LoadIntervalElement(iter->first.first);
      }
      else if (iter->first.second == G4TrajectoryDrawByAttribute::SingleValue) {
	filter->LoadSingleValueElement(iter->first.first);
      }
      iter++;
    }
  }
 
  // Get attribute value
  G4AttValue attVal;

  // Expect value to exist
  if (!G4AttUtils::ExtractAttValue(object, fAttName, attVal)) {
    static G4bool warnedUnableToExtract = false;
    if (!warnedUnableToExtract) {
      G4ExceptionDescription ed;
      ed <<"Unable to extract attribute definition named "<<fAttName << '\n'
      << "Available attributes:\n"
      << *object.GetAttDefs();
      G4Exception
	("G4TrajectoryDrawByAttribute::Draw",
	 "modeling0118", JustWarning, ed, ". Invalid attribute name");
      warnedUnableToExtract = true;
    }
      return;
  }
  
  G4VisTrajContext myContext(GetContext());
  G4String key;

  // If attribute value passes filter, get corresponding interval/single value 
  // key loaded into G4AttValue filter.
  if (filter->GetValidElement(attVal, key)) {

    // Extract context corresponding to valid key.
    // Single value match should have overriden interval match.
    ContextMap::const_iterator iter = fContextMap.begin();

    G4bool gotContext(false);

    while (!gotContext && (iter != fContextMap.end())) {
      if (iter->first.first == key) {
	myContext = *(iter->second);
	gotContext = true;
      }
      iter++;
    }

    assert (gotContext);
  }
  
  if (GetVerbose()) {
    G4cout<<"G4TrajectoryDrawByAttribute drawer named "<<Name();
    G4cout<<", drawing style selected according to value of attribute "<<fAttName;
    G4cout<<" : "<<attVal.GetValue()<<".  Selected context:"<<G4endl;
    myContext.Print(G4cout);
  }

  // Draw the trajectory
  G4TrajectoryDrawerUtils::DrawLineAndPoints(object, myContext);
}

void
G4TrajectoryDrawByAttribute::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryDrawByAttribute, dumping configuration for model named "<< Name() <<":"<<std::endl;;

  ostr<<"Default configuration:"<<G4endl;
  GetContext().Print(ostr);

  ostr<<"\nAttribute name "<<fAttName<<std::endl;
  ostr<<"\nKey<->Context map dump:"<<std::endl;

  ContextMap::const_iterator iter = fContextMap.begin();

  while (iter != fContextMap.end()) {
    ostr<<"Context for key "<<iter->first.first<<":"<<std::endl;
    iter->second->Print(ostr);
    
    iter++;
  }
}

void
G4TrajectoryDrawByAttribute::Set(const G4String& name)
{
  fAttName = name;
}

void
G4TrajectoryDrawByAttribute::AddIntervalContext(const G4String& name, G4VisTrajContext* context)
{
  // Takes ownership of context
  std::pair<G4String, Config> myPair(name, G4TrajectoryDrawByAttribute::Interval);

  ContextMap::iterator iter = fContextMap.find(myPair);
  
  if (iter != fContextMap.end()) {
    G4ExceptionDescription ed;
    ed <<"Interval "<< name <<" already exists";
    G4Exception
      ("G4TrajectoryDrawByAttribute::AddIntervalContext",
       "modeling0119", FatalErrorInArgument, ed, ". Invalid interval");
  }

  fContextMap[myPair] = context;
}

void
G4TrajectoryDrawByAttribute::AddValueContext(const G4String& name, G4VisTrajContext* context)
{
  // Takes ownership of context
  std::pair<G4String, Config> myPair(name, G4TrajectoryDrawByAttribute::SingleValue);

  ContextMap::iterator iter = fContextMap.find(myPair);
  
  if (iter != fContextMap.end()) {
    G4ExceptionDescription ed;
    ed <<"Single value "<< name <<" already exists";
    G4Exception
      ("G4TrajectoryDrawByAttribute::AddSingleValueContext",
       "modeling0120", FatalErrorInArgument, ed, ". Invalid value");
  }

  fContextMap[myPair] = context;
}
