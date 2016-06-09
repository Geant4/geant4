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
// $Id: OlapSteppingAction.cc,v 1.3 2006-06-29 17:23:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// OlapSteppingAction
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "OlapSteppingAction.hh"
#include "OlapEventAction.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4NavigationHistory.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "SolidAnalyser.hh"

std::ostream&
operator << (std::ostream& os, std::vector<OlapStepInfo*>& vec) 
{
   std::vector<OlapStepInfo*>::iterator it = vec.begin();
   while (it!=vec.end())
   {
      G4VPhysicalVolume * pv = (*it)->theHist.GetTopVolume();
      os << (*it)->thePoint << " [" << pv->GetName() << ":" << pv->GetCopyNo() 
         << "]" << G4endl;
      it++;
   }
   return os;
}

std::ostream&
operator << (std::ostream& os, const OlapStepInfo& aStepInfo)
{
  //G4cout << "History depth="<<aStepInfo.theHist.GetDepth()<< G4endl;
  for (G4int i=0;i<=aStepInfo.theHist.GetDepth();i++)
  {
      G4VSolid * solid =
        aStepInfo.theHist.GetVolume(i)->GetLogicalVolume()->GetSolid();
      SolidAnalyser * solAna = SolidAnalyser::GetSolidAnalyser();
      std::vector< std::pair<G4String,G4double> > param;
      solAna->GetParam(solid,param);

      os << "["<<i<<"]: " ;
      const G4AffineTransform & transP = aStepInfo.theHist.GetTransform(i);
      G4ThreeVector locP(transP.TransformPoint(aStepInfo.thePoint));
      char c = '?';
      if (solid->Inside(locP)==kInside)
        c='i';
      if (solid->Inside(locP)==kOutside)
        c='o';
      if (solid->Inside(locP)==kSurface)
        c='s';
      os << " ins=[" << c << "] " ; // << "lok.pt.=" << locP;
      if( aStepInfo.theHist.GetVolume(i) != 0 ) {
	 os   << " PVName=["<< aStepInfo.theHist.GetVolume(i)->GetName() << ":"
	      << aStepInfo.theHist.GetVolume(i)->GetCopyNo() 
	      << "] Type=[";
	 switch(aStepInfo.theHist.GetVolumeType(i))
	   {
	   case kNormal:
	     os <<"N";
	     break;
	   case kReplica:
	     os <<"R" << aStepInfo.theHist.GetReplicaNo(i);
	     break;
	   case kParameterised:
	     os <<"P" << aStepInfo.theHist.GetReplicaNo(i);
	     break;
	   }
	 os << "] ";
      }else{
	 os << "Phys = <Null>";
      }
      
      const G4AffineTransform & trans =
         aStepInfo.theHist.GetTransform(i).Inverse();
      
      os << "P=" << trans.NetTranslation() 
         << " R=[" << trans.NetRotation().phiX()/degree << "," 
	            << trans.NetRotation().thetaX()/degree << ","  
	            << trans.NetRotation().phiY()/degree << "," 
		    << trans.NetRotation().thetaY()/degree << ","
		    << trans.NetRotation().phiZ()/degree << "," 
		    << trans.NetRotation().thetaZ()/degree << "] "; 
      	    
      // solid specification
      
      os << " " << solid->GetEntityType() << ": "; 		  
      std::vector< std::pair<G4String,G4double> >::iterator it =
         param.begin();
      while ( it != param.end() )
      {
        os << (*it).first << '=' << (*it).second << " ";
        it++;
      }
      os << G4endl;
  }
  return os;
}


OlapSteppingAction::OlapSteppingAction(OlapEventAction * aEvAct)
  : theEventAction(aEvAct)
{
}


OlapSteppingAction::~OlapSteppingAction()
{
} 


void OlapSteppingAction::UserSteppingAction(const G4Step * aStep)
{
  const G4Track * aTrack = aStep->GetTrack();
  const G4NavigationHistory & aHist = *(aTrack->GetTouchable()->GetHistory()); 
  G4int aTrackID = aTrack->GetTrackID();
  G4int aStepNo = aTrack->GetCurrentStepNumber();
  //G4cout << "Stepping: " << aTrack->GetVolume()->GetName() 
  //       << " CpNr: " << aTrack->GetVolume()->GetCopyNo() << G4endl;
	 
   if (aStepNo>999) {
      G4cerr << "OlapSteppingAction(): to many steps, killing track" << G4endl
             << "      pv=[" << aTrack->GetVolume()->GetName() 
	     << "] cpnr=[" << aTrack->GetVolume()->GetCopyNo() 
	     << "] pos=" << aTrack->GetPosition()
             << G4endl;
      //G4EventManager::GetEventManager()->AbortCurrentEvent();
      //G4RunManager::GetRunManager()->AbortRun();
      G4Track * aNonConstTrack = const_cast<G4Track*>(aTrack);
      aNonConstTrack->SetTrackStatus(fStopAndKill);
      return;
  }
  OlapStepInfo * aStepInfo = new OlapStepInfo;
  
  aStepInfo->thePoint = aTrack->GetPosition();
  aStepInfo->theHist = aHist;
  
  if (aTrackID % 2)
  {
     theEventAction->ABSteps.push_back(aStepInfo);
  }
  else 
  {
     theEventAction->BASteps.push_back(aStepInfo);
  }   
  //G4cout << *aStepInfo << G4endl;
}
