#ifndef MLStackingAction_h
#define MLStackingAction_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              
//
// Version:		0.A
// Date:		
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd
// Customer:		UK Ministry of Defence
// Contract:		CRP MAT Domain
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DISCLAIMER
// ----------
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
////////////////////////////////////////////////////////////////////////////////
//

#include "G4ThreeVector.hh"
#include "G4UserStackingAction.hh"

#include <vector>

class G4Track;
////////////////////////////////////////////////////////////////////////////////
//
class MLStackingAction : public G4UserStackingAction
{
public:
  MLStackingAction();
  virtual ~MLStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

};
////////////////////////////////////////////////////////////////////////////////
//
#endif
