#ifndef MLSteppingAction_h
#define MLSteppingAction_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4UserSteppingAction.hh"

class MLGeometryConstruction;
////////////////////////////////////////////////////////////////////////////////
//
class MLSteppingAction : public G4UserSteppingAction
{
  public:
    MLSteppingAction (MLGeometryConstruction*);
    virtual ~MLSteppingAction ();

    virtual void UserSteppingAction (const G4Step*);
  
  private:

    MLGeometryConstruction *geometry;
};
////////////////////////////////////////////////////////////////////////////////
#endif
