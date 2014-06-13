// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFRunAction.hh
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Definition of the FFRunAction class
//!
//  ================ End Documentation Comments ================
//
//  Modified: 
//
// -------------------------------------------------------------

#ifndef FFRUNACTION
#define FFRUNACTION

#include "G4UserRunAction.hh"

class FFRunAction
:   public G4UserRunAction
{
public:
// Constructor
    FFRunAction();
    
// Function
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    
// Destructor
    virtual ~FFRunAction();
};

#endif //FFRUNACTION


