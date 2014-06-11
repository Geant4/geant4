// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFActionInitialization.hh
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Definition of the FFActionInitialization class
//!
//  ================ End Documentation Comments ================
//
//  Modified: 
//
// -------------------------------------------------------------

#ifndef FFACTIONINITIALIZATION
#define FFACTIONINITIALIZATION

#include "G4VUserActionInitialization.hh"

#include "FFRunAction.hh"


class FFActionInitialization
:   public G4VUserActionInitialization
{
public:
// Constructor
    FFActionInitialization();

// Functions
    virtual void BuildForMaster() const;
    virtual void Build() const;

// Destructor
    virtual ~FFActionInitialization();
    
private:
// Fields
    FFRunAction* const masterRunAction;
};

#endif // FFACTIONINITIALIZATION


