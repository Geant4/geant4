//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalVOrganization_h
#define CCalVOrganization_h

#include "G4Step.hh"
#include "CCaloOrganization.hh"

class CCalVOrganization {

public:
  CCalVOrganization(){};
  virtual ~CCalVOrganization(){};
	 
  virtual unsigned int GetUnitID(const G4Step* aStep) const = 0;
  virtual int  Levels(const G4Step*) const;
  virtual void DetectorLevel(const G4Step*, int&, int*, G4String*) const;
     
protected:
  CCaloOrganization theOrg;
};

#endif
