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
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalEcalOrganization.hh
// Description: Defines numbering schema for the Electromagnetic Calorimeter
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalEcalOrganization_h
#define CCalEcalOrganization_h

#include "CCalVOrganization.hh"

class CCalEcalOrganization: public CCalVOrganization {

public:
  CCalEcalOrganization(){};
  ~CCalEcalOrganization();
         
  virtual unsigned int GetUnitID(const G4Step* aStep) const ;
      
};

#endif
