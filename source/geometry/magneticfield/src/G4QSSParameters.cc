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
// G4QSSParameters
//
// Hold parameters for QSS Integrator driver
// 
// Author: John Apostolakis (CERN), 19.08.2025
// --------------------------------------------------------------------

#include "G4QSSParameters.hh"
#include "G4Exception.hh"

const G4double dQMin_smallest = 1.0e-20;   // Minimum limit for deviation
const G4double dQMin_biggest  = 1.0e+20;   // Maximum limit for deviation
// These can be problem dependent - wide values

const G4double dQRel_smallest = 1.0e-08;  // Minimum limit for relative deviation
const G4double dQRel_biggest  = 0.01;     // Maximum limit for relative deviation
// Achievable relative accuracy needs tight bounds 

G4QSSParameters* G4QSSParameters::Instance()
{
  static G4QSSParameters theQSSParametersSingleton;
  return &theQSSParametersSingleton;
}

// --------------------------------------------------------------------------------

G4bool  G4QSSParameters::Set_dQMin( G4double value )
{ 
  G4bool good_value= dQMin_smallest < value && value < dQMin_biggest;
  if( good_value )
  { 
    fdQMin = value; 
    //-------------
  } 
  else 
  { 
    G4ExceptionDescription errmsg;
    errmsg << " OUT of RANGE: Proposed value of dQMin is not within expected limits: " 
           << " Smallest = " << dQMin_smallest 
           << " Biggest  = " << dQMin_biggest
           << " Proposed value = " << value;
    G4Exception("G4QSSParameters::Set_dQMin", "QSS-Parameters-0002", 
                FatalException, errmsg );
  }
  return good_value;
}

// --------------------------------------------------------------------------------

G4bool  G4QSSParameters::Set_dQRel( G4double value )
{ 
  G4bool good_value= dQRel_smallest < value && value < dQRel_biggest ;
  if( good_value )
  { 
    fdQRel = value; 
    //-------------
  } 
  else 
  { 
    G4ExceptionDescription errmsg;
    errmsg << " OUT of RANGE: Proposed value of dQRel is not within expected limits: " 
           << " Smallest = " << dQRel_smallest 
           << " Biggest  = " << dQRel_biggest
           << " Proposed value = " << value;
    G4Exception("G4QSSParameters::Set_dQRel", "QSS-Parameters-0003", 
                FatalException, errmsg );
  }
  return good_value;
}

// --------------------------------------------------------------------------------

G4bool  G4QSSParameters::SetQssOrder( G4int order_val,  G4bool onlyWarn )
{
  G4bool good_value= false;
  if( 2 <= order_val && order_val <= 3 )
  {
    good_value= true;
    fQssOrder = order_val;
  }
  else
  {
    G4ExceptionDescription err_msg;     
    err_msg << "Requested order for G4QSStepper= "
            << order_val << "  is neither 2 nor 3.";
    G4Exception("G4QSSParameters::SetQssOrder", "QSS-Parameters-0001", 
                onlyWarn ? JustWarning : FatalException, err_msg );
  }
  return good_value;
}

// --------------------------------------------------------------------------------

G4bool G4QSSParameters::SetMaxSubsteps( G4int maxSubSteps )
{
  const G4int maxSubstep_smallest= 2;     // Min acceptable value
  const G4int maxSubstep_biggest = 10000; // Maxacceptable value
   
  G4bool good_value= false;
  if( maxSubstep_smallest <= maxSubSteps && maxSubSteps <= maxSubstep_biggest )
  {
    good_value= true;
    fMaxSubsteps = maxSubSteps;
  }
  else
  {
    G4ExceptionDescription err_msg;     
    err_msg << "Requested max-substeps for G4QSStepper= " << maxSubSteps
            << "  is not between "
            << maxSubstep_smallest << "  and " << maxSubstep_biggest;
    G4Exception("G4QSSParameters::SetMaxSubsteps", "QSS-Parameters-0001",
                JustWarning, err_msg );
  }
  return good_value;
}
