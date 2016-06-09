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
// $Id: G4VIntraNuclearTransportModel.cc,v 1.5 2007/01/11 05:28:56 dennis Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// $Id: G4VIntraNuclearTransportModel.cc,v 1.0 1998/06/30
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, A. Feliciello, 30th June 1998
//               Removed delete of DeExcitation model, deleted elsewhere.
//                  F.W.Jones, 06-JUL-99
// -----------------------------------------------------------------------------

#include "G4VIntraNuclearTransportModel.hh"


// Class G4VIntraNuclearTransportModel 


G4VIntraNuclearTransportModel::G4VIntraNuclearTransportModel(const G4String& modelName) :
  G4HadronicInteraction(modelName),
  theTransportModelName(modelName), the3DNucleus(NULL), theDeExcitation(NULL)
{}


G4VIntraNuclearTransportModel::
G4VIntraNuclearTransportModel(const G4VIntraNuclearTransportModel& right) : 
  G4HadronicInteraction(right.GetModelName() )
{
  theTransportModelName = right.GetModelName();
  the3DNucleus = right.Get3DNucleus();
  theDeExcitation = right.GetDeExcitation();
}


G4VIntraNuclearTransportModel::~G4VIntraNuclearTransportModel()
{
//  if(the3DNucleus!=NULL) delete the3DNucleus;
  // This is deleted by ~G4HadronicInteractionRegistry
  // if(theDeExcitation!=NULL) delete theDeExcitation;
}


const G4VIntraNuclearTransportModel& 
G4VIntraNuclearTransportModel::
operator=(const G4VIntraNuclearTransportModel& right)
{
 if (this != &right)
   {
     theTransportModelName = right.GetModelName();
     the3DNucleus = right.Get3DNucleus();
     theDeExcitation = right.GetDeExcitation();
   }
 return *this;
}


int G4VIntraNuclearTransportModel::
operator==(const G4VIntraNuclearTransportModel& right) const
{
 return (this == (G4VIntraNuclearTransportModel *) & right);
}

int G4VIntraNuclearTransportModel::
operator!=(const G4VIntraNuclearTransportModel& right) const
{
 return (this != (G4VIntraNuclearTransportModel *) & right);
}
