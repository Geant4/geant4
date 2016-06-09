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
//
// $Id: G4VIntraNuclearTransportModel.cc,v 1.3 2005/06/04 13:40:04 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
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



G4VIntraNuclearTransportModel::G4VIntraNuclearTransportModel() :
                              the3DNucleus(NULL),
                              theDeExcitation(NULL)
{
}

G4VIntraNuclearTransportModel::
G4VIntraNuclearTransportModel(const G4VIntraNuclearTransportModel& right) : G4HadronicInteraction()
{
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
