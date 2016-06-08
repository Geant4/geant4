// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VIntraNuclearTransportModel.cc,v 1.2.8.1 1999/12/07 20:51:46 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// $Id: G4VIntraNuclearTransportModel.cc,v 1.0 1998/06/30
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD Group
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
G4VIntraNuclearTransportModel(const G4VIntraNuclearTransportModel& right)
{
 the3DNucleus = right.Get3DNucleus();
 theDeExcitation = right.GetDeExcitation();
}


G4VIntraNuclearTransportModel::~G4VIntraNuclearTransportModel()
{
  if(the3DNucleus!=NULL) delete the3DNucleus;
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
