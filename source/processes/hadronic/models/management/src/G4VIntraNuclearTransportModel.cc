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
// $Id: G4VIntraNuclearTransportModel.cc 69717 2013-05-13 09:47:57Z gcosmo $
//
// $Id: G4VIntraNuclearTransportModel.cc,v 1.0 1998/06/30
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, A. Feliciello, 30th June 1998
//               Removed delete of DeExcitation model, deleted elsewhere.
//                  F.W.Jones, 06-JUL-99
//      V.Ivanchenko 03.01.2012
//          Added G4VPreCompoundModel pointer to the constructor and cleanup
//      V. Uzhinsky Nov. 2012
//          Added method PropagateNuclNucl for simulation of nucleus-nucleus inter. 
// -----------------------------------------------------------------------------

#include "G4VIntraNuclearTransportModel.hh"

G4VIntraNuclearTransportModel::G4VIntraNuclearTransportModel(
        const G4String& modName, G4VPreCompoundModel* ptr)
  : G4HadronicInteraction(modName),theTransportModelName(modName),
    the3DNucleus(0),theDeExcitation(ptr),thePrimaryProjectile(0)
{}

G4VIntraNuclearTransportModel::~G4VIntraNuclearTransportModel()
{
  //  if(the3DNucleus!=NULL) delete the3DNucleus;
  // This is deleted by ~G4HadronicInteractionRegistry
  // if(theDeExcitation!=NULL) delete theDeExcitation;
}

void G4VIntraNuclearTransportModel::ModelDescription(std::ostream& outFile) const
{
	outFile << "G4VIntraNuclearTransportModel is abstract class" << G4endl;
	G4Exception("G4VIntraNuclearTransportModel::ModelDescription()","G4VINT01",FatalException,
			"G4VIntraNuclearTransportModel is abstract class, no description available");
}

void G4VIntraNuclearTransportModel::PropagateModelDescription(std::ostream& outFile) const
{
	outFile << "G4VIntraNuclearTransportModel is abstract class, missing description" << G4endl;
//	G4Exception("G4VIntraNuclearTransportModel::ModelDescription()","G4VINT01",FatalException,
//			"G4VIntraNuclearTransportModel is abstract class, no description available");
}

G4ReactionProductVector* G4VIntraNuclearTransportModel::PropagateNuclNucl(G4KineticTrackVector* ,
               G4V3DNucleus* , G4V3DNucleus* )
{
   G4Exception("G4VIntraNuclearTransportModel::Propagate()","G4VINT02",FatalException,
         "Propagate method for nucleus-nucleus interactions not implemented");
   return 0;
}
