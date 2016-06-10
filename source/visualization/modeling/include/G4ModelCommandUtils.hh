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
// $Id: G4ModelCommandUtils.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay September 2006
//
// Command utilities
//
#ifndef G4MODELCOMMANDUTILS_HH
#define G4MODELCOMMANDUTILS_HH

#include "G4String.hh"
#include "G4ModelCommandsT.hh"
#include "G4UImessenger.hh"
#include "G4VisTrajContext.hh"

namespace G4ModelCommandUtils {

  void AddContextMsgrs(G4VisTrajContext* context, std::vector<G4UImessenger*>& messengers,
		       const G4String& placement)
  {
    messengers.push_back(new G4ModelCmdCreateContextDir<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetDrawLine<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetLineVisible<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetLineColour<G4VisTrajContext>(context, placement));
    
    messengers.push_back(new G4ModelCmdSetDrawStepPts<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsVisible<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsColour<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsSize<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsSizeType<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsType<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetStepPtsFillStyle<G4VisTrajContext>(context, placement));
    
    messengers.push_back(new G4ModelCmdSetDrawAuxPts<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsVisible<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsColour<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsSize<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsSizeType<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsType<G4VisTrajContext>(context, placement));
    messengers.push_back(new G4ModelCmdSetAuxPtsFillStyle<G4VisTrajContext>(context, placement));

    messengers.push_back(new G4ModelCmdSetTimeSliceInterval<G4VisTrajContext>(context, placement));
  }  
}

#endif
