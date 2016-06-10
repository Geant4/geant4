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
// $Id: G4TrajectoryGenericDrawer.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay May 2006
//
#include "G4TrajectoryGenericDrawer.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VisTrajContext.hh"
#include "G4VTrajectory.hh"
#include <sstream>

G4TrajectoryGenericDrawer::G4TrajectoryGenericDrawer(const G4String& name, G4VisTrajContext* context)
  :G4VTrajectoryModel(name, context)
{}

G4TrajectoryGenericDrawer::~G4TrajectoryGenericDrawer() {}

void
G4TrajectoryGenericDrawer::Draw(const G4VTrajectory& traj, const G4bool& visible) const
{
  G4VisTrajContext myContext(GetContext());
  myContext.SetVisible(visible);

  if (GetVerbose()) {
    G4cout<<"G4TrajectoryGenericDrawer named "<<Name();
    G4cout<<", drawing trajectory with configuration: "<<G4endl;
    myContext.Print(G4cout);
  }
  
  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, myContext);
}

void
G4TrajectoryGenericDrawer::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryGenericDrawer model "<< Name()<< ", default configuration :"<<G4endl;
  GetContext().Print(G4cout);
}
