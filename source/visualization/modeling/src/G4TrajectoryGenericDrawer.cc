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
// $Id: G4TrajectoryGenericDrawer.cc,v 1.1 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
G4TrajectoryGenericDrawer::Draw(const G4VTrajectory& traj, const G4int& i_mode, const G4bool& visible) const
{
  G4VisTrajContext myContext(GetContext());
  myContext.SetVisible(visible);

  if (GetVerbose()) {
    G4cout<<"G4TrajectoryGenericDrawer named "<<Name();
    G4cout<<", drawing trajectory with configuration: "<<G4endl;
    myContext.Print(G4cout);
  }
  
  G4TrajectoryDrawerUtils::DrawLineAndPoints(traj, myContext, i_mode);
}

void
G4TrajectoryGenericDrawer::Print(std::ostream& ostr) const
{
  ostr<<"G4TrajectoryGenericDrawer model "<< Name()<< ", default configuration :"<<G4endl;
  GetContext().Print(G4cout);
}
