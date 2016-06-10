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
// $Id: G4TrajectoryDrawByAttribute.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// Jane Tinslay September 2006
//
// Draw trajectories according to attribute.
//
#ifndef G4TRAJECTORYDRAWBYATTRIBUTE_HH
#define G4TRAJECTORYDRAWBYATTRIBUTE_HH

#include "globals.hh"
#include "G4VTrajectoryModel.hh"
#include <map>

class G4VAttValueFilter;
class G4VisTrajContext;

class G4TrajectoryDrawByAttribute : public G4VTrajectoryModel {

public:

  // Constructor
  G4TrajectoryDrawByAttribute(const G4String& name = "Unspecified", G4VisTrajContext* context=0);

  // Destructor
  virtual ~G4TrajectoryDrawByAttribute();

  // Draw the trajectory
  virtual void Draw(const G4VTrajectory& trajectory, 
		    const G4bool& visible = true) const;

  // Print configuration
  virtual void Print(std::ostream& ostr) const;

  // Configuration methods
  void Set(const G4String& attribute);
  void AddIntervalContext(const G4String& name, G4VisTrajContext* context);
  void AddValueContext(const G4String& name, G4VisTrajContext* context);

private:

  enum Config {Interval, SingleValue};

  typedef std::pair<G4String, Config> Pair;
  typedef std::map<Pair, G4VisTrajContext*> ContextMap;

  // Data members
  G4String fAttName;
  ContextMap fContextMap;

  // Caching
  mutable G4bool fFirst;
  mutable G4bool fWarnedMissingAttribute;
  mutable G4VAttValueFilter* filter;

};

#endif
