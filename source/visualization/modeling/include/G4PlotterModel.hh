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
// Guy Barrand 25th September 2021
//

#ifndef G4PLOTTERMODEL_HH
#define G4PLOTTERMODEL_HH

#include "G4VModel.hh"

#include "G4Plotter.hh"

class G4PlotterModel: public G4VModel {
public:
  G4PlotterModel(G4Plotter& plotter,const G4String& global_description = "",const G4Transform3D& transform = G4Transform3D());
  virtual ~G4PlotterModel () = default;
  void DescribeYourselfTo (G4VGraphicsScene&) override;
  G4Plotter& plotter() {return fPlotter;}

private:
  G4PlotterModel (const G4PlotterModel&);
  G4PlotterModel& operator = (const G4PlotterModel&);

  G4Plotter& fPlotter;
  G4Transform3D fTransform;
};

#endif
