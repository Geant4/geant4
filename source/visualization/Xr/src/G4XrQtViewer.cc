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

#include "G4XrQtViewer.hh"

#include "G4LogicalVolume.hh"
#include "G4Qt.hh"
#include "G4UIQt.hh"
#include "G4UImanager.hh"
#include "G4VSceneHandler.hh"

G4XrQtViewer::G4XrQtViewer(G4VSceneHandler& sceneHandler, const G4String& name)
  : G4XrViewer(sceneHandler, name)
{
  G4Qt::getInstance();
}

G4XrQtViewer::~G4XrQtViewer()
{
}

void G4XrQtViewer::Initialise()
{
}

#ifdef G4MULTITHREADED
void G4XrQtViewer::DoneWithMasterThread()
{
}

void G4XrQtViewer::MovingToVisSubThread()
{
}

void G4XrQtViewer::SwitchToVisSubThread()
{
}

void G4XrQtViewer::DoneWithVisSubThread()
{
}

void G4XrQtViewer::SwitchToMasterThread()
{
}

#endif

void G4XrQtViewer::FinishView()
{
}

void G4XrQtViewer::createSceneTreeWidget() {}

void G4XrQtViewer::createSceneTreeComponent() {}

QTreeWidgetItem*
G4XrQtViewer::createTreeWidgetItem(const PVPath& /*fullPath*/, const QString& /*name*/,
                                    int /*copyNb*/, int /*POIndex*/, const QString& /*logicalName*/,
                                    Qt::CheckState /*state*/, QTreeWidgetItem* /*parentTreeNode*/,
                                    const G4Colour& /*color*/)
{
  QTreeWidgetItem* newItem = nullptr;
  return newItem;
}

void G4XrQtViewer::addNonPVSceneTreeElement(const G4String& /*model*/, G4Visible& /*visible*/,
                                             int /*currentPOIndex*/)
{}

void G4XrQtViewer::addPVSceneTreeElement(const G4String& /*model*/,
                                          G4PhysicalVolumeModel* /*pPVModel*/,
                                          int /*currentPOIndex*/)
{}

QString G4XrQtViewer::getModelShortName(const G4String& /*model*/)
{
  QString modelShortName;
  return modelShortName;
}

bool G4XrQtViewer::parseAndInsertInSceneTree(QTreeWidgetItem* /*parentItem*/,
                                              G4PhysicalVolumeModel* /*pPVModel*/,
                                              unsigned int /*fullPathIndex*/,
                                              const QString& /*parentRoot*/,
                                              unsigned int /*currentIndexInTreeSceneHandler*/,
                                              int /*currentPVPOIndex*/)
{
  return false;
}