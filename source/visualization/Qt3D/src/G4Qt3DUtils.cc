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
// John Allison, 18th July 2020

#include "G4Qt3DUtils.hh"

#include "G4Qt3DQEntity.hh"
#include "G4PhysicalVolumeModel.hh"

Qt3DCore::QTransform* G4Qt3DUtils::CreateQTransformFrom(const G4Transform3D& g)
{
  auto* q = new Qt3DCore::QTransform;
  q->setMatrix
  (QMatrix4x4
  (g.xx(),g.xy(),g.xz(),g.dx(),
   g.yx(),g.yy(),g.yz(),g.dy(),
   g.zx(),g.zy(),g.zz(),g.dz(),
   0,0,0,1));
  q->setObjectName("transform");
  return q;
}

QColor G4Qt3DUtils::ConvertToQColor(const G4Colour& c) {
  QColor qColor;
  qColor.setRgbF(c.GetRed(),c.GetGreen(),c.GetBlue(),c.GetAlpha());
  return qColor;
}

QVector3D G4Qt3DUtils::ConvertToQVector3D(const G4ThreeVector& v) {
  return QVector3D(v.x(),v.y(),v.z());
}

// https://stackoverflow.com/questions/45759274/how-can-i-delete-all-nodes-recursively-in-the-root-entity-of-qt3dwindow
void G4Qt3DUtils::delete_entity_recursively(Qt3DCore::QNode *node){
#ifdef G4QT3DDEBUG
  G4Qt3DUtils::LogFile << "node " << node->objectName().toStdString() << std::endl;
#endif
  Qt3DCore::QEntity* entity = dynamic_cast<Qt3DCore::QEntity*>(node);
  if(entity == nullptr){
#ifdef G4QT3DDEBUG
    G4String name = node->objectName().toStdString();
    if (name == "") name = "X";
    G4Qt3DUtils::LogFile << (void*)node << ": "
    << "Deleting non-entity node " << name << std::endl;
#endif
    delete node;
    node = nullptr;
    return;
  }
  for (auto component: entity->components()) {
#ifdef G4QT3DDEBUG
    G4String name = component->objectName().toStdString();
    if (name == "") name = "X";
    G4Qt3DUtils::LogFile << (void*)node << ": " << "Deleting component " << name
    << " of " << entity->objectName().toStdString() << std::endl;
#endif
    entity->removeComponent(component);
    delete component;
    component = nullptr;
  }
  for (auto child_node: entity->childNodes()) {
    G4String name = child_node->objectName().toStdString();
    if (name == "") name = "X";
#ifdef G4QT3DDEBUG
    G4Qt3DUtils::LogFile << (void*)child_node << ": " << "Child node " << name
    << " of " << entity->objectName().toStdString() << std::endl;
#endif
    delete_entity_recursively(child_node);
  }
  G4String name = entity->objectName().toStdString();
  if (name == "") name = "X";
#ifdef G4QT3DDEBUG
  G4Qt3DUtils::LogFile << (void*)entity << ": " << "Deleting entity " << name << std::endl;
#endif
  delete entity;
  entity = nullptr;
}

void G4Qt3DUtils::delete_components_and_children_of_entity_recursively(Qt3DCore::QNode *node){
  Qt3DCore::QEntity* entity = dynamic_cast<Qt3DCore::QEntity*>(node);
  if(entity == nullptr){
#ifdef G4QT3DDEBUG
    G4String name = node->objectName().toStdString();
    if (name == "") name = "X";
    G4Qt3DUtils::LogFile << (void*)node << ": " << "Found non-entity node " << name << std::endl;
#endif
    return;
  }
  for (auto component: entity->components()){
#ifdef G4QT3DDEBUG
    G4String name = component->objectName().toStdString();
    if (name == "") name = "X";
    G4Qt3DUtils::LogFile << (void*)entity << ": " << "Deleting component " << name
    << " of " << entity->objectName().toStdString() << std::endl;
#endif
    entity->removeComponent(component);
    delete(component);
    component = nullptr;
  }
  auto child_nodes = entity->childNodes();
  for (auto child_node: child_nodes) {
    G4String name = child_node->objectName().toStdString();
    if (name == "") name = "X";
#ifdef G4QT3DDEBUG
    G4Qt3DUtils::LogFile << (void*)child_node << ": " << "Child node " << name
    << " of " << entity->objectName().toStdString() << std::endl;
#endif
    delete_entity_recursively(child_node);
  }
  G4String name = entity->objectName().toStdString();
  if (name == "") name = "X";
#ifdef G4QT3DDEBUG
  G4Qt3DUtils::LogFile << (void*)entity << ": " << "Clearing child nodes of " << name << std::endl;
#endif
  child_nodes.clear();
}

#ifdef G4QT3DDEBUG
std::ofstream G4Qt3DUtils::LogFile("LogFile.txt");
void G4Qt3DUtils::PrintQObjectTree
 (const QObject* node,
  const G4String& where)
{
  auto& logFile = G4Qt3DUtils::LogFile;
  if (where.length()) logFile << "\n===== QObjectTree at " << where << std::endl;
  static G4int iDep = -1;
  ++iDep;
  G4String nodeName = node->objectName().toStdString();
  if (nodeName == "") nodeName = "X";
  for (G4int i = 0; i < iDep; ++i) logFile << "  ";
  logFile << (void*)node << ": "
  << "Node at depth " << iDep << ": " << nodeName << ": "
  << "thread: " << node->thread() << ": "
  << "parent: " << node->parent() << ": ";
  const auto* g4node = dynamic_cast<const G4Qt3DQEntity*>(node);
  if (g4node) {
    logFile << g4node->GetPVNodeID() << std::endl;
  } else {
    logFile << typeid(node).name() << std::endl;
  }
  if (g4node) {
    for (const auto& component: g4node->components()) {
      G4String name = component->objectName().toStdString();
      if (name == "") name = "X";
      for (G4int i = 0; i < iDep; ++i) logFile << "  ";
      logFile << (void*)component << ": "<< "Component at depth " << iDep << " "
      << name << " of " << nodeName << std::endl;
    }
  }
  for (const auto& child: node->children()) {
    PrintQObjectTree(child);
  }
  --iDep;
  if (where.length()) logFile << "===== End: QObjectTree at " << where << std::endl;
  return;
}
#endif
