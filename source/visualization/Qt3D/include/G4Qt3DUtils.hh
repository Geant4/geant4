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

#if defined (G4VIS_BUILD_QT3D_DRIVER) || defined (G4VIS_USE_QT3D)

#ifndef G4QT3DUTILS_HH
#define G4QT3DUTILS_HH

#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"

#include <Qt3DCore>
#include <fstream>

class G4Colour;

namespace G4Qt3DUtils {

Qt3DCore::QTransform* CreateQTransformFrom(const G4Transform3D&);
QColor ConvertToQColor(const G4Colour& c);
QVector3D ConvertToQVector3D(const G4ThreeVector& v);

void delete_entity_recursively(Qt3DCore::QNode *node);
void delete_components_and_children_of_entity_recursively(Qt3DCore::QNode *node);

#ifdef G4QT3DDEBUG
extern std::ofstream LogFile;
void PrintQObjectTree
(const QObject* node,
 const G4String& where = "");
#endif

}  // namespace G4Qt3DUtils

#endif

#endif  // #if defined (G4VIS_BUILD_QT3D_DRIVER) || defined (G4VIS_USE_QT3D)
