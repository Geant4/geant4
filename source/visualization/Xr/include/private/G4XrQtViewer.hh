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

#ifndef G4XRKQTVIEWER_HH
#define G4XRKQTVIEWER_HH


#include "G4VViewer.hh"

class G4UIManager;
class G4UIQt;

class QString;
class QWidget;
class QLineEdit;
class QSlider;
class QTreeWidget;
class QTreeWidgetItem;

#include "G4XrViewer.hh"

#include <qcoreapplication.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qopenglcontext.h>
#include <qpushbutton.h>
#include <qthread.h>
#include <qtreewidget.h>

#include <map>
#include <vector>

class G4XrQtViewer :  public G4XrViewer
{
  public:
    using PVNodeID = G4PhysicalVolumeModel::G4PhysicalVolumeNodeID;
    using PVPath = std::vector<PVNodeID>;

  public:
    G4XrQtViewer(G4VSceneHandler&, const G4String& name);
    ~G4XrQtViewer() override;
    void Initialise() override;
#ifdef G4MULTITHREADED
    // For switching threads in MT mode
    // Note: the order of calling of MovingToVisSubThread and SwitchToVisSubThread
    // is undefined, so we have to use mutexes to ensure required information,
    // namely the vis sub-thread address, is available before moving objects.
    // To summarise, the order of calling is
    //   DoneWithMasterThread
    //   MovingToVisSubThread ) or ( SwitchToVisSubThread
    //   SwitchToVisSubThread )    ( MovingToVisSubThread
    //   DoneWithVisSubThread
    //   MovingToMasterThread
    //   SwitchToMasterThread
    // Called on the master thread before starting the vis sub-thread.
    void DoneWithMasterThread() override;
    // Called on the master thread after starting the vis sub-thread.
    void MovingToVisSubThread() override;
    // Called on the vis sub-thread when waiting for events.
    void SwitchToVisSubThread() override;
    // Called on the vis sub-thread when all events have been processed.
    void DoneWithVisSubThread() override;
    // Called on the vis sub-thread when all events have been processed.
    // virtual void MovingToMasterThread ();  Not used in G4VtkQtViewer.
    // Called on the master thread after the vis sub-thread has terminated.
    void SwitchToMasterThread() override;
#endif

    void FinishView() override;
    void createSceneTreeWidget();
    void createSceneTreeComponent();
    QTreeWidgetItem* createTreeWidgetItem(const PVPath& fullPath, const QString& name, int copyNb,
                                          int POIndex, const QString& logicalName,
                                          Qt::CheckState state, QTreeWidgetItem* parentTreeNode,
                                          const G4Colour& color);
    void addNonPVSceneTreeElement(const G4String& model, G4Visible& visible, int currentPOIndex);
    void addPVSceneTreeElement(const G4String& model, G4PhysicalVolumeModel* pPVModel,
                               int currentPOIndex);

    QString getModelShortName(const G4String& model);
    bool parseAndInsertInSceneTree(QTreeWidgetItem* parentItem, G4PhysicalVolumeModel* pPVModel,
                                   unsigned int fullPathIndex, const QString& parentRoot,
                                   unsigned int currentIndexInTreeSceneHandler,
                                   int currentPVPOIndex);

  private:
    //G4UIQt* fUiQt;
    //QWidget* fGLWidget;

#ifdef G4MULTITHREADED
    QThread* fQVtkContextVisSubThread;
    QThread* fQVtkContextMainThread;
#endif
};

#endif