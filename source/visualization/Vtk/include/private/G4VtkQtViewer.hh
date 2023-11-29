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

#ifndef G4VTKQTVIEWER_HH
#define G4VTKQTVIEWER_HH

// #define G4VTKDEBUG  // Comment this out to suppress debug code.

#include "G4VViewer.hh"

class G4UIManager;
class G4UIQt;

class QString;
class QWidget;
class QLineEdit;
class QSlider;
class QTreeWidget;
class QTreeWidgetItem;

#include "G4VtkViewer.hh"

#include "QVTKOpenGLNativeWidget.h"

#include <map>
#include <vector>

class G4VtkQtViewer : public QVTKOpenGLNativeWidget, public G4VtkViewer
{
  public:
    using PVNodeID = G4PhysicalVolumeModel::G4PhysicalVolumeNodeID;
    using PVPath = std::vector<PVNodeID>;

  public:
    G4VtkQtViewer(G4VSceneHandler&, const G4String& name);
    ~G4VtkQtViewer() override;
    void Initialise() override;
    virtual void CreateMainWindow(QVTKOpenGLNativeWidget*, const QString&);
#ifdef G4MULTITHREADED
    // In MT mode these functions are called in the following order for each run:
    // Called on the master thread before starting the vis sub-thread.
    void DoneWithMasterThread() override;
    // Called on the master thread after starting the vis sub-thread.
    void MovingToVisSubThread() override;
    // Called on the vis sub-thread when waiting for events.
    void SwitchToVisSubThread() override;
    // Called on the vis sub-thread when all events have been processed.
    void DoneWithVisSubThread() override;
    // Called on the vis sub-thread when all events have been processed.
    // virtual void MovingToMasterThread ();  Not used in G4OpenGLQtViewer.
    // Called on the master thread after the vis sub-thread has terminated.
    void SwitchToMasterThread() override;

    inline void SetQGLContextVisSubThread(QThread* th) { fQGLContextVisSubThread = th; }
    inline void SetQGLContextMainThread(QThread* th) { fQGLContextMainThread = th; }
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

    void EnableClipperWidget() override;

    void SetWidgetInteractor(vtkAbstractWidget* widget) override;

  private:
    G4UIQt* fUiQt;
    QWidget* fGLWidget;

#ifdef G4MULTITHREADED
    QThread* fQGLContextVisSubThread;
    QThread* fQGLContextMainThread;
#endif

    // safe to use in serial mode
    G4AutoLock* lWaitForVisSubThreadQtOpenGLContextInitialized;
    G4AutoLock* lWaitForVisSubThreadQtOpenGLContextMoved;
};

#endif
