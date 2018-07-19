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
#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TGFrame.h>

class TGWindow;
class TGMainFrame;
class TRootEmbeddedCanvas;
class TGTab;
class TCanvas;

class CanvasInTab : public TGMainFrame
{
private:
  TGTab* fpTab;
  std::vector<TRootEmbeddedCanvas*> fEcanvas;
  TGLayoutHints *fHintPlots;
  
public:
  CanvasInTab(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~CanvasInTab();
  
  size_t AddCanvas(const char* name = "New tab");
  TCanvas* GetCanvas(int i);
  size_t GetNCanvas() const
  {
    return fEcanvas.size();
  }
  
  void SaveCanvas();
  
  ClassDef(CanvasInTab,0)
};
