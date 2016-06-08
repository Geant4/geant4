//
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
//
// $Id: G4LizardSystem.cc,v 1.8.4.1 2001/06/28 19:07:46 gunter Exp $
// GEANT4 tag $Name:  $
//
// 
// Guy Barrand 14 September 2000

#ifdef G4ANALYSIS_BUILD_LIZARD

#include <IHistogram1D.h>
#include <IHistogram2D.h>
//#include <IHistogramFactory.h>

#ifdef AIDA_DONT_USE_STD
#define AIDA_STD
#else
#define AIDA_STD std
#endif

// Lizard :
#include <AIDAHistogramFactory.h>
#include "AIDA_Plotter/AIDAPlotter.h"
#include "Interfaces/IVector.h"
#include "Interfaces/IVectorFactory.h"

#include "G4ios.hh"

#include "G4LizardSystem.hh"

G4LizardSystem::G4LizardSystem (
 const G4String& aName
)
:fName(aName)
,fHistogramFactory(0)
,fVectorFactory(0)
,fPlotter(0)
{
  fHistogramFactory = createIHistogramFactory();
  if (fHistogramFactory == 0) {
    G4cout << "Can't find the histogram factory." << G4std::endl;
  }

  fVectorFactory = createIVectorFactory();
  if (fVectorFactory == 0) {
    G4cout << "Can't find the vector factory." << G4std::endl;
  }

  fPlotter = createIPlotter(); // create plotter with 1 zone
  if (fPlotter == 0) {
    G4cout << "Can't create a plotter." << G4std::endl;
  }
}
G4LizardSystem::~G4LizardSystem (){
  delete fHistogramFactory;
}
const G4String& G4LizardSystem::GetName() const {
  return fName;
}
IHistogramFactory* G4LizardSystem::GetHistogramFactory() {
  return fHistogramFactory;
}
void G4LizardSystem::Store(IHistogram* aHistogram,const G4String& aStorageID){
  if (!aHistogram) return;
  if (aStorageID == "") return;

  if (aHistogram->dimensions()==1) {
    IVector* v = 
      fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(aHistogram));
    if (v == 0) {
      G4cerr << "G4LizardSystem::Store> ERROR: could not create vector !" <<
G4std::endl;
      return;
    }
    v->toAscii(aStorageID.c_str());
    delete v;			// don't forget this !
  }

  if(aHistogram->dimensions()==2) {
    IVector* v = 
      fVectorFactory->from2D(dynamic_cast<IHistogram2D*>(aHistogram));
    if (v == 0) {
      G4cerr << "G4LizardSystem::Store> ERROR: could not create vector !" <<
G4std::endl;
      return;
    }
    v->toAscii(aStorageID.c_str());
    delete v;			// don't forget this !
  }
}
void G4LizardSystem::Plot(IHistogram* aHistogram){
  if(!fVectorFactory) return;
  if(!fPlotter) return;
  if(!aHistogram) return;

  if(aHistogram->dimensions()==1) {
    IVector* v = 
      fVectorFactory->from1D(dynamic_cast<IHistogram1D*>(aHistogram));
    if (v == 0) {
      G4cerr << "ERROR: could not create vector !" << G4std::endl;
      return;
    }
    fPlotter->plot(v);
    fPlotter->refresh();
    delete v;			// don't forget this !
  }

  if(aHistogram->dimensions()==2) {
    IVector* v = 
      fVectorFactory->from2D(dynamic_cast<IHistogram2D*>(aHistogram));
    if (v == 0) {
      G4cerr << "ERROR: could not create vector !" << G4std::endl;
      return;
    }
    fPlotter->plot(v);
    fPlotter->refresh();
    delete v;			// don't forget this !
  }

  //fPlotter->psPrint("posplot.ps");
}

#endif
