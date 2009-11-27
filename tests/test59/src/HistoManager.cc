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
// $Id: HistoManager.cc,v 1.1 2009-11-27 16:06:28 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------

//#define debug

#include "HistoManager.hh"
#include "HistoMessenger.hh"

HistoManager::HistoManager()
{
  histoMessenger = new HistoMessenger(this);
}

HistoManager::~HistoManager()
{
  for(unsigned i=0; i < hisBins.size(); ++i)
  {
    delete hisBins[i];
    delete hisBinX[i];
  }
  delete histoMessenger;
}

void HistoManager::BookHisto()
{
  //for(G4int i=0; i < 2; ++i)
  //{
    G4int nBins = 50;
    nBinHis.push_back(nBins);
    G4double xMin=sqr(std::log(.001));
    minValH.push_back(xMin);
    G4double xMax=sqr(std::log(.09));
    G4double step=(xMax-xMin)/nBins;
    stepHis.push_back(step);
    std::vector<G4double>* hist = new std::vector<G4double>;
    hist->resize(nBins);
    std::vector<G4double>* binx = new std::vector<G4double>;
    binx->resize(nBins);
    for(G4int i=0; i<nBins; ++i) // Not necessary after resizing?
    {
      (*hist)[i]=0.;
      (*binx)[i]=0.;
    }
    hisBins.push_back(hist);
    hisBinX.push_back(binx);
  //}
}

void HistoManager::SaveHisto()
{
  G4double lgx=minValH[0];
  G4double low=std::exp(lgx);
  G4double dlg=stepHis[0];
  std::vector<G4double>* hb=hisBins[0];
  std::vector<G4double>* xb=hisBinX[0];
  for(G4int i=0; i<nBinHis[0]; ++i)
  {
    lgx+=dlg;
    G4double hx=std::exp(-std::sqrt(lgx));
    G4double dx=hx-low;
    G4double z=(*hb)[i];
    G4double x=(hx+low)/2.;
    G4double y=0.;
    if(z>0.)
    {
      x=(*xb)[i]/z;
      y=z/dx;
    }
    G4double d=std::sqrt(z+1.)/dx;
    G4cout<<x<<" "<<y<<" "<<d<<G4endl;
    low=hx;
  }
}

// Called from TrackingManager
void HistoManager::FillHisto(G4int ih, G4double e, G4double w)
{
  G4int nh = minValH.size();
  if(ih >= nh) G4cout<<"-Warning-HistoManager::FillHisto: Hist# "<<ih<<" >= "<<nh<<G4endl;
  else if(!ih)
  {
    G4double ae=std::abs(e);
    if(e > std::exp(-std::sqrt(minValH[0])))
    {
      G4double dx=stepHis[ih];
#ifdef debug
      G4cout<<"HistoManager::FillHisto: ae = "<<ae<<", e = "<<e<<G4endl;
#endif
      G4int ix=(G4int)((sqr(std::log(ae))-minValH[ih])/dx);
#ifdef debug
      G4cout<<"HistoManager::FillHisto: ix = "<<ix<<", dx = "<<dx<<", ih="<<ih<<G4endl;
#endif
      if(ix < nBinHis[0])
      {
        ((*(hisBins[ih]))[ix])+=w;
        ((*(hisBinX[ih]))[ix])+=ae*w;
      }
    }
  }
}
