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
// G4VRangeToEnergyConverter class implementation
//
// Author: H.Kurashige, 05 October 2002 - First implementation
// --------------------------------------------------------------------

#include "G4VRangeToEnergyConverter.hh"
#include "G4ParticleTable.hh"
#include "G4Element.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

#ifdef G4MULTITHREADED
G4Mutex G4VRangeToEnergyConverter::theMutex = G4MUTEX_INITIALIZER;
#endif

G4double  G4VRangeToEnergyConverter::Emin = 0.0;
G4double  G4VRangeToEnergyConverter::Emax = 0.0;

std::vector<G4double>* G4VRangeToEnergyConverter::Energy = nullptr;

G4int G4VRangeToEnergyConverter::NbinPerDecade = 50;
G4int G4VRangeToEnergyConverter::Nbin = 350;

// --------------------------------------------------------------------
G4VRangeToEnergyConverter::G4VRangeToEnergyConverter()
{
  if(nullptr == Energy)
  {
    Energy = new std::vector<G4double>(Nbin + 1);
    FillEnergyVector(1*CLHEP::keV, 10.0*CLHEP::GeV);
  }
}

// --------------------------------------------------------------------
G4VRangeToEnergyConverter::~G4VRangeToEnergyConverter()
{
  if(nullptr != Energy) 
  { 
    delete Energy;
    Energy = nullptr; 
  }
}

// --------------------------------------------------------------------
G4double G4VRangeToEnergyConverter::Convert(const G4double rangeCut, 
                                            const G4Material* material) 
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>3)
  {
    G4cout << "G4VRangeToEnergyConverter::Convert() - ";
    G4cout << "Convert for " << material->GetName() 
	   << " with Range Cut " << rangeCut/mm << "[mm]" << G4endl;
  }
#endif

  G4double cut = 0.0;
  if(fPDG == 22) 
  {  
    cut = ConvertForGamma(rangeCut, material);
  }
  else 
  {
    cut = ConvertForElectron(rangeCut, material);

    const G4double tune = 0.025*CLHEP::mm*CLHEP::g/CLHEP::cm3;
    const G4double lowen = 30.*CLHEP::keV; 
    if(cut < lowen)
    {
      //  corr. should be switched on smoothly   
      cut /= (1.+(1.-cut/lowen)*tune/(rangeCut*material->GetDensity())); 
    }
  }

  cut = std::max(Emin, std::min(cut, Emax));
  return cut;
}

// --------------------------------------------------------------------
void G4VRangeToEnergyConverter::SetEnergyRange(const G4double lowedge,
                                               const G4double highedge)
{
  G4double ehigh = std::min(Emax, highedge);
  if(ehigh > lowedge)
  {
    FillEnergyVector(lowedge, ehigh);
  }  
}

// --------------------------------------------------------------------
G4double G4VRangeToEnergyConverter::GetLowEdgeEnergy()
{
  return Emin;
}
    
// --------------------------------------------------------------------
G4double G4VRangeToEnergyConverter::GetHighEdgeEnergy()
{
  return Emax;
}

// --------------------------------------------------------------------

G4double G4VRangeToEnergyConverter::GetMaxEnergyCut()
{
  return Emax;
}

// --------------------------------------------------------------------
void G4VRangeToEnergyConverter::SetMaxEnergyCut(const G4double value)
{
  if(value > Emin)
  {
    FillEnergyVector(Emin, value);
  }
}

// --------------------------------------------------------------------
void G4VRangeToEnergyConverter::FillEnergyVector(const G4double emin, 
                                                 const G4double emax)
{
  if(emin == Emin && emax == Emax) { return; }
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&theMutex);
  if(emin == Emin && emax == Emax) { return; }
#endif
  Emin = emin;
  Emax = emax;
  Nbin = NbinPerDecade*static_cast<G4int>(std::log10(emax/emin));
  Energy->resize(Nbin + 1);
  (*Energy)[0] = emin;
  (*Energy)[Nbin] = emax;
  G4double fact = G4Log(emax/emin)/Nbin;
  for(G4int i=1; i<Nbin; ++i)
  {
    (*Energy)[i] = emin*G4Exp(i * fact);
  }
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&theMutex);
#endif
}

// --------------------------------------------------------------------
G4double 
G4VRangeToEnergyConverter::ConvertForGamma(const G4double rangeCut, 
                                           const G4Material* material)
{
  const G4ElementVector* elm = material->GetElementVector();
  const G4double* dens = material->GetAtomicNumDensityVector();

  // fill absorption length vector
  G4int nelm = material->GetNumberOfElements();
  G4double range1 = 0.0;
  G4double range2 = 0.0;
  G4double e1 = 0.0;
  G4double e2 = 0.0;
  for (G4int i=0; i<Nbin; ++i)
  {
    e2 = (*Energy)[i];
    G4double sig = 0.;
    
    for (G4int j=0; j<nelm; ++j)
    {
      sig += dens[j]*ComputeValue((*elm)[j]->GetZasInt(), e2); 
    }
    range2 = (sig > 0.0) ? 5./sig : DBL_MAX;
    if(i == 0 || range2 < rangeCut)
    {
      e1 = e2;
      range1 = range2;
    }
    else
    {
      break;
    }
  }
  return LiniearInterpolation(e1, e2, range1, range2, rangeCut);
}

// --------------------------------------------------------------------
G4double 
G4VRangeToEnergyConverter::ConvertForElectron(const G4double rangeCut, 
                                              const G4Material* material)
{
  const G4ElementVector* elm = material->GetElementVector();
  const G4double* dens = material->GetAtomicNumDensityVector();

  // fill absorption length vector
  G4int nelm = material->GetNumberOfElements();
  G4double dedx1 = 0.0;
  G4double dedx2 = 0.0;
  G4double range1 = 0.0;
  G4double range2 = 0.0;
  G4double e1 = 0.0;
  G4double e2 = 0.0;
  G4double range = 0.;
  for (G4int i=0; i<Nbin; ++i)
  {
    e2 = (*Energy)[i];
    dedx2 = 0.0;
    for (G4int j=0; j<nelm; ++j)
    {
      dedx2 += dens[j]*ComputeValue((*elm)[j]->GetZasInt(), e2); 
    }
    range += (dedx1 + dedx2 > 0.0) ? 2*(e2 - e1)/(dedx1 + dedx2) : 0.0;
    range2 = range;
    if(range2 < rangeCut)
    {
      e1 = e2;
      dedx1 = dedx2;
      range1 = range2;
    }
    else
    {
      break;
    }
  }
  return LiniearInterpolation(e1, e2, range1, range2, rangeCut);
}

// --------------------------------------------------------------------
