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
//---------------------------------------------------------------------------
//
// ClassName:   EmAnalysis
//
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "EmAnalysis.hh"
#include "EmAnalysisMessenger.hh"
#include "Histo.hh"
#include "G4EmCalculator.hh"
#include "G4MuonPlus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EmAnalysis::EmAnalysis()
{
  verbose = 0;
  nHisto  = 0;
  nbins   = 120;
  cut     = 1000.0*MeV;
  xmin    = -3.0;
  xmax    = 3.0;
  histo   = Histo::GetInstance();
  messenger  = new EmAnalysisMessenger(this);
  calculator = new G4EmCalculator();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EmAnalysis::~EmAnalysis()
{
  delete messenger;
  delete calculator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::clear()
{
  histo->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int EmAnalysis::AddHistOnCrossSection(const G4String& part_name, 
					const G4String& mat_name,
					const G4String& proc_name, 
					const G4String& proc_type,
					const G4String& h_id)
{
  G4String name = proc_type + " " + part_name + " " + mat_name + " " + proc_name;
  G4String hid  = h_id;
  if(hid == "") {
    char buffer [4];
    sprintf(buffer,"%d",nHisto);
    hid = buffer;
    if(nHisto>0) {
      G4bool end = false;
      do {
	end = true;
	for(G4int i=0; i<nHisto-1; i++) {
	  if(hid == idstr[i]) {
	    hid += "000";
	    end = false;
	    break;
	  }
	}
      } while (end);
    }
  }
  particle.push_back(part_name);
  material.push_back(mat_name);
  process.push_back(proc_name);
  idstr.push_back(hid);
  atype.push_back(proc_type);
  G4int id = histo->add1D(hid,name,nbins,xmin,xmax,1.0);
  nHisto++;
  histid.push_back(id);
  if(0 < verbose) {
    G4cout << "EmAnalysis: Add 1D histo ID= " << id << " <" << name << "> " << G4endl;
  }
  return id;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::setHisto1D(G4int id, G4int nbins, G4double emin, G4double emax)
{
  for(G4int i=0; i<nHisto; i++) {
    if(id == histid[i]) {
      histo->setHisto1D(id,nbins,log10(emin),log10(emax),1.0);
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::setHisto1D(const G4String& h_id, G4int nbins, 
                                  G4double emin, G4double emax)
{
  for(G4int i=0; i<nHisto; i++) {
    if(h_id == idstr[i]) {
      setHisto1D(histid[i],nbins,emin,emax);
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::activate1D(G4int id, G4bool val)
{
  for(G4int i=0; i<nHisto; i++) {
    if(id == histid[i]) {
      histo->activate(id,val);
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::activate1D(const G4String& h_id, G4bool val)
{
  for(G4int i=0; i<nHisto; i++) {
    if(h_id == idstr[i]) {
      histo->activate(histid[i],val);
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::setVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::PrintHist(G4int val)
{
  histo->ListHistogram(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EmAnalysis::saveToFile()
{
  G4EmCalculator  calc;
  histo->book();

  G4double mu_mass = G4MuonPlus::MuonPlus()->GetPDGMass();
  G4double eth     = 2.0*mu_mass*mu_mass/electron_mass_c2 - electron_mass_c2;

  G4double r_muon  = elm_coupling/mu_mass;  //classical particle radius
  G4double mufac   = pi*r_muon*r_muon/3.;   //constant in crossSection


  for(G4int ii=0; ii<nHisto; ii++) {
    G4int i = histid[ii];
    if(histo->IsActive(i)) {
      const G4Material* mat = calc.FindMaterial(material[ii]);
      G4int nb = histo->NumberOfBins(i);
      G4double x0    = histo->MinBin(i);
      G4double step0 = (histo->MaxBin(i) - x0)/((G4double)nb);
      G4double x     = pow(10.0, x0 + step0*0.5);
      G4double step  = pow(10.0, step0);
      G4bool isdedx = false;
      G4bool iscspe = false;
      G4bool ismumu = false;
      G4double unit = microbarn;

      if(atype[ii] == "dedx") isdedx = true;
      else if(atype[ii] == "cs_per_e") iscspe = true;
      else if(atype[ii] == "mumu") ismumu = true;

      for(G4int j=0; j<nb; j++) {
        G4double s = 0.0;
        if(isdedx) {
          s = calc.ComputeDEDX(particle[ii],material[ii],process[ii],x,cut);
          s *= mm/MeV;
        } else if(iscspe) { 
          s = calc.ComputeCrossSection(particle[ii],material[ii],process[ii],x,cut);
  	  s /= (unit*mat->GetElectronDensity());
        } else if(ismumu) {
          G4double y = eth/x;
          if(y < 1.0) s = mufac*y*(1.0 + 0.5*y)*sqrt(1.0 - y)/unit; 
        } else { 
          s = calc.ComputeCrossSection(particle[ii],material[ii],process[ii],x,cut);
  	  s /= (unit*mat->GetTotNbOfAtomsPerVolume());
        }
	histo->fill(i,log10(x),s); 
	x *= step;
      }
    }
  }
  histo->save();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
