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
// ClassName:   test31Histo
//
//
// Author:      V.Ivanchenko 30/01/01
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "test31Histo.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmCalculator.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo* test31Histo::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo* test31Histo::GetPointer()
{
  if(!fManager) {
    fManager = new test31Histo();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::test31Histo()
{
  verbose = 0;
  histName = G4String("histo.hbook");
  ntup = 0;
  nHisto = 1;
  maxEnergy = 0.0;
  nTuple = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::~test31Histo()
{
  histo.clear();
  G4cout << "test31Histo: Histograms are deleted for " << theName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::BeginOfHisto(G4int num)
{
  if(0 < verbose) G4cout << "test31Histo # " << num << " started " << G4endl;
  zend     = 0.0;
  zend2    = 0.0;
  zEvt     = 0.0;
  etot     = 0.0;

  n_evt  = 0;
  n_elec = 0;
  n_posit= 0;
  n_gam  = 0;
  n_step = 0;

  n_charged_leak = 0;
  n_gam_leak = 0;
  n_charged_back = 0;
  n_gam_back = 0;

  if(0 < nHisto) bookHisto();

  if(verbose > 0) {
    G4cout << "test31Histo: Histograms are booked and run has been started"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::EndOfHisto()
{

  G4cout << "test31Histo: End of run actions are started" << G4endl;

  // Zend average

  G4cout<<"===================================================================="<<G4endl;
  if(zEvt > 0.0) {
    zend  /= zEvt;
    zend2 /= zEvt;
    zend2 -= zend*zend;
    G4double sig = 0.0;
    if(zend2 > 0.) sig = sqrt(zend2);
    zend2 = sig / sqrt(zEvt);
    G4cout << std::setprecision(5) << "Range(mm)= " << zend/mm
           << "; Stragling(mm)= " << sig/mm
           << std::setprecision(2) << " +- " << zend2/mm
           << "    " << zEvt << " events for range" << G4endl;
  }
  G4double x = (G4double)n_evt;
  if(n_evt > 0) x = 1.0/x;
  etot *= x;
  G4double xe = x*(G4double)n_elec;
  G4double xg = x*(G4double)n_gam;
  G4double xp = x*(G4double)n_posit;
  G4double xs = x*(G4double)n_step;
  G4double xcl = x*(G4double)n_charged_leak;
  G4double xgl = x*(G4double)n_gam_leak;
  G4double xcb = x*(G4double)n_charged_back;
  G4double xgb = x*(G4double)n_gam_back;
  G4cout                    << "Number of events               " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average energy deposit         " << etot/MeV << " MeV" << G4endl;
  G4cout << std::setprecision(4) << "Average number of e-           " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma        " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+           " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps        " << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of leak changed " << xcl << G4endl;
  G4cout << std::setprecision(4) << "Average number of leak gamma   " << xgl << G4endl;
  G4cout << std::setprecision(4) << "Average number of back changed " << xcb << G4endl;
  G4cout << std::setprecision(4) << "Average number of back gamma   " << xgb << G4endl;
  G4cout<<"===================================================================="<<G4endl;

  //TableControl();

   // Write histogram file
  if(0 < nHisto) {
    for(G4int i=0; i<nHisto; i++) {(histo[i])->scale(x);}
    tree->commit();
    std::cout << "Closing the tree..." << std::endl;
    tree->close();
    G4cout << "Histograms and Ntuples are saved" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveEvent()
{
  if(ntup) {
    ntup->addRow();
  }                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& parname, G4double val)
{
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& parname,G4double val, G4double)
{
  if(ntup) ntup->fill( ntup->findColumn(parname), (float)val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::bookHisto()
{
  G4double zmax = (AbsorberThickness + gap) * NumberOfAbsorbers / mm;
  G4cout << "test31Histo: Histograms will be saved to the file <" 
         << histName << ">"
         << " AbsThick(mm)= " << AbsorberThickness/mm
         << " Nabs= " << NumberOfAbsorbers
         << " zmax= " << zmax
         << " nHisto= " << nHisto
         << G4endl;

  // Creating the analysis factory
  std::auto_ptr< AIDA::IAnalysisFactory > af( AIDA_createAnalysisFactory() );

  // Creating the tree factory
  std::auto_ptr< AIDA::ITreeFactory > tf( af->createTreeFactory() );

  // Creating a tree mapped to a new hbook file.
  tree = tf->create(histName,"hbook",false,false);
  G4cout << "Tree store : " << tree->storeName() << G4endl;
 
  histo.resize(nHisto);

  // Creating a histogram factory, whose histograms will be handled by the tree
  std::auto_ptr< AIDA::IHistogramFactory > hf(af->createHistogramFactory( *tree ));

  // Creating an 1-dimensional histograms in the root directory of the tree

  if(0 < nHisto) histo[0] = hf->createHistogram1D("10",
    "Energy deposit (MeV) in absorber (mm)",NumberOfAbsorbers,0.0,zmax);

  if(1 < nHisto) histo[1] = hf->createHistogram1D("11",
    "Energy (MeV) of secondary electrons",50,0.0,maxEnergy/MeV);

  if(2 < nHisto) histo[2] = hf->createHistogram1D("12",
    "Theta (degrees) of delta-electrons",36,0.0,180.);

  if(3 < nHisto) histo[3] = hf->createHistogram1D("13",
    "Energy (MeV) of secondary gamma",50,0.0,maxEnergy/MeV);

  if(4 < nHisto) histo[4] = hf->createHistogram1D("14",
    "Theta (degrees) of secondary gamma",36,0.0,180.);

  // Creating a tuple factory, whose tuples will be handled by the tree
  std::auto_ptr< AIDA::ITupleFactory > tpf( af->createTupleFactory( *tree ) );

  // If using Anaphe HBOOK implementation, there is a limitation on the
  // length of the variable names in a ntuple
  if(nTuple) ntup = tpf->create( "100", "Range/Energy",
  "float tkin mass beta xend, yend, zend, ltpk, tend, teta, loss, dedx, back, leak, edep" );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEnergy(G4double edep, G4double z)
{
  etot += edep;
  if(0 < nHisto) histo[0]->fill((float)z/mm, (float)edep/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEndPoint(G4double z)
{
  zend  += z;
  zend2 += z*z;
  zEvt  += 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddDeltaElectron(const G4DynamicParticle* elec)
{
  n_elec++;
  if(1 < nHisto) histo[1]->fill((float)elec->GetKineticEnergy()/MeV,1.0);
  if(2 < nHisto)
     histo[2]->fill((float)(elec->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddPhoton(const G4DynamicParticle* ph)
{
  n_gam++;
  if(3 < nHisto) histo[3]->fill((float)ph->GetKineticEnergy()/MeV,1.0);
  if(4 < nHisto)
     histo[4]->fill((float)(ph->GetMomentumDirection()).theta()/deg,1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddParticleLeak(const G4DynamicParticle* dp)
{
  if(dp->GetDefinition() == G4Gamma::Gamma()) {
    n_gam_leak++;
  } else if (dp->GetCharge() != 0.0) {
    n_charged_leak++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddParticleBack(const G4DynamicParticle* dp)
{
  if(dp->GetDefinition() == G4Gamma::Gamma()) {
    n_gam_back++;
  } else if (dp->GetCharge() != 0.0) {
    n_charged_back++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::TableControl()
{
  G4EmCalculator cal;

// parameters
  G4double tmin = 0.001*keV;
  G4double tmax = 100.*GeV;
  G4int    nbin = 1100;
  G4int   index = 1;
  G4ParticleDefinition* part = G4Proton::Proton();
  //  G4ParticleDefinition* part = G4Electron::Electron();
//  G4ParticleDefinition* part = G4Alpha::Alpha();
  cal.PrintDEDXTable(part);
  cal.PrintRangeTable(part);
  cal.PrintInverseRangeTable(part);

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4LossTableManager* theManager = G4LossTableManager::Instance();
  const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(index);

  G4double xmin = log10(tmin);
  G4double xmax = log10(tmax);
  G4double step = (xmax - xmin)/(G4double)nbin;
  G4double x    = xmin;
  G4cout << "====================================================================" << G4endl;
  G4cout << "   Tables control for  " << part->GetParticleName() << G4endl;
  G4cout << "   Material            " << couple->GetMaterial()->GetName() << G4endl;
  G4cout << "====================================================================" << G4endl;

  for(G4int i=0; i<=nbin; i++) {
    G4double e  = pow(10.,x);
    G4double dedx = theManager->GetDEDX(part,e,couple);
    G4double r = theManager->GetRange(part,e,couple);
    G4double de = (theManager->GetEnergy(part,r,couple)/e - 1.)*100.;
    G4cout << i << ".   e(MeV)= " << e/MeV << ";    dedx(MeV/mm)= " << dedx*mm/MeV
           << ";  r(mm)= " << r/mm << ";  deltaE(%)= " << de << G4endl;
    x += step;
  }

  G4cout << "====================================================================" << G4endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
