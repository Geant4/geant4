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
#include "Histo.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4LossTableManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EmCalculator.hh"
#include "EmAnalysis.hh"
#include "G4EmCorrections.hh"
#include "G4NistManager.hh"
#include "G4BraggModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4EnergyLossForExtrapolator.hh"
#include "G4NucleiProperties.hh"
#include "G4WaterStopping.hh"
#include <iomanip>
#include <fstream>

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
  nHisto  = 0;
  maxEnergy = 0.0;
  nTuple  = false;
  tables  = true;
  histo   = Histo::GetInstance();
  //  ema = new EmAnalysis();
  histoID.resize(7);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31Histo::~test31Histo()
{
  delete extra;
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

  n_mumu = 0;
  n_pipi = 0;

  if(0 < nHisto) {
    if(num == 0) bookHisto();
    histo->book();

    if(verbose > 0) {
      G4cout << "test31Histo: Histograms are booked and run has been started"
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::EndOfHisto()
{

  G4cout << "test31Histo: End of run actions" << G4endl;

  // Zend average

  G4cout<<"===================================================================="<<G4endl;
  G4cout                         << "Initial particle               " 
				 << beamParticle->GetParticleName() 
                                 << "   Ekin(GeV)= " << beamEnergy/GeV << G4endl;
  if(zEvt > 0.0) {
    zend  /= zEvt;
    zend2 /= zEvt;
    zend2 -= zend*zend;
    G4double sig = 0.0;
    if(zend2 > 0.) sig = std::sqrt(zend2);
    zend2 = sig / std::sqrt(zEvt);
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
  G4double xmu = x*(G4double)n_mumu;
  G4double xpi = x*(G4double)n_pipi;
  G4cout                         << "Number of events               " << n_evt <<G4endl;
  G4cout << std::setprecision(4) << "Average energy deposit         " << etot/MeV 
	 << " MeV" << G4endl;
  G4cout << std::setprecision(4) << "Average number of e-           " << xe << G4endl;
  G4cout << std::setprecision(4) << "Average number of gamma        " << xg << G4endl;
  G4cout << std::setprecision(4) << "Average number of e+           " << xp << G4endl;
  G4cout << std::setprecision(4) << "Average number of steps        " << xs << G4endl;
  G4cout << std::setprecision(4) << "Average number of leak charged " << xcl << G4endl;
  G4cout << std::setprecision(4) << "Average number of leak gamma   " << xgl << G4endl;
  G4cout << std::setprecision(4) << "Average number of back charged " << xcb << G4endl;
  G4cout << std::setprecision(4) << "Average number of back gamma   " << xgb << G4endl;
  G4cout << std::setprecision(4) << "Average number of mu+mu-       " << xmu << G4endl;
  G4cout << std::setprecision(4) << "Average number of pi+pi-       " << xpi << G4endl;
  G4cout<<"===================================================================="<<G4endl;

  if(tables) {
    TableControl();
    tables = false;
  }
  //  MuonTest();
  //  ElectronTest();
  if(0 < nHisto) {

    // normalise histograms
    for(G4int i=0; i<nHisto; i++) {
      histo->scale(i,x);
    }
    histo->save();
  }
  G4cout<<"=========   End of tets31Histo  ============================"<<G4endl;  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveEvent()
{
  if(nTuple) histo->addRow(0);        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& parname, G4double val)
{
  if(nTuple) histo->fillTuple(0, parname, val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::SaveToTuple(const G4String& parname,G4double val, G4double)
{
  if(nTuple) histo->fillTuple(0, parname, val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::bookHisto()
{
  zmax = (AbsorberThickness + gap) * NumberOfAbsorbers / mm;
  G4cout << "test31Histo: "
         << " AbsThick(mm)= " << AbsorberThickness/mm
         << " Nabs= " << NumberOfAbsorbers
         << " zmax= " << zmax
         << " nHisto= " << nHisto
         << G4endl;

  extra = new G4EnergyLossForExtrapolator(1);

  // Creating an 1-dimensional histograms in the root directory of the tree
  
  if(nHisto >= 0) {
    G4double em = maxEnergy/MeV;
    histoID[0] = 
      histo->add1D("10","Energy deposit (MeV) in absorber (mm)",
		   NumberOfAbsorbers,0.0,zmax/mm,mm);

    histoID[1] =  histo->add1D("11",
      "Energy (MeV) of secondary electrons",50,0.0,em,MeV);

    histoID[2] =  histo->add1D("12",
      "Theta (degrees) of delta-electrons",36,0.0,180.,degree);

    histoID[3] =  histo->add1D("13",
      "Energy (MeV) of secondary gamma",50,0.0,em,MeV);

    histoID[4] =  histo->add1D("14",
      "Theta (degrees) of secondary gamma",36,0.0,180.,degree);

    histoID[5] =  histo->add1D("15",
      "Theta (degrees) of primary",50,0.0,10.,degree);

    histoID[6] =  histo->add1D("16",
      "Delta Energy (MeV) of Reconstruction",100,-em,em,MeV);

    for(G4int i=0; i<nHisto; i++) {histo->activate(histoID[i], true);}
  }

  if(nTuple){ 
    histo->addTuple( "100", "Range/Energy",
     "float tkin mass beta xend, yend, zend, ltpk, tend, teta, loss, dedx, back, leak, edep" );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddEnergy(G4double edep, G4double z)
{
  etot += edep;
  if(0 < nHisto) histo->fill(histoID[0], z, edep/MeV);
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
  if(1 < nHisto) histo->fill(histoID[1],elec->GetKineticEnergy(),1.0);
  if(2 < nHisto) histo->fill(histoID[2],elec->GetMomentumDirection().theta(),1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddPhoton(const G4DynamicParticle* ph)
{
  n_gam++;
  if(3 < nHisto) histo->fill(histoID[3],ph->GetKineticEnergy(),1.0);
  if(4 < nHisto) histo->fill(histoID[4],(ph->GetMomentumDirection()).theta(),1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddParticleLeak(const G4Track* track)
{
  const G4DynamicParticle* dp = track->GetDynamicParticle();
  if(dp->GetDefinition() == G4Gamma::Gamma()) {
    n_gam_leak++;
  } else if (dp->GetCharge() != 0.0) {
    n_charged_leak++;
    if(track->GetTrackID() == 1) {
      G4double tet = (dp->GetMomentumDirection()).theta();
      if(5 < nHisto) histo->fill(histoID[5],tet,1.0);
      G4double e0 = track->GetVertexKineticEnergy();
      G4double e1 = dp->GetKineticEnergy();
      G4double e2 = 
        extra->EnergyBeforeStep(e1,zmax,absMaterial,dp->GetDefinition());
      if(n_evt < 10)
        G4cout << "Extrapolation of primary E0(MeV)= " << e0/MeV 
	       << " E1(MeV)= " << e1/MeV 
	       << " Erec-E0(MeV)= " << (e2-e0)/MeV << G4endl;
      if(6 < nHisto) histo->fill(histoID[6],e2 - e0,1.0);      
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::AddParticleBack(const G4Track* track)
{
  if(track->GetDynamicParticle()->GetDefinition() == G4Gamma::Gamma()) {
    n_gam_back++;
  } else if (track->GetDynamicParticle()->GetCharge() != 0.0) {
    n_charged_back++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::TableControl()
{
  G4NistManager*  mman = G4NistManager::Instance();
  G4EmCorrections* emc = G4LossTableManager::Instance()->EmCorrections();
  G4EmCalculator cal;
  cal.SetVerbose(0);

  // parameters
  // G4double tmin = 1.*keV;
  // G4double tmax = 1.*GeV;
  // G4int    nbin = 60;
  G4String ion_name  = "ionIoni";
  G4String h_name    = "hIoni";
  G4String mu_name   = "muIoni";
  G4String proc_name = "eIoni";
  G4String part_name = "proton";
 
  const G4ParticleDefinition* part = cal.FindParticle(part_name);
  if(!part) return;
  // cal.PrintDEDXTable(part);
  // cal.PrintRangeTable(part);
  // cal.PrintInverseRangeTable(part);

  G4String mat_name  = "G4_WATER";
  G4Material* mat = mman->FindOrBuildMaterial(mat_name);
  //mat->SetChemicalFormula("H_2O");
  G4double fact = 0.001*gram/(MeV*cm2*mat->GetDensity());

  //  G4double xmin = std::log10(tmin);
  // G4double xmax = std::log10(tmax);
  // G4double step = (xmax - xmin)/(G4double)nbin;
  // G4double x    = xmin;

  const G4int ne = 42;
  G4double e0[ne]={0.025, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 
		   0.5,   1.0,   1.1,  1.2,  1.3, 1.4, 1.5, 
		   1.6,   1.7,   1.8,  1.9,  2.0, 2.1, 2.2,
                   2.3,   2.4,   2.5,  3.0,  5.0, 10., 15.0, 
                   20.0,   30.,  50.,  100., 200., 300., 400., 
		   500., 1000., 2000., 3000., 5000., 10000., 20000.};
  //                  30000., 100000., 300000., 1000000., 10000000.,};
  const G4int np = 8;
  G4String namep[np] = {"e-","mu+","pi-","proton","alpha", "C12[0.0]", 
			"Ar40[0.0]", "Pb208[0.0]"};

  G4int ii1 = 5;
  G4int ii2 = 6;

  G4WaterStopping wst;

  for(G4int ii=ii1; ii<ii2; ii++) {
    mat_name  = "G4_WATER";
    mat = mman->FindOrBuildMaterial(mat_name);
    fact = 0.001*gram/(MeV*cm2*mat->GetDensity());

    const G4ParticleDefinition* part = cal.FindParticle(namep[ii]);
    if(!part) break;
      
    if(ii == 1) proc_name = mu_name;
    else if(ii == 2) proc_name = h_name;
    else if(ii == 4) {proc_name = ion_name;}

    G4double AA = 1.0;
    if(ii >= 5) {
      G4int ZZ = part->GetAtomicNumber();
      AA = mman->GetAtomicMassAmu(ZZ);
    }

    G4cout << "================================================================" << G4endl;
    G4cout << "   Tables control for  " << namep[ii] << G4endl;
    G4cout << "   Material            " << mat_name << "  AA= "<< AA << G4endl;
    G4cout << "================================================================" << G4endl;

    G4cout << "  N   E(MeV/N)  Esc(MeV) dEdx_T/NIST "
	   << "dEdx_C/NIST  dEdx_G/NIST  dedx(MeV*cm^2/mg)" 
	   << std::setprecision(6) << G4endl;

    for(G4int ij=0; ij<ne; ij++) {
    
      G4double e = e0[ij];
      G4double e1 = e;
      if(ii >= 5) e1 *= AA;
      G4double dedx0 = cal.ComputeTotalDEDX(e1,part,mat,e1);
      G4double dedx  = cal.ComputeElectronicDEDX(e1,part,mat,e1);
      G4double dedx1 = cal.GetDEDX(e1,part,mat);
      G4double dedx2 = wst.GetElectronicDEDX(6, e1);
      G4cout << std::setw(3) << ij << "." 
             << std::setw(10) << e/MeV 
	     << std::setw(10) << e1/MeV
	     << std::setw(11) << dedx0/dedx2
	     << std::setw(11) << dedx/dedx2
	     << std::setw(11) << dedx1/dedx2
	     << std::setw(11) << dedx0*fact
	     << G4endl; 
    }
    /*
    G4double ee0 = 0.0;
    G4cout << "######## Ranges ########" << G4endl;
    for(G4int ik=0; ik<1000; ik++) {
      G4double e = std::pow(10.,ee0);
      G4double r = cal.GetRange(e,part,mat);
      G4cout << "E(MeV)= " << e << "  R(mm)= " << r << G4endl;
      ee0 += 0.003;
    }
    */
    //    G4bool icorr = true;
    G4bool icorr = false;
    if(icorr) {
      G4cout << "================================================================" << G4endl;
      G4cout << "             Ionisation Corrections" << G4endl;
      G4cout << "================================================================" << G4endl;

      const G4int nm = 7;
      G4String nmat[nm] = {"G4_H", "G4_C", "G4_Al", "G4_Cu", "G4_Ag", "G4_Au", "G4_Pb"};
      const G4int kkk = 25;
      G4double ek[kkk] = {0.3, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5,  
		          4.0, 4.5, 5.0, 6.5, 8.0,  12.5, 20.0, 30., 100., 
			  300., 1000., 3000., 10000., 100000., 1000000., 10000000.};

      G4double mass = part->GetPDGMass();
  
      G4double L, L0, L1, L2, Spin, KS, LS, S, S0, del, mk, dedx, fac, fact, dc, fs(0), nuc;
  
      for(G4int i=0; i<nm; i++) {
	dc = 0.0;
	mat = mman->FindOrBuildMaterial(nmat[i]);
	fact = gram/(MeV*cm2*mat->GetDensity());
	fac = 2.0*twopi_mc2_rcl2*(mat->GetElectronDensity())*fact;
	G4cout << "   New Material  " << mat->GetName() << G4endl;
	for(G4int j=0; j<kkk; j++) {
	  G4double e = ek[j]*MeV;
	  G4double tau = e/mass;
	  G4double gamma = 1.0 + tau;
	  G4double beta2 = tau*(tau + 2.0)/(gamma*gamma);
	  //        G4double dedx0 = cal.ComputeDEDX(e,part,"hIoni",mat)*fact;
    
	  L0   = emc->Bethe(part,mat,e);
	  Spin = emc->SpinCorrection(part,mat,e);
	  KS   = emc->KShellCorrection(part,mat,e);
	  LS   = emc->LShellCorrection(part,mat,e);
	  S    = emc->ShellCorrection(part,mat,e);
	  S0   = emc->ShellCorrectionSTD(part,mat,e);
	  L1   = emc->BarkasCorrection(part,mat,e);
	  L2   = emc->BlochCorrection(part,mat,e);
	  del  = -0.5*emc->DensityCorrection(part,mat,e);
	  mk   = 0.5*emc->MottCorrection(part,mat,e);
	  nuc  = fact*emc->NuclearDEDX(part,mat,e,false);
	  L = L0 + L1 + L2 +fs + del -S;
	  dedx = L*fac/beta2;
	  // if(0 == j) dc = (dedx0 - dedx)*MeV*MeV;
	  //	dedx += dc/(e*e);
	  //	  G4double x = S + del - L1 - L2;
	  G4cout << j+1 << ". " << ek[j] << " MeV "
		 << " L0= " << L0
	    //		 << " Spin= " << Spin
		 << " KSh= " << KS
		 << " LSh= " << LS
		 << " Sh= " << S 
	    //		 << " Sh0= " << S0 
		 << " L1= " << L1
		 << " L2= " << L2
		 << " dn= " << del
		 << " mott= " << mk
	    //	    	 << " fs= " << fs
		 << " L= " << L
		 << " dedx= " << dedx
	    //     << " dedx0= " << dedx0
		 << " nuc= " << nuc
		 << G4endl;
	}
	G4cout << "==============================================================" << G4endl;
      }
    }
  }

  G4cout << "==============  End of Table Control ===============================" << G4endl;

  G4bool ish = false;
  if (ish) {

    const G4ParticleDefinition* proton = cal.FindParticle("proton");
    G4BraggModel        bragg;
    G4BetheBlochModel   bethe;
    G4DataVector        empty;
    bragg.Initialise(proton, empty);
    bethe.Initialise(proton, empty);
    const G4Element* elm;
    const G4Material* ma;
   
    G4double e     = 2.0*MeV;
    G4double tau   = e/proton_mass_c2;
    G4double gam   = tau + 1.0;
    G4double bg2   = tau * (tau+2.0);
    G4double beta2 = bg2/(gam*gam);
    G4double eta   = beta2/(fine_structure_const*fine_structure_const);
    G4double fact  = 0.5 * beta2*eta*eta / twopi_mc2_rcl2;
    G4double dedx0, dedx1;

    for(G4int z=1; z<93; z++) {
      elm = mman->FindOrBuildElement(z, false);
      G4String nam = "G4_"+elm->GetSymbol();
      ma = mman->FindOrBuildMaterial(nam, false);
      //      G4cout << "Elm " << elm << "  mat " << ma << "   " << nam << G4endl;
    }
    for(G4int z=1; z<93; z++) {
      elm = mman->FindOrBuildElement(z, false);
      G4String nam = "G4_"+elm->GetSymbol();
      ma = mman->FindOrBuildMaterial(nam, false);
      dedx0 = bragg.ComputeDEDXPerVolume(ma,proton,e,GeV);
      dedx1 = bethe.ComputeDEDXPerVolume(ma,proton,e,GeV);
      G4cout << " " << (dedx1 - dedx0)*fact/(ma->GetElectronDensity()) << ",";
      if(z/10*10 == z) G4cout << G4endl;
    }
    G4cout << G4endl;
    G4double fe     = 8.0*MeV;
    G4double ftau   = fe/proton_mass_c2;
    G4double fgam   = ftau + 1.0;
    G4double fbg2   = ftau * (ftau+2.0);
    G4double fbeta2 = fbg2/(fgam*fgam);
    G4double ffact  = 0.5 * fbeta2 / twopi_mc2_rcl2;
    G4double fdedx0, fdedx1, s0, s1;
    G4cout << G4endl;
    for(G4int z=1; z<93; z++) {
      elm = mman->FindOrBuildElement(z, false);
      G4String nam = "G4_"+elm->GetSymbol();
      ma = mman->FindOrBuildMaterial(nam, false);
      fdedx0 = bragg.ComputeDEDXPerVolume(ma,proton,e,GeV)*ffact/ma->GetElectronDensity();
      fdedx1 = bethe.ComputeDEDXPerVolume(ma,proton,e,GeV)*ffact/ma->GetElectronDensity();
      s0 = fdedx1 - fdedx0;
      s1 = emc->ShellCorrectionSTD(proton,ma,fe);

      G4cout << " " << (s1*std::log(tau) - s0*std::log(ftau))/(ftau - tau) << ",";
      // << " s0= " << s0 << " s1= " << s1 << G4endl;;
      if(z/10*10 == z) G4cout << G4endl;
    }
    G4cout << G4endl;
  }

  G4bool ihist = false;
  if (ihist) {
    G4cout << "=================================================================" << G4endl;
    G4cout << "             Stopping Powers" << G4endl;
    G4cout << "=================================================================" << G4endl;

    G4String nm[7] = {"G4_Be", "G4_Al", "G4_Si", "G4_Ge", "G4_Fe", "G4_Ag", "G4_Au"};
    const G4String partc[2] = {"proton", "alpha"};
    const G4String proc[2] = {"hIoni", "ionIoni"};
    
    G4double e, se, sn, st, mce, mcn, mct; 
    char line[200];

    for(G4int ii=0; ii<2; ii++) {
    
      const G4ParticleDefinition* part = cal.FindParticle(partc[ii]);

      for(G4int i=0; i<7; i++) {
	mat = mman->FindOrBuildMaterial(nm[i]);
        G4cout << "  Particle  " << partc[ii] << " in  Material  " 
	       << mat->GetName() << G4endl;
	G4double fact = gram/(MeV*cm2*mat->GetDensity());
        std::ifstream* fin = new std::ifstream();
        std::string fname = "stopping/" + partc[ii] + "_" + nm[i];
        std::string fnamef= fname + ".txt";
        fin->open(fnamef.c_str());
        if( !fin->is_open()) {
          G4cout << "Input file <" << fname << "> does not exist! Exit" << G4endl;
          exit(1);
        }
        //G4int i1 = histo->addCloud1D(fname);
        for(G4int j=0; j<8; j++) {fin->getline( line, 200);}
        G4int i2 =0;
        do {
          i2++;
          (*fin) >> e >> se >> sn >> st;
          e *= MeV;
          mce = fact*(cal.ComputeDEDX(e,part,proc[ii],mat));
          mcn = fact*emc->NuclearDEDX(part,mat,e,false);
          mct = mce + mcn;
          G4double diff = 100.*(mct/st - 1.0);
	  /*
          G4cout << e/MeV 
                 << " NIST:  dedx= " << se
                 << " nuc= " << sn
                 << " tot= " << st
                 << " G4:  dedx= " << mce
                 << " nuc= " << mcn
                 << " tot= " << mct
                 << " diff= " << diff << " %"
                 << G4endl;
	  */
          if(ii==0 && i==0) G4cout << " " << e;
	  else              G4cout << " " << -diff;
          //histo->fill(i1,e,diff);
              
        } while ( std::fabs(e - 1000.) > MeV);
	  G4cout << G4endl;
	  G4cout << "========= n= " << i2 
		 << " ===========================================================" 
		 << G4endl;
        fin->close();
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::MuonTest()
{
  G4NistManager*  mman = G4NistManager::Instance();
  G4EmCalculator cal;
  cal.SetVerbose(0);

  const G4ParticleDefinition* part = cal.FindParticle("mu+");

  G4bool imu = false;
  if (imu) {
    G4cout << "====================================================================" << G4endl;
    G4cout << "             Stopping Powers" << G4endl;
    G4cout << "====================================================================" << G4endl;

    G4String nm[3] = {"G4_WATER", "G4_Al", "G4_Fe"};
    
    G4double energy[43] = {10.,  14.,  20., 30., 40., 80.,
                           100., 140., 200., 300., 400., 800., 
                           1000., 1400., 2000., 3000., 4000., 8000., 
                           10000., 14000., 20000., 30000., 40000., 80000., 
                           100000., 140000., 200000., 300000., 400000., 800000., 
                           1000000., 1400000., 2000000., 3000000., 4000000., 8000000., 
                           10000000., 14000000., 20000000., 30000000., 40000000., 80000000., 100000000.}; 
    G4double dedx[3][43] = {
      { 7.965, 6.213, 4.852, 3.764, 3.214, 2.413,
        2.270, 2.116, 2.026, 1.992, 1.999, 2.075,
        2.109, 2.166, 2.229, 2.300, 2.351, 2.470,
        2.507, 2.564, 2.625, 2.701, 2.760, 2.942,
        3.020, 3.166, 3.376, 3.709, 4.039, 5.352,
        6.014, 7.330, 9.327,12.654,16.018,29.603,
	36.462, 50.192, 70.964, 105.740, 140.768, 282.138, 353.358},
      { 6.188, 4.849, 3.802, 2.961, 2.533, 1.908,
        1.798, 1.688, 1.630, 1.616, 1.630, 1.711,
        1.745, 1.799, 1.858, 1.925, 1.971, 2.082,
        2.117, 2.172, 2.233, 2.312, 2.377, 2.594,
        2.694, 2.885, 3.167, 3.628, 4.093, 5.964,
        6.914, 8.805, 11.628,16.474,21.32,40.865,
	50.72, 70.429, 100.206, 149.938, 199.946, 401.196, 502.351},
      { 5.494, 4.321, 3.399, 2.654, 2.274, 1.717,
        1.616, 1.516, 1.463, 1.453, 1.467, 1.548,
        1.582, 1.637, 1.697, 1.767, 1.816, 1.936,
        1.975, 2.039, 2.113, 2.214, 2.303, 2.623,
        2.777, 3.082, 3.543, 4.304, 5.079, 8.221,
        9.820, 13.013, 17.877,25.974,34.162,67.147,
	83.761, 116.947, 167.024, 250.537, 334.408, 671.133, 840.063}
    }; 
    G4double dedxn[3][43] = {
      { 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0.,
	0., 0.001, 0.001, 0.001, 0.002, 0.004,
	0.005, 0.007, 0.009, 0.013, 0.018, 0.034,
	0.042, 0.059, 0.084, 0.125, 0.167, 0.337,
	0.423, 0.601, 0.870, 1.332, 1.803, 3.763,
	4.773, 6.854, 10.051, 15.6, 21.296, 45.199, 57.59},
      { 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0.,
	0., 0.001, 0.001, 0.001, 0.002, 0.004,
	0.005, 0.006, 0.009, 0.013, 0.017, 0.033,
	0.040, 0.056, 0.080, 0.120, 0.160, 0.323,
	0.405, 0.575, 0.832, 1.274, 1.723, 3.59,
	4.551, 6.528, 9.562, 14.821, 20.213, 42.793, 54.48},
      { 0., 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0., 0.,
	0., 0.001, 0.001, 0.001, 0.002, 0.003,
	0.004, 0.006, 0.008, 0.012, 0.016, 0.031,
	0.038, 0.054, 0.076, 0.114, 0.152, 0.307, 
        0.386, 0.547, 0.791, 1.211, 1.637, 3.406,
	4.316, 6.184, 9.05, 14.009, 19.089, 40.329, 51.31}
    };

    G4cout << "###  Energy (MeV) ### n= 43" << G4endl;  
    G4int i, j;
    for(i=0; i<43; i++) {G4cout << energy[i] << " ";}
    G4cout << G4endl;  
    G4cout << G4endl;  

    for(i=0; i<3; i++) {
      const G4Material* mat = mman->FindOrBuildMaterial(nm[i]);
      G4double fact = gram/(MeV*cm2*mat->GetDensity());
      G4cout << "###  Material ### " << mat->GetName() << " Data" << G4endl;  
      for(j=0; j<43; j++) {
        dedx[i][j]   -= dedxn[i][j];
        dedxn[i][j]  = fact*(cal.ComputeDEDX(energy[j],part,"muIoni",mat));
        dedxn[i][j] += fact*(cal.ComputeDEDX(energy[j],part,"muBrems",mat));
        dedxn[i][j] += fact*(cal.ComputeDEDX(energy[j],part,"muPairProd",mat));
	G4cout << dedx[i][j] << " ";
      }
      G4cout << G4endl;  
      G4cout << "### Geant4 " << G4endl;  
      for(j=0; j<43; j++) {G4cout << dedxn[i][j] << " ";}
      G4cout << G4endl;  
      G4cout << "### 1 - Geant4/Data " << G4endl;  
      for(j=0; j<43; j++) {G4cout << 1.0 - dedxn[i][j]/dedx[i][j] << " ";}
      G4cout << G4endl;  
    }
    G4cout << G4endl;  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::ElectronTest()
{
  G4NistManager*  mman = G4NistManager::Instance();
  G4EmCalculator cal;
  cal.SetVerbose(0);

  //const G4ParticleDefinition* part = cal.FindParticle("e-");
  const G4ParticleDefinition* part = cal.FindParticle("e+");

  G4bool imu = true;
  if (imu) {
    G4cout << "====================================================================" << G4endl;
    G4cout << "             Stopping Powers" << G4endl;
    G4cout << "====================================================================" << G4endl;

    G4String nm[3] = {"G4_WATER", "G4_Si", "G4_W"};
    const G4int nmax = 30;
    
    G4double e = MeV;    

    G4cout << "###  Electron test for Energy (MeV) = " << e << G4endl;  
    G4int i, j;
    G4double cs;
    G4double cut[nmax];

    for(j=0; j<nmax; j++) {
      cut[j] = std::pow(10.0, -3.0 + 0.1*G4double(j));
      G4cout << cut[j] << " ";
    }
    G4cout << G4endl;  
    for(i=0; i<3; i++) {
      const G4Material* mat = mman->FindOrBuildMaterial(nm[i]);
      //      G4double fact = gram/(barn*cm3*mat->GetDensity());
      G4double fact = 1.0;
      G4cout << "###  Material ### " << mat->GetName() << "   Cross Sections: " <<G4endl;  
      for(j=0; j<30; j++) {
        cs  = fact*(cal.ComputeCrossSectionPerVolume(e,part,"eIoni",mat,cut[j]));
	G4cout << cs << " ";
      }
      G4cout << G4endl;  
    }
    G4cout << G4endl;  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31Histo::CountProcess(const G4String& name)
{
  //  G4cout << "### " << name << G4endl;
  if(name == "AnnihiToMuPair") n_mumu++;  
  else if(name == "ee2hadr") n_pipi++;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
