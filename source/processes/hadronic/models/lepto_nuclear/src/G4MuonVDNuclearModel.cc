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
//
// Author:      D.H. Wright
// Date:        2 February 2011
//
// Description: model of muon nuclear interaction in which a gamma from
//              the virtual photon spectrum interacts in the nucleus as
//              a real gamma at low energies and as a pi0 at high energies.
//              Kokoulin's muon cross section and equivalent gamma spectrum
//              are used.
//

#include "G4MuonVDNuclearModel.hh"

#include "Randomize.hh"
#include "G4Log.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4CascadeInterface.hh"
#include "G4TheoFSGenerator.hh" 
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4PreCompoundModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4FTFModel.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4KokoulinMuonNuclearXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ElementData.hh" 
#include "G4Physics2DVector.hh" 
#include "G4Pow.hh" 
#include "G4PhysicsModelCatalog.hh"

const G4int G4MuonVDNuclearModel::zdat[] = {1, 4, 13, 29, 92};
const G4double G4MuonVDNuclearModel::adat[] = {1.01,9.01,26.98,63.55,238.03};
const G4double G4MuonVDNuclearModel::tdat[] = {
  1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, 
  1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, 
  1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, 
  1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, 
  1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7, 
  1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8, 
  1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9, 
  1.e10,2.e10,3.e10,4.e10,5.e10,6.e10,7.e10,8.e10,9.e10,1.e11}; 

G4ElementData* G4MuonVDNuclearModel::fElementData = nullptr;             

G4MuonVDNuclearModel::G4MuonVDNuclearModel()
  : G4HadronicInteraction("G4MuonVDNuclearModel"),isMaster(false)
{
  muNucXS = (G4KokoulinMuonNuclearXS*)G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4KokoulinMuonNuclearXS::Default_Name());

  SetMinEnergy(0.0);
  SetMaxEnergy(1*CLHEP::PeV);
  CutFixed = 0.2*CLHEP::GeV;

  if(!fElementData && G4Threading::IsMasterThread()) { 
    fElementData = new G4ElementData();
    MakeSamplingTable(); 
    isMaster = true;
  }
                    
  // reuse existing pre-compound model
  G4GeneratorPrecompoundInterface* precoInterface 
    = new G4GeneratorPrecompoundInterface();
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  G4VPreCompoundModel* pre = static_cast<G4VPreCompoundModel*>(p);
  if(!pre) { pre = new G4PreCompoundModel(); }
  precoInterface->SetDeExcitation(pre);

  // Build FTFP model
  ftfp = new G4TheoFSGenerator();
  ftfp->SetTransport(precoInterface);
  theFragmentation = new G4LundStringFragmentation();
  theStringDecay = new G4ExcitedStringDecay(theFragmentation);    
  G4FTFModel* theStringModel = new G4FTFModel;
  theStringModel->SetFragmentationModel(theStringDecay);
  ftfp->SetHighEnergyGenerator(theStringModel);

  // Build Bertini cascade
  bert = new G4CascadeInterface();

  // Creator model ID
  secID = G4PhysicsModelCatalog::GetModelID( "model_" + GetModelName() );
}

G4MuonVDNuclearModel::~G4MuonVDNuclearModel()
{
  delete theFragmentation;
  delete theStringDecay;

  if(isMaster) { 
    delete fElementData;
    fElementData = nullptr;
  } 
}
  
G4HadFinalState*
G4MuonVDNuclearModel::ApplyYourself(const G4HadProjectile& aTrack,
                                    G4Nucleus& targetNucleus)
{
  theParticleChange.Clear();

  // For very low energy, return initial track
  G4double epmax = aTrack.GetTotalEnergy() - 0.5*proton_mass_c2;
  if (epmax <= CutFixed) {
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }

  // Produce recoil muon and transferred photon
  G4DynamicParticle* transferredPhoton = CalculateEMVertex(aTrack, targetNucleus);

  // Interact the gamma with the nucleus
  CalculateHadronicVertex(transferredPhoton, targetNucleus);
  return &theParticleChange;
}

G4DynamicParticle*
G4MuonVDNuclearModel::CalculateEMVertex(const G4HadProjectile& aTrack,
                                        G4Nucleus& targetNucleus)
{
  // Select sampling table
  G4double KineticEnergy = aTrack.GetKineticEnergy();
  G4double TotalEnergy = aTrack.GetTotalEnergy();
  G4double Mass = G4MuonMinus::MuonMinus()->GetPDGMass();
  G4Pow* g4calc = G4Pow::GetInstance();
  G4double lnZ = g4calc->logZ(targetNucleus.GetZ_asInt());

  G4double epmin = CutFixed;
  G4double epmax = TotalEnergy - 0.5*proton_mass_c2;
  G4double m0 = CutFixed;

  G4double delmin = 1.e10;
  G4double del;
  G4int izz = 0;
  G4int itt = 0;

  for (G4int iz = 0; iz < nzdat; ++iz) {
    del = std::abs(lnZ - g4calc->logZ(zdat[iz]));
    if (del < delmin) {
      delmin = del;
      izz = iz;
    }
  }
 
  delmin = 1.e10;
  for (G4int it = 0; it < ntdat; ++it) {
    del = std::abs(G4Log(KineticEnergy)-G4Log(tdat[it]) );
    if (del < delmin) {
      delmin = del;
      itt = it;
    }
  }

  // Sample the energy transfer according to the probability table
  G4double r = G4UniformRand();

  G4int iy;

  G4int Z = zdat[izz];

  for(iy = 0; iy<NBIN; ++iy)  {

    G4double pvv = fElementData->GetElement2DData(Z)->GetValue(iy, itt); 
    if(pvv >= r) { break; }
  }       

  // Sampling is done uniformly in y in the bin
  G4double pvx = fElementData->GetElement2DData(Z)->GetX(iy); 
  G4double pvx1 = fElementData->GetElement2DData(Z)->GetX(iy+1); 
  G4double y = pvx + G4UniformRand() * (pvx1 - pvx);

  G4double x = G4Exp(y);
  G4double ep = epmin*G4Exp(x*G4Log(epmax/epmin) );

  // Sample scattering angle of mu, but first t should be sampled.
  G4double yy = ep/TotalEnergy;
  G4double tmin = Mass*Mass*yy*yy/(1.-yy);
  G4double tmax = 2.*proton_mass_c2*ep;
  G4double t1;
  G4double t2;
  if (m0 < ep) {
    t1 = m0*m0;
    t2 = ep*ep;
  } else {
    t1 = ep*ep;
    t2 = m0*m0;
  }

  G4double w1 = tmax*t1;
  G4double w2 = tmax+t1;
  G4double w3 = tmax*(tmin+t1)/(tmin*w2);
  G4double y1 = 1.-yy;
  G4double y2 = 0.5*yy*yy;
  G4double y3 = y1+y2;

  G4double t = 0.0;
  G4double rej = 0.0;

  // Now sample t
  G4int ntry = 0; 
  do
  {
    ntry += 1;
    if (ntry > 10000) {
      G4ExceptionDescription eda;
      eda << " While count exceeded " << G4endl;
      G4Exception("G4MuonVDNuclearModel::CalculateEMVertex()", "HAD_RPG_100", JustWarning, eda);
      break;
    }
  
    t = w1/(w2*G4Exp(G4UniformRand()*G4Log(w3))-tmax);
    rej = (1.-t/tmax)*(y1*(1.-tmin/t)+y2)/(y3*(1.-t/t2)); 
  } while (G4UniformRand() > rej) ;   /* Loop checking, 01.09.2015, D.Wright */

  // compute angle from t
  G4double sinth2 =
             0.5*(t-tmin)/(2.*(TotalEnergy*(TotalEnergy-ep)-Mass*Mass)-tmin);
  G4double theta = std::acos(1. - 2.*sinth2);

  G4double phi = twopi*G4UniformRand();
  G4double sinth = std::sin(theta);
  G4double dirx = sinth*std::cos(phi);
  G4double diry = sinth*std::sin(phi);
  G4double dirz = std::cos(theta);
  G4ThreeVector finalDirection(dirx,diry,dirz);
  G4ThreeVector ParticleDirection(aTrack.Get4Momentum().vect().unit() );
  finalDirection.rotateUz(ParticleDirection);

  G4double NewKinEnergy = KineticEnergy - ep;
  G4double finalMomentum = std::sqrt(NewKinEnergy*(NewKinEnergy+2.*Mass) );
  G4double Ef = NewKinEnergy + Mass;
  G4double initMomentum = std::sqrt(KineticEnergy*(TotalEnergy+Mass) );

  // Set energy and direction of scattered primary in theParticleChange
  theParticleChange.SetStatusChange(isAlive);
  theParticleChange.SetEnergyChange(NewKinEnergy);
  theParticleChange.SetMomentumChange(finalDirection);

  // Now create the emitted gamma 
  G4LorentzVector primaryMomentum(initMomentum*ParticleDirection, TotalEnergy);
  G4LorentzVector fsMomentum(finalMomentum*finalDirection, Ef);
  G4LorentzVector momentumTransfer = primaryMomentum - fsMomentum;

  G4DynamicParticle* gamma = 
           new G4DynamicParticle(G4Gamma::Gamma(), momentumTransfer);
 
  return gamma;
}

void
G4MuonVDNuclearModel::CalculateHadronicVertex(G4DynamicParticle* incident,
                                              G4Nucleus& target)
{
  G4HadFinalState* hfs = 0;
  G4double gammaE = incident->GetTotalEnergy();

  if (gammaE < 10*GeV) {
    G4HadProjectile projectile(*incident);
    hfs = bert->ApplyYourself(projectile, target);
  } else {
    // convert incident gamma to a pi0
    G4double piMass = G4PionZero::PionZero()->GetPDGMass();
    G4double piKE = incident->GetTotalEnergy() - piMass;
    G4double piMom = std::sqrt(piKE*(piKE + 2*piMass) );
    G4ThreeVector piMomentum(incident->GetMomentumDirection() );
    piMomentum *= piMom;
    G4DynamicParticle theHadron(G4PionZero::PionZero(), piMomentum);
    G4HadProjectile projectile(theHadron);
    hfs = ftfp->ApplyYourself(projectile, target);
  }

  delete incident;

  // Assign the creator model ID to the secondaries
  for ( size_t i = 0; i < hfs->GetNumberOfSecondaries(); ++i ) {
    hfs->GetSecondary( i )->SetCreatorModelID( secID );
  }
  
  // Copy secondaries from sub-model to model
  theParticleChange.AddSecondaries(hfs);
} 


void G4MuonVDNuclearModel::MakeSamplingTable()
{
  G4int nbin;
  G4double KineticEnergy;
  G4double TotalEnergy;
  G4double Maxep;
  G4double CrossSection;

  G4double c;
  G4double y;
  G4double ymin,ymax;
  G4double dy,yy;
  G4double dx,x;
  G4double ep;

  G4double AtomicNumber;
  G4double AtomicWeight;
 
  G4double mumass = G4MuonMinus::MuonMinus()->GetPDGMass();

  for (G4int iz = 0; iz < nzdat; ++iz) {
    AtomicNumber = zdat[iz];
    AtomicWeight = adat[iz]*(g/mole);

    G4Physics2DVector* pv = new G4Physics2DVector(NBIN+1,ntdat+1); 
    G4double pvv; 

    for (G4int it = 0; it < ntdat; ++it) {
      KineticEnergy = tdat[it];
      TotalEnergy = KineticEnergy + mumass;
      Maxep = TotalEnergy - 0.5*proton_mass_c2;

      CrossSection = 0.0;

      // Calculate the differential cross section
      // numerical integration in log .........
      c = G4Log(Maxep/CutFixed);
      ymin = -5.0;
      ymax = 0.0;
      dy = (ymax-ymin)/NBIN; 

      nbin=-1;              

      y = ymin - 0.5*dy;
      yy = ymin - dy;
      for (G4int i = 0; i < NBIN; ++i) {
        y += dy;
        x = G4Exp(y);
        yy += dy;
        dx = G4Exp(yy+dy)-G4Exp(yy);
      
        ep = CutFixed*G4Exp(c*x);

        CrossSection +=
           ep*dx*muNucXS->ComputeDDMicroscopicCrossSection(KineticEnergy,
							   AtomicNumber,
							   AtomicWeight, ep);
        if (nbin < NBIN) {
          ++nbin;
          pv->PutValue(nbin, it, CrossSection); 
          pv->PutX(nbin, y); 
        }
      }
      pv->PutX(NBIN, 0.); 
 
      if (CrossSection > 0.0) {
        for (G4int ib = 0; ib <= nbin; ++ib) { 
          pvv = pv->GetValue(ib, it); 
          pvv = pvv/CrossSection; 
          pv->PutValue(ib, it, pvv); 
        } 
      }
    } // loop on it

    fElementData->InitialiseForElement(zdat[iz], pv);
  } // loop on iz

  // G4cout << " Kokoulin XS = "
  //       << muNucXS->ComputeDDMicroscopicCrossSection(1*GeV, 20.0, 
  //         40.0*g/mole, 0.3*GeV)/millibarn
  //       << G4endl; 
}

void G4MuonVDNuclearModel::ModelDescription(std::ostream& outFile) const 
{
  outFile << "G4MuonVDNuclearModel handles the inelastic scattering\n"
          << "of mu- and mu+ from nuclei using the equivalent photon\n"
          << "approximation in which the incoming lepton generates a\n"
          << "virtual photon at the electromagnetic vertex, and the\n"
          << "virtual photon is converted to a real photon.  At low\n"
          << "energies, the photon interacts directly with the nucleus\n"
          << "using the Bertini cascade.  At high energies the photon\n"
          << "is converted to a pi0 which interacts using the FTFP\n"
          << "model.  The muon-nuclear cross sections of R. Kokoulin \n"
          << "are used to generate the virtual photon spectrum\n";
}

