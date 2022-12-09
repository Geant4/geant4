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
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#include "G4FermiFragmentsPoolVI.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4DeexPrecoParameters.hh"
#include <iomanip>

G4FermiFragmentsPoolVI::G4FermiFragmentsPoolVI()
{
  //  G4cout << "### G4FermiFragmentsPoolVI is constructed" << G4endl;
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  tolerance = param->GetMinExcitation();
  timelim   = (G4float)param->GetMaxLifeTime();

  elim = param->GetFBUEnergyLimit();
  elimf= (G4float)elim;
  /*
  G4cout << "G4FermiFragmentsPoolVI: tolerance= " << tolerance
         << " timelim= " << timelim << " elim= " << elim << G4endl;
  */
  fragment_pool.reserve(991);
  Initialise();
}

G4FermiFragmentsPoolVI::~G4FermiFragmentsPoolVI()
{
  for(G4int i=0; i<maxA; ++i) {
    for(auto & ptr : list_p[i]) { delete ptr; ptr = nullptr; }
    for(auto & ptr : list_c[i]) { delete ptr; ptr = nullptr; }
  }
  for(auto & ptr : fragment_pool) { delete ptr; ptr = nullptr; }
}

const G4FermiChannels* 
G4FermiFragmentsPoolVI::ClosestChannels(G4int Z, G4int A, G4double e) const
{
  const G4FermiChannels* res = nullptr;
  G4double demax = 1.e+9;

  // stable channels
  for(std::size_t j=0; j<(list_c[A]).size(); ++j) {
    const G4FermiFragment* frag = (list_f[A])[j];
    if(frag->GetZ() != Z) { continue; }
    G4double de = e - frag->GetTotalEnergy();
    //G4cout << "  Stab check " << j << " channel de= " << de 
    // << " tol= " << tolerance << G4endl;
    // an excitation coincide with a level
    if(std::abs(de) <= tolerance) { 
      res = (list_c[A])[j];
      break; 
    } else {
      // closest level selected
      de += tolerance;
      if(de >= 0.0 && de <= demax) {
        res = (list_c[A])[j];
        demax = de;
      }
      //G4cout << "  Stab chan: " << j << " N= " 
      //<< res->GetNumberOfChannels() << G4endl;
    }
  }
  return res;
}

G4bool G4FermiFragmentsPoolVI::IsPhysical(G4int Z, G4int A) const
{
  for(auto const& ptr : list_f[A]) {
    if(ptr->GetZ() == Z) { return true; }
  }
  return false;
}

G4bool G4FermiFragmentsPoolVI::IsInThePool(G4int Z, G4int A,  
                                           G4double exc) const
{
  for(auto const& fr : fragment_pool) {
    if(fr->GetZ() == Z && fr->GetA() == A &&  
       std::abs(exc - fr->GetExcitationEnergy()) < tolerance) 
      { return true; }
  }
  return false;
}

G4bool 
G4FermiFragmentsPoolVI::HasChannels(G4int Z, G4int A, G4double exc) const
{
  // stable fragment
  for(std::size_t j=0; j<(list_f[A]).size(); ++j) {
    const G4FermiFragment* frag = (list_f[A])[j];
    if(frag->GetZ() == Z) {
      if(exc > frag->GetExcitationEnergy() && 
         (list_c[A])[j]->GetNumberOfChannels() > 0) { return true; } 
    }
  }
  return false;
}

G4bool G4FermiFragmentsPoolVI::IsInPhysPairs(
       const G4FermiFragment* f1, const G4FermiFragment* f2) const
{
  const G4int A = f1->GetA() + f2->GetA();
  for(auto const& ptr : list_p[A]) {
    if(f1 == ptr->GetFragment1() && f2 == ptr->GetFragment2()) {
      return true;
    }
  }  
  return false;
}

void G4FermiFragmentsPoolVI::Initialise()
{
  //G4cout << "G4FermiFragmentsPoolVI::Initialise main loop @@@@@@" << G4endl;

  // stable particles
  fragment_pool.push_back(new G4FermiFragment(1, 0, 1, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(1, 1, 1, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(2, 1, 2, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(3, 1, 1, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(3, 2, 1, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(4, 2, 0, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(5, 2, 3, 0.0)); 
  fragment_pool.push_back(new G4FermiFragment(5, 3, 3, 0.0));

  // use level data and construct the pool
  G4NuclearLevelData* ndata = G4NuclearLevelData::GetInstance();
  for(G4int Z=1; Z<maxZ; ++Z) {
    G4int Amin = ndata->GetMinA(Z);
    G4int Amax = std::min(maxA, ndata->GetMaxA(Z)+1);
    for(G4int A=Amin; A<Amax; ++A) {
      const G4LevelManager* man = ndata->GetLevelManager(Z, A);
      if(man) {
        std::size_t nn = man->NumberOfTransitions();
        // very unstable state
        if(ndata->MaxLevelEnergy(Z, A) == 0.0f && man->LifeTime(0) == 0.0f) { 
          continue;
        }
        for(std::size_t i=0; i<=nn; ++i) {
          G4float exc = man->LevelEnergy(i);
          /*
          G4cout << "Z= " << Z << " A= " << A << " Eex= " << exc 
                 << " elimf= " << elimf << " toler= " << tolerance 
                 << " time= " << man->LifeTime(i) << " i= " << i << G4endl;
          */ 
          // only levels below limit are consided 
          if(exc >= elimf) { continue; }
          G4double excd = (G4double)exc;
          // only new are considered
          if(IsInThePool(Z, A, excd)) { continue; }
          fragment_pool.push_back(new G4FermiFragment(A,Z,man->SpinTwo(i),excd));
        }
      }
    }
  }
  G4int nfrag = (G4int)fragment_pool.size();
  // prepare structures per A for normal fragments
  const std::size_t lfmax[maxA] = {
    0, 2, 1, 2, 1, 2, 8, 19, 28, 56, 70, 104, 74, 109, 143, 212, 160};
  for(G4int A=1; A<maxA; ++A) {
    list_f[A].reserve(lfmax[A]);
    list_c[A].reserve(lfmax[A]);
  }
  const std::size_t lfch[maxA] = {
    0, 0, 0, 0, 0, 1, 4, 8, 6, 13, 27, 40, 29, 21, 31, 32, 30};

  for(auto const& f : fragment_pool) {
    G4int A = f->GetA();
    G4double exc = f->GetExcitationEnergy();
    list_f[A].push_back(f);
    list_c[A].push_back(new G4FermiChannels(lfch[A], exc, f->GetTotalEnergy()));
  }
  /*  
  G4cout << "Defined fragments @@@@@@" 
         << " PhysicalFrag= " << nfrag 
         << " UnphysicalFrag= " << funstable.size() << G4endl;
  */
  // list of fragment pairs ordered by A
  for(G4int i=0; i<nfrag; ++i) {
    const G4FermiFragment* f1 = fragment_pool[i];
    G4int Z1 = f1->GetZ();
    G4int A1 = f1->GetA();
    G4double e1 = f1->GetTotalEnergy();
    for(G4int j=0; j<nfrag; ++j) {
      const G4FermiFragment* f2 = fragment_pool[j];
      G4int Z2 = f2->GetZ();
      G4int A2 = f2->GetA();
      if(A2 < A1 || (A2 == A1 && Z2 < Z1)) { continue; }
      G4int Z = Z1 + Z2;
      G4int A = A1 + A2;

      if(Z >= maxZ || A >= maxA || IsInPhysPairs(f1, f2)) { continue; } 

      G4double e2 = f2->GetTotalEnergy();
      G4double minE = e1 + e2;
      G4double exc = 0.0;
      if(IsPhysical(Z, A)) {
        minE += f1->GetCoulombBarrier(A2, Z2, 0.0);
        exc = minE - G4NucleiProperties::GetNuclearMass(A, Z);
      }
      /*
        G4cout << "Z= " << Z << " A= " << A 
             << " Z1= " << Z1 << " A1= " << A1
             << " Z2= " << Z2 << " A2= " << A2 << " Eex= " << exc 
             << "  Qb= " << f1->GetCoulombBarrier(A2, Z2, 0.0) 
             << "  " << e1
             << "  " << e2
             << "  " << G4NucleiProperties::GetNuclearMass(A, Z)
             << G4endl; 
      */
      // ignore very excited case
      if(exc >= elim) { continue; }
      G4FermiPair* fpair = nullptr;
      std::size_t kmax = list_f[A].size();
      for(std::size_t k=0; k<kmax; ++k) {
        const G4FermiFragment* f3 = (list_f[A])[k];
        if(Z == f3->GetZ() &&  
           f3->GetTotalEnergy() - minE + tolerance >= 0.0) {
          if(!fpair) {
            fpair = new G4FermiPair(f1, f2);
            list_p[A].push_back(fpair);
          }
          (list_c[A])[k]->AddChannel(fpair); 
        }
      }
    }
  }
  // compute static probabilities
  for(G4int A=1; A<maxA; ++A) {
    for(std::size_t j=0; j<list_c[A].size(); ++j) {
      G4FermiChannels* ch = (list_c[A])[j];
      const G4FermiFragment* frag = (list_f[A])[j];
      std::size_t nch = ch->GetNumberOfChannels();
      if(1 < nch) {
        std::vector<G4double>& prob = ch->GetProbabilities();
        const std::vector<const G4FermiPair*>& pairs = ch->GetChannels();
        G4double ptot = 0.0;
        for(std::size_t i=0; i<nch; ++i) {
          ptot += theDecay.ComputeProbability(frag->GetZ(), frag->GetA(),
                                              frag->GetSpin(), 
                                              frag->GetTotalEnergy(), 
                                              pairs[i]->GetFragment1(),
                                              pairs[i]->GetFragment2());
          prob[i] = ptot;
        }
        if(0.0 == ptot) {
          prob[0] = 1.0;
        } else {
          ptot = 1./ptot;
          for(std::size_t i=0; i<nch-1; ++i) { prob[i] *= ptot; }
          prob[nch-1] = 1.0;
        }
      }
    }
  }
}

void G4FermiFragmentsPoolVI::DumpFragment(const G4FermiFragment* f) const
{
  if(f) {
    G4long prec = G4cout.precision(6);
    G4cout << "   Z= " << f->GetZ() << " A= " << std::setw(2) << f->GetA() 
           << " Mass(GeV)= " << std::setw(8) << f->GetFragmentMass()/GeV
           << " Eexc(MeV)= " << std::setw(7) << f->GetExcitationEnergy()
           << " 2s= " << f->GetSpin() << " IsStable: " 
           << HasChannels(f->GetZ(), f->GetA(), f->GetExcitationEnergy())
           << G4endl;
    G4cout.precision(prec);
  }
}

void G4FermiFragmentsPoolVI::Dump() const
{
  G4cout <<"----------------------------------------------------------------"
         <<G4endl;
  G4cout << "##### List of Fragments in the Fermi Fragment Pool #####" 
         << G4endl;
  std::size_t nfrag = fragment_pool.size();
  G4cout << "      For stable " << nfrag << " Elim(MeV) = " 
         << elim/CLHEP::MeV << G4endl; 
  for(std::size_t i=0; i<nfrag; ++i) {
    DumpFragment(fragment_pool[i]);
  }
  G4cout << G4endl;


  G4cout << "----------------------------------------------------------------"
         << G4endl;
  G4cout << "### G4FermiFragmentPoolVI: fragments sorted by A" << G4endl; 

  G4long prec = G4cout.precision(6);
  std::size_t ama[maxA];
  ama[0] =  0;
  for(G4int A=1; A<maxA; ++A) {
    G4cout << " # A= " << A << G4endl;
    std::size_t am(0); 
    for(std::size_t j=0; j<list_f[A].size(); ++j) {
      const G4FermiFragment* f = (list_f[A])[j];
      G4int a1 = f->GetA();
      G4int z1 = f->GetZ();
      std::size_t nch = (list_c[A])[j]->GetNumberOfChannels();
      am = std::max(am, nch);
      G4cout << "   ("<<a1<<","<<z1<<");  Eex(MeV)= "
             << f->GetExcitationEnergy() 
             << " 2S= " << f->GetSpin()
             << "; Nchannels= " << nch
             << " MassExcess= " << f->GetTotalEnergy() - 
        (z1*proton_mass_c2 + (a1 - z1)*neutron_mass_c2)
             << G4endl; 
      for(std::size_t k=0; k<nch; ++k) {
        const G4FermiPair* fpair = ((list_c[A])[j]->GetChannels())[k];
        G4cout << "         (" << fpair->GetFragment1()->GetZ()
               << ", "  << fpair->GetFragment1()->GetA() 
               << ",  " << fpair->GetFragment1()->GetExcitationEnergy()
               << ")  ("<< fpair->GetFragment2()->GetZ()
               << ", "  << std::setw(3)<< fpair->GetFragment2()->GetA()  
               << ",  " << std::setw(8)<< fpair->GetFragment2()->GetExcitationEnergy()
               << ")  prob= " << ((list_c[A])[j]->GetProbabilities())[k]
               << G4endl; 
      }
    }
    ama[A] = am; 
  }
  G4cout.precision(prec);
  G4cout << G4endl;

  G4cout << "    Number of fragments per A:" << G4endl;
  for(G4int j=0; j<maxA; ++j) { G4cout << list_f[j].size() << ", "; }
  G4cout << G4endl;

  G4cout << "    Max number of channels per A:" << G4endl; 
  for (std::size_t j=0; j<maxA; ++j) { G4cout << ama[j] << ", "; }
  G4cout << G4endl;

  G4cout << "    Number of fragment pairs per A:" << G4endl;
  for(G4int j=0; j<maxA; ++j) { G4cout << list_p[j].size() << ", "; }
  G4cout << G4endl;

  G4cout << "----------------------------------------------------------------"
         << G4endl;
  G4cout << "### Pairs of stable fragments: " << G4endl;

  prec = G4cout.precision(6);
  for(G4int A=2; A<maxA; ++A) {
    G4cout << "  A= " << A<<G4endl; 
    for(std::size_t j=0; j<list_p[A].size(); ++j) {
      const G4FermiFragment* f1 = (list_p[A])[j]->GetFragment1();
      const G4FermiFragment* f2 = (list_p[A])[j]->GetFragment2();
      G4int a1 = f1->GetA();
      G4int z1 = f1->GetZ();
      G4int a2 = f2->GetA();
      G4int z2 = f2->GetZ();
      G4cout << "("<<a1<<","<<z1<<")("<<a2<<","<<z2<<") % Eex(MeV)= "
             << std::setw(8)<< (list_p[A])[j]->GetExcitationEnergy() 
             << " Eex1= " << std::setw(8)<< f1->GetExcitationEnergy()
             << " Eex2= " << std::setw(8)<< f2->GetExcitationEnergy()
             << G4endl; 
    }
    G4cout << G4endl;
    G4cout <<"----------------------------------------------------------------"
           << G4endl;
  }
  G4cout.precision(prec);
}

