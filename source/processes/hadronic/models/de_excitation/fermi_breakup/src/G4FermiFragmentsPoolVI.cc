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
// $Id: G4FermiFragmentsPoolVI.cc,v 1.5 2006-06-29 20:13:13 gunter Exp $
//
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#include "G4FermiFragmentsPoolVI.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4StableFermiFragment.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include <iomanip>

G4FermiFragmentsPoolVI::G4FermiFragmentsPoolVI()
{
  maxZ = 9;
  maxA = 17;
  tolerance = 0.1*CLHEP::keV;
  elim = 10*CLHEP::MeV;
  elim_unstable = 10*CLHEP::MeV;
  Initialise();
}

G4FermiFragmentsPoolVI::~G4FermiFragmentsPoolVI()
{
  size_t nn;
  for(G4int i=0; i<maxA; ++i) {
    nn = list_p[i].size();
    for(size_t j=0; j<nn; ++j) { delete (list_p[i])[j]; }
    nn = list_c[i].size();
    for(size_t j=0; j<nn; ++j) { delete (list_c[i])[j]; }
    nn = list_d[i].size();
    for(size_t j=0; j<nn; ++j) { delete (list_d[i])[j]; }
    nn = list_u[i].size();
    for(size_t j=0; j<nn; ++j) { delete (list_u[i])[j]; }
  }
  nn = fragment_pool.size();
  for(size_t j=0; j<nn; ++j) { delete fragment_pool[j]; }
  nn = funstable.size();
  for(size_t j=0; j<nn; ++j) { delete funstable[j]; }
}

G4bool 
G4FermiFragmentsPoolVI::IsApplicable(G4int Z, G4int A, G4double etot) const
{
  G4bool isInList = false;
  size_t nn = list_f[A].size();
  for(size_t i=0; i<nn; ++i) {
    if(Z == (list_f[A])[i]->GetZ()) { 
      isInList = true;
      if(etot <= (list_f[A])[i]->GetFragmentMass() + elim) { return true; }
    }
  }
  if(isInList) { return false; }
  nn = list_g[A].size();
  for(size_t i=0; i<nn; ++i) {
    if(Z == (list_g[A])[i]->GetZ() &&
       etot <= (list_g[A])[i]->GetFragmentMass() + elim) { return true; }
  }
  return false;
}

const G4FermiChannels* 
G4FermiFragmentsPoolVI::ClosestChannels(G4int Z, G4int A, G4double e) const
{
  const G4FermiChannels* res = nullptr; 
  G4double demax = e;

  // stable channels;
  for(size_t j=0; j<list_c[A].size(); ++j) {
    const G4FermiFragment* frag = (list_f[A])[j];
    if(frag->GetZ() != Z) { continue; }
    G4double de = e - frag->GetTotalEnergy();
    if(std::abs(de) <= tolerance) { 
      res = (list_c[A])[j];
      break;
    } else if(de > 0.0 && de < demax) {
      res = (list_c[A])[j];
      demax = de;
    } else if(de > 0.0 && de >= demax) {
      break;
    }
  }
  // stable + unstable channels
  if(!res) {
    for(size_t j=0; j<list_d[A].size(); ++j) {
      const G4FermiFragment* frag = (list_g[A])[j];
      if(frag->GetZ() != Z) { continue; }
      G4double de = e - frag->GetTotalEnergy();
      if(std::abs(de) <= tolerance || de > 0.0) { 
	res = (list_d[A])[j];
	break;
      }
    }  
  }
  return res;
}

G4bool G4FermiFragmentsPoolVI::IsInThePool(G4int Z, G4int A,  
					   G4double exc) const
{
  G4bool res = false;
  G4int nfrag = fragment_pool.size();
  for(G4int i=0; i<nfrag; ++i) {
    const G4FermiFragment* fr = fragment_pool[i];
    if(fr->GetZ() == Z && fr->GetA() == A && 
       std::abs(exc - fr->GetExcitationEnergy()) < tolerance) {
      res = true;
      break;
    }
  }
  return res;
}

G4bool G4FermiFragmentsPoolVI::IsInPhysPairs(const G4FermiFragment* f1, 
					     const G4FermiFragment* f2,
					     G4double exc) const
{
  G4bool res = false;
  G4int A1 = f1->GetA();
  G4int A2 = f2->GetA();
  G4int A = A1 + A2;
  G4int nn = list_p[A].size();
  for(G4int i=0; i<nn; ++i) {
    if(f1 == (list_p[A])[i]->GetFragment1() && f2 == (list_p[A])[i]->GetFragment2()
       && std::abs((list_p[A])[i]->GetExcitationEnergy() - exc) < tolerance) {
      res = true;
      break;
    }
  }  
  return res;
}

G4bool G4FermiFragmentsPoolVI::IsInUnphysPairs(const G4FermiFragment* f1, 
					       const G4FermiFragment* f2,
					       G4double exc) const
{
  G4bool res = false;
  G4int A1 = f1->GetA();
  G4int A2 = f2->GetA();
  G4int A = A1 + A2;
  G4int nn = list_u[A].size();
  for(G4int i=0; i<nn; ++i) {
    if(f1 == (list_u[A])[i]->GetFragment1() && f2 == (list_u[A])[i]->GetFragment2()
       && std::abs((list_u[A])[i]->GetExcitationEnergy() - exc) < tolerance) {
      res = true;
      break;
    }
  }  
  return res;
}

void G4FermiFragmentsPoolVI::Initialise()
{
  static const G4int nmin = 8;

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
  for(G4int Z=3; Z<maxZ; ++Z) {
    G4int Amin = ndata->GetMinA(Z);
    G4int Amax = std::min(maxA, ndata->GetMaxA(Z)+1);
    for(G4int A=Amin; A<Amax; ++A) {
      const G4LevelManager* man = ndata->GetLevelManager(Z, A);
      if(man) {
	size_t nn = man->NumberOfTransitions();
        for(size_t i=0; i<=nn; ++i) {
          G4double exc = (G4double)man->LevelEnergy(i);
	  // only levels below limit are consided 
          if(exc < elim) {
	    fragment_pool.push_back(new G4FermiFragment(A, Z, 
							std::abs(man->SpinParity(i)), exc));
	  }
	}
      }
    }
  }
  // prepare structures per A
  G4int nfrag = fragment_pool.size();
  for(G4int i=0; i<nfrag; ++i) {
    const G4FermiFragment* f = fragment_pool[i];
    G4int A = f->GetA();
    G4double exc = f->GetExcitationEnergy();
    list_f[A].push_back(f);
    list_c[A].push_back(new G4FermiChannels(6, exc, f->GetTotalEnergy()));
  }
  // list of unphysical fragments
  for(G4int Z=0; Z<maxZ; ++Z) {
    G4int A0 = std::max(Z, 1);
    for(G4int A=A0; A<maxA; ++A) {
      if(IsInThePool(Z, A, 0.0)) { continue; }
      const G4FermiFragment* f = new G4FermiFragment(A, Z, -1, 0.0, false);
      funstable.push_back(f);
      list_g[A].push_back(f);
      list_d[A].push_back(new G4FermiChannels(nmin+1,0.0,f->GetTotalEnergy()));
    }
  }
  // G4cout << "@@@@@@ main loop @@@@@@" << G4endl;
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
      if(Z >= maxZ || A >= maxA) { continue; }
      G4double e2 = f2->GetTotalEnergy();
      G4double minE = e1 + e2 + f1->GetCoulombBarrier(A2, Z2, 0.0);
      G4double exc = minE - G4NucleiProperties::GetNuclearMass(A, Z);
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
                  
      G4bool isFound = false;
      G4int kmax = list_f[A].size();
      for(G4int k=0; k<kmax; ++k) {
        if(Z != (list_f[A])[k]->GetZ()) { continue; }
        G4double exc1 = (list_f[A])[k]->GetExcitationEnergy();
        //G4cout << "## Z= " << Z << " A= " << A 
	//       << " Exc1= " << exc1 << " Exc= " << exc << G4endl; 
        if(exc1 + tolerance < exc || IsInPhysPairs(f1, f2, exc1)) { continue; }
        G4FermiPair* fpair = new G4FermiPair(f1, f2, exc1);
	list_p[A].push_back(fpair);
        (list_c[A])[k]->AddChannel(fpair); 
        isFound = true;
      }
      if(isFound) { continue; }
      kmax = list_g[A].size();
      for(G4int k=0; k<kmax; ++k) {
	if(Z != (list_g[A])[k]->GetZ()) { continue; }
        G4double exc1 = (list_g[A])[k]->GetExcitationEnergy();
	if(exc1 + tolerance < exc || IsInUnphysPairs(f1, f2, exc) 
	   || (list_d[A])[k]->GetNumberOfChannels() > nmin) { continue; }
	G4FermiPair* fpair = new G4FermiPair(f1, f2, exc);
	list_u[A].push_back(fpair);
	(list_d[A])[k]->AddChannel(fpair); 
      }
    }
  }

  // G4cout << "@@@@@@ sec loop @@@@@@" << G4endl;
  // list of fragment pairs (stable+unstable) ordered by A
  G4int unphys = funstable.size();
  for(G4int i=0; i<nfrag; ++i) {
    const G4FermiFragment* f1 = fragment_pool[i];
    G4int Z1 = f1->GetZ();
    G4int A1 = f1->GetA();
    G4double e1 = f1->GetTotalEnergy();
    for(G4int j=0; j<unphys; ++j) {
      const G4FermiFragment* f2 = funstable[j];
      G4int Z2 = f2->GetZ();
      G4int A2 = f2->GetA();
      G4int Z = Z1 + Z2;
      G4int A = A1 + A2;
      if(Z >= maxZ || A >= maxA) { continue; }
      G4double e2 = f2->GetTotalEnergy();
      G4double minE = e1 + e2;
      /*
      G4cout << "Z= " << Z << " A= " << A << " Z1= " << Z1 << " A1= " << A1
	     << " Z2= " << Z2 << " A2= " << A2 << G4endl;
      */

      // check if this is the list of stable pairs
      G4bool isFound = false;
      G4int kmax = list_f[A].size();
      for(G4int k=0; k<kmax; ++k) {
        if(Z != (list_f[A])[k]->GetZ() || 
	   (list_c[A])[k]->GetNumberOfChannels() > 0) { continue; } 
	G4double etot = (list_f[A])[k]->GetTotalEnergy();
	if(etot + tolerance >= minE) {
	  G4FermiPair* fpair = 
	    new G4FermiPair(f1,f2,(list_f[A])[k]->GetExcitationEnergy());
	  list_p[A].push_back(fpair);
	  (list_c[A])[k]->AddChannel(fpair); 
	  isFound = true;
	  break; 
	}
      }
      if(isFound) { continue; }
      kmax = list_g[A].size();
      for(G4int k=0; k<kmax; ++k) {
	if(Z != (list_g[A])[k]->GetZ()) { continue; }
	G4double etot = (list_g[A])[k]->GetTotalEnergy();
        G4double exc  = etot - minE;
	//G4cout << "@@ Exc= " << exc << G4endl; 
	if(exc + tolerance < 0.0 || IsInUnphysPairs(f1, f2, exc)
	   || (list_d[A])[k]->GetNumberOfChannels() > nmin) { continue; }
	G4FermiPair* fpair = new G4FermiPair(f1, f2, std::max(exc, 0.0));
	list_u[A].push_back(fpair);
	(list_d[A])[k]->AddChannel(fpair); 
      }
    }
  }

  // compute static probabilities
  for(G4int A=1; A<maxA; ++A) {
    for(size_t j=0; j<list_c[A].size(); ++j) {
      G4FermiChannels* ch = (list_c[A])[j];
      const G4FermiFragment* frag = (list_f[A])[j];
      size_t nch = ch->GetNumberOfChannels();
      if(1 < nch) {
	std::vector<G4double>& prob = ch->GetProbabilities();
	const std::vector<const G4FermiPair*>& pairs = ch->GetChannels();
        G4double ptot = 0.0;
        for(size_t i=0; i<nch; ++i) {
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
	  for(size_t i=0; i<nch-1; ++i) { prob[i] *= ptot; }
	  prob[nch-1] = 1.0;
	}
      }
    }
  }
  for(G4int A=1; A<maxA; ++A) {
    for(size_t j=0; j<list_d[A].size(); ++j) {
      G4FermiChannels* ch = (list_d[A])[j];
      const G4FermiFragment* frag = (list_g[A])[j];
      size_t nch = ch->GetNumberOfChannels();
      if(1 < nch) {
	std::vector<G4double>& prob = ch->GetProbabilities();
	const std::vector<const G4FermiPair*>& pairs = ch->GetChannels();
        G4double ptot = 0.0;
        for(size_t i=0; i<nch; ++i) {
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
	  for(size_t i=0; i<nch-1; ++i) { prob[i] *= ptot; }
	  prob[nch-1] = 1.0;
	}
      }
    }
  }
}

void G4FermiFragmentsPoolVI::DumpFragment(const G4FermiFragment* f) const
{
  if(f) {
    G4int prec = G4cout.precision(6);
    G4cout << "   Z= " << f->GetZ() << " A= " << std::setw(2) << f->GetA() 
	   << " Mass(GeV)= " << std::setw(8) << f->GetFragmentMass()/GeV
	   << " Eexc(MeV)= " << std::setw(7) << f->GetExcitationEnergy()
	   << " 2s= " << f->GetSpin() 
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
  G4cout << "      For stable Elim(MeV) = " << elim/CLHEP::MeV 
	 << "   for unstable Elim(MeV) = " << elim_unstable/CLHEP::MeV << G4endl; 
  G4int nfrag = fragment_pool.size();
  for(G4int i=0; i<nfrag; ++i) {
    DumpFragment(fragment_pool[i]);
  }
  G4cout << G4endl;

  G4int prec = G4cout.precision(6);

  G4cout << "----------------------------------------------------------------"
	 << G4endl;
  G4cout << "### G4FermiFragmentPoolVI: fragments sorted by A" << G4endl; 
  for(G4int A=1; A<maxA; ++A) {
    G4cout << " # A= " << A << G4endl; 
    for(size_t j=0; j<list_f[A].size(); ++j) {
      const G4FermiFragment* f = (list_f[A])[j];
      G4int a1 = f->GetA();
      G4int z1 = f->GetZ();
      size_t nch = (list_c[A])[j]->GetNumberOfChannels();
      G4cout << "   ("<<a1<<","<<z1<<");  Eex(MeV)= "
	     << f->GetExcitationEnergy() 
	     << " 2S= " << f->GetSpin()
	     << "; Nchannels= " << nch
	     << " MassExcess= " << f->GetTotalEnergy() - 
	(z1*proton_mass_c2 + (a1 - z1)*neutron_mass_c2)
	     << G4endl; 
      for(size_t k=0; k<nch; ++k) {
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
  }
  G4cout << G4endl;
  G4cout << "----------------------------------------------------------------"
	 << G4endl;
 
  G4cout << "### G4FermiFragmentPoolVI: " << funstable.size() 
	 << " unphysical fragments" << G4endl;
  for(G4int A=1; A<maxA; ++A) {
    G4cout << " # A= " << A << G4endl; 
    for(size_t j=0; j<list_g[A].size(); ++j) {
      const G4FermiFragment* f = (list_g[A])[j];
      G4int a1 = f->GetA();
      G4int z1 = f->GetZ();
      size_t nch = (list_d[A])[j]->GetNumberOfChannels();
      G4cout << "("<<a1<<","<<z1<<");  Eex(MeV)= "
	     << std::setw(8) << f->GetExcitationEnergy()
	     << "; Nchannels= " << nch
	     << " MassExcess= " << f->GetTotalEnergy() -
	(z1*proton_mass_c2 + (a1 - z1)*neutron_mass_c2)
	     << G4endl; 
      for(size_t k=0; k<nch; ++k) {
        const G4FermiPair* fpair = ((list_d[A])[j]->GetChannels())[k];
        G4cout << "         (" << fpair->GetFragment1()->GetZ()
	       << ", "  << fpair->GetFragment1()->GetA() 
	       << ",  " << std::setw(8)<< fpair->GetFragment1()->GetExcitationEnergy()
	       << ")  ("<< fpair->GetFragment2()->GetZ()
	       << ", "  << std::setw(3)<< fpair->GetFragment2()->GetA()  
	       << ",  " << std::setw(8)<< fpair->GetFragment2()->GetExcitationEnergy()
	       << ")  prob= " << ((list_d[A])[j]->GetProbabilities())[k]
	       << G4endl; 
      }
    }
    G4cout << G4endl; 
  }
  G4cout << G4endl;
  G4cout << "----------------------------------------------------------------"
	 << G4endl;
  G4cout << G4endl;
  G4cout << "### Pairs of stable fragments: " << G4endl;
  for(G4int A=2; A<maxA; ++A) {
    G4cout << "  A= " << A<<G4endl; 
    for(size_t j=0; j<list_p[A].size(); ++j) {
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
  G4cout << "### Pairs of stable+unstable fragments: " << G4endl;
  for(G4int A=2; A<maxA; ++A) {
    G4cout << "  A= " << A << G4endl; 
    for(size_t j=0; j<list_u[A].size(); ++j) {
      const G4FermiFragment* f1 = (list_u[A])[j]->GetFragment1();
      const G4FermiFragment* f2 = (list_u[A])[j]->GetFragment2();
      G4int a1 = f1->GetA();
      G4int z1 = f1->GetZ();
      G4int a2 = f2->GetA();
      G4int z2 = f2->GetZ();
      G4cout << "("<<a1<<","<<z1<<")("<<a2<<","<<z2<<") % Eex(MeV)= "
	     << std::setw(8)<< (list_u[A])[j]->GetExcitationEnergy() 
             << " Eex1= " << std::setw(8)<< f1->GetExcitationEnergy()
             << " Eex2= " << std::setw(8)<< f2->GetExcitationEnergy()
	     << G4endl; 
    }
    G4cout << G4endl;
    G4cout << "----------------------------------------------------------------"
	   << G4endl;
  }
  G4cout.precision(prec);
}

