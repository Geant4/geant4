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
// $Id: G4QHadronInelasticDataSet.cc,v 1.2 2010-05-26 12:19:06 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// GEANT4 physics class: G4QHadronInelasticDataSet -- header file
// Created by M. Kosov (Mikhail.Kossov@cern.ch) 11.11.09
//
// ------------------------------------------------------------------------
// Short description: G4hadr wrapper for CHIPS inelastic hA cross-sections.
// ------------------------------------------------------------------------
// 

#include "G4QHadronInelasticDataSet.hh"

// Initialization of static vectors
std::vector<G4int> G4QHadronInelasticDataSet::ElementZ; // Z of the element(i) in LastCalc
std::vector<std::vector<G4int>*> G4QHadronInelasticDataSet::ElIsoN; // N of iso(j) of El(i)
std::vector<std::vector<G4double>*> G4QHadronInelasticDataSet::IsoProbInEl;//SumProbIsoInEl

G4QHadronInelasticDataSet::G4QHadronInelasticDataSet()
{
  //CHIPSpAin    = G4QProtonNuclearCrossSection::GetPointer();
  //CHIPSnAin    = G4QNeutronNuclearCrossSection::GetPointer();
  //CHIPSpimAin  = G4QPionMinusNuclearCrossSection::GetPointer();
  //CHIPSpipAin  = G4QPionPlusNuclearCrossSection::GetPointer();
  //CHIPSkpAin   = G4QKaonPlusNuclearCrossSection::GetPointer();
  //CHIPSkmAin   = G4QKaonMinusNuclearCrossSection::GetPointer();
  //CHIPSk0Ain   = G4QKaonZeroNuclearCrossSection::GetPointer();
  //CHIPShAin    = G4QHyperonNuclearCrossSection::GetPointer();
  //CHIPShpAin   = G4QHyperonPlusNuclearCrossSection::GetPointer();
  //CHIPSabpAin  = G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
  //CHIPSabAin   = G4QAntiBaryonNuclearCrossSection::GetPointer();
  ////CHIPSphAin   = G4QPhotonNuclearCrossSection::GetPointer();
  ////CHIPSeAin    = G4QElectronNuclearCrossSection::GetPointer();
  ////CHIPSmuAin   = G4QMuonNuclearCrossSection::GetPointer();
  ////CHIPStauAin  = G4QTauNuclearCrossSection::GetPointer();
  ////CHIPSnumAin  = G4QNuMuNuclearCrossSection::GetPointer();
  ////CHIPSanumAin = G4QANuMuNuclearCrossSection::GetPointer();
  ////CHIPSnueAin  = G4QNuENuclearCrossSection::GetPointer();
  ////CHIPSanueAin = G4QANuENuclearCrossSection::GetPointer();
  ////CHIPSnunuAin = G4QNuNuNuclearCrossSection::GetPointer();
  ////CHIPSananAin = G4QANuANuNuclearCrossSection::GetPointer();

  Isotopes = G4QIsotope::Get(); // Pointer to the G4QIsotopes singleton
}

G4bool G4QHadronInelasticDataSet::IsApplicable(const G4DynamicParticle* P,const G4Element*)
{
  G4ParticleDefinition* particle = P->GetDefinition();
  if      (particle ==         G4Neutron::Neutron()        ) return true; 
  else if (particle ==          G4Proton::Proton()         ) return true;
  else if (particle ==       G4PionMinus::PionMinus()      ) return true;
  else if (particle ==        G4PionPlus::PionPlus()       ) return true;
  else if (particle ==        G4KaonPlus::KaonPlus()       ) return true;
  else if (particle ==       G4KaonMinus::KaonMinus()      ) return true;
  else if (particle ==    G4KaonZeroLong::KaonZeroLong()   ) return true;
  else if (particle ==   G4KaonZeroShort::KaonZeroShort()  ) return true;
  else if (particle ==          G4Lambda::Lambda()         ) return true;
  else if (particle ==       G4SigmaPlus::SigmaPlus()      ) return true;
  else if (particle ==      G4SigmaMinus::SigmaMinus()     ) return true;
  else if (particle ==       G4SigmaZero::SigmaZero()      ) return true;
  else if (particle ==         G4XiMinus::XiMinus()        ) return true;
  else if (particle ==          G4XiZero::XiZero()         ) return true;
  else if (particle ==      G4OmegaMinus::OmegaMinus()     ) return true;
  else if (particle ==     G4AntiNeutron::AntiNeutron()    ) return true;
  else if (particle ==      G4AntiProton::AntiProton()     ) return true;
  else if (particle ==      G4AntiLambda::AntiLambda()     ) return true;
  else if (particle ==   G4AntiSigmaPlus::AntiSigmaPlus()  ) return true;
  else if (particle ==  G4AntiSigmaMinus::AntiSigmaMinus() ) return true;
  else if (particle ==   G4AntiSigmaZero::AntiSigmaZero()  ) return true;
  else if (particle ==     G4AntiXiMinus::AntiXiMinus()    ) return true;
  else if (particle ==      G4AntiXiZero::AntiXiZero()     ) return true;
  else if (particle ==  G4AntiOmegaMinus::AntiOmegaMinus() ) return true;
  //else if (particle ==        G4MuonPlus::MuonPlus()       ) return true;
  //else if (particle ==       G4MuonMinus::MuonMinus()      ) return true;
  //else if (particle ==           G4Gamma::Gamma()          ) return true;
  //else if (particle ==        G4Electron::Electron()       ) return true;
  //else if (particle ==        G4Positron::Positron()       ) return true;
  //else if (particle ==         G4TauPlus::TauPlus()        ) return true;
  //else if (particle ==        G4TauMinus::TauMinus()       ) return true;
  //else if (particle ==   G4AntiNeutrinoE::AntiNeutrinoE()  ) return true;
  //else if (particle ==       G4NeutrinoE::NeutrinoE()      ) return true;
  //else if (particle ==  G4AntiNeutrinoMu::AntiNeutrinoMu() ) return true;
  //else if (particle ==      G4NeutrinoMu::NeutrinoMu()     ) return true;
  //else if (particle == G4AntiNeutrinoTau::AntiNeutrinoTau()) return true;
  //else if (particle ==     G4NeutrinoTau::NeutrinoTau()    ) return true;
  return false;
}

G4bool G4QHadronInelasticDataSet::IsZAApplicable(const G4DynamicParticle* Pt,
                                                 G4double, G4double)
{
  G4ParticleDefinition* particle = Pt->GetDefinition();
  if      (particle ==         G4Neutron::Neutron()        ) return true; // @@ isotopes?
  else if (particle ==          G4Proton::Proton()         ) return true;
  else if (particle ==       G4PionMinus::PionMinus()      ) return true;
  else if (particle ==        G4PionPlus::PionPlus()       ) return true;
  else if (particle ==        G4KaonPlus::KaonPlus()       ) return true;
  else if (particle ==       G4KaonMinus::KaonMinus()      ) return true;
  else if (particle ==    G4KaonZeroLong::KaonZeroLong()   ) return true;
  else if (particle ==   G4KaonZeroShort::KaonZeroShort()  ) return true;
  else if (particle ==          G4Lambda::Lambda()         ) return true;
  else if (particle ==       G4SigmaPlus::SigmaPlus()      ) return true;
  else if (particle ==      G4SigmaMinus::SigmaMinus()     ) return true;
  else if (particle ==       G4SigmaZero::SigmaZero()      ) return true;
  else if (particle ==         G4XiMinus::XiMinus()        ) return true;
  else if (particle ==          G4XiZero::XiZero()         ) return true;
  else if (particle ==      G4OmegaMinus::OmegaMinus()     ) return true;
  else if (particle ==     G4AntiNeutron::AntiNeutron()    ) return true;
  else if (particle ==      G4AntiProton::AntiProton()     ) return true;
  else if (particle ==      G4AntiLambda::AntiLambda()     ) return true;
  else if (particle ==   G4AntiSigmaPlus::AntiSigmaPlus()  ) return true;
  else if (particle ==  G4AntiSigmaMinus::AntiSigmaMinus() ) return true;
  else if (particle ==   G4AntiSigmaZero::AntiSigmaZero()  ) return true;
  else if (particle ==     G4AntiXiMinus::AntiXiMinus()    ) return true;
  else if (particle ==      G4AntiXiZero::AntiXiZero()     ) return true;
  else if (particle ==  G4AntiOmegaMinus::AntiOmegaMinus() ) return true;
  //else if (particle ==        G4MuonPlus::MuonPlus()       ) return true;
  //else if (particle ==       G4MuonMinus::MuonMinus()      ) return true;
  //else if (particle ==           G4Gamma::Gamma()          ) return true;
  //else if (particle ==        G4Electron::Electron()       ) return true;
  //else if (particle ==        G4Positron::Positron()       ) return true;
  //else if (particle ==         G4TauPlus::TauPlus()        ) return true;
  //else if (particle ==        G4TauMinus::TauMinus()       ) return true;
  //else if (particle ==   G4AntiNeutrinoE::AntiNeutrinoE()  ) return true;
  //else if (particle ==       G4NeutrinoE::NeutrinoE()      ) return true;
  //else if (particle ==  G4AntiNeutrinoMu::AntiNeutrinoMu() ) return true;
  //else if (particle ==      G4NeutrinoMu::NeutrinoMu()     ) return true;
  //else if (particle == G4AntiNeutrinoTau::AntiNeutrinoTau()) return true;
  //else if (particle ==     G4NeutrinoTau::NeutrinoTau()    ) return true;
  return false;
}

G4double G4QHadronInelasticDataSet::GetCrossSection(const G4DynamicParticle* Pt,
                                                    const G4Element* pElement,
                                                    G4double)
{
  G4int IPIE=IsoProbInEl.size();          // How many old elements?
  if(IPIE) for(G4int ip=0; ip<IPIE; ++ip) // Clean up the SumProb's of Isotopes (SPI)
  {
    std::vector<G4double>* SPI=IsoProbInEl[ip]; // Pointer to the SPI vector
    SPI->clear();
    delete SPI;
    std::vector<G4int>* IsN=ElIsoN[ip];   // Pointer to the N vector
    IsN->clear();
    delete IsN;
  }
  ElementZ.clear();                       // Clear the body vector for Z of Elements
  IsoProbInEl.clear();                    // Clear the body vector for SPI
  ElIsoN.clear();                         // Clear the body vector for N of Isotopes
  G4int Z = static_cast<G4int>(pElement->GetZ()); // Z of the Element
  ElementZ.push_back(Z);                  // Remember Z of the Element
  G4int isoSize=0;                        // The default for the isoVectorLength is 0
  G4int indEl=0;                          // Index of non-trivial element or 0(default)
  G4IsotopeVector* isoVector=pElement->GetIsotopeVector(); // Get the predefined IsoVect
  if(isoVector) isoSize=isoVector->size();// Get size of the existing isotopeVector
  if(isoSize)                             // The Element has non-trivial abundance set
  {
    indEl=pElement->GetIndex()+1;         // Index of the non-trivial element
    if(!Isotopes->IsDefined(Z,indEl))     // This index is not defined for this Z: define
    {
      std::vector<std::pair<G4int,G4double>*>* newAbund =
                                             new std::vector<std::pair<G4int,G4double>*>;
      G4double* abuVector=pElement->GetRelativeAbundanceVector();
      for(G4int j=0; j<isoSize; j++)      // Calculation of abundance vector for isotopes
      {
        G4int N=pElement->GetIsotope(j)->GetN()-Z; // N means A=N+Z !
        if(pElement->GetIsotope(j)->GetZ()!=Z)
          G4cerr<<"G4QHadronInelasticDataSet::GetCrossSection"<<": Z="
                <<pElement->GetIsotope(j)->GetZ()<<" # "<<Z<<G4endl;
        G4double abund=abuVector[j];
        std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
        newAbund->push_back(pr);
      }
      indEl=G4QIsotope::Get()->InitElement(Z,indEl,newAbund); // definition of the newInd
      for(G4int k=0; k<isoSize; k++) delete (*newAbund)[k];   // Cleaning temporary
      delete newAbund; // Was "new" in the beginning of the name space
    }
  }
  std::vector<std::pair<G4int,G4double>*>* cs= Isotopes->GetCSVector(Z,indEl);//CSPointer
  std::vector<G4double>* SPI = new std::vector<G4double>; // Pointer to the SPI vector
  IsoProbInEl.push_back(SPI);
  std::vector<G4int>* IsN = new std::vector<G4int>; // Pointer to the N vector
  ElIsoN.push_back(IsN);
  G4int nIs=cs->size();                   // A#Of Isotopes in the Element
  G4double susi=0.;                       // sum of CS over isotopes
  if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
  {
    std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
    G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
    IsN->push_back(N);                    // Remember Min N for the Element
    G4double CSI=GetIsoZACrossSection(Pt,Z,Z+N,0.);//CrossSection(j,i) for the isotope
    curIs->second = CSI;
    susi+=CSI;                            // Make a sum per isotopes
    SPI->push_back(susi);                 // Remember summed cross-section
  } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
  return Isotopes->GetMeanCrossSection(Z,indEl); // MeanCS over isotopes
}

G4double G4QHadronInelasticDataSet::GetIsoZACrossSection(const G4DynamicParticle* Pt,
                                                         G4double Z, G4double A, G4double)
{
  G4ParticleDefinition* particle = Pt->GetDefinition();
  G4double Momentum=Pt->GetTotalMomentum();
  G4VQCrossSection* CSmanager=0;
  //G4VQCrossSection* CSmanager2=0;
  G4int pPDG=0;
  if(particle == G4Neutron::Neutron())
  {
    CSmanager=G4QNeutronNuclearCrossSection::GetPointer();
    pPDG=2112;
  }
  else if(particle == G4Proton::Proton())
  {
    CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    pPDG=2212;
  }
  else if(particle == G4PionMinus::PionMinus())
  {
    CSmanager=G4QPionMinusNuclearCrossSection::GetPointer();
    pPDG=-211;
  }
  else if(particle == G4PionPlus::PionPlus())
  {
    CSmanager=G4QPionPlusNuclearCrossSection::GetPointer();
    pPDG=211;
  }
  else if(particle == G4KaonMinus::KaonMinus())
  {
    CSmanager=G4QKaonMinusNuclearCrossSection::GetPointer();
    pPDG=-321;
  }
  else if(particle == G4KaonPlus::KaonPlus())
  {
    CSmanager=G4QKaonPlusNuclearCrossSection::GetPointer();
    pPDG=321;
  }
  else if(particle == G4KaonZeroLong::KaonZeroLong()   ||
          particle == G4KaonZeroShort::KaonZeroShort() ||
          particle == G4KaonZero::KaonZero()           ||
          particle == G4AntiKaonZero::AntiKaonZero()   )
  {
    CSmanager=G4QKaonZeroNuclearCrossSection::GetPointer();
    if(G4UniformRand() > 0.5) pPDG= 311;
    else                      pPDG=-311;
  }
  else if(particle == G4Lambda::Lambda())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3122;
  }
  else if(particle == G4SigmaPlus::SigmaPlus())
  {
    CSmanager=G4QHyperonPlusNuclearCrossSection::GetPointer();
    pPDG=3222;
  }
  else if(particle == G4SigmaMinus::SigmaMinus())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3112;
  }
  else if(particle == G4SigmaZero::SigmaZero())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3212;
  }
  else if(particle == G4XiMinus::XiMinus())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3312;
  }
  else if(particle == G4XiZero::XiZero())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3322;
  }
  else if(particle == G4OmegaMinus::OmegaMinus())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3334;
  }
  else if(particle == G4AntiNeutron::AntiNeutron())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-2112;
  }
  else if(particle == G4AntiProton::AntiProton())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-2212;
  }
  else if(particle == G4AntiLambda::AntiLambda())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3122;
  }
  else if(particle == G4AntiSigmaPlus::AntiSigmaPlus())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3222;
  }
  else if(particle == G4AntiSigmaMinus::AntiSigmaMinus())
  {
    CSmanager=G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
    pPDG=-3112;
  }
  else if(particle == G4AntiSigmaZero::AntiSigmaZero())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3212;
  }
  else if(particle == G4AntiXiMinus::AntiXiMinus())
  {
    CSmanager=G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
    pPDG=-3312;
  }
  else if(particle == G4AntiXiZero::AntiXiZero())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3322;
  }
  else if(particle == G4AntiOmegaMinus::AntiOmegaMinus())
  {
    CSmanager=G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
    pPDG=-3334;
  }
  //else if(particle == G4Gamma::Gamma())
  //{
  //  CSmanager=G4QPhotonNuclearCrossSection::GetPointer();
  //  pPDG=22;
  //}
  //else if(particle == G4Electron::Electron() ||
  //        particle == G4Positron::Positron())
  //{
  //  CSmanager=G4QElectronNuclearCrossSection::GetPointer();
  //  pPDG=11;
  //}
  //else if(particle == G4MuonPlus::MuonPlus() ||
  //        particle == G4MuonMinus::MuonMinus())
  //{
  //  CSmanager=G4QMuonNuclearCrossSection::GetPointer();
  //  pPDG=13;
  //}
  //else if(particle == G4TauPlus::TauPlus() ||
  //        particle == G4TauMinus::TauMinus())
  //{
  //  CSmanager=G4QTauNuclearCrossSection::GetPointer();
  //  pPDG=15;
  //}
  //else if(particle == G4NeutrinoMu::NeutrinoMu() )
  //{
  //  CSmanager=G4QNuMuNuclearCrossSection::GetPointer();
  //  CSmanager2=G4QNuNuNuclearCrossSection::GetPointer();
  //  pPDG=14;
  //}
  //else if(particle == G4AntiNeutrinoMu::AntiNeutrinoMu() )
  //{
  //  CSmanager=G4QANuMuNuclearCrossSection::GetPointer();
  //  CSmanager2=G4QANuANuNuclearCrossSection::GetPointer();
  //  pPDG=-14;
  //}
  //else if(particle == G4NeutrinoE::NeutrinoE() )
  //{
  //  CSmanager=G4QNuENuclearCrossSection::GetPointer();
  //  CSmanager2=G4QNuNuNuclearCrossSection::GetPointer();
  //  pPDG=12;
  //}
  //else if(particle == G4AntiNeutrinoE::AntiNeutrinoE() )
  //{
  //  CSmanager=G4QANuENuclearCrossSection::GetPointer();
  //  CSmanager2=G4QANuANuNuclearCrossSection::GetPointer();
  //  pPDG=-12;
  //}
  else G4cout<<"-Warning-G4QHadronInelasticDataSet::GetIsoZACrossSection: PDG="
             <<particle->GetPDGEncoding()<<" isn't supported by CHIPS"<<G4endl;
  G4int tZ=(G4int)(Z);
  G4int tN=(G4int)(A-Z);
  G4double CSI=CSmanager->GetCrossSection(true, Momentum, tZ, tN, pPDG); // CS(j,i) basic
  //if(CSmanager2) CSI+=CSmanager2->GetCrossSection(true,Momentum,Z,N,pPDG); // e.g.(nu,nu)
  return CSI;
}
