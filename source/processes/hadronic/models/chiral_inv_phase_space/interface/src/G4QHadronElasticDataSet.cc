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
// $Id: G4QHadronElasticDataSet.cc,v 1.3 2010-05-26 12:19:06 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// GEANT4 physics class: G4QHadronElasticDataSet -- header file
// Created by M. Kosov (Mikhail.Kossov@cern.ch) 21.01.10
//
// ----------------------------------------------------------------------
// Short description: G4hadr wrapper for CHIPS elastic hA cross-sections.
// ----------------------------------------------------------------------
// 

#include "G4QHadronElasticDataSet.hh"

// Initialization of static vectors
std::vector<G4int> G4QHadronElasticDataSet::ElementZ; // Z of the element(i) in LastCalc
std::vector<std::vector<G4int>*> G4QHadronElasticDataSet::ElIsoN; // N of iso(j) of El(i)
std::vector<std::vector<G4double>*> G4QHadronElasticDataSet::IsoProbInEl;//SumProbIsoInEl

G4QHadronElasticDataSet::G4QHadronElasticDataSet()
{
  Isotopes = G4QIsotope::Get(); // Pointer to the G4QIsotopes singleton
}

G4bool G4QHadronElasticDataSet::IsApplicable(const G4DynamicParticle* P,const G4Element*)
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
  return false;
}

G4bool G4QHadronElasticDataSet::IsZAApplicable(const G4DynamicParticle* Pt,
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
  return false;
}

G4double G4QHadronElasticDataSet::GetCrossSection(const G4DynamicParticle* Pt,
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
          G4cerr<<"G4QHadronElasticDataSet::GetCrossSection"<<": Z="
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

G4double G4QHadronElasticDataSet::GetIsoZACrossSection(const G4DynamicParticle* Pt,
                                                         G4double Z, G4double A, G4double)
{
  G4ParticleDefinition* particle = Pt->GetDefinition();
  G4double Momentum=Pt->GetTotalMomentum();
  G4VQCrossSection* CSmanager=0;
  G4VQCrossSection* CSmanager2=0;
  //G4VQCrossSection* CSmanager2=0;
  G4int pPDG=0;
  if(particle == G4Neutron::Neutron())
  {
    CSmanager=G4QNeutronElasticCrossSection::GetPointer();
    pPDG=2112;
  }
  else if(particle == G4Proton::Proton())
  {
    CSmanager=G4QProtonElasticCrossSection::GetPointer();
    pPDG=2212;
  }
  else if(particle == G4PionMinus::PionMinus())
  {
    CSmanager=G4QPionMinusElasticCrossSection::GetPointer();
    pPDG=-211;
  }
  else if(particle == G4PionPlus::PionPlus())
  {
    CSmanager=G4QPionPlusElasticCrossSection::GetPointer();
    pPDG=211;
  }
  else if(particle == G4KaonMinus::KaonMinus())
  {
    CSmanager=G4QKaonMinusElasticCrossSection::GetPointer();
    pPDG=-321;
  }
  else if(particle == G4KaonPlus::KaonPlus())
  {
    CSmanager=G4QKaonPlusElasticCrossSection::GetPointer();
    pPDG=321;
  }
  else if(particle == G4KaonZeroLong::KaonZeroLong()   ||
          particle == G4KaonZeroShort::KaonZeroShort() ||
          particle == G4KaonZero::KaonZero()           ||
          particle == G4AntiKaonZero::AntiKaonZero()   )
  {
    CSmanager=G4QKaonMinusElasticCrossSection::GetPointer();
    CSmanager2=G4QKaonPlusElasticCrossSection::GetPointer();
    if(G4UniformRand() > 0.5) pPDG= 321;
    else                      pPDG=-321;
  }
  else if(particle == G4Lambda::Lambda())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3122;
  }
  else if(particle == G4SigmaPlus::SigmaPlus())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3222;
  }
  else if(particle == G4SigmaMinus::SigmaMinus())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3112;
  }
  else if(particle == G4SigmaZero::SigmaZero())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3212;
  }
  else if(particle == G4XiMinus::XiMinus())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3312;
  }
  else if(particle == G4XiZero::XiZero())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3322;
  }
  else if(particle == G4OmegaMinus::OmegaMinus())
  {
    CSmanager=G4QHyperonElasticCrossSection::GetPointer();
    pPDG=3334;
  }
  else if(particle == G4AntiNeutron::AntiNeutron())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-2112;
  }
  else if(particle == G4AntiProton::AntiProton())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-2212;
  }
  else if(particle == G4AntiLambda::AntiLambda())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3122;
  }
  else if(particle == G4AntiSigmaPlus::AntiSigmaPlus())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3222;
  }
  else if(particle == G4AntiSigmaMinus::AntiSigmaMinus())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3112;
  }
  else if(particle == G4AntiSigmaZero::AntiSigmaZero())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3212;
  }
  else if(particle == G4AntiXiMinus::AntiXiMinus())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3312;
  }
  else if(particle == G4AntiXiZero::AntiXiZero())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3322;
  }
  else if(particle == G4AntiOmegaMinus::AntiOmegaMinus())
  {
    CSmanager=G4QAntiBaryonElasticCrossSection::GetPointer();
    pPDG=-3334;
  }
  else G4cout<<"-Warning-G4QHadronElasticDataSet::GetIsoZACrossSection: PDG="
             <<particle->GetPDGEncoding()<<" isn't supported by CHIPS"<<G4endl;
  G4int tZ=(G4int)(Z);
  G4int tN=(G4int)(A-Z);
  G4double CSI=CSmanager->GetCrossSection(true, Momentum, tZ, tN, pPDG); // CS(j,i) basic
  if(CSmanager2) CSI= (CSI  +CSmanager2->GetCrossSection(true, Momentum, tZ, tN, pPDG))/2.;
  return CSI;
}
