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
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TParticlesDAL.h"


G4TParticlesDAL *gParticlesDAL = new G4TParticlesDAL();


ClassImp(G4TParticlesDAL)

using namespace std;
using namespace ROOT;
using namespace TMath;


//______________________________________________________________________________
TString G4TParticlesDAL::GetCut(Int_t PDG)
{
  if     (PDG == 1000000010) return "!ND && PDG == 2112 || (!NS && !NZ && NN==1)";//neutron
  else if(PDG == 1000010010) return "!ND && PDG == 2112 || (!NS && NZ==1 && !NN)";//proton
  else if(PDG == 1000010020) return "!ND && NS == 0 && NZ == 1 && NN == 1";      //deuteron
  else if(PDG == 1000010030) return "!ND && NS == 0 && NZ == 1 && NN == 2";      //triton
  else if(PDG == 1000020030) return "!ND && NS == 0 && NZ == 2 && NN == 1";      //he3
  else if(PDG == 1000020040) return "!ND && NS == 0 && NZ == 2 && NN == 2";      //he4
  else if(PDG == 211) return "!ND && PDG == 211 || (!NS && NZ==1 && NN==-1)";    //pi+
  else if(PDG ==-211) return "!ND && PDG ==-211 || (!NS && NZ==-1 && NN==1)";    //pi-
  else if(PDG == 321) return "!ND && PDG == 321 || (NS==-1 && NZ==1 && !NN)";    //K+
  else if(PDG ==-321) return "!ND && PDG == 321 || (NS==1 && NZ==-1 && !NN)";    //K-

  cout << "Warning>G4TParticlesDAL::GetCut: Cut for PDG=" << PDG << " isn't found" << endl;
  return "";
}


//____________________________________________________@@_keV?____________________
Double_t G4TParticlesDAL::ComputeNuclearMass(Double_t AtomicMass, Double_t A, Int_t Z)
{
  Double_t const u  = 931.494043;   // Unified atomic mass unit
  Double_t const me  = 0.510998918; // Electron mass
  return AtomicMass/1.e3 + u*A - Z*me +
         ( 14.4381 * std::pow(Z, 2.39) + 1.55468e-6 * std::pow(Z, 5.35) ) / 1.e6;
}


//______________________________________________________________________________
TGeoElementRN* G4TParticlesDAL::GetParticle(Int_t PDG)
{
  if(!fNuclearMassesRead) ReadNuclearMasses();

  Int_t A = GetA(PDG);                    // @@ improve GetA
  Int_t Z = GetZ(PDG);                    // @@ improve GetZ
  return fTable->GetElementRN(A,Z);
}

//______________________________________________________________________________
Int_t G4TParticlesDAL::GetPDG(TString const& particleName)
{
  cout<<"G4TParticlesDAL::GetPDG(in): string="<<particleName<<endl;
  Int_t PdgCode = 0;
  if(!fNuclearMassesRead) ReadNuclearMasses();
  if     (particleName == "p"  ) PdgCode = 1000010010;
  else if(particleName == "n"  ) PdgCode = 1000000010;
  else if(particleName == "d"  ) PdgCode = 1000010020;
  else if(particleName == "t"  ) PdgCode = 1000010030;
  else if(particleName == "pi+") PdgCode =  211;
  else if(particleName == "pi-") PdgCode = -211;
  else if(particleName == "K+" ) PdgCode =  321;
  else if(particleName == "K-" ) PdgCode = -321;
  else
  {
    Int_t A = GetA(particleName);
    Int_t Z = GetZ(particleName);
    PdgCode = 1000000000;
    if(A > 0 && Z > 0)
    {
      PdgCode += Z * 10000;
      PdgCode += A * 10;
    }
    else cout<< "Warning>G4TParticlesDAL::GetPDG: No PDG for Name=" << particleName <<endl;
  }
  cout<<"G4TParticlesDAL::GetPDG(out): PDG="<<PdgCode<<endl;
  return PdgCode;
}


//______________________________________________________________________________
Int_t G4TParticlesDAL::GetPDG(Int_t Z, Int_t A, Int_t S)
{
  Bool_t aB=false;
  Int_t PdgCode = 0;
  if(!A)                               // Mesons
  {
    if(!S)                             // Not strange mesons
    {
      if     (Z== 1) return 211;
      else if(Z==-1) return -211;
      else cout<<"Warning>G4TParticlesDAL::GetPDG: No PDG for A=0, S=0, Z="<< Z <<endl;
    }
    else                               // Strange mesons
    {
      if     (S==-1 && Z== 1) return  321;
      else if(S== 1 && Z==-1) return -321;
      else cout<<"Warning>G4TParticlesDAL::GetPDG: No PDG for A=0, S="<<S<<", Z="<<Z<<endl;
    }
  }
  else if(A < 0)
  {
    A=-A;
    Z=-Z;
    aB=true;
  }
  if(Z < 0) cout<<"Warning>G4TParticlesDAL::GetPDG: No PDG for A="<< A <<", Z="<< Z <<endl;
  else
  {
    PdgCode = 1000000000;
    PdgCode += Z * 10000;
    PdgCode += A * 10;
    if(aB) PdgCode=-PdgCode;
  }
  return PdgCode;
}


//______________________________________________________________________________
TString G4TParticlesDAL::GetParticleName(Int_t PDG, Bool_t inLatex)
{
  TString result;
  if     (PDG == 1000010010) return "p";
  else if(PDG == 1000000010) return "n";
  else if(PDG == 211)        return "pi+";
  else if(PDG ==-211)        return "pi-";
  else if(PDG == 321)        return "k+";
  else if(PDG ==-321)        return "k-";
  else if(PDG == 0)          return "";
  else
  {    
    Int_t Z = GetZ(PDG);
    Int_t A = GetA(PDG);
    if(!A)    cout<< "Warning>G4TParticlesDAL::GetParticleName: No Name for A = 0" <<endl;
    if(A < 0) cout<< "Warning>G4TParticlesDAL::GetParticleName: No Name for A < 0" <<endl;
    if(Z < 0) cout<< "Warning>G4TParticlesDAL::GetParticleName: No Name for Z < 0" <<endl;
    TString name(GetElementName(Z));

    if(!inLatex) result = TString::Format("%s%d", name.Data(), A);
    else         result = TString::Format("^{%d}%s",  A, name.Data());
  }
  return result;
}


//______________________________________________________________________________
TString G4TParticlesDAL::GetFileName(Int_t PDG)
{
  TString result = "";
  if     (PDG == 1000010010) return "p";
  else if(PDG == 1000000010) return "n";
  else if(PDG == 211) return "pip";
  else if(PDG ==-211) return "pim";
  else if(PDG == 321) return "kp";
  else if(PDG ==-321) return "km";
  else
  {
    Int_t Z = GetZ(PDG);
    Int_t A = GetA(PDG);
    if(!A)    cout<< "Warning>G4TParticlesDAL::GetFileName: No Name for A = 0" <<endl;
    if(A < 0) cout<< "Warning>G4TParticlesDAL::GetFileName: No Name for A < 0" <<endl;
    if(Z < 0) cout<< "Warning>G4TParticlesDAL::GetFileName: No Name for Z < 0" <<endl;

    TString name(GetElementName(Z));
    name.ToLower();

     // getting name for old ascii files
    result = TString::Format("%s%d%d", name.Data(), Z, A-Z);
  }
  return result;
}


//______________________________________________________________________________
Double_t G4TParticlesDAL::GetParticleMass(Int_t PDG)
{
  Double_t result = 0.;  
  if(PDG < 0) PDG=-PDG;
  if(PDG == 2212) PDG=1000010010;
  if(PDG >= 1000010010) // Charged nuclei/NuclearFragments
  {
    TGeoElementRN* particle = GetParticle(PDG);
    Double_t A = particle->A();
    Int_t Z = particle->Z();
    result = ComputeNuclearMass(particle->MassEx(), A, Z );
  }
  else if(PDG == 2112) result = 939.56536;
  else if(PDG ==  211) result = 139.57018;
  else if(PDG ==  321) result = 493.677;
  else cout << "Warning>G4TParticlesDAL::GetParticleMass: No mass for PDG=" << PDG << endl;
  return result;
}

//______________________________________________________________________________
Int_t G4TParticlesDAL::GetA(Int_t PDG)
{
  PDG /= 10;
  return PDG % 1000;
}

//______________________________________________________________________________
Int_t G4TParticlesDAL::GetZ(Int_t PDG)
{
  PDG /= 10000;
  return PDG % 1000;
}

//______________________________________________________________________________
Int_t G4TParticlesDAL::GetN(Int_t PDG)
{
  PDG /= 10;
  Int_t A = PDG % 1000;
  PDG /= 1000;
  Int_t Z = PDG % 1000;
  return A - Z;
}


//______________________________________________________________________________
Int_t G4TParticlesDAL::GetA(TString const& elementName)
{
  Int_t result = 0;
  if     (elementName == "p" || elementName == "n") return 1;
  else if(elementName == "pi+" || elementName == "pi-" || elementName == "K+" ||
          elementName == "K-") return 0;
  else
  {
    if(!fNuclearMassesRead) ReadNuclearMasses();

    TString numb(elementName);
    while(!numb.IsDigit() && numb.Length() != 0) numb = numb.Replace(0, 1, "");

    result = atoi(numb);
  }
  return result;
}

//______________________________________________________________________________
Int_t G4TParticlesDAL::GetZ(TString const& elementName)
{
  Int_t result = 0;
  if     (elementName == "n") return 0;
  else if(elementName == "p" || elementName == "pi+" || elementName == "K+") return 1;
  else if(elementName == "pi-"  || elementName == "K-") return -1;
  else
  {
    if(!fNuclearMassesRead) ReadNuclearMasses();

    TString name(elementName);
    TString numb(elementName);

    while(!numb.IsDigit() && numb.Length() != 0) numb = numb.Replace(0,1,"");

    name = name.ReplaceAll(numb,""); // pure name

    for(UInt_t i = 0; i < fMendeleevTable.size(); ++i)
    {
      MTableEntry_t item = fMendeleevTable[i];

      if(item.fName == name)
      {
        result = item.fZ;
        break;
      }
    }
  }
  return result;
}

//______________________________________________________________________________
TString G4TParticlesDAL::GetElementName(Int_t Z)
{
  if(!fNuclearMassesRead) ReadNuclearMasses();
  TString result;

  for(UInt_t i = 0; i < fMendeleevTable.size(); ++i)
  {
    MTableEntry_t item = fMendeleevTable[i];

    if(item.fZ == Z)
    {
      result = item.fName;
      break;
    }
  }
  return result;
}

//______________________________________________________________________________
void G4TParticlesDAL::ReadNuclearMasses()
{
  if(fTable == 0)
  {
    fGeoManager = new TGeoManager("particlesDB","Particles DB Manager");
    fTable = fGeoManager->GetElementTable();
  }

  fMendeleevTable.push_back( MTableEntry_t(  1, "H" ) );
  fMendeleevTable.push_back( MTableEntry_t(  2, "He") );
  fMendeleevTable.push_back( MTableEntry_t(  3, "Li") );
  fMendeleevTable.push_back( MTableEntry_t(  4, "Be") );
  fMendeleevTable.push_back( MTableEntry_t(  5, "B" ) );
  fMendeleevTable.push_back( MTableEntry_t(  6, "C" ) );
  fMendeleevTable.push_back( MTableEntry_t(  7, "N" ) );
  fMendeleevTable.push_back( MTableEntry_t(  8, "O" ) );
  fMendeleevTable.push_back( MTableEntry_t(  9, "F" ) );
  fMendeleevTable.push_back( MTableEntry_t( 10, "Ne") );
  fMendeleevTable.push_back( MTableEntry_t( 11, "Na") );
  fMendeleevTable.push_back( MTableEntry_t( 12, "Mg") );
  fMendeleevTable.push_back( MTableEntry_t( 13, "Al") );
  fMendeleevTable.push_back( MTableEntry_t( 14, "Si") );
  fMendeleevTable.push_back( MTableEntry_t( 15, "P" ) );
  fMendeleevTable.push_back( MTableEntry_t( 16, "S" ) );
  fMendeleevTable.push_back( MTableEntry_t( 17, "Cl") );
  fMendeleevTable.push_back( MTableEntry_t( 18, "Ar") );
  fMendeleevTable.push_back( MTableEntry_t( 19, "K" ) );
  fMendeleevTable.push_back( MTableEntry_t( 20, "Ca") );
  fMendeleevTable.push_back( MTableEntry_t( 21, "Sc") );
  fMendeleevTable.push_back( MTableEntry_t( 22, "Ti") );
  fMendeleevTable.push_back( MTableEntry_t( 23, "V" ) );
  fMendeleevTable.push_back( MTableEntry_t( 24, "Cr") );
  fMendeleevTable.push_back( MTableEntry_t( 25, "Mn") );
  fMendeleevTable.push_back( MTableEntry_t( 26, "Fe") );
  fMendeleevTable.push_back( MTableEntry_t( 27, "Co") );
  fMendeleevTable.push_back( MTableEntry_t( 28, "Ni") );
  fMendeleevTable.push_back( MTableEntry_t( 29, "Cu") );
  fMendeleevTable.push_back( MTableEntry_t( 30, "Zn") );
  fMendeleevTable.push_back( MTableEntry_t( 31, "Ga") );
  fMendeleevTable.push_back( MTableEntry_t( 32, "Ge") );
  fMendeleevTable.push_back( MTableEntry_t( 33, "As") );
  fMendeleevTable.push_back( MTableEntry_t( 34, "Se") );
  fMendeleevTable.push_back( MTableEntry_t( 35, "Br") );
  fMendeleevTable.push_back( MTableEntry_t( 36, "Kr") );
  fMendeleevTable.push_back( MTableEntry_t( 37, "Rb") );
  fMendeleevTable.push_back( MTableEntry_t( 38, "Sr") );
  fMendeleevTable.push_back( MTableEntry_t( 39, "Y" ) );
  fMendeleevTable.push_back( MTableEntry_t( 40, "Zr") );
  fMendeleevTable.push_back( MTableEntry_t( 41, "Nb") );
  fMendeleevTable.push_back( MTableEntry_t( 42, "Mo") );
  fMendeleevTable.push_back( MTableEntry_t( 43, "Tc") );
  fMendeleevTable.push_back( MTableEntry_t( 44, "Ru") );
  fMendeleevTable.push_back( MTableEntry_t( 45, "Rh") );
  fMendeleevTable.push_back( MTableEntry_t( 46, "Pd") );
  fMendeleevTable.push_back( MTableEntry_t( 47, "Ag") );
  fMendeleevTable.push_back( MTableEntry_t( 48, "Cd") );
  fMendeleevTable.push_back( MTableEntry_t( 49, "In") );
  fMendeleevTable.push_back( MTableEntry_t( 50, "Sn") );
  fMendeleevTable.push_back( MTableEntry_t( 51, "Sb") );
  fMendeleevTable.push_back( MTableEntry_t( 52, "Te") );
  fMendeleevTable.push_back( MTableEntry_t( 53, "I" ) );
  fMendeleevTable.push_back( MTableEntry_t( 54, "Xe") );
  fMendeleevTable.push_back( MTableEntry_t( 55, "Cs") );
  fMendeleevTable.push_back( MTableEntry_t( 56, "Ba") );
  fMendeleevTable.push_back( MTableEntry_t( 57, "La") );
  fMendeleevTable.push_back( MTableEntry_t( 58, "Ce") );
  fMendeleevTable.push_back( MTableEntry_t( 59, "Pr") );
  fMendeleevTable.push_back( MTableEntry_t( 60, "Nd") );
  fMendeleevTable.push_back( MTableEntry_t( 61, "Pm") );
  fMendeleevTable.push_back( MTableEntry_t( 62, "Sm") );
  fMendeleevTable.push_back( MTableEntry_t( 63, "Eu") );
  fMendeleevTable.push_back( MTableEntry_t( 64, "Gd") );
  fMendeleevTable.push_back( MTableEntry_t( 65, "Tb") );
  fMendeleevTable.push_back( MTableEntry_t( 66, "Dy") );
  fMendeleevTable.push_back( MTableEntry_t( 67, "Ho") );
  fMendeleevTable.push_back( MTableEntry_t( 68, "Er") );
  fMendeleevTable.push_back( MTableEntry_t( 69, "Tm") );
  fMendeleevTable.push_back( MTableEntry_t( 70, "Yb") );
  fMendeleevTable.push_back( MTableEntry_t( 71, "Lu") );
  fMendeleevTable.push_back( MTableEntry_t( 72, "Hf") );
  fMendeleevTable.push_back( MTableEntry_t( 73, "Ta") );
  fMendeleevTable.push_back( MTableEntry_t( 74, "W" ) );
  fMendeleevTable.push_back( MTableEntry_t( 75, "Re") );
  fMendeleevTable.push_back( MTableEntry_t( 76, "Os") );
  fMendeleevTable.push_back( MTableEntry_t( 77, "Ir") );
  fMendeleevTable.push_back( MTableEntry_t( 78, "Pt") );
  fMendeleevTable.push_back( MTableEntry_t( 79, "Au") );
  fMendeleevTable.push_back( MTableEntry_t( 80, "Hg") );
  fMendeleevTable.push_back( MTableEntry_t( 81, "Tl") );
  fMendeleevTable.push_back( MTableEntry_t( 82, "Pb") );
  fMendeleevTable.push_back( MTableEntry_t( 83, "Bi") );
  fMendeleevTable.push_back( MTableEntry_t( 84, "Po") );
  fMendeleevTable.push_back( MTableEntry_t( 85, "At") );
  fMendeleevTable.push_back( MTableEntry_t( 86, "Rn") );
  fMendeleevTable.push_back( MTableEntry_t( 87, "Fr") );
  fMendeleevTable.push_back( MTableEntry_t( 88, "Ra") );
  fMendeleevTable.push_back( MTableEntry_t( 89, "Ac") );
  fMendeleevTable.push_back( MTableEntry_t( 90, "Th") );
  fMendeleevTable.push_back( MTableEntry_t( 91, "Pa") );
  fMendeleevTable.push_back( MTableEntry_t( 92, "U" ) );
  fMendeleevTable.push_back( MTableEntry_t( 93, "Np") );
  fMendeleevTable.push_back( MTableEntry_t( 94, "Pu") );
  fMendeleevTable.push_back( MTableEntry_t( 95, "Am") );
  fMendeleevTable.push_back( MTableEntry_t( 96, "Cm") );
  fMendeleevTable.push_back( MTableEntry_t( 97, "Bk") );
  fMendeleevTable.push_back( MTableEntry_t( 98, "Cf") );
  fMendeleevTable.push_back( MTableEntry_t( 99, "Es") );
  fMendeleevTable.push_back( MTableEntry_t(100, "Fm") );
  fMendeleevTable.push_back( MTableEntry_t(101, "Md") );
  fMendeleevTable.push_back( MTableEntry_t(102, "No") );
  fMendeleevTable.push_back( MTableEntry_t(103, "Lr") );
  fMendeleevTable.push_back( MTableEntry_t(104, "Rf") );
  fMendeleevTable.push_back( MTableEntry_t(105, "Db") );
  fMendeleevTable.push_back( MTableEntry_t(106, "Sg") );
  fMendeleevTable.push_back( MTableEntry_t(107, "Bh") );
  fMendeleevTable.push_back( MTableEntry_t(108, "Hs") );
  fMendeleevTable.push_back( MTableEntry_t(109, "Mt") );
  fMendeleevTable.push_back( MTableEntry_t(110, "Ds") );
  fMendeleevTable.push_back( MTableEntry_t(111, "Rg") );

  fNuclearMassesRead = true;
}



