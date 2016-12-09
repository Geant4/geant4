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
// $Id: G4LevelReader.cc 88407 2015-02-18 09:18:44Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NucLevel
//
//      Author:        V.Ivanchenko (M.Kelsey reading method is used)
// 
//      Creation date: 4 January 2012
//
//      Modifications:
//
// -------------------------------------------------------------------

#include "G4LevelReader.hh"
#include "G4NucleiProperties.hh"
#include "G4NucLevel.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4SystemOfUnits.hh"
#include "G4Pow.hh"
#include <vector>
#include <fstream>
#include <sstream>

G4String G4LevelReader::fFloatingLevels[] = {
  "-", "+X", "+Y", "+Z", "+U", "+V", "+W", "+R", "+S", "+T", "+A", "+B", "+C"};

G4LevelReader::G4LevelReader(G4NuclearLevelData* ptr) 
  : fData(ptr),fMinProbability(1.e-8),fVerbose(0),fLevelMax(632),fTransMax(30)
{
  fParam = fData->GetParameters();
  fTimeFactor = CLHEP::second/G4Pow::GetInstance()->logZ(2);
  char* directory = getenv("G4LEVELGAMMADATA");
  if(directory) {
    fDirectory = directory;
  } else {
    G4Exception("G4LevelReader()","had0707",FatalException,
		"Environment variable G4LEVELGAMMADATA is not defined");
    fDirectory = "";
  } 
  fFile = fDirectory + "/z100.a200";
  fPol = "  ";
  for(G4int i=0; i<10; ++i) { fICC[i] = 0.0; }
  for(G4int i=0; i<nbufmax; ++i) { buffer[i] = ' '; }
  for(G4int i=0; i<nbuf2; ++i)   { buff2[i] = ' '; }
  bufp[0] = bufp[1] = ' ';

  fEnergy = fCurrEnergy = fTrEnergy = fProb = fTime = 
    fSpin = fAlpha = fRatio = 0.0;
  fNorm1 = fNorm2 = 0.0f;

  vIndex.resize(fTransMax,0);
  vTrans.resize(fTransMax,0);
  vRatio.resize(fTransMax,0.0f);
  vGammaCumProbability.resize(fTransMax,0.0f);
  vGammaECumProbability.resize(fTransMax,0.0f);
  vGammaProbability.resize(fTransMax,0.0f);
  vShellProbability.resize(fTransMax,nullptr);
  vMpRatio.resize(fTransMax,0.0f);

  vEnergy.resize(fLevelMax,0.0f);
  vTime.resize(fLevelMax,0.0f);
  vTimeg.resize(fLevelMax,0.0f);
  vSpin.resize(fLevelMax,0);
  vLevel.resize(fLevelMax,nullptr);
  vMeta.resize(fLevelMax,0);
  vIndexDB.resize(fLevelMax,-1);
}

const G4LevelManager* 
G4LevelReader::CreateLevelManager(G4int Z, G4int A)
{
  std::ostringstream ss;
  ss << "/z" << Z << ".a" << A;
  G4String st =  G4String(ss.str());
  fFile = fDirectory + st;
  return MakeLevelManager(Z, A, fFile);
}

const G4LevelManager* 
G4LevelReader::MakeLevelManager(G4int Z, G4int A, const G4String& filename)
{
  vEnergy.resize(1,0.0f);
  vTime.resize(1,FLT_MAX);
  vTimeg.resize(1,FLT_MAX);
  vSpin.resize(1,0);
  vMeta.resize(1,0);
  vLevel.resize(1,nullptr);

  std::ifstream infile(filename, std::ios::in);

  // file is not opened
  if (!infile.is_open()) {
    if (fVerbose > 0) {
      G4cout << " G4LevelReader: fail open file for Z= " 
	     << Z << " A= " << A 
	     << " <" << filename << ">" <<  G4endl;
    }

  } else {

    if (fVerbose > 0) {
      G4cout << "G4LevelReader: open file for Z= " 
	     << Z << " A= " << A 
	     << " <" << filename << ">" <<  G4endl;
    }
    // read line by line
    G4bool end = false;
    G4bool next = true;
    G4int nline = -1;
    G4String xl = "- ";
    fCurrEnergy = DBL_MAX;
    do {
      fPol = "  ";
      ++nline;
      G4bool isOK = (ReadDataItem(infile,fEnergy)   &&
		     ReadDataItem(infile,fTrEnergy) &&
		     ReadDataItem(infile,fProb)     &&
		     ReadDataItem(infile,fPol)      &&
		     ReadDataItem(infile,fTime)     &&
		     ReadDataItem(infile,fSpin)     &&
		     ReadDataItem(infile,fAlpha));

      fEnergy   *= CLHEP::keV;
      fTrEnergy *= CLHEP::keV;

      if(isOK) {
	for(G4int i=0; i<10; ++i) { 
	  isOK = (isOK && (ReadDataItem(infile,fICC[i])));
	}
      }
      if(!isOK) {
	end = true; 
	next = false; 
      }
      if(fVerbose > 1) {
	G4cout << "#Line " << nline << " " << fEnergy << " " << fTrEnergy 
	       << " " << fProb << " " <<  fPol << " " << fTime << " " 
	       << fSpin << " " << fAlpha << G4endl;
	G4cout << "      ";
	for(G4int i=0; i<10; ++i) { G4cout << fICC[i] << " "; }
	G4cout << G4endl;
      }
      // end of nuclear level data
      if(end || fEnergy > fCurrEnergy) {
	size_t nn = vTrans.size();
	if(nn > 0) {
	  --nn;
	  if(fVerbose > 1) {
	    G4cout << "Reader: new level E= " << fEnergy 
		   << " Ntransitions= " << nn+1 << " fNorm1= " << fNorm1 
		   << " fNorm2= " << fNorm2 << G4endl;
	  } 
	  if(fNorm1 > 0.0f) {
	    fNorm1 = 1.0f/fNorm1; 
	    vTimeg.push_back(((G4float)(fTime*fTimeFactor))*fNorm2*fNorm1);
	    fNorm2 = 1.0f/fNorm2; 
	    for(size_t i=0; i<nn; ++i) {
	      vGammaCumProbability[i]  *= fNorm1;
	      vGammaECumProbability[i] *= fNorm2;
	      if(fVerbose > 2) {
		G4cout << "Probabilities[" << i 
		       << "]= " << vGammaCumProbability[i]
		       << "  " << vGammaECumProbability[i] 
		       << " idxTrans= " << vIndex[i]
		       << G4endl;
	      }
	    }
	    vGammaCumProbability[nn]  = 1.0f;
	    vGammaECumProbability[nn] = 1.0f;
	    if(fVerbose > 2) {
	      G4cout << "Probabilities[" << nn << "]= " << vGammaCumProbability[nn]
		     << "  " << vGammaECumProbability[nn] 
		     << " IdxTrans= " << vIndex[nn]
		     << G4endl;
	    }
            vMeta.push_back(0);
	    //case of X-level
	  } else {
            vMeta.push_back(1);
	    vTimeg.push_back(0.0f);
	    vGammaCumProbability[0]  = 0.0f;
	    vGammaECumProbability[0] = 0.0f;
	  }
      
	  vLevel.push_back(new G4NucLevel(vIndex.size(),
					  vIndex,
					  vTrans,
					  vGammaCumProbability,
					  vGammaECumProbability,
					  vGammaProbability,
					  vMpRatio,
					  vShellProbability));
          vIndex.clear();
          vTrans.clear();
          vGammaCumProbability.clear();
          vGammaECumProbability.clear();
          vGammaProbability.clear();
          vShellProbability.clear();
          vMpRatio.clear();
	}
        if(!end) { next = true; }
      }
      fCurrEnergy = fEnergy;
      // begin nuclear level data
      if(next) {
	if(fVerbose > 2) {
	  G4cout << "== Reader: begin of new level E= " << fEnergy << G4endl;
	} 
	// protection for bad level energy
	size_t nn = vEnergy.size();
        G4float ener = (G4float)fEnergy;
	if(0 < nn && vEnergy[nn-1] > ener) { ener = vEnergy[nn-1]; } 
	vEnergy.push_back(ener);
	vTime.push_back((G4float)(fTime*fTimeFactor));
	if(fSpin > 20.0) { fSpin = 0.0; }
	fProb = std::max(fProb, fMinProbability);  
	vSpin.push_back((G4int)(fSpin+fSpin));
	fNorm1 = 0.0f;
	fNorm2 = 0.0f;
	next = false;
      } 
      // continue filling level data
      if(!end) {
        if(fProb > 0.0) {
	  // by default transition to a ground state
          G4float efinal = std::max((G4float)(fEnergy - fTrEnergy),0.0f);
          G4float elevel = 0.0f;
          size_t  idxLevel = 0;
	  G4int tnum = 0;
	  // do not check initial energy
	  size_t nn = vEnergy.size();
	  static const G4float x_energy = (G4float)(0.1*CLHEP::eV);
	  if(1 < nn) {
	    G4float ediffMin = fEnergy;
	    for(size_t i=0; i<nn-1; ++i) {
              G4float ediff = std::abs(efinal - vEnergy[i]);
	      /*
	      G4cout << "Elevel[" << i << "]= " << vEnergy[i] 
		     << " Efinal= " << efinal
		     << " Ediff= " << ediff
		     << " EdiffMin= " << ediffMin << G4endl;
	      */
	      if(ediff < ediffMin) {
		ediffMin = ediff;
		elevel   = vEnergy[i];
		idxLevel = i;
                if(ediff <= x_energy) { break; }
	      }
	    }
	    if(std::abs(vEnergy[nn-1] - elevel) < x_energy) { tnum = 1; }
	  }
	  G4double x = 1.0 + fAlpha;
	  fNorm1 += (G4float)fProb;
	  fNorm2 += (G4float)(fProb*x);
	  vIndex.push_back(idxLevel);
	  vGammaCumProbability.push_back(fNorm1);
	  vGammaECumProbability.push_back(fNorm2);
	  vGammaProbability.push_back((G4float)(1.0/x));
	  vMpRatio.push_back(0.0f);
	  vTrans.push_back(tnum);
	  if(fAlpha > 0.0) {
	    vShellProbability.push_back(NormalizedICCProbability(Z));
	  } else {
	    vShellProbability.push_back(nullptr);
	  }
	}
      }
      if(nline > 10000) {
	G4cout << "G4LevelReader: Line #" << nline << " Z= " << Z << " A= "
	       << " this file is too long - stop loading" << G4endl;
	end = true;
      }
      // Loop checking, 05-Aug-2015, Vladimir Ivanchenko
    } while (!end);
    infile.close();
  }

  G4LevelManager* man = nullptr;
  if(vEnergy.size() >= 2) {
    man = new G4LevelManager(vEnergy.size(),vEnergy,vTime,vTimeg,vSpin,vMeta,vLevel); 
    if(fVerbose > 0) {
      G4cout << "=== Reader: new manager for Z= " << Z << " A= " << A 
	     << " Nlevels= " << vEnergy.size() << " E[0]= " 
	     << vEnergy[0]/CLHEP::MeV << " MeV   Emax= " 
	     << man->MaxLevelEnergy()/CLHEP::MeV << " MeV "
	     << " S: " <<  vEnergy.size() << " " << vTime.size()
	     << " " << vTimeg.size() << " " << vSpin.size() << " " << vLevel.size() 
	     << G4endl;
    }
  } 
  return man;
}

G4bool G4LevelReader::ReadData(std::istringstream& stream, G4double& x)
{
  stream >> x;
  return stream.fail() ? false : true;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4double& x)
{
  x = 0.0;
  for(G4int i=0; i<nbufmax; ++i) { buffer[i] = ' '; }
  G4bool okay = true;
  dataFile >> buffer;
  if(dataFile.fail()) { okay = false; }
  else { x = strtod(buffer, 0); }

  return okay;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4int& ix)
{
  ix = 0;
  for(G4int i=0; i<nbuf2; ++i) { buff2[i] = ' '; }
  G4bool okay = true;
  dataFile >> buff2;
  if(dataFile.fail()) { okay = false; }
  else { ix = atoi(buff2); }

  return okay;
}

G4bool G4LevelReader::ReadDataItem(std::istream& dataFile, G4String& x)
{
  G4bool okay = true;
  dataFile >> bufp;
  if(dataFile.fail()) { okay = false; }
  else { x = G4String(bufp, 2); }
  return okay;
}

const std::vector<G4float>* G4LevelReader::NormalizedICCProbability(G4int Z)
{
  std::vector<G4float>* vec = nullptr;
  G4int LL = 3;
  G4int M = 5;
  G4int N = 1;
  if(Z <= 27) {
    M = N = 0;
    if(Z <= 4) {
      LL = 1;
    } else if(Z <= 6) {
      LL = 2;
    } else if(Z <= 10) {
    } else if(Z <= 12) {
      M = 1;
    } else if(Z <= 17) {
      M = 2;
    } else if(Z == 18) {
      M = 3;
    } else if(Z <= 20) {
      M = 3;
      N = 1;
    } else {
      M = 4;
      N = 1;
    }
    if(LL < 3) { for(G4int i=LL+1; i<=4; ++i) { fICC[i] = 0.0; } }
    if(M < 5)  { for(G4int i=M+4;  i<=8; ++i) { fICC[i] = 0.0; } }
    if(N < 1)  { fICC[9] = 0.0; }
  }
  G4float norm = 0.0;
  for(G4int i=0; i<10; ++i) {
    norm += fICC[i]; 
    fICC[i] = norm;
  }
  if(norm > 0.0f) {
    norm = 1.0f/norm;
    vec = new std::vector<G4float>;
    G4float x;
    for(G4int i=0; i<10; ++i) {
      x = (G4float)(fICC[i]*norm);
      if(x > 0.995f) {
	vec->push_back(1.0f);
	break;
      } 
      vec->push_back(x); 
    }
    if (fVerbose > 2) {
      G4int prec = G4cout.precision(3);
      G4cout << "# InternalConv: ";
      G4int nn = vec->size();
      for(G4int i=0; i<nn; ++i) { G4cout << " " << (*vec)[i]; }
      G4cout << G4endl;
      G4cout.precision(prec);
    }
  }
  return vec;
}

const G4LevelManager* 
G4LevelReader::CreateLevelManagerNEW(G4int Z, G4int A)
{
  std::ostringstream ss;
  ss << "/correlated_gamma/z" << Z << ".a" << A;
  G4String st = G4String(ss.str());
  fFile = fDirectory + st;
  std::ifstream infile(fFile, std::ios::in);

  // file is not opened
  if (!infile.is_open()) {
    if (fVerbose > 0) {
      G4cout << " G4LevelReader: fail open file for Z= " 
	     << Z << " A= " << A 
	     << " <" << fFile << ">" <<  G4endl;
    }
    return nullptr;
  }
  if (fVerbose > 0) {
    G4cout << "G4LevelReader: open file for Z= " 
	   << Z << " A= " << A 
	   << " <" << fFile << ">" <<  G4endl;
  }
  return LevelManager(Z, A, 0, infile);
}

const G4LevelManager* 
G4LevelReader::MakeLevelManagerNEW(G4int Z, G4int A,
				   const G4String& filename)
{
  std::ifstream infile(filename, std::ios::in);

  // file is not opened
  if (!infile.is_open()) {
    if (fVerbose > 0) {
      G4cout << " G4LevelReader: fail open file for Z= " 
	     << Z << " A= " << A 
	     << " <" << filename << ">" <<  G4endl;
    }
    return nullptr;
  }
  if (fVerbose > 0) {
    G4cout << "G4LevelReader: open file for Z= " 
	   << Z << " A= " << A 
	   << " <" << filename << ">" <<  G4endl;
  }
  return LevelManager(Z, A, 0, infile);
}

const G4LevelManager* 
G4LevelReader::LevelManager(G4int Z, G4int A, G4int nlev,
			    std::ifstream& infile)
{
  G4bool allLevels = fParam->StoreAllLevels();
  G4float emax = fData->GetMaxLevelEnergy(Z, A);

  static const G4double fkev = CLHEP::keV;
  G4int nlevels = (0 == nlev) ? fLevelMax : nlev;
  if(fVerbose > 0) {
    G4cout << "## New isotope Z= " << Z << "  A= " << A;
    if(nlevels < fLevelMax) { G4cout << " Nlevels= " << nlevels; }
    G4cout << G4endl;
  }
  if(nlevels > fLevelMax) {
    fLevelMax = nlevels;
    vEnergy.resize(fLevelMax,0.0f);
    vTime.resize(fLevelMax,0.0f);
    vTimeg.resize(fLevelMax,0.0f);
    vSpin.resize(fLevelMax,0);
    vLevel.resize(fLevelMax,0);
    vMeta.resize(fLevelMax,0);
    vIndexDB.resize(fLevelMax,-1);
  }
  G4int ntrans(0), i(0), i1, i2, i3, j, k;
  G4String xf("  ");
  G4float x, x1;

  for(G4int ii=0; ii<nlevels; ++ii) {

    infile >> i1 >> xf;
    if(infile.eof()) {
      if(fVerbose > 1) { 
	G4cout << "### End of file Z= " << Z << " A= " << A 
	       << " Nlevels= " << ii << G4endl;
      }
      break;
    }
    if(!(ReadDataItem(infile,fEnergy) &&
	 ReadDataItem(infile,fTime)   &&
	 ReadDataItem(infile,fSpin)   &&
	 ReadDataItem(infile,ntrans))) {
      if(fVerbose > 1) { 
	G4cout << "### End of file Z= " << Z << " A= " << A 
	       << " Nlevels= " << ii << G4endl;
      }
      break;
    }
    fTime = std::max(fTime, 0.0); 
    fEnergy *= fkev;
    for(k=0; k<nfloting; ++k) {
      if(xf == fFloatingLevels[k]) {
	break;
      }
    }

    // if a previous level has not transitions it may be ignored
    if(0 < ii) {
      // do not store level without transitions
      if(!allLevels && 0 == k && 0 == ntrans) { continue; } 

      // protection
      if(fEnergy < vEnergy[i-1]) {
	G4cout << "### G4LevelReader: broken level " << ii 
	       << " E(MeV)= " << fEnergy << " < " << vEnergy[i-1]
	       << " for isotope Z= " << Z << " A= " 
	       << A << " level energy increased" << G4endl; 
	fEnergy = vEnergy[i-1];
      }
      // upper limit
      if(fEnergy > emax) { break; }
    }
    vEnergy[i] = (G4float)fEnergy;
    vTime[i]   = (G4float)(fTime*fTimeFactor);
    vTimeg[i]  = vTime[i];
    if(fSpin > 20.0) { fSpin = 0.0; }
    vSpin[i]   = (G4int)(fSpin + fSpin);
    vMeta[i]   = k; 
    vIndexDB[ii] = i;
    if(fVerbose > 1) {
      G4cout << "   Level #" << i1 << " E(MeV)= " << fEnergy/CLHEP::MeV
	     << "  LTime(s)= " << fTime << " 2S= " << vSpin[i]
	     << "  meta= " << vMeta[i] << " idx= " << i << " ii= " << ii 
	     << " ntr= " << ntrans << G4endl;
    }
    vLevel[i] = nullptr;
    if(ntrans > 0) {

      // there are transitions
      if(ntrans > fTransMax) {
	fTransMax = ntrans;
	vIndex.resize(fTransMax);
	vTrans.resize(fTransMax);
	vRatio.resize(fTransMax);
	vGammaCumProbability.resize(fTransMax);
	vGammaECumProbability.resize(fTransMax);
	vGammaProbability.resize(fTransMax);
	vShellProbability.resize(fTransMax);
	vMpRatio.resize(fTransMax);
      }
      fNorm1 = fNorm2 = 0.0f;
      j = 0; 
      for(G4int jj=0; jj<ntrans; ++jj) {
       
	if(!(ReadDataItem(infile,i2)        &&
	     ReadDataItem(infile,fTrEnergy) &&
	     ReadDataItem(infile,fProb)     &&
	     ReadDataItem(infile,vTrans[j]) &&
	     ReadDataItem(infile,fRatio)    &&
	     ReadDataItem(infile,fAlpha))) {
	  //infile >>i2 >> fTrEnergy >> fProb >> vTrans[j] >> fRatio >> fAlpha; 
	  //if(infile.fail()) { 
	  if(fVerbose > 0) { 
	    G4cout << "### Fail to read transition j= " << j 
		   << "  Z= " << Z << " A= " << A << G4endl; 
	  }
	  break;
	}
        if(i2 >= ii) {
	  G4cout << "### G4LevelReader: broken transition " << j 
	     << " from level " << ii << " to " << i2
	     << " for isotope Z= " << Z << " A= " 
	     << A << " - use ground level" << G4endl; 
          i2 = 0;
	}
	i3 = vIndexDB[std::abs(i2)];
        
        if(i3 >= 0) {
	  vIndex[j] = i3;
	  x = 1.0f + std::max((G4float)fAlpha,0.0f);
          x1= (G4float)fProb;
	  fNorm1 += x1;
	  fNorm2 += x*x1;
	  vGammaCumProbability[j] = fNorm1;
	  vGammaECumProbability[j]= fNorm2;
	  vGammaProbability[j] = 1.0f/x;
	  vRatio[j] = (G4float)fRatio;
	  vShellProbability[j] = nullptr;
	  if(fVerbose > 1) { 
	    fTrEnergy *= fkev;
	    G4int prec = G4cout.precision(4);
	    G4cout << "### Transition #" << j << " to level " << vIndex[j] 
                   << " i2= " << i2 <<  " Etrans(MeV)= " << fTrEnergy
		   << "  fProb= " << fProb << " MultiP= " << vTrans[j]
		   << "  fMpRatio= " << fRatio << " fAlpha= " << fAlpha 
                   << G4endl;
	    G4cout.precision(prec);
	  }
	}
	if(fAlpha > 0.0) {
	  for(k=0; k<10; ++k) {
	    //infile >> fICC[k];
	    if(!ReadDataItem(infile,fICC[k])) {
	      //if(infile.fail()) { 
	      if(fVerbose > 0) { 
		G4cout << "### Fail to read convertion coeff k= " << k 
		       << " for transition j= " << j 
		       << "  Z= " << Z << " A= " << A << G4endl; 
	      }
	      break;
	    }
	  }
	  if(i3 >= 0) { 
	    vShellProbability[j] = NormalizedICCProbability(Z);
            if(!vShellProbability[j]) { vGammaProbability[j] = 1.0f; }
	  }
	}
	if(i3 >= 0) { ++j; }
      }
      if(j > 0) {
	if(0.0f < fNorm1) { 
	  fNorm1 = 1.0f/fNorm1; 
	  vTimeg[i] *= fNorm2*fNorm1;
	  fNorm2 = 1.0f/fNorm2;
	}
	G4int nt = j - 1;      
	for(k=0; k<nt; ++k) {
	  vGammaCumProbability[k]  *= fNorm1;
	  vGammaECumProbability[k] *= fNorm2;
	  if(fVerbose > 2) {
	    G4cout << "Probabilities[" << k 
		   << "]= " << vGammaCumProbability[k]
		   << "  " << vGammaECumProbability[k] 
		   << " idxTrans= " << vIndex[k]
		   << G4endl;
	  }
	}
	vGammaCumProbability[nt]  = 1.0f;
	vGammaECumProbability[nt] = 1.0f;
	if(fVerbose > 2) {
	  G4cout << "Probabilities[" << nt << "]= " 
		 << vGammaCumProbability[nt]
		 << "  " << vGammaECumProbability[nt] 
		 << " IdxTrans= " << vIndex[nt]
		 << G4endl;
	}
	vLevel[i] = new G4NucLevel((size_t)j,
				   vIndex,
				   vTrans,
				   vGammaCumProbability,
				   vGammaECumProbability,
				   vGammaProbability,
				   vMpRatio,
				   vShellProbability);
      }
    }
    ++i;
  }
  G4LevelManager* lman = nullptr;
  if(1 < i) { 
    lman = new G4LevelManager((size_t)i,vEnergy,vTime,vTimeg,
			      vSpin,vMeta,vLevel);
    if(fVerbose > 0) {
      G4cout << "=== Reader: new manager for Z= " << Z << " A= " << A 
	     << " Nlevels= " << i << " E[0]= " 
	     << vEnergy[0]/CLHEP::MeV << " MeV   Emax= " 
	     << vEnergy[i-1]/CLHEP::MeV << " MeV "
	     << G4endl;
    }
  }
  for(G4int ii=0; ii<nlevels; ++ii) {
    vIndexDB[ii] = -1;
  }
  return lman;
}
