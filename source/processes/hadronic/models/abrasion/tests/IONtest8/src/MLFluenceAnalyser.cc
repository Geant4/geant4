////////////////////////////////////////////////////////////////////////////////
//
#include "MLFluenceAnalyser.hh"

#include <iomanip>
#include "G4RunManager.hh"
#include "MLGeometryConstruction.hh"
#include "MLAnalysisManager.hh"
#ifndef USEHBOOK
#include "MLVersion.hh"
#endif

#include <strstream>
////////////////////////////////////////////////////////////////////////////////
//
MLFluenceAnalyser::MLFluenceAnalyser (MLAnalysisManager* anMan,
  MLGeometryConstruction *det)
{
//
//
// Set default conditions.
//
  Nb          = -9999;
  fluenceUnit = "cm2";
  DefaultEEdg();
  DefaultAEdg();
  sPart.clear();
  //  sPart.push_back("proton");
  //  sPart.push_back("neutron");
  divideByCosT = true;
  factor.clear();
//
//
// Set Analysis Manager.
//
  analysisManager = anMan;
//
//
// Get the geometry.
//
  geometry = det;
}
////////////////////////////////////////////////////////////////////////////////
//
MLFluenceAnalyser::~MLFluenceAnalyser ()
{}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetEType (G4String sval)
{
  if (sval == "lin" || sval == "LIN" || sval == "linear" || sval == "LINEAR") {
    eType = LIN;
    eEdg.clear();

  } else if (sval == "log" || sval == "LOG" || sval == "logarithmic" ||
    sval == "LOGARITHMIC") {
    eType = LOG;
    eEdg.clear();

  } else if (sval == "arb" || sval == "ARB" || sval == "arbitrary" ||
    sval == "ARBITRARY") {
    eType = ARB;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetEngMax (G4double max)
{
  if (max < mineng) {
    G4cerr <<G4endl;
    G4cerr <<"SetEngMax : maximum energy " <<max/keV <<" keV less than minimum "
           <<mineng/keV <<" keV." <<G4endl;
    G4cerr <<"(Try changing minimum first?)" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    maxeng = max;
//
//
// If an arbitrary binning scheme is used, check for and delete any which are
// above maxeng.
//
    if (eType == ARB && eEdg.size() > 0) {
      size_t i = 0;
      while (i <= eEdg.size()) {
        if (eEdg[i] >= maxeng) {
          eEdg.erase(eEdg.begin()+i);
        } else {
          i++;
        }
      }
      eNBin = eEdg.size()+1;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetEngMin (G4double min)
{
  if (min > maxeng) {
    G4cerr <<G4endl;
    G4cerr <<"SetEngMin : minimum energy " <<min/keV
           <<" keV greater than maximum " <<maxeng/keV <<" keV." <<G4endl;
    G4cerr <<"(Try changing maximum first?)" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    mineng = min;
//
//
// If an arbitrary binning scheme is used, check for and delete any which are
// below mineng.
//
    if (eType == ARB && eEdg.size() > 0) {
      size_t i = 0;
      while (i < eEdg.size()) {
        if (eEdg[i] <= mineng) {
          eEdg.erase(eEdg.begin()+i);
        } else {
          i++;
        }
      eNBin = eEdg.size()+1;
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::DefaultEEdg ()
{
  eType  = LOG;
  maxeng = 1.0*GeV;
  mineng = 1.0*keV;
  eNBin  = 60;
  eEdg.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetENBin (G4int ival)
{
  if (eType == ARB) {
    G4cerr <<G4endl;
    G4cerr <<"SetENbin : number of bins cannot be set for arbitrary scale."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    eNBin = ival;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::AddEEdg (G4double dval)
{
//
//
// A bin edge can only be added if the energy-scale is arbitrary.
//
  if (eType != ARB) {
    G4cerr <<"AddEEdg: an bin-edge cannot be specified unless the energy "
           <<"scale is arbitrary" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
//
//
// Check that the new bin edge is within the maximum and minimum limits
// defined by the user for the fluence energy scale.
//
  } else if (dval > mineng && dval < maxeng) {
//
//
// Check that the edge does not correspond to an existing edge.
//
    size_t i;
    for (i = 0; i < eEdg.size(); i++) {
      if (eEdg[i] == dval) break;
    }
    if (i == eEdg.size()) {
      eEdg.push_back(dval);
      sort(eEdg.begin(),eEdg.end());
      eNBin = eEdg.size()+1;
    } else {
      G4cerr <<"AddEEdg: a bin edge already exists with that energy!"
             <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    }
  } else {
    G4cerr <<"AddEEdg: the energy specified is outside the energy range specified!"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::DeleteEEdg (G4double dval)
{
//
//
// A bin edge can only be deleted if the energy-scale is arbitrary.
//
  if (eType != ARB) {
    G4cerr <<"DeleteEEdg: an bin-edge cannot be deleted unless the energy "
           <<"scale is arbitrary" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    if (find(eEdg.begin(),eEdg.end(),dval) != eEdg.end() ) 
      eEdg.erase(find(eEdg.begin(),eEdg.end(),dval));
    eNBin = eEdg.size()+1;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::ListEEdg ()
{
  if (eType == LIN) {
    G4cout <<"The energy range is linear over:"
           <<G4endl;
    G4cout <<mineng/keV <<" keV  to " <<maxeng/keV <<" keV" << G4endl;
    G4cout <<"and comprises " <<eNBin <<" bins."  <<G4endl;
    G4cout <<G4endl;

  } else if (eType == LOG) {
    G4cout <<"The energy range is logarithmic over:"
           <<G4endl;
    G4cout <<mineng/keV <<" keV  to " <<maxeng/keV <<" keV" << G4endl;
    G4cout <<"and comprises " <<eNBin <<" bins."  <<G4endl;
    G4cout <<G4endl;

  } else if (eType == ARB) {
    G4cout <<"The energy binning scheme is arbitrary between "
           <<mineng/keV <<" keV to " <<maxeng/keV <<" keV" <<G4endl;
    G4cout <<"with bin edges at:"  <<G4endl;
    for (size_t i = 0; i < eEdg.size(); i++)
      G4cout <<"     " <<eEdg[i]/keV <<" keV " <<G4endl;
    G4cout <<G4endl;

  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetAType (G4String sval)
{
  if (sval == "lin" || sval == "LIN" || sval == "linear" || sval == "LINEAR") {
    aType = LIN;
    aEdg.clear();

  } else if (sval == "arb" || sval == "ARB" || sval == "arbitrary" ||
    sval == "ARBITRARY") {
    aType = ARB;

  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetAngMax (G4double max)
{
  if (max < minang) {
    G4cerr <<G4endl;
    G4cerr <<"SetAngMax : maximum angle " <<max/deg <<" deg less than minimum "
           <<minang/deg <<" deg." <<G4endl;
    G4cerr <<"(Try changing minimum first?)" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    maxang = max;
//
//
// If an arbitrary binning scheme is used, check for and delete any which are
// above maxang.
//
    if (aType == ARB && aEdg.size() > 0) {
      size_t i = 0;
      while (i < aEdg.size()) {
        if (aEdg[i] >= maxang) {
          aEdg.erase(aEdg.begin()+i);
        } else {
          i++;
        }
      }
      aNBin = aEdg.size()+1;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetAngMin (G4double min)
{
  if (min > maxang) {
    G4cerr <<G4endl;
    G4cerr <<"SetAngMin : minimum angle " <<min/deg
           <<" deg greater than maximum " <<maxang/deg <<" deg." <<G4endl;
    G4cerr <<"(Try changing maximum first?)" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    minang = min;
//
//
// If an arbitrary binning scheme is used, check for and delete any which are
// below minang.
//
    if (aType == ARB && aEdg.size() > 0 ) {
      size_t i = 0;
      while (i < aEdg.size()) {
        if (aEdg[i] <= minang) {
          aEdg.erase(aEdg.begin()+i);
        } else {
          i++;
        }
      }
      aNBin = aEdg.size()+1;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::DefaultAEdg ()
{
  aType  = LIN;
  maxang = 180.*deg;
  minang = 0.*deg;
  aNBin  = 2 ;
  aEdg.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetANBin (G4int ival)
{
  if (aType == ARB) {
    G4cerr <<G4endl;
    G4cerr <<"SetANbin : number of bins cannot be set for arbitrary scale."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    aNBin = ival;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::AddAEdg (G4double dval)
{
//
//
// A bin edge can only be added if the angle-scale is arbitrary.
//
  if (aType != ARB) {
    G4cerr <<"AddAEdg: an bin-edge cannot be specified unless the angle "
           <<"scale is arbitrary" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
//
//
// Check that the new bin edge is within the maximum and minimum limits
// defined by the user for the fluence angle scale.
//
  } else if (dval > minang && dval < maxang) {
//
//
// Check that the edge does not correspond to an existing edge.
//
    size_t i;
    for (i = 0; i < aEdg.size(); i++) {
      if (aEdg[i] == dval) break;
    }
    if (i == aEdg.size()) {
      aEdg.push_back(dval);
      sort (aEdg.begin(),aEdg.end());
      aNBin = aEdg.size()+1;
    } else {
      G4cerr <<"AddAEdg: a bin edge already exists with that angle!"
             <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    }
  } else {
    G4cerr <<"AddEEdg: the angle specified is outside the angle range specified!"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::DeleteAEdg (G4double dval)
{
//
//
// A bin edge can only be deleted if the angle-scale is arbitrary.
//
  if (aType != ARB) {
    G4cerr <<"DeleteAEdg: an bin-edge cannot be deleted unless the angle "
           <<"scale is arbitrary" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    if (find(aEdg.begin(),aEdg.end(),dval) != aEdg.end() )  
      aEdg.erase(find(aEdg.begin(),aEdg.end(),dval));
    aNBin = aEdg.size()+1;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::ListAEdg ()
{
  if (aType == LIN) {
    G4cout <<"The angle range is linear over:"
           <<G4endl;
    G4cout <<minang/deg <<" deg  to " <<maxang/deg <<" deg" << G4endl;
    G4cout <<"and comprises " <<aNBin <<" bins."  <<G4endl;
    G4cout <<G4endl;

  } else if (aType == ARB) {
    G4cout <<"The angle binning scheme is arbitrary between "
           <<minang/deg <<" deg to " <<maxang/deg <<" deg" <<G4endl;
    G4cout <<"with bin edges at:"  <<G4endl;
    for (size_t i = 0; i < aEdg.size(); i++)
      G4cout <<"     " <<aEdg[i]/deg <<" deg " <<G4endl;
    G4cout <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::AddSPart (G4String sval)
{
  if (find(sPart.begin(),sPart.end(),sval) == sPart.end()) {
    sPart.push_back(sval);
  } else {
    G4cerr <<"Particle " <<sval <<" already selected for fluence analysis"
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::DeleteSPart (G4String sval)
{
  if (find(sPart.begin(),sPart.end(),sval) != sPart.end() )
    sPart.erase(find(sPart.begin(),sPart.end(),sval));
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::ListSPart ()
{
  G4cout <<"The following particles are selected for fluence analysis: "
         <<G4endl;
  for (size_t i = 0; i < sPart.size(); i++)
    G4cout <<"     " <<sPart[i] <<G4endl;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::SetFluenceType (G4String type)
{
  if (type == "OMNI" || type == "omni" || type == "OMNIDIRECTIONAL" ||
    type == "omnidirectional") {
    divideByCosT = true;

  } else if (type == "PLANAR" || type == "planar" || type == "BOUNDARY" ||
    type == "boundary") {
    divideByCosT = false;

  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4String MLFluenceAnalyser::GetFluenceType ()
{
  if (divideByCosT) return "OMNIDIRECTIONAL";
  else return "PLANAR";
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLFluenceAnalyser::GetNBlocks ()
{
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLFluenceAnalyser::GetAddedNBlocks (G4int i)
{
  Nb = i + sPart.size() * (aNBin+1) * (geometry->GetNbOfFLayers());
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLFluenceAnalyser::GetHistIdx (G4int i, G4int j, G4int k)
{
//
//
// Warning - there are no checks done to make sure this arguments supplied are
// valid.  It's therefore efficient, but watchout for errors.
//
  return (i*6 + j) * (1+aNBin) + k;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::BeginOfRunAction (G4double nFactor)
{
  factor.clear();
  G4int iLayer      = 0;
  G4int NbOfFLayers = geometry->GetNbOfFLayers();
  for (G4int i = 0; i < NbOfFLayers; i++) {
    iLayer = geometry->GetFLayerIdx(i);
//
//
// Normalise the results depending upon the units of fluence used in the output.
//
    if (fluenceUnit == "cm2") {
      factor.push_back (nFactor * cm2);
    } else if (fluenceUnit == "m2") {
      factor.push_back (nFactor * m2);
    }
//
//
// If the geometry is spherical, normalise according to the surface area of the
// outer-most sphere compared with the surface of at which the fluence
// measurement is required.
//
    if (geometry->GetShape() == SPHERE) {
      G4double rtot = geometry->GetLayerRadius(0); // the outer layer radius
      G4double rout = geometry->GetLayerRadius(iLayer);
      factor[i]    *= rtot*rtot/rout/rout;
    }
  }

  G4double xbin[1001];
  G4int nbin = 0;
#ifdef USEHBOOK
  HbookHistogram* aHisto = NULL;
#else
  MLHisto1D* aHisto = NULL;
#endif

  sort(sPart.begin(),sPart.end());

  if (eType == LIN) {
    G4double dx = (maxeng - mineng) / G4double(eNBin);
    for (G4int i = 0; i < eNBin+1; i++) xbin[i] = (mineng + i*dx) / keV;
    nbin = eNBin;
  } else if (eType == LOG) {
    G4double dx = (log10(maxeng) - log10(mineng)) / G4double(eNBin);
    for (G4int i = 0; i < eNBin+1; i++) {
      G4double dd = log10(mineng) + i*dx;
      xbin[i]     = pow(10,dd) / keV;
    }
    nbin = eNBin;
  } else if (eType == ARB) {
    if (!(binary_search(eEdg.begin(),eEdg.end(),mineng))) eEdg.push_back(mineng);
    if (!(binary_search(eEdg.begin(),eEdg.end(),maxeng))) eEdg.push_back(maxeng);
    sort(eEdg.begin(),eEdg.end());
    nbin = 0;
    for (size_t i = 0; i < eEdg.size(); i++) {
      xbin[nbin] = eEdg[i] / keV;
      nbin++;
    }
  }
  // now angular distribution
  G4int i;
  if (aType == LIN) {
    aEdg.clear();
    G4double dx = (maxang - minang) / G4double(aNBin);
    for ( i = 0; i < aNBin+1; i++) aEdg.push_back(minang + i*dx);
  } else if (aType == ARB) {
    if (!(binary_search(aEdg.begin(),aEdg.end(),minang))) aEdg.push_back(minang);
    if (!(binary_search(aEdg.begin(),aEdg.end(),maxang))) aEdg.push_back(maxang);
    sort (aEdg.begin(),aEdg.end());
    aNBin = aEdg.size()-1;
  }

#ifdef USEHBOOK
  theHbookManager.SetID(1000);
#endif
  G4String particleTypes[6] = {"Proton", "Neutron", "Electron", "Gamma",
    "Charged muon", "Charged pion"};
  G4int j        = 0;
  G4int k        = 0;
  G4int l        = 0;
  G4String unitx = "keV";
  G4String unity = "particles/bin/" + fluenceUnit;

  for (i = 0; i < NbOfFLayers; i++) {
    k = geometry->GetFLayerIdx(i);

    for (l = 0; l < 6; l++) {
      char title[80] = {' '};
      std::ostrstream os(title,80);
      os <<particleTypes[l]
         <<" spectrum crossing the boundary between layers " <<k+1
         <<" and " <<k+2 <<".";
#ifdef USEHBOOK
      aHisto = new HbookHistogram(title,nbin,xbin);
      Histo.push_back(aHisto);
#else
      aHisto = new MLHisto1D(title,unitx,unity,xbin,nbin+1,RIGHT);
      Histo.push_back(aHisto);
#endif
      for (j = 0; j < aNBin; j++) {
        char title[80] = {' '};
        std::ostrstream os(title,80);
        os <<particleTypes[l]
           <<"s crossing at angles between " <<aEdg[j]/deg
           <<" - " <<aEdg[j+1]/deg <<" degrees at boundary "
           <<k+1 <<" - " <<k+2 <<".";
#ifdef USEHBOOK
        aHisto = new HbookHistogram(title,nbin,xbin);
        Histo.push_back(aHisto);
#else
        aHisto = new MLHisto1D(title,unitx,unity,xbin,nbin+1,RIGHT);
        Histo.push_back(aHisto);
#endif
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::TallyFluenceEvent (G4int layer, G4int part,
  G4double theta, G4double energy, G4double weight)
{
  G4int histd = 0; //proton
  if (part == 2112) { // neutron
    histd = 1;
  } else if ( part == 11) { // electron
    histd = 2;
  } else if ( part == 22) { // gamma
    histd = 3;
  } else if ( part == 13 || part == -13) { // muon
    histd = 4;
  } else if ( part == 211 || part == -211) { // pion
    histd = 5;
  }

  size_t i = 0;
  while (i < aEdg.size()) {
    if (theta < aEdg[i]) break;
    i++;
  }
  G4int hista = G4int(i);

  G4int j = 0;
  while (j < geometry->GetNbOfFLayers()) {
    if (layer == geometry->GetFLayerIdx(j)) break;
    j++;
  }

  G4double w = weight;
  if (divideByCosT) {
    G4double cosT = fabs(cos(theta*rad));
    if (cosT < 1.0e-4) cosT = 1.0e-4;
    w /= cosT;
  }
  G4int hist = GetHistIdx(j,histd,0);
  FillHbook(hist, energy, 0., w);
  hist = GetHistIdx(j,histd,hista);
  FillHbook(hist, energy, 0., w);

}
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::FillHbook (G4int id, G4double x, G4double , G4double w)
{
#ifdef USEHBOOK
  Histo[id]->accumulate(x/keV,0.,w);
#else
  Histo[id]->fill(x/keV,w);
#endif
}
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
void MLFluenceAnalyser::NormaliseOutput (G4double nevent)
{
  G4int hist    = 0;
//
//
// Loop through layers in geometry, particles and angles.
//
  for (G4int i = 0; i < geometry->GetNbOfFLayers(); i++) {
    for (G4int j = 0; j < 6; j++) {
      for (G4int l = 0; l < aNBin+1; l++) {
        hist = GetHistIdx(i,j,l);
        Histo[hist]->SetNormalisation(factor[i]/nevent);
      }
    }
  }
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
void MLFluenceAnalyser::EndOfRunAction ()
{
  G4int hist = 0;
  for (G4int i = 0; i < geometry->GetNbOfFLayers(); i++) {
    for (G4int j = 0; j < 6; j++) {
      for (G4int l = 0; l < aNBin+1; l++) {
        hist = GetHistIdx(i,j,l);
        delete Histo[hist];
      }
    }
  }
  Histo.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
CSVofstream & operator << (CSVofstream &CSVFile, MLFluenceAnalyser &q)
{
//
//
// Initialise some variables.
//
  MLVersion codeVersion = MLVersion();
  G4RunManager *runManager = G4RunManager::GetRunManager();
  MLGeometryConstruction *geometry =
    (MLGeometryConstruction*)(runManager->GetUserDetectorConstruction());
//
//
// Check whether any layers were selected for PHS analysis - if not don't
// output anything.
//
  if (geometry->GetNbOfFLayers() == 0) return CSVFile;
  G4int hist  = 0;
  G4int block = 0;

#ifdef SPENVIS
  G4int Nc =  4; // Number of comment lines.
#else
  G4int Nc =  2; // Number of comment lines.
#endif
  G4int Nm =  4; // Number of meta-variable lines.
  G4int Na =  0; // Number of annotation lines.
  G4int Nv =  5; // Number of variable lines.
  G4int Nd =  5; // Number of data columns.
  G4int Nl = q.GetENBin(); // Number of data lines.
  G4int Nh = Nc + Nm + Nv + 1;
                  // Total number of header records(incl. current one.
//
//
// Loop through layers in geometry.
//
  for (G4int i = 0; i < geometry->GetNbOfFLayers(); i++) {
//
//
// Loop through particles.
//
    G4String part = "";
    for (G4int j = 0; j < 6; j++) {
      G4String part;
      switch (j) {
      case 0 :
        part = "proton";
        break;
      case 1 :
        part = "neutron";
        break;
      case 2 :
        part = "e-";
        break;
      case 3 :
        part = "gamma";
        break;
      case 4 :
        part = "muon";
        break;
      case 5 :
        part = "pion";
        break;
      }
//
//
// Loop through angle bins.
//
      for (G4int l = 0; l < q.GetANBin()+1; l++) { // loop over angle bins
        if (q.IsParticleOutput(part)) {
//
//
// Output header record for CSV file.
//
          CSVFile <<"'*'"
                  <<std::setw(6) <<Nh <<","
                  <<std::setw(6) <<Nc <<","
                  <<std::setw(6) <<Nm <<","
                  <<std::setw(6) <<Na <<","
                  <<std::setw(6) <<Nv <<","
                  <<std::setw(6) <<Nd <<","
                  <<std::setw(6) <<Nl <<","
                  <<std::setw(6) <<q.GetNBlocks()-block-1 <<G4endl;
//
//
// Output comment lines - basically the version number of MULASSIS and of
// SPENVIS if applicable.
//
          if (l == 0) {
            CSVFile <<"'" <<q.GetFluenceType()
                    <<" FLUENCE ANALYSIS'"<<G4endl;
          } else {
            CSVFile <<"'" <<q.GetFluenceType()
                    <<" FLUENCE ANALYSIS AS A FUNCTION OF ANGLE'" <<G4endl;
          }
          CSVFile <<"'MULASSIS " <<codeVersion.version <<"'" <<G4endl;
#ifdef SPENVIS
          CSVFile <<"'SPENVIS " <<codeVersion.SPENVISver <<"'" <<G4endl;
          CSVFile <<codeVersion.implementation <<G4endl;
#endif
//
//
// Output meta-variable lines.
//
          CSVFile <<"'FLU_BDY',  1, " <<geometry->GetFLayerIdx(i)+1
                  <<", ' '" <<G4endl;
          CSVFile <<"'FLU_PAR', -1, '" <<part <<"', ' '" <<G4endl;
          if (l == 0) {
            CSVFile <<"'FLU_AMN',  1,   0., 'deg'" <<G4endl;
            CSVFile <<"'FLU_AMX',  1, 180., 'deg'" <<G4endl;
          } else {
            CSVFile <<std::setiosflags(std::ios::fixed);
            CSVFile <<std::setprecision(5);
            CSVFile <<"'FLU_AMN',  1, " <<q.GetAEdg(l-1)/deg
                    <<", 'deg'" <<G4endl;
            CSVFile <<"'FLU_AMX',  1, " <<q.GetAEdg(l)/deg
                    <<", 'deg'" <<G4endl;
            CSVFile <<std::resetiosflags(std::ios::fixed);
          }
//
//
// Output annotation lines.  (CURRENTLY NONE)
//

//
//
// Output variable lines.
//
          CSVFile <<"'Elo','keV',  1,'Lower edge of energy bin'"
                  <<G4endl;
          CSVFile <<"'Eup','keV',  1,'Upper edge of energy bin'"
                  <<G4endl;
          CSVFile <<"'Emean','keV',  1,'Mean energy of bin'"
                  <<G4endl;
          CSVFile <<"'Value','particles/" <<q.GetFluenceUnit()
                  <<"/bin',  1,'Fluence/flux'" <<G4endl;
          CSVFile <<"'Error','particles/" <<q.GetFluenceUnit()
                  <<"/bin',  1,'Error in fluence/flux'" <<G4endl;

          hist = q.GetHistIdx(i,j,l);
          CSVFile << *(q.GetHisto(hist));
          if (q.GetNBlocks()-block > 1) {
            CSVFile <<"'End of Block'" <<G4endl;
          } else {
            CSVFile <<"'End of File'" <<G4endl;
          }
          block++;
        }
      }
    }
  }
  return CSVFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
RPTofstream & operator << (RPTofstream &RPTFile, MLFluenceAnalyser &)
//
//
// Definition of the insertion operator << to provide the histogram output to
// RPTofstream.
//
{
  RPTFile << "-------------------------------------------------------------" << G4endl;
  RPTFile << "Fluence Analysis:" << G4endl;
  RPTFile << "-------------------------------------------------------------" << G4endl;
  RPTFile << " Outputed in the .csv file " <<G4endl;
  return RPTFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////
