////////////////////////////////////////////////////////////////////////////////
//
#include "MLPHSAnalyser.hh"

#include "globals.hh"
#include <iomanip>
#include <vector>

#include "MLGeometryConstruction.hh"
#include "G4RunManager.hh"

#ifndef USEHBOOK
#include "MLVersion.hh"
#endif
////////////////////////////////////////////////////////////////////////////////
//
MLPHSAnalyser::MLPHSAnalyser (MLAnalysisManager *anMan,
  MLGeometryConstruction *det)
{
//
//
// Set default conditions.
//
  Nb       = -9999;
  DefaultPEdg();
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
MLPHSAnalyser::~MLPHSAnalyser ()
{
  for (size_t i = 0; i < Histo.size(); i++) {
    delete Histo[i];
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::SetPType (G4String sval)
{
  if (sval == "lin" || sval == "LIN" || sval == "linear" || sval == "LINEAR") {
    pType = LIN;
    pEdg.clear();

  } else if (sval == "log" || sval == "LOG" || sval == "logarithmic" ||
    sval == "LOGARITHMIC") {
    pType = LOG;
    pEdg.clear();

  } else if (sval == "arb" || sval == "ARB" || sval == "arbitrary" ||
    sval == "ARBITRARY") {
    pType = ARB;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::SetPHSMax (G4double max)
{
  if (max < minphs) {
    G4cerr <<G4endl;
    G4cerr <<"SetPHSMax : maximum energy " <<max/keV <<" keV less than minimum "
           <<minphs/keV <<" keV." <<G4endl;
    G4cerr <<"(Try changing minimum first?)" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    maxphs = max;
//
//
// If an arbitrary binning scheme is used, check for and delete any which are
// above maxphs.
//
    if (pType == ARB) {
      size_t i = 0;
      while (i < pEdg.size()) {
        if (pEdg[i] >= maxphs) {
          pEdg.erase(pEdg.begin()+i);
        } else {
          i++;
        }
      }
      pNBin = pEdg.size()+1;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::SetPHSMin (G4double min)
{
  if (min > maxphs) {
    G4cerr <<G4endl;
    G4cerr <<"SetPHSMin : minimum energy " <<min/keV
           <<" keV greater than maximum " <<maxphs/keV <<" keV." <<G4endl;
    G4cerr <<"(Try changing maximum first?)" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;

  } else {
    minphs = min;
//
//
// If an arbitrary binning scheme is used, check for and delete any which are
// below mineng.
//
    if (pType == ARB) {
      size_t i = 0;
      while (i < pEdg.size()) {
        if (pEdg[i] <= minphs) {
          pEdg.erase(pEdg.begin()+i);
        } else {
          i++;
        }
      }
      pNBin = pEdg.size()+1;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::SetPNBin (G4int ival)
{
  if (pType == ARB) {
    G4cerr <<G4endl;
    G4cerr <<"SetPNbin : number of bins cannot be set for arbitrary scale."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    pNBin = ival;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::DefaultPEdg ()
{
  pType  = LIN;
  maxphs = 1000.0*keV;
  minphs = 0.0*keV;
  pNBin  = 100;
  pEdg.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::AddPEdg (G4double dval)
{
//
//
// A bin edge can only be added if the energy-scale is arbitrary.
//
  if (pType != ARB) {
    G4cerr <<"AddPEdg: an bin-edge cannot be specified unless the energy "
           <<"scale is arbitrary" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
//
//
// Check that the new bin edge is within the maximum and minimum limits
// defined by the user for the PHS scale.
//
  } else if (dval > minphs && dval < maxphs) {
//
//
// Check that the edge does not correspond to an existing edge.
//
    size_t i;
    for (i = 0; i < pEdg.size(); i++) {
      if (pEdg[i] == dval) break;
    }
    if (i == pEdg.size()) {
      pEdg.push_back(dval);
      sort (pEdg.begin(),pEdg.end());
      pNBin = pEdg.size()+1;
    } else {
      G4cerr <<"AddPEdg: a bin edge already exists with that energy."
             <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    }
  } else {
    G4cerr <<"AddPEdg: the energy specified is outside the PHS range specified."
           <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::DeletePEdg (G4double dval)
{
//
//
// A bin edge can only be deleted if the energy-scale is arbitrary.
//
  if (pType != ARB) {
    G4cerr <<"DeletePEdg: an bin-edge cannot be deleted unless the energy "
           <<"scale is arbitrary" <<G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  } else {
    if (find(pEdg.begin(),pEdg.end(),dval) != pEdg.end() )  
      pEdg.erase(find(pEdg.begin(),pEdg.end(),dval));
    pNBin = pEdg.size()+1;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::ListPEdg ()
{
  if (pType == LIN) {
    G4cout <<"The PHS energy range is linear over:"
           <<G4endl;
    G4cout <<minphs/keV <<" keV  to " <<maxphs/keV <<" keV" << G4endl;
    G4cout <<"and comprises " <<pNBin <<" bins."  <<G4endl;
    G4cout <<G4endl;

  } else if (pType == LOG) {
    G4cout <<"The PHS energy range is logarithmic over:"
           <<G4endl;
    G4cout <<minphs/keV <<" keV  to " <<maxphs/keV <<" keV" << G4endl;
    G4cout <<"and comprises " <<pNBin <<" bins."  <<G4endl;
    G4cout <<G4endl;

  } else if (pType == ARB) {
    G4cout <<"The PHS binning scheme is arbitrary between "
           <<minphs/keV <<" keV to " <<maxphs/keV <<" keV" <<G4endl;
    G4cout <<"with bin edges at:"  <<G4endl;
    for (size_t i = 0; i < pEdg.size(); i++)
      G4cout <<"     " <<pEdg[i]/keV <<" keV " <<G4endl;
    G4cout <<G4endl;

  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLPHSAnalyser::GetNBlocks ()
{
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
G4int MLPHSAnalyser::GetAddedNBlocks (G4int i)
{
  Nb = i + geometry->GetNbOfELayers();
  return Nb;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::BeginOfRunAction (G4double nFactor)
{
  factor = nFactor * cm2;

  G4double xbin[1001] = {0.};
  G4int nbin = 0;
  char title[80] ;
#ifdef USEHBOOK
  HbookHistogram* aHisto = NULL;
#else
  MLHisto1D* aHisto = NULL;
#endif

  if (pType == LIN) {
    G4double dx = (maxphs - minphs) / G4double(pNBin);
    for ( G4int i = 0; i < pNBin+1; i++) {
      xbin[i] = (minphs + i*dx) / keV;
    }
    nbin = pNBin;
  } else if (pType == LOG) {
    G4double dx = (log10(maxphs) - log10(minphs)) / G4double(pNBin);
    for (G4int i = 0; i < pNBin+1; i++) {
      G4double dd = log10(minphs) + i*dx;
      xbin[i] = pow(10,dd) / keV;
    }
    nbin = pNBin;
  } else if (pType == ARB) {
    xbin[0] = minphs;
    nbin    = 0;
    for (size_t i = 0; i < pEdg.size(); i++) {
      nbin++;
      xbin[nbin] = pEdg[i] / keV;
    }
    nbin++;
    xbin[nbin] = maxphs / keV;
    pEdg.push_back(minphs);
    pEdg.push_back(maxphs);
    sort (pEdg.begin(),pEdg.end());
  }

  G4int NbOfELayers = geometry->GetNbOfELayers();

  for (G4int i = 0; i < NbOfELayers; i++) {
    G4int k = geometry->GetELayerIdx(i);
    sprintf(title,"Pulse Height Spectrum in layer %i",k+1);

#ifdef USEHBOOK
    aHisto = new HbookHistogram(title,nbin,xbin);
    Histo.push_back(aHisto);
#else
    G4String unitx = "keV";
    G4String unity = "counts/cm2/bin";
    aHisto = new MLHisto1D(title,unitx, unity,xbin,nbin+1,RIGHT);
    Histo.push_back(aHisto);

#endif
  }
}
////////////////////////////////////////////////////////////////////////////////
//
void MLPHSAnalyser::FillHbook (G4int id, G4double x, G4double, G4double w)
{
  //
#ifdef USEHBOOK
  Histo[id]->accumulate(x/keV,0.,w);
#else
  Histo[id]->fill(x/keV,w);
#endif
}
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
void MLPHSAnalyser::NormaliseOutput (G4double nevent)
{
//
//
// Loop through layers in geometry, particles and angles.
//
  for (size_t i = 0; i < (size_t) geometry->GetNbOfELayers(); i++) {
    Histo[i]->SetNormalisation(factor/nevent);
  }
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
CSVofstream & operator << (CSVofstream &CSVFile, MLPHSAnalyser &q)
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
  if (geometry->GetNbOfELayers() == 0) return CSVFile;
#ifdef SPENVIS
  G4int Nc =  4;            // Number of comment lines.
#else
  G4int Nc =  2;            // Number of comment lines.
#endif
  G4int Nm =  1;            // Number of meta-variable lines.
  G4int Na =  0;            // Number of annotation lines.
  G4int Nv =  5;            // Number of variable lines.
  G4int Nd =  5;            // Number of data columns.
  G4int Nl = q.GetPNBin();  // Number of data lines.
  G4int Nh = Nc + Nm + Nv + 1;
                         // Total number of header records(incl. current one.
//
//
// Loop through layers in geometry for which dose information is to be output.
//
  for (G4int i = 0; i < geometry->GetNbOfELayers(); i++) {
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
            <<std::setw(6) <<q.GetNBlocks()-i-1 <<G4endl;
//
//
// Output comment lines - basically the version number of MULASSIS and of
// SPENVIS if applicable.
//
    CSVFile <<"'PULSE-HEIGHT ANALYSIS'" <<G4endl;
    CSVFile <<"'MULASSIS " <<codeVersion.version <<"'" <<G4endl;
#ifdef SPENVIS
    CSVFile <<"'SPENVIS " <<codeVersion.SPENVISver <<"'" <<G4endl;
    CSVFile <<"'" <<codeVersion.implementation <<"'" <<G4endl;
#endif
//
//
// Output meta-variable lines.
//
    CSVFile <<"'PHS_LYR',  1, " <<geometry->GetELayerIdx(i)+1 <<G4endl;
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
    CSVFile <<"'Value','events/cm2/bin',  1,'Number of events per unit detector area'"
            <<G4endl;
    CSVFile <<"'Error','events/cm2/bin',  1,'Error in number of events'"
            <<G4endl;


    CSVFile << *(q.GetHisto(i));

    if (q.GetNBlocks()-i > 1) {
      CSVFile <<"'End of Block'" <<G4endl;
    } else {
      CSVFile <<"'End of File'" <<G4endl;
    }
  }

  return CSVFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////
//
#ifndef USEHBOOK
RPTofstream & operator << (RPTofstream &RPTFile, MLPHSAnalyser &)
{
  RPTFile << "-------------------------------------------------------------" << G4endl;
  RPTFile << "PHS Analysis:" << G4endl;
  RPTFile << "-------------------------------------------------------------" << G4endl;
  RPTFile << " Outputed in the .csv file " <<G4endl;
  return RPTFile;
}
#endif
////////////////////////////////////////////////////////////////////////////////
