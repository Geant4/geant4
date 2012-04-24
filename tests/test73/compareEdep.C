#include <vector>
#include <algorithm>
#include "Math/GoFTest.h"
#include <TFile.h>
#include <TH1F.h>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using namespace ROOT::Math;

void compareEdep( TString fname1 , TString fname2 , TString hExt="0.33", Int_t testtype = 1)
{
 //Test the Edep (binned) distributions from two files.
//Apply statistic from testtype
// 1=Anderson-Darling
  TFile * f1 = TFile::Open( fname1 );
  TFile * f2 = TFile::Open( fname2 );
    if ( f1->IsZombie() ) { cerr<<"Cannot open file: "<<fname1<<endl; return; }
    if ( f2->IsZombie() ) { cerr<<"Cannot open file: "<<fname2<<endl; return; }
    TH1F* h1 = (TH1F*)f1->Get("EDepPerTrack"+hExt);
    TH1F* h2 = (TH1F*)f2->Get("EDepPerTrack"+hExt);
    if ( !h1 || !h2 ) { cerr<<"Cannot get histogram: EDepPerTrack"<<hExt<<endl; return; }
    
    //Create vector of entries for binned distribution
    //Get histos as vectors (of doubles) and sort them
    std::vector<Double_t> v1;
    std::vector<Double_t> v2;
    for ( Int_t i = h1->GetXaxis()->GetFirst(); i<h1->GetXaxis()->GetLast() ; ++i) {
      Double_t bcenter = h1->GetXaxis()->GetBinCenter( i );
      Double_t bcont   = h1->GetBinContent( i );
      if ( bcont > 0 ) {
	for ( Int_t ctr = 0 ; ctr < std::floor(bcont) ; ++ctr ) {
	  v1.push_back( bcenter );
	}
      }
    }
    for ( Int_t i = h2->GetXaxis()->GetFirst(); i<h2->GetXaxis()->GetLast() ; ++i) {
      Double_t bcenter = h2->GetXaxis()->GetBinCenter( i );
      Double_t bcont   = h2->GetBinContent( i );
      if ( bcont > 0 ) {
	for ( Int_t ctr = 0 ; ctr < std::floor(bcont) ; ++ctr ) {
	  v2.push_back( bcenter );
	}
      }
    }
    std::sort( v1.begin() , v1.end() );
    std::sort( v2.begin() , v2.end() );
    //Get back an array of doubles
    Double_t array1[v1.size()];
    Double_t array2[v2.size()];
    for ( size_t i=0;i<v1.size();i++) array1[i]=v1[i];
    for ( size_t i=0;i<v2.size();i++) array2[i]=v2[i];
    //for ( Int_t i = 0 ; i < v1.size() ; ++i ) cout<<i<<" "<<array1[i]<<endl;
    //for ( Int_t i = 0 ; i < v2.size() ; ++i ) cout<<array2[i]<<endl;
    GoFTest gof( v1.size() , array1 , v2.size(), array2 );
    Double_t pvalue , tstat;
    if ( testtype == 1 ) //AndersonDaling
      {
	gof.AndersonDarling2SamplesTest(pvalue,tstat);
      }
    
    cout<<"For histograms with Extension: "<<hExt<<" t-value="<<tstat<<" p-value="<<pvalue<<endl;
    f1->Close();
    f2->Close();
    return;
}
