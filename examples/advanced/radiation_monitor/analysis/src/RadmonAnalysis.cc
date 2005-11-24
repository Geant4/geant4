//
// File name:     RadmonAnalysis.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysis.cc,v 1.1 2005-11-24 02:33:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonAnalysis.hh"

#include "RadmonVDataAnalysisFactory.hh"
#include "RadmonVDataAnalysis.hh"
#include "RadmonVAnalysisLayout.hh"

                                                RadmonAnalysis :: RadmonAnalysis(RadmonVAnalysisLayout * layout, RadmonVDataAnalysisFactory * factory, AIDA::IAnalysisFactory * analysis)
:
 analysisLayout(layout),
 dataFactory(factory),
 analysisFactory(analysis),
 treeFactory(0),
 tree(0),
 changed(true)
{
 if (analysisLayout==0)
  G4Exception("RadmonAnalysis::RadmonAnalysis: layout==0.");

 if (dataFactory==0)
  G4Exception("RadmonAnalysis::RadmonAnalysis: factory==0.");
 
 if (analysisFactory==0)
  G4Exception("RadmonAnalysis::RadmonAnalysis: analysis==0.");
 
 analysisLayout->AttachObserver(this);
 
 treeFactory=analysisFactory->createTreeFactory();
}



                                                RadmonAnalysis :: ~RadmonAnalysis()
{
 analysisLayout->DetachObserver(this);
 
 Destruct();
 
 delete treeFactory;
 delete analysisFactory;
 delete dataFactory;
}





void                                            RadmonAnalysis :: OnLayoutChange(void)
{
 changed=true;
}





void                                            RadmonAnalysis :: OnBeginOfEvent(const G4Event * /* event */)
{
 if (changed)
 {
  Destruct();
 
  // TO DO
 
  changed=false;
 }
}





void                                            RadmonAnalysis :: OnEndOfEvent(const G4Event * /* event */)
{
}





void                                            RadmonAnalysis :: Destruct(void)
{
}
