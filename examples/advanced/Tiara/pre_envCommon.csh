# $Id: pre_envCommon.csh,v 1.1 2004-06-09 15:04:34 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: not supported by cvs2svn $
# -------------------------------------------------------------------
# Before sourcing this script make sure you have set the 
# environment variables according to the description in README.
# -------------------------------------------------------------------
# The following environmental variables should be set before sourcing 
# envCommon.csh, in order to fix the base dirs for some packages in afs
# CLHEP_BASE_DIR should be set before sourcing envCommon.csh
#
setenv PI_BASE_DIR /afs/cern.ch/sw/geant4/dev/PI
setenv SWIG_BASE_DIR /afs/cern.ch/sw/lhcxx/specific/rh73_gcc32/PublicDomainPackages/2.0.0
setenv G4ANALYSIS_USE 1 
