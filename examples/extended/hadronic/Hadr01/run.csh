#/bin/csh

setenv HISTODIR 20070516
mkdir -p $HISTODIR

setenv PHYSLIST LHEP
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QGSP
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QGSP_EMV
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QGSP_EMX
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QGSC
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QGSP_BERT
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QGSP_BIC
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

setenv PHYSLIST QBBC
$G4MY/hadr01 p_pb.in >& $HISTODIR/$PHYSLIST.out

echo "Done!"
#
