
echo "Start of run for " $TARGET

$G4MY/cascade10gev $TARGET/run.mac  >& /dev/null

mv bic.paw   $TARGET/
mv ber.paw   $TARGET/
mv qgsp.paw   $TARGET/
mv qgsc.paw   $TARGET/

echo $TARGET " is done!"

