
echo "Start of run for " $TARGET

rm $TARGET/res8.0.log

$G4MY/cascade_test $TARGET/run_fast.mac  >& $TARGET/res8.0.log

mv *.paw   $TARGET/

echo $TARGET " is done!"

