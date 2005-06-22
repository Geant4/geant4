
echo "Start of run for " $TARGET

#rm $TARGET/res.log

$G4MY/cascade_test $TARGET/run_fast.mac  >& $TARGET/res7.1.log

mv *.paw   $TARGET/

echo $TARGET " is done!"

