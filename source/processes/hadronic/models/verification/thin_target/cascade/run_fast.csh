
echo "Start of run for " $TARGET

rm $TARGET/res.log

$G4MY/test30 $TARGET/run_fast.mac  >& $TARGET/res.log

mv *.paw   $TARGET/

echo $TARGET " is done!"

