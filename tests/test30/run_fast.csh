
echo "Start of run for " $TARGET

$G4MY/test30 $TARGET/run_fast.mac  > /dev/null

mv kin.paw   $TARGET/

echo $TARGET " is done!"

