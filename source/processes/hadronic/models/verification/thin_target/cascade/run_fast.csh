
echo "Start of run for " $TARGET

rm kin.paw
rm ber.paw

$G4MY/test30 $TARGET/run_fast.mac  >& /dev/null

mv kin*.paw   $TARGET/
mv ber.paw   $TARGET/

echo $TARGET " is done!"

