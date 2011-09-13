#
#
foreach fl ( *newout* *newerr* )
  set oldfl=`echo $fl | perl -pe 's#new##;' `
  cp -p $fl $oldfl
end
