BEGIN{start=0;}
{
  if ($1=="Compiling") {
    s=$2;
    if(start==1) start = 0;
  }
  if ($1=="Error") {
    if(start==0) {
      print "******************************";
      print "Error_in : " s;
      print "******************************";
      start = 1;
    }
  }
  if ($1=="Warning") {
    if(start==0) {
      print "******************************";
      print "Warning_in : " s;
      print "******************************";
      start = 1;
    }
  }
  if(start==1) print $0;
}
