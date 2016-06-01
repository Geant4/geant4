// GGE (Geant4 Geometry Editor)
// Tetsuyay Yamada


import java.io.Serializable;
import java.util.*;
import java.lang.*;
import java.lang.Math;

class BREPItem extends SolidItem {
  static final int PCONE = 0, PGON = 1;
  static final String solidName[] = {"PolyConeSegment","PolyGonSegment"};
  static final String brepName[] = {"PCone","Polyhedra"};
  static final String solidConstr[] = {"PolyConeSegment","PolyGonSegment"};
  static final String paraName[][] = {
          {"StartPhiAngle", "DeltaPhiAngle","Z of a Section", "OuterRadius",
                                                              "InnerRadius"},
          {"StartPhiAngle", "DeltaPhiAngle", "NumberOfSides","Z of a Section",
                                                "OuterRadius","InnerRadius"}
     };

  double phi[];
  private double boundZ, boundRmax, boundBox; 
  private double values[];
  int nSides, nZ;
  String sAngUnit, dAngUnit;
//, lenUnit;
  Vector zVector;
  String dzA, rminA, rmaxA;
  String sphi, dphi;
  StringBuffer dz, rmin, rmax;
  BREPItem(int solidType){
    super(solidType);
    values = new double[paraName[solidType].length];
    phi = new double[2];
    zVector = new Vector();
    isEmpty = true;
  }
  public String toString(){
    return solidName[solidType];
  }
  int getParaCount(){
    return paraName[solidType].length;
  }
  double getValues(int pos){
    return values[pos];
  }

  String getParaName(int pos){
    return paraName[solidType][pos];
  }

  String getCPP(String name){
    nZ = zVector.size();
    dz = new StringBuffer("");
    rmin = new StringBuffer("");
    rmax = new StringBuffer("");
    for(int s=0; s<zVector.size(); s++){
      if(s!=zVector.size()-1){
        Vector array = (Vector)zVector.elementAt(s);
        dz.append(array.elementAt(0)+"*"+array.elementAt(3)+", ");
        rmin.append(array.elementAt(2)+"*"+array.elementAt(3)+", ");
        rmax.append(array.elementAt(1)+"*"+array.elementAt(3)+", ");
      }else if(s==zVector.size()-1){
        Vector lastArray = (Vector)zVector.elementAt(s);
        dz.append(lastArray.elementAt(0)+"*"+lastArray.elementAt(3));
        rmin.append(lastArray.elementAt(2)+"*"+lastArray.elementAt(3));
        rmax.append(lastArray.elementAt(1)+"*"+lastArray.elementAt(3));
      }
    }
    dzA = dz.toString();
    rminA = rmin.toString();
    rmaxA = rmax.toString();
    StringBuffer cpp = new StringBuffer("G4double DzArray"+name+"   [] = {"+
             dzA+"};\nG4double RminArray"+name+" [] = {"+rminA+
             "};\nG4double RmaxArray"+name+" [] = {"+rmaxA+"};\n\n");

         cpp.append("G4BREPSolid"+brepName[solidType]+
             " *solid"+name+" = new G4BREPSolid"+
             brepName[solidType]+"(\"solid"+name+"\", \t //its name\n");

         cpp.append("\t\t "+phi[0]+"*"+sAngUnit+", \t\t //its start angle\n");
         cpp.append("\t\t "+phi[1]+"*"+dAngUnit+",\t\t //its opening angle\n");
         if (solidType == PGON) cpp.append("\t\t "+nSides+",\t\t //its sides\n");
         cpp.append("\t\t "+nZ + ", \t\t //its nZ\n\t\t DzArray"+
                    name+"[0], \t\t //z start \n\t\t DzArray"+name+
                    ", \t\t //z value \n\t\t RminArray"+name+
                    ", \t\t //rmin \n\t\t RmaxArray"+name+" ); \t\t //rmax\n");
    return cpp.toString();
  }

  public String getPrim(){
// BoundingBox dZ of SUM, Rmax of MAX
   boundZ = 0.0;
   boundRmax = 0.0;

// sphi
   if(sAngUnit.equals("rad")){
     sphi = Double.toString(phi[0]);
   }else if(sAngUnit.equals("mrad")){
     sphi = Double.toString(phi[0]/1000);
   }else if(sAngUnit.equals("deg")){
     sphi = Double.toString(phi[0]*2*3.14159265358979323/360);
   }
// dphi
   if(dAngUnit.equals("rad")){
     dphi = Double.toString(phi[1]);
   }else if(dAngUnit.equals("mrad")){
     dphi = Double.toString(phi[1]/1000);
   }else if(dAngUnit.equals("deg")){
     dphi = Double.toString(phi[1]*2*3.14159265358979323/360);
   }
// nZ
   nZ = zVector.size();
// dz, rmin, rmax --start--
   dz = new StringBuffer("");
   rmin = new StringBuffer("");
   rmax = new StringBuffer("");
    for(int s=0; s<zVector.size(); s++){
     Vector array = (Vector)zVector.elementAt(s);
     if(array.elementAt(3).equals("km")){        // km
      try{
        Double dzk = Double.valueOf(array.elementAt(0).toString());
        double dk = dzk.doubleValue();
        Double rMik = Double.valueOf(array.elementAt(2).toString());
        double rmik = rMik.doubleValue();
        Double rMak = Double.valueOf(array.elementAt(1).toString());
        double rmak = rMak.doubleValue();
        Double kM = Double.valueOf("1000000");        
        double km = kM.doubleValue();
        dz.append(dk*km+" ");
        boundZ += dk*km;
        rmin.append(rmik*km+" ");
        rmax.append(rmak*km+" ");
        boundRmax = Math.max(boundRmax, rmak*km);        
      }catch(NumberFormatException e){
        System.out.println("");
      }
     }else if(array.elementAt(3).equals("m")){  // m
      try{
        Double dzm = Double.valueOf(array.elementAt(0).toString());
        double dm = dzm.doubleValue();
        Double rMim = Double.valueOf(array.elementAt(2).toString());
        double rmim = rMim.doubleValue();
        Double rMam = Double.valueOf(array.elementAt(1).toString());
        double rmam = rMam.doubleValue();
        double m = Double.valueOf("1000").doubleValue();
        dz.append(dm*m+" ");
        boundZ += dm*m;
        rmin.append(rmim*m+" ");
        rmax.append(rmam*m+" ");
        boundRmax = Math.max(boundRmax, rmam*m);
      }catch(NumberFormatException e){
        System.out.println("");
      }       
     }else if(array.elementAt(3).equals("cm")){  // cm
      try{
        Double dzc = Double.valueOf(array.elementAt(0).toString());
        double dc = dzc.doubleValue();
        Double rMic = Double.valueOf(array.elementAt(2).toString());
        double rmic = rMic.doubleValue();
        Double rMac = Double.valueOf(array.elementAt(1).toString());
        double rmac = rMac.doubleValue();
        double cm = Double.valueOf("10").doubleValue();
        dz.append(dc*cm+" ");
        boundZ += dc*cm;
        rmin.append(rmic*cm+" ");
        rmax.append(rmac*cm+" ");
        boundRmax = Math.max(boundRmax, rmac*cm);
      }catch(NumberFormatException e){
        System.out.println("");
      } 
     }else if(array.elementAt(3).equals("mm")){  // mm
      try{
        Double dzmm = Double.valueOf(array.elementAt(0).toString());
        double dmm = dzmm.doubleValue();
        Double rMimm = Double.valueOf(array.elementAt(2).toString());
        double rmimm = rMimm.doubleValue();
        Double rMamm = Double.valueOf(array.elementAt(1).toString());
        double rmamm = rMamm.doubleValue();
        dz.append(dmm+" ");
        boundZ += dmm;
        rmin.append(rmimm+" ");
        rmax.append(rmamm+" ");
        boundRmax = Math.max(boundRmax, rmamm);
      }catch(NumberFormatException e){
        System.out.println("");
      }
     }else if(array.elementAt(3).equals("micrometer")){  // micrometer
      try{
         Double dzz = Double.valueOf(array.elementAt(0).toString());
         double d = dzz.doubleValue();
         Double rMi = Double.valueOf(array.elementAt(2).toString());
         double rmi = rMi.doubleValue();
         Double rMa = Double.valueOf(array.elementAt(1).toString());
         double rma = rMa.doubleValue();
         double mic = Double.valueOf("1000").doubleValue();
         dz.append(d/mic+" ");
         boundZ += d/mic;
         rmin.append(rmi/mic+" ");
         rmax.append(rma/mic+" ");
         boundRmax = Math.max(boundRmax, rma/mic);
      }catch(NumberFormatException e){
         System.out.println("");
      }
     }else if(array.elementAt(3).equals("namoometer")){  // nanoometer
      try{
         Double dzn = Double.valueOf(array.elementAt(0).toString());
         double dn = dzn.doubleValue();
         Double rMina = Double.valueOf(array.elementAt(2).toString());
         double rmina = rMina.doubleValue();
         Double rMana = Double.valueOf(array.elementAt(1).toString());
         double rmana = rMana.doubleValue();
         double nan = Double.valueOf("1000000").doubleValue();
         dz.append(dn/nan+" ");
         boundZ += dn/nan;
         rmin.append(rmina/nan+" ");
         rmax.append(rmana/nan+" ");
         boundRmax = Math.max(boundRmax, rmana/nan);
      }catch(NumberFormatException e){
         System.out.println("");
      }
     }else if(array.elementAt(3).equals("fermi")){  // fermi
      try{
         Double dzf = Double.valueOf(array.elementAt(0).toString());
         double df = dzf.doubleValue();
         Double rMif = Double.valueOf(array.elementAt(2).toString());
         double rmif = rMif.doubleValue();
         Double rMaf = Double.valueOf(array.elementAt(1).toString());
         double rmaf = rMaf.doubleValue();
         double fer = Double.valueOf("1000000000").doubleValue();
         dz.append(df/fer+" ");
         boundZ += df/fer;
         rmin.append(rmif/fer+" ");
         rmax.append(rmaf/fer+" ");
         boundRmax = Math.max(boundRmax, rmaf/fer);
      }catch(NumberFormatException e){
         System.out.println("");
      }
     }
   }
// BoundingBox


//   System.out.println(boundZ);
//   System.out.println(boundRmax);
//   System.out.println(boundBox);
   boundBox = Math.sqrt(Math.pow(boundZ, 2)+Math.pow(boundRmax, 2))*2;
   dzA = dz.toString();
   rminA = rmin.toString();
   rmaxA = rmax.toString();




   StringBuffer prim = new StringBuffer("#################################\n");
    prim.append("###         GGESolid.prim      ###\n");
    prim.append("#################################\n");
    prim.append("#G4.PRIM-FORMAT-2.4\n\n");
    prim.append("/BoundingBox  -"+boundBox+" -"+boundBox+" -"+boundBox+
                " "+boundBox+" "+boundBox+" "+boundBox+"\n");
    prim.append("!SetCamera\n!OpenDevice\n!BeginModeling\n\n");
    prim.append("#"+solidName[solidType]+"\n");
       if(solidName[solidType].equals("PolyConeSegment")){
            prim.append("/PolyCone    ");
       }else if(solidName[solidType].equals("PolyGonSegment")){
            prim.append("/PolyGon    ");
       }

       prim.append(sphi+"  "+dphi+"  ");   
       if(solidName[solidType].equals("PolyGonSegment")){
            prim.append(nSides+"  ");
       }
    
       prim.append(nZ+"  "+dzA+"  "+rminA+"  "+rmaxA+"  ");
       

    prim.append("\n\n!EndModeling\n!DrawAll\n!CloseDevice");
    return prim.toString();
  }
}












