/*
 *7.5
 */

import com.sun.java.swing.*;
import com.sun.java.swing.border.*;
import com.sun.java.swing.event.*;

import com.sun.java.swing.preview.filechooser.*;
import ExampleFileFilter;
import com.sun.java.swing.preview.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

class VolumesFrame extends AbstractGGEFrame implements ActionListener {
  MaterialFrame matFrame; 
  SourceFrame sf;
  DAWNSourceFrame dawnSFrame;
  LogicalPanel logPanel;
  VisAtbPanel visPanel;
  PhysicalPanel phyPanel;
  private ReplicasFrame repFrame;   //8.20 Replicas 
  private RepeatedVolume repeatVol; //8.20 RepeatedVolume
  private RepeatedRotateVolume repeatRotVol; //9.14 Rotate
  private SingleVolume singleVol;   //8.21 Single Positioned Volume
  private SlicingAxialFrame slicingFrame; //1998 10.26
  private JMenuItem load, append, save, clear, exit;

//, mat;
//Makes C++ code
  private JMenuItem CPP, closeCPP;
//Preview
//  private JMenuItem dawn , david;
//Elements Table
  private JMenuItem tOpen ,tClose; 
//  private JCheckBoxMenuItem sourceWin, materialWin;
  JLabel label;
  private JFileChooser loadGGE, saveGGE;
  private String fileName;
  static String fileDir;
  VolumesFrame(){
    super("Volumes");
    getContentPane().setLayout(new BorderLayout());
    matFrame = new MaterialFrame();
    sf = new SourceFrame();
    dawnSFrame = new DAWNSourceFrame();
    repFrame = new ReplicasFrame();  //8.21 Replicas 
    slicingFrame = new SlicingAxialFrame(); //1998 10.26 SlicingAxial
    repeatVol = new RepeatedVolume(); // 8.21 Repeated Volume  
    repeatRotVol = new RepeatedRotateVolume(); // 9.14 Rotate
    singleVol = new SingleVolume();   // 8.21 SinglePositionedVolume
    logPanel = new LogicalPanel(this);
    visPanel = new VisAtbPanel(this);
    phyPanel = new PhysicalPanel(this);
    JPanel panel1 = new JPanel(new BorderLayout());

    JPanel panel2 = new JPanel(new BorderLayout());

    panel1.add("Center", logPanel);
    panel1.add("East", visPanel);
    panel2.add("Center", phyPanel);
    JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, panel1, panel2);
    getContentPane().add("Center", splitPane);
    setJMenuBar( createMenubar() );
    setSize(780, 450);
//    pack();
    setVisible(true);
//    setBackground(Color.green);
//    setForeground(Color.black);
  }
  private JMenuBar createMenubar(){
    JMenu file ,mkCode;     
    JMenuBar mb = new JMenuBar();
    file = new JMenu("File");
    file.add(load = new JMenuItem("Load a Volumes File"));
    file.add(append = new JMenuItem("Append a Volumes File"));
    file.add(save = new JMenuItem("Save Volumes"));
    file.add(clear = new JMenuItem("Clear Volumes"));
    file.addSeparator();
    file.add(exit = new JMenuItem("Exit"));
    load.addActionListener(this);
    append.addActionListener(this);
    save.addActionListener(this);
    clear.addActionListener(this);
    exit.addActionListener(this);
    mb.add(file);
    mkCode = new JMenu("MakeSource");
    mkCode.add(CPP = new JMenuItem("Make C++ code"));
    CPP.addActionListener(this);
    mb.add(mkCode);
    return mb;
  }
  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if( o == load ){GGEfileLoad();}
    if( o == append ){GGEfileAppend();}
    if( o == save ){GGEfileSave();}
    if( o == clear ){ volClear();}
    if( o == exit ){ System.exit(0); }
    if( o == CPP){
//for debug      System.out.println(makeCPP());
        sf.setVisible(true);
        sf.editor.setText(makeCPP());
    }
  }
//------single, replica, repeated Volume -----open method---------
  public void singleOpen(){
    singleVol.setVisible(true);
  }
  public void replicaOpen(){
    repFrame.setVisible(true);
  }
  public void slicingOpen(){  // 1998 10.26 SlicingAxial
    slicingFrame.setVisible(true);
  }
  public void repeatOpen(){
    repeatVol.setVisible(true);
  }
  public void repeatRotateOpen(){
    repeatRotVol.setVisible(true);
  }
  public void parameOpen(){;}

//Create C++ -------------------------------------------------

  String makeCPP(){
	StringBuffer smatcpp = new StringBuffer();
        smatcpp.append(headerCPP());
	smatcpp.append(matFrame.getCPP());
        smatcpp.append(visPanel.getCPP());
	smatcpp.append(logPanel.getCPP());
        smatcpp.append(phyPanel.getCPP());
        smatcpp.append(singleVol.getCPP());
        smatcpp.append(repeatVol.getCPP());
        smatcpp.append(repeatRotVol.getCPP());
        smatcpp.append(repFrame.getCPP());
        smatcpp.append(slicingFrame.getCPP()); //1998 10.26 
        smatcpp.append(getCPP());
        return smatcpp.toString();
  }
  private synchronized void resume(){
    notify();
  }

  String getCPP(){
    StringBuffer cpp = new StringBuffer("\n// return the physical World\n\n");
    if(singleVol.singleTable.getRowCount()==0){
      return cpp.toString();
    }else{
      String mother;
      for (int i=0;i<singleVol.singleDataModel.getRowCount();i++){
        mother = ((String)singleVol.singleDataModel.getValueAt(i,3));
        if(!mother.equals("NULL")){
           //next single data
        }else if(mother.equals("NULL")){
          cpp.append("\n return physical"
            +singleVol.singleDataModel.getValueAt(i,1)+";\n}");            
        }
      }
    }
    return cpp.toString();
  }

//Create HeaderFile ----------------------------------------------------

  String headerCPP(){
    StringBuffer cpp = new StringBuffer("//------HeaderFile-\n ");
    cpp.append("#include \"MyDetectorConstruction.hh\"\n\n");

    cpp.append("#include \"G4UnitsTable.hh\"\n\n");   // 8.17

//    cpp.append("#include \"MyCalorimeterSD.cc\"\n");
//    cpp.append("#include \"MyCalorimeterHit.cc\"\n");
//    cpp.append("#include \"MyCalorimeterHitsCollection.cc\"\n\n");
    cpp.append("#include \"G4VUserDetectorConstruction.hh\"\n\n");
    cpp.append("#include \"globals.hh\"\n");
    cpp.append("#include \"G4Material.hh\"\n");
    cpp.append("#include \"G4MaterialTable.hh\"\n");
    cpp.append("#include \"G4Element.hh\"\n");
    cpp.append("#include \"G4ElementTable.hh\"\n");
//   G4Solid

//   solid selected

    SolidItem solidItem;
    Object solidData;
    int a=0,b=0,c=0,d=0,e=0,f=0,g=0,h=0,i=0,j=0;
    for (int s=0; s<logPanel.logDataModel.getRowCount();s++){
      solidItem = ((SolidItem)logPanel.logDataModel.getValueAt(s,1));
      if(solidItem instanceof CSGItem){
        solidData = new CSGItem(solidItem.getSolidType());

         if(solidData.toString()=="Box"){
            a++;
         }else if(solidData.toString()=="TubeSegment"){
            b++;
         }else if(solidData.toString()=="ConeSegment"){
            c++;
         }else if(solidData.toString()=="SymTrapezoid"){
            d++;
         }else if(solidData.toString()=="SphereSegment"){
            e++;
         }else if(solidData.toString()=="Parallelepiped"){
            f++;
         }else if(solidData.toString()=="TorusSegment"){
            g++;
         }else if(solidData.toString()=="Hype"){
            h++;
         }
      }else if(solidItem instanceof BREPItem){
        solidData = new BREPItem(solidItem.getSolidType());
        if(solidData.toString()=="PolyConeSegment"){
         i++;
        }else if(solidData.toString()=="PolyGonSegment"){
         j++;
        }      

      }
    }
    if(a!=0)cpp.append("#include \"G4Box.hh\"\n");
    if(b!=0)cpp.append("#include \"G4Tubs.hh\"\n");
    if(c!=0)cpp.append("#include \"G4Cons.hh\"\n");
    if(d!=0)cpp.append("#include \"G4Trd.hh\"\n");
    if(e!=0)cpp.append("#include \"G4Sphere.hh\"\n");
    if(f!=0)cpp.append("#include \"G4Para.hh\"\n");
    if(g!=0)cpp.append("#include \"G4Torus.hh\"\n");
    if(h!=0)cpp.append("#include \"G4Hype.hh\"\n");    
    if(i!=0)cpp.append("#include \"G4BREPSolidPCone.hh\"\n");
    if(j!=0)cpp.append("#include \"G4BREPSolidPolyhedra.hh\"\n");
    cpp.append("#include \"G4LogicalVolume.hh\"\n");
    cpp.append("#include \"G4ThreeVector.hh\"\n");
    cpp.append("#include \"G4PVPlacement.hh\"\n");
    cpp.append("#include \"G4PVReplica.hh\"\n");
    cpp.append("#include \"G4SDManager.hh\"\n");
    cpp.append("#include \"G4VisAttributes.hh\"\n");
    cpp.append("#include \"G4Colour.hh\"\n\n");
    cpp.append("MyDetectorConstruction::MyDetectorConstruction( )\n{ ; }\n");
    cpp.append("MyDetectorConstruction::~MyDetectorConstruction( )\n{ ; }\n");
    cpp.append("G4VPhysicalVolume* MyDetectorConstruction::Construct( )\n{\n");
    return cpp.toString();
  }
//Detector Clear-------------------------------------------------

  void volClear(){
      logPanel.logClear();
      visPanel.visClear();
      singleVol.singleClear();
      repFrame.replicasClear();
      slicingFrame.slicingClear(); //1998 10.26 SlicingAxial
      repeatVol.repeatedClear();
      repeatRotVol.repeatedRotClear();
  }

//Detector load--------------------------------------------------
  void GGEfileLoad(){
      JFileChooser loadGGE = new JFileChooser(".");
      ExampleFileFilter dtFile = new ExampleFileFilter("g4dt","DetectorSource");
      loadGGE.setFileFilter(dtFile);
      loadGGE.setDialogTitle("Load Detector");
      loadGGE.showOpenDialog(this);
      loadGGE.setMultiSelectionEnabled(true);
      File fl = loadGGE.getSelectedFile();

      if (fl.isFile()) {
        fileName = fl.getPath();

        try{
          FileInputStream filein = new FileInputStream(fileName);
          ObjectInputStream objin =new ObjectInputStream(filein);
          Vector detectorData = (Vector)objin.readObject();
          Vector logData =  (Vector)detectorData.elementAt(0);
          Vector visData =  (Vector)detectorData.elementAt(1);
          Vector replicaData = (Vector)detectorData.elementAt(2);
          Vector repeatData = (Vector)detectorData.elementAt(3);
          Vector repeatRotData = (Vector)detectorData.elementAt(4);
          Vector singleData = (Vector)detectorData.elementAt(5);
          Vector matSData = (Vector)detectorData.elementAt(6);
          Vector matCData = (Vector)detectorData.elementAt(7);
          Vector slicingData = (Vector)detectorData.elementAt(8);//1998 10.26  
// Logical
          logPanel.logClear();
          int logVectorCount = logData.size();
          for(int i=0; i<logVectorCount; i++){
           Vector lt = (Vector)logData.elementAt(i);
           logPanel.logDataModel.addRow(lt);
          }
          logPanel.repaint();
// VisAtb
          visPanel.visClear();
          int visVectorCount = visData.size();
          for(int i=0; i<visVectorCount; i++){
           Object vis = (Object)visData.elementAt(i);
           visPanel.listModel.addElement(vis);
          }
          visPanel.repaint();
// Replicas
          repFrame.replicasClear();
          int replicaVectorCount = replicaData.size();
          for(int i=0; i<replicaVectorCount; i++){
           Vector replica = (Vector)replicaData.elementAt(i);
           repFrame.repDataModel.addRow(replica);
          }
          repFrame.repaint();
// Repeated Volume
          repeatVol.repeatedClear();
          int repeatVectorCount = repeatData.size();
          for(int i=0; i<repeatVectorCount; i++){
           Vector repeat = (Vector)repeatData.elementAt(i);
           repeatVol.pDataModel.addRow(repeat);
          }
          repeatVol.repaint();

// RepeatedRotate Volume
          repeatRotVol.repeatedRotClear();
          int repeatRotVectorCount = repeatRotData.size();
          for(int i=0; i<repeatRotVectorCount; i++){
           Vector repeatRot = (Vector)repeatRotData.elementAt(i);
           repeatRotVol.pDataModel.addRow(repeatRot);
          }
          repeatRotVol.repaint();

// Single Positioned Volume
          singleVol.singleClear();
          int singleVectorCount = singleData.size();
          for(int i=0; i<singleVectorCount; i++){
           Vector single = (Vector)singleData.elementAt(i);
           singleVol.singleDataModel.addRow(single);
          }
          singleVol.repaint();
// Used Scratch Materials
          matFrame.MaterialfileClear();
          int msVectorCount = matSData.size();
          for(int i=0; i<msVectorCount; i++){
           Vector scratch = (Vector)matSData.elementAt(i);
           matFrame.msDataModel.addRow(scratch);
          }
// Used Combi Materials
          int mcVectorCount = matCData.size();
          for(int i=0; i<mcVectorCount; i++){
           Vector combi = (Vector)matCData.elementAt(i);
           matFrame.mcDataModel.addRow(combi);
          }
// SlicingAxial
          slicingFrame.slicingClear();
          int slicingVectorCount = slicingData.size();
          for(int i=0; i<slicingVectorCount; i++){
           Vector slicing = (Vector)slicingData.elementAt(i);
           slicingFrame.slicingDataModel.addRow(slicing);
          }
          slicingFrame.repaint();


          objin.close();
        }catch(Exception e){
          System.out.println(e.toString());
        }
      }else{
        System.out.println("error="+fl.toString()+"isNotaFile");
      }
      return;
  }

//Detector Append--------------------------------------------------
  void GGEfileAppend(){
      JFileChooser loadGGE = new JFileChooser(".");
      ExampleFileFilter dtFile = 
                         new ExampleFileFilter("g4dt","DetectorSource");
      loadGGE.setFileFilter(dtFile);
      loadGGE.setDialogTitle("Load Detector");
      loadGGE.showOpenDialog(this);
      loadGGE.setMultiSelectionEnabled(true);
      File fl = loadGGE.getSelectedFile();

      if (fl.isFile()) {
        fileName = fl.getPath();

        try{
          FileInputStream filein = new FileInputStream(fileName);
          ObjectInputStream objin =new ObjectInputStream(filein);
          Vector detectorData = (Vector)objin.readObject();
          Vector logData =  (Vector)detectorData.elementAt(0);
          Vector visData =  (Vector)detectorData.elementAt(1);
          Vector replicaData = (Vector)detectorData.elementAt(2);
          Vector repeatData = (Vector)detectorData.elementAt(3);
          Vector repeatRotData = (Vector)detectorData.elementAt(4);
          Vector singleData = (Vector)detectorData.elementAt(5);
          Vector matSData = (Vector)detectorData.elementAt(6);
          Vector matCData = (Vector)detectorData.elementAt(7);
          Vector slicingData = (Vector)detectorData.elementAt(8);     

// Logical
          int logVectorCount = logData.size();
          for(int i=0; i<logVectorCount; i++){
           Vector lt = (Vector)logData.elementAt(i);
           logPanel.logDataModel.addRow(lt);
          }
          logPanel.repaint();
// VisAtb
          int visVectorCount = visData.size();
          for(int i=0; i<visVectorCount; i++){
           Object vis = (Object)visData.elementAt(i);
           visPanel.listModel.addElement(vis);
          }
          visPanel.repaint();
// Replicas
          int replicaVectorCount = replicaData.size();
          for(int i=0; i<replicaVectorCount; i++){
           Vector replica = (Vector)replicaData.elementAt(i);
           repFrame.repDataModel.addRow(replica);
          }
          repFrame.repaint();
// Repeated Volume
          int repeatVectorCount = repeatData.size();
          for(int i=0; i<repeatVectorCount; i++){
           Vector repeat = (Vector)repeatData.elementAt(i);
           repeatVol.pDataModel.addRow(repeat);
          }
          repeatVol.repaint();
// RepeatedRotate Volume
          int repeatRotVectorCount = repeatRotData.size();
          for(int i=0; i<repeatRotVectorCount; i++){
           Vector repeatRot = (Vector)repeatRotData.elementAt(i);
           repeatVol.pDataModel.addRow(repeatRot);
          }
          repeatRotVol.repaint();

// Single Positioned Volume
          int singleVectorCount = singleData.size();
          for(int i=0; i<singleVectorCount; i++){
           Vector single = (Vector)singleData.elementAt(i);
           singleVol.singleDataModel.addRow(single);
          }
          singleVol.repaint();
// Used Scratch Materials
          int msVectorCount = matSData.size();
          for(int i=0; i<msVectorCount; i++){
           Vector scratch = (Vector)matSData.elementAt(i);
           matFrame.msDataModel.addRow(scratch);
          }
// Used Combi Materials
          int mcVectorCount = matCData.size();
          for(int i=0; i<mcVectorCount; i++){
           Vector combi = (Vector)matCData.elementAt(i);
           matFrame.mcDataModel.addRow(combi);
          }
// SlicingAxial
          int slicingVectorCount = slicingData.size();
          for(int i=0; i<slicingVectorCount; i++){
           Vector slicing = (Vector)slicingData.elementAt(i);
           slicingFrame.slicingDataModel.addRow(slicing);
          }
          slicingFrame.repaint();


          objin.close();
        }catch(Exception e){
          System.out.println(e.toString());
        }
      }else{
        System.out.println("error="+fl.toString()+"isNotaFile");
      }
      return;
  }


//Detector Save----------------------------------------------------
  void GGEfileSave(){
    JFileChooser saveGGE = new JFileChooser(".");
    ExampleFileFilter dtFile = new ExampleFileFilter("g4dt","DetectorSource");
    saveGGE.setFileFilter(dtFile);
    saveGGE.setDialogTitle("Save Detector");
    if ( saveGGE.showSaveDialog(this) == -1 ) return;
    File fs = saveGGE.getSelectedFile();
    fileName = fs.getPath();
   if (fileName != null){
     try{
      FileOutputStream fileout = new FileOutputStream(fileName);
      ObjectOutputStream objout = new ObjectOutputStream(fileout);
      Vector detectorVector = new Vector();
// Logical Save
      detectorVector.addElement((Object)logPanel.logDataModel.getDataVector());

// VisAtb Save
      Vector visVector = new Vector();
      for (int i=0; i<visPanel.listModel.getSize(); i++){
        visVector.addElement((Object)visPanel.listModel.elementAt(i));
      }
      detectorVector.addElement((Object)visVector);
// Replicas Save
      detectorVector.addElement((Object)repFrame.repDataModel.getDataVector());
// Repeated Volume Save
      detectorVector.addElement((Object)repeatVol.pDataModel.getDataVector());
// RepeatedRotate Volume Save
    detectorVector.addElement((Object)repeatRotVol.pDataModel.getDataVector());
// Single Positioned Volume
  detectorVector.addElement((Object)singleVol.singleDataModel.getDataVector());

// Used ScratchMaterials ----- Start --------------------------------
      Vector msdata = matFrame.msDataModel.getDataVector();
      Vector msRowData;
      Vector msUsed = new Vector();
// Density must be checked before generating C++ code. 8.21 not yet!!

      for (int i=0; i<msdata.size(); i++){
        msRowData = (Vector)msdata.elementAt(i);
        if(!msRowData.elementAt(0).equals("Used")){
//Use G4 Units
        }else if(msRowData.elementAt(0).equals("Used")){
           msUsed.addElement((Vector)msdata.elementAt(i));
        }
      }
      detectorVector.addElement((Object)msUsed);
// Used ScratchMaterials ----- End --------------------------------

// Used CombiMaterial -------- Start ------------------------------
      Vector mcdata = matFrame.mcDataModel.getDataVector();
      Vector mcRowData; 
      Vector mcUsed = new Vector();
// Density must be checked before generating C++ code. 8.21 not yet!!

      for (int i=0; i<mcdata.size(); i++){
        mcRowData = (Vector)mcdata.elementAt(i);
        if(!mcRowData.elementAt(0).equals("Used")){
//Use G4 Units
        }else if(mcRowData.elementAt(0).equals("Used")){
           mcUsed.addElement((Vector)mcdata.elementAt(i));
        }
      }
      detectorVector.addElement((Object)mcUsed);
// Used ScratchMaterials ----- End --------------------------------

// SlicingAxial Save
      detectorVector.addElement((Object)slicingFrame.slicingDataModel.getDataVector());


      objout.writeObject(detectorVector);
      objout.flush();
      fileout.close();
     }catch(IOException e){ System.out.println(e.getMessage()); }
      return;
   }
   if (fileName == null){   return;}
   
  }
//
//Used Material----------------------------------------------------
//
  public void okCommand(){
    Object mtclear[] = new Object[1];
    mtclear[0] = "";
    for (int i=0; i<matFrame.msTable.getRowCount(); i++){
       matFrame.msTable.setValueAt(mtclear[0],i,0);
    }    
    for (int i=0; i<matFrame.mcTable.getRowCount(); i++){
       matFrame.mcTable.setValueAt(mtclear[0],i,0);
    }
    String logMaterial, msMaterial, mcMaterial;
    for (int i=0;i<logPanel.logTable.getRowCount();i++){
     logMaterial = ((String)logPanel.logDataModel.getValueAt(i,2));
      if(matFrame.msDataModel.getRowCount()!=-1){
       for (int j=0;j<matFrame.msDataModel.getRowCount();j++){
           msMaterial = ((String)matFrame.msDataModel.getValueAt(j,1));
           if(!logMaterial.equals(msMaterial)){
           }else if(logMaterial.equals(msMaterial)){
            matFrame.msDataModel.setValueAt("Used",j,0);
           }
       }
      }
      if(matFrame.mcDataModel.getRowCount()!=-1){
        for (int h=0;h<matFrame.mcDataModel.getRowCount();h++){
          mcMaterial = ((String)matFrame.mcDataModel.getValueAt(h,1));
          if(!logMaterial.equals(mcMaterial)){
          }else if(logMaterial.equals(mcMaterial)){
           matFrame.mcDataModel.setValueAt("Used",h,0);
          }
        }
      }
    }  
  }
  
}











