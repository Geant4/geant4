/* GGE
 * Axially symmetric arrangement of repeated volumes by rotation
 * Active rotation of Bodies or rotation of frames are available
 *   See source/geometry/management/include/G4PVPlacement.hh
 */

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

class RepeatedRotateVolume extends JFrame implements ActionListener {
  PhysicalItem phyItem;
  JTable pTable;
  private JFrame frame1, frame2;
  private JMenuItem close;
  private JButton create, ins, dell;
  DefaultTableModel pDataModel;
    //  private CreateRepeatRotateDialog createRepeatDialog;
  private DellRepeatRotateDialog dellRepeat;
  private InsertRepeatRotateDialog insertRepeatDialog;
  private String move, name, log, type, many, x0, y0, z0;
  private String radius, rotAxis, startAng, incAng, elem;
  private String pValue;
  private String tableIdent[] = { "Move", "pName", "pLogic", "MomType",
              "pMother", "pMany", "X0", "Y0", "Z0", 
              "Radius", "lengthUnit","RotAxis", "Phi_0", "dPhi",
               "angleUnit", "nElem" };
  RepeatedRotateVolume(){
    super("Repeated Placement Axial Symmetric");
    //    createRepeatDialog =new CreateRepeatRotateDialog(this);
    dellRepeat =new DellRepeatRotateDialog(this);
    insertRepeatDialog =new InsertRepeatRotateDialog(this);
    
    JPanel panel = new JPanel();
//    phyPanel = new PhysicalPanel(parent);
    panel.setLayout(new BorderLayout(2, 1));
    pTable = new JTable(pDataModel = new DefaultTableModel(tableIdent,0));
    pTable.setAutoCreateColumnsFromModel(false);
    pTable.setColumnSelectionAllowed(false);
    pTable.setRowSelectionAllowed(true);
    pTable.getTableHeader().setReorderingAllowed(false);
    pTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    JComboBox moveCombo = new JComboBox();  // frame or body
       moveCombo.addItem("Frame");
       moveCombo.addItem("Body");
    JComboBox manyCombo = new JComboBox();
       manyCombo.addItem("true");
       manyCombo.addItem("false");
    JComboBox vectorCombo = new JComboBox();
       vectorCombo.addItem("X");
       vectorCombo.addItem("Y");
       vectorCombo.addItem("Z");
    JComboBox motherCombo = new JComboBox();
       motherCombo.addItem("logical");
       motherCombo.addItem("physical");
//       motherCombo.addItem("NULL");
    LenUnitCombo widthCombo = new LenUnitCombo();
    AngUnitCombo angCombo = new AngUnitCombo();
    TableColumn columns[] = new TableColumn[tableIdent.length];
    for (int i=0; i<tableIdent.length; i++){
       columns[i] = pTable.getColumn(tableIdent[i]);
    }
    DefaultTableCellRenderer nameCellRenderer = new DefaultTableCellRenderer();
      nameCellRenderer.setBackground(Color.pink);

    columns[0].setCellEditor(new DefaultCellEditor(moveCombo));
    columns[0].setMinWidth(40);

    columns[1].setCellEditor(new NameCellEditor());
    columns[1].setCellRenderer(nameCellRenderer);
    columns[1].setMinWidth(60);
    columns[1].setResizable(true);

    columns[2].setCellEditor(new NameCellEditor());
    columns[2].setMinWidth(60);

    columns[3].setCellEditor(new DefaultCellEditor(motherCombo));
    columns[3].setMinWidth(60);

    columns[4].setCellEditor(new NameCellEditor());
    columns[4].setMinWidth(50);

    columns[5].setCellEditor(new DefaultCellEditor(manyCombo));
    columns[5].setMinWidth(60);

    columns[6].setCellEditor(new DoubleCellEditor());
    columns[6].setMinWidth(40);


    columns[7].setCellEditor(new DoubleCellEditor());
    columns[7].sizeWidthToFit();
    columns[7].setMinWidth(40);


    columns[8].setCellEditor(new DoubleCellEditor());
    columns[8].sizeWidthToFit();
    columns[8].setMinWidth(40);


    columns[9].setCellEditor(new DoubleCellEditor());
    columns[9].setMinWidth(45);


    columns[10].setCellEditor(new DefaultCellEditor(widthCombo));
    columns[10].setMinWidth(70);

    columns[11].setCellEditor(new DefaultCellEditor(vectorCombo));
    columns[11].setMinWidth(55);

    columns[12].setCellEditor(new DoubleCellEditor());
    columns[12].setMinWidth(50);


    columns[13].setCellEditor(new DoubleCellEditor());
    columns[13].setMinWidth(50);

    columns[14].setCellEditor(new DefaultCellEditor(angCombo));
    columns[14].setMinWidth(60);

    columns[15].setCellEditor(new NameCellEditor());
    columns[15].setMinWidth(50);

    JScrollPane scrollpane = new JScrollPane(pTable);
    scrollpane.setPreferredSize(new Dimension(850, 180));
    scrollpane.setMinimumSize(new Dimension(50, 25));
    JPanel panel2 = new JPanel();
      panel2.setLayout(new GridLayout(1,3));
      panel2.add(create = new JButton(" Append "));
      panel2.add(ins = new JButton(" Insert "));
      panel2.add(dell = new JButton(" Delete "));
    panel.add("Center",scrollpane);
    panel.add("North",panel2);
    getContentPane().add(panel, BorderLayout.CENTER);
    create.addActionListener(this);
    dell.addActionListener(this);
    ins.addActionListener(this);
    PhysicalItem phyItem = new PhysicalItem();
    setJMenuBar(createMenubar());
    pack();
    setResizable(false);
  }
  private JMenuBar createMenubar(){
    JMenu file;
    JMenuBar mb = new JMenuBar();
    file = new JMenu("File");
//    file.add(save = new JMenuItem("Save"));
    file.add(close = new JMenuItem("Close"));
    mb.add(file);
    close.addActionListener(this);
    return mb;
  }
  
  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if(o==dell){
       dellRepeat.setVisible(true);
       dell.setForeground(Color.red);
    }else if(o==create){
      appendRepeat();
//      create.setForeground(Color.red);
    }else if(o==ins){
      insertRepeatDialog.setVisible(true);
      ins.setForeground(Color.red);
    }else if(o==close){
      setVisible(false);
    }
  }
  public void dellRepeatCloseAct(){dell.setForeground(Color.black);repaint();}
  public void createRepeatCloseAct(){create.setForeground(Color.black);repaint();}
  public void insertRepeatCloseAct(){ins.setForeground(Color.black);repaint();}

  public void repeatedRotClear(){
     pTable.setModel(pDataModel = new DefaultTableModel(tableIdent,0));
     repaint();
  }

  public void appendRepeat(){
     Object tmp[] = new Object[20];
      tmp[0] = "Body";
      tmp[1] = "";
      tmp[2] = "";
      tmp[3] ="physical";
      tmp[4] = "";
      tmp[5] = "false";
      tmp[6] = "0.0";    // X0
      tmp[7] = "0.0";
      tmp[8] = "0.0";
      tmp[9] = "0.0";    
      tmp[10] = "mm";    // radius
      tmp[11] = "Z";      // rot axis
      tmp[12] = "0.0";    // start angle
      tmp[13] = "0.0";    // inc angle
      tmp[14] = "deg";
      tmp[15] = "0";
     pDataModel.addRow(tmp);
  }

  public void insertRepeat(){
     Object tmp[] = new Object[20];
      tmp[0] = "Body";
      tmp[1] = "";
      tmp[2] = "";
      tmp[3] ="physical";
      tmp[4] = "";
      tmp[5] = "false";
      tmp[6] = "0.0";
      tmp[7] = "0.0";
      tmp[8] = "0.0";
      tmp[9] = "0.0";
      tmp[10] = "mm";
      tmp[11] = "Z";
      tmp[12] = "0.0";
      tmp[13] = "0.0";
      tmp[14] = "deg";
      tmp[15] = "0";


     if(pTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the pName","Warning Dialog"
, JOptionPane.WARNING_MESSAGE);
     }else{
      int pRowCount = pTable.getSelectedRow()+1;
      pDataModel.insertRow(pRowCount,tmp);
     }
  }

  public void dellRepeat(){
     if(pTable.getSelectedRow()==-1){
      frame2 = new JFrame();
      JOptionPane.showMessageDialog(frame2, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      pDataModel.removeRow(pTable.getSelectedRow());
      dellRepeat.setVisible(false);
      dell.setForeground(Color.black);
      repaint();
     }
  }

  //C++ source code
  String getCPP(){

   StringBuffer cpp = new StringBuffer("\n// Repeated Placement AxialSymmetoric\n\n");

   if(pDataModel.getRowCount()==0){
      return cpp.toString();
   }else{

    for(int i=0;i<pDataModel.getRowCount();i++){
       String pMother = (String)pDataModel.getValueAt(i,3)+
                             (String)pDataModel.getValueAt(i,4);      
// move 
      move = pDataModel.getValueAt(i,0).toString();
// pName
      name = pDataModel.getValueAt(i,1).toString();
// pLogical
      log = pDataModel.getValueAt(i,2).toString();
// Mother Type
      type = pDataModel.getValueAt(i,3).toString();
// pMany
      many = pDataModel.getValueAt(i,5).toString();
// X0
      x0 = pDataModel.getValueAt(i,6).toString()+"*"+
                             pDataModel.getValueAt(i,10).toString();
// Y0
      y0 = pDataModel.getValueAt(i,7).toString()+"*"+
                             pDataModel.getValueAt(i,10).toString();
// Z0
      z0 = pDataModel.getValueAt(i,8).toString()+"*"+
                             pDataModel.getValueAt(i,10).toString();
// Radius
      radius = pDataModel.getValueAt(i,9).toString()+"*"+
                             pDataModel.getValueAt(i,10).toString();
// Rotation Axis
      rotAxis = pDataModel.getValueAt(i,11).toString();
// StartAngle
      startAng = pDataModel.getValueAt(i,12).toString()+"*"+
                             pDataModel.getValueAt(i,14).toString();
// IncAngle
      incAng = pDataModel.getValueAt(i,13).toString()+"*"+
                             pDataModel.getValueAt(i,14).toString();
// nElement
      elem = pDataModel.getValueAt(i,15).toString();
// 
 if(move.equals("Body")){   // active rotation
     cpp.append("// Active Rotation of Bodies\n");
     cpp.append("G4int copy"+name+";\ncopy"+name+"=0;\n");
     cpp.append("G4RotationMatrix rotationMatrix"+name+";");
     cpp.append("G4double start"+name+" = "+startAng+";\n");
     cpp.append("G4double inc"+name+" = "+incAng+";\n");


       cpp.append("rotationMatrix"+name+".rotate"+rotAxis+"(start"+name+");\n");
       cpp.append("G4double x0"+name+" = "+x0+";\n");
       cpp.append("G4double y0"+name+" = "+y0+";\n");
       cpp.append("G4double z0"+name+" = "+z0+";\n");
       cpp.append("G4double radius"+name+" = "+radius+";\n\n");
       cpp.append("for (G4int "+name+"=1; "+name+"<="+elem+"; "+name+"++){\n");
       cpp.append("  G4double trans"+name+" = start"+name+"+"+
                                                "inc"+name+"*("+name+"-1);\n\n");
       cpp.append("  G4double x"+name+", y"+name+", z"+name+";\n");

     if(rotAxis.equals("X")){
       cpp.append("  x"+name+" = x0"+name+";\n");
       cpp.append("  y"+name+" = y0"+name+" + radius"+name+
                                          " * cos(trans"+name+");\n");
       cpp.append("  z"+name+" = z0"+name+" + radius"+name+            
                                          " * sin(trans"+name+");\n\n");
     }else if(rotAxis.equals("Y")){
       cpp.append("  x"+name+" = x0"+name+" + radius"+name+
                                         " * sin(trans"+name+");\n");
       cpp.append("  y"+name+" = y0"+name+";\n");
       cpp.append("  z"+name+" = z0"+name+" + radius"+name+
                                        " * cos(trans"+name+");\n");
     }else if(rotAxis.equals("Z")){
       cpp.append("  x"+name+" = x0"+name+" + radius"+name+
                                          " * cos(trans"+name+");\n");
       cpp.append("  y"+name+" = y0"+name+" + radius"+name+
                                          " * sin(trans"+name+");\n");
       cpp.append("  z"+name+" = z0"+name+";\n");
     }

       cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                " = new G4PVPlacement(G4Transform3D(rotationMatrix"+name+",\t//rotate\n"+
                "\t\t G4ThreeVector("+"x"+name +", y"+name+", z"+name+")),\n");

		   if(type.equals("physical")){
                     cpp.append("\t\t \"physical"+name+"\",   //its name \n"+
                     "\t\t logical"+log+",         //its logical volume\n");
	           }else if(type.equals("logical")){//mother logical
                     cpp.append("\t\t logical"+log+",\t//its logical volume\n"+
                     "\t\t \"physical"+name+"\",   //its name \n");
		   }
                   cpp.append(
                   "\t\t "+pMother+",            //its mother volume\n"+
                   "\t\t "+many+",             //no boolean operation\n"+
                   "\t\t copy"+name+"++);      //copy number \n"+
                   "  rotationMatrix"+name+".rotate"+rotAxis+"(inc"+name+");\n"+
                   "}\n\n");
                     
     
 }else if(move.equals("Frame")){
   cpp.append("// Rotation of Frames\n\n");

     cpp.append("G4int copy"+name+";\ncopy"+name+"=0;\n");
     cpp.append("G4RotationMatrix rotationMatrix"+name+"; //rotation matrix \n");
     cpp.append("G4double start"+name+" = "+startAng+";  // starting angle\n");
     cpp.append("G4double inc"+name+" = "+incAng+";  // incremental angle\n");


       cpp.append("rotationMatrix"+name+".rotate"+rotAxis+"(start"+name+");  //rotate to the starting angle\n");
       cpp.append("G4double x0"+name+" = "+x0+";  //X0 position\n");
       cpp.append("G4double y0"+name+" = "+y0+";  //Y0 position\n");
       cpp.append("G4double z0"+name+" = "+z0+";  //Z0 position\n");
       cpp.append("G4double radius"+name+" = "+radius+";  //radius\n\n");
       cpp.append("for (G4int "+name+"=1; "+name+"<="+elem+"; "+name+"++){  // loop for incremental rotation\n");
       cpp.append("  G4double trans"+name+" = -(start"+name+") + "+
                                                "-(inc"+name+")*("+name+"-1);  // angle for every increment\n\n");
       cpp.append("  G4double x"+name+", y"+name+", z"+name+";  // coordinates for every incremental rotation\n");

     if(rotAxis.equals("X")){
       cpp.append("  x"+name+" = x0"+name+";  // rotation around the X axisi\n");
       cpp.append("  y"+name+" = y0"+name+" + radius"+name+
                                          " * cos(trans"+name+");\n");
       cpp.append("  z"+name+" = z0"+name+" + radius"+name+            
                                          " * sin(trans"+name+");\n\n");
     }else if(rotAxis.equals("Y")){
       cpp.append("  x"+name+" = x0"+name+" + radius"+name+
                                         " * sin(trans"+name+");  // rotation around the Y axis\n");
       cpp.append("  y"+name+" = y0"+name+";\n");
       cpp.append("  z"+name+" = z0"+name+" + radius"+name+
                                        " * cos(trans"+name+");\n");
     }else if(rotAxis.equals("Z")){
       cpp.append("  x"+name+" = x0"+name+" + radius"+name+
                                          " * cos(trans"+name+");  // rotation around the Z axis\n");
       cpp.append("  y"+name+" = y0"+name+" + radius"+name+
                                          " * sin(trans"+name+");\n");
       cpp.append("  z"+name+" = z0"+name+";\n");
     }

       cpp.append("  G4VPhysicalVolume * "+"physical"+name+
                " = new G4PVPlacement( new G4RotationMatrix(rotationMatrix"+name+"),\t//rotate\n"+
                "\t\t G4ThreeVector("+"x"+name +", y"+name+", z"+name+"),  // translate\n");

		   if(type.equals("physical")){
                     cpp.append("\t\t \"physical"+name+"\",   //mother is physical. its name \n"+
                     "\t\t logical"+log+",         //its logical volume\n");
	           }else if(type.equals("logical")){//mother logical
                     cpp.append("\t\t logical"+log+",\t// mother is logical. Current logical volume\n"+
                     "\t\t \"physical"+name+"\",   //its name \n");
		   }
                   cpp.append(
                   "\t\t "+pMother+",            //its mother volume\n"+
                   "\t\t "+many+",             //no boolean operation\n"+
                   "\t\t copy"+name+"++);      //copy number \n"+
                   "  rotationMatrix"+name+".rotate"+rotAxis+"(inc"+name+");\n"+
                   "}\n\n");
                     

 } else { // do nothing. never happens.
   cpp.append("// No Rotation \n\n");
 }
    }
    return cpp.toString();
   }
  }
}














