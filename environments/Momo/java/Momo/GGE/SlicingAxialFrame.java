//1998 10 26
//GGE 
//Tetsuya Yamada

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

class SlicingAxialFrame extends JFrame implements ActionListener{
  private JButton create, ins, dell;
  private JFrame frame1, frame2;
  private JMenuItem close;
  private JComboBox motherCombo;
  JTable slicingTable;
  DefaultTableModel slicingDataModel;
  private CreateSlicingDialog createSlicingDialog;
  private DellSlicingDialog dellSlicing;
  private InsertSlicingDialog insertSlicingDialog;
  private String pValue;
  private String mother;
  private String tableIdent[] = { "pName", "pLogic", "MomType", "pMother", "EAxis", "nReplicas", "Width", "WidthUnit", "Offset", "OffsetUnit" };  
   
  SlicingAxialFrame(){
    super("Slicing Axial Symmetric");
    createSlicingDialog = new CreateSlicingDialog(this);
    dellSlicing = new DellSlicingDialog(this);
    insertSlicingDialog = new InsertSlicingDialog(this);

    JPanel panel = new JPanel();
//    phyPanel = new PhysicalPanel(parent);
    panel.setLayout(new BorderLayout(2, 1));
    slicingTable = new JTable(slicingDataModel = 
                                  new DefaultTableModel(tableIdent,0));
    slicingTable.setAutoCreateColumnsFromModel(false);
    slicingTable.setColumnSelectionAllowed(false);
    slicingTable.setRowSelectionAllowed(true);
    slicingTable.getTableHeader().setReorderingAllowed(false);
    slicingTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

    JComboBox motherCombo = new JComboBox();
       motherCombo.addItem("logical");
       motherCombo.addItem("physical");
//       motherCombo.addItem("NULL");
//       motherCombo.addActionListener(this);
    JComboBox eAxisCombo = new JComboBox();
//       eAxisCombo.addItem("kXAxis");
//       eAxisCombo.addItem("kYAxis");
//       eAxisCombo.addItem("kZAxis");
//       eAxisCombo.addItem("kRho");
       eAxisCombo.addItem("kPhi");

    JComboBox angleCombo = new JComboBox();
       angleCombo.addItem("rad");
       angleCombo.addItem("mrad");
       angleCombo.addItem("deg");

    JComboBox widCombo = new JComboBox();
       widCombo.addItem("rad");
       widCombo.addItem("mrad");
       widCombo.addItem("deg");

    JComboBox offsetCombo = new JComboBox();
       offsetCombo.addItem("rad");
       offsetCombo.addItem("mrad");
       offsetCombo.addItem("deg");

    TableColumn column[] = new TableColumn[tableIdent.length];
    for (int i=0; i<tableIdent.length; i++){
       column[i] = slicingTable.getColumn(tableIdent[i]);
    }
    DefaultTableCellRenderer nameCellRenderer = new DefaultTableCellRenderer();
      nameCellRenderer.setBackground(Color.pink);
              column[0].setCellEditor(new NameCellEditor());
              column[0].setCellRenderer(nameCellRenderer);
              column[0].sizeWidthToFit();
              column[0].setResizable(true);

              column[1].setCellEditor(new NameCellEditor());

              column[2].setCellEditor( new DefaultCellEditor(motherCombo));

              column[3].setCellEditor(new NameCellEditor());

              column[4].setCellEditor( new DefaultCellEditor(eAxisCombo));

              column[5].setCellEditor(new IntCellEditor());

              column[6].setCellEditor( new DoubleCellEditor());

              column[7].setCellEditor( new DefaultCellEditor(widCombo));

              column[8].setCellEditor( new DoubleCellEditor());

              column[9].setCellEditor(new DefaultCellEditor(offsetCombo));
    JScrollPane scrollpane = new JScrollPane(slicingTable);
    scrollpane.setPreferredSize(new Dimension(750, 180));
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
    ins.addActionListener(this);
    dell.addActionListener(this);

    setJMenuBar(createMenubar());
    pack();
    setResizable(false);
//    setBackground(Color.green);
    setVisible(false);
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
       dellSlicing.setVisible(true);
       dell.setForeground(Color.red);
    }else if(o==create){
       appendSlicing();
//      createSlicingDialog.setVisible(true);
//      create.setForeground(Color.red);
    }else if(o==ins){
      insertSlicingDialog.setVisible(true);
      ins.setForeground(Color.red);
    }else if(o==close){
      setVisible(false);
    }
  }
  public void dellSlicingCloseAct(){
     dell.setForeground(Color.black);
     dell.repaint();
  }
  public void createSlicingCloseAct(){
     create.setForeground(Color.black);
     create.repaint();
  }
  public void insertSlicingCloseAct(){
     ins.setForeground(Color.black);
     ins.repaint();
  }
  public void slicingClear(){
     slicingTable.setModel(slicingDataModel =
                             new DefaultTableModel(tableIdent,0));
     repaint();
  } 
  public void appendSlicing(){
     Object slicingTmp[] = new Object[10];
      slicingTmp[0] = "";
      slicingTmp[1] = "";
      slicingTmp[2] = "logical";
      slicingTmp[3] = "";
      slicingTmp[4] = "kRhi";
      slicingTmp[5] = "";
      slicingTmp[6] = "";
      slicingTmp[7] = "rad";
      slicingTmp[8] = "0.0";
      slicingTmp[9] = "deg";
     slicingDataModel.addRow(slicingTmp);
  }
  public void insertSlicing(){
     Object tmp[] = new Object[10];
      tmp[0] = "";
      tmp[1] = "";
      tmp[2] = "logical";
      tmp[3] = "";
      tmp[4] = "kRhi";
      tmp[5] = "";
      tmp[6] = "";
      tmp[7] = "rad";
      tmp[8] = "0.0";
      tmp[9] = "deg";
     if(slicingTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      int slicingRowCount = slicingTable.getSelectedRow()+1;
      slicingDataModel.insertRow(slicingRowCount,tmp);
     }
  }
  
  public void dellSlicing(){
     if(slicingTable.getSelectedRow()==-1){
      frame2 = new JFrame();
      JOptionPane.showMessageDialog(frame2, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      slicingDataModel.removeRow(slicingTable.getSelectedRow());
      dellSlicing.setVisible(false);
      dell.setForeground(Color.black);
      repaint();
     }
  }

//C++ source code
  String getCPP(){
    StringBuffer cpp = new StringBuffer("\n// Slicing AxialSymmetric \n\n");
    if(slicingDataModel.getRowCount()==0){
      return cpp.toString();
    }else{
     for(int i=0;i<slicingDataModel.getRowCount();i++){
       String pMother = (String)slicingDataModel.getValueAt(i,2)+
                                   (String)slicingDataModel.getValueAt(i,3);

       cpp.append( "G4PVReplica *  physical"+slicingDataModel.getValueAt(i,0)+
	  "= new G4PVReplica(\"physical"+slicingDataModel.getValueAt(i,0)+"\",  //name\n" 
	  +"\t\tlogical"+slicingDataModel.getValueAt(i,1)+",  // its logical\n"  
	  + "\t\t" + pMother+",  // its mother\n"
	  + "\t\t" + slicingDataModel.getValueAt(i,4)+", // along Axis\n"
	  + "\t\t" + slicingDataModel.getValueAt(i,5)+", // No of replicas\n"
	  + "\t\t" + slicingDataModel.getValueAt(i,6)+"*"+
                         slicingDataModel.getValueAt(i,7)+", //  width\n"
	  + "\t\t" + slicingDataModel.getValueAt(i,8)+"*"+
                         slicingDataModel.getValueAt(i,9)+"); // offset\n\n");
    }
   }
   return cpp.toString();
  }
  


}











