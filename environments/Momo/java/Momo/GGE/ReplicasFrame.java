//1998 8 20
//GGE 
//Tetsuya Yamada

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

class ReplicasFrame extends JFrame implements ActionListener{
  private JButton create, ins, dell;
  private JFrame frame1, frame2;
  private JMenuItem close;
  private JComboBox motherCombo;
  JTable repTable;
  DefaultTableModel repDataModel;
  private CreateReplicaDialog createRepDialog;
  private DellReplicaDialog dellRep;
  private InsertReplicaDialog insertRepDialog;
  private String pValue;
  private String mother;
  private String tableIdent[] = { "pName", "pLogic", "MomType", "pMother", "EAxis", "nReplicas", "Width", "WidthUnit", "Offset", "OffsetUnit" };  
   
  ReplicasFrame(){
    super("Slicing Translation ");
    createRepDialog = new CreateReplicaDialog(this);
    dellRep = new DellReplicaDialog(this);
    insertRepDialog = new InsertReplicaDialog(this);

    JPanel panel = new JPanel();
//    phyPanel = new PhysicalPanel(parent);
    panel.setLayout(new BorderLayout(2, 1));
    repTable = new JTable(repDataModel = new DefaultTableModel(tableIdent,0));
    repTable.setAutoCreateColumnsFromModel(false);
    repTable.setColumnSelectionAllowed(false);
    repTable.setRowSelectionAllowed(true);
    repTable.getTableHeader().setReorderingAllowed(false);
    repTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

    JComboBox motherCombo = new JComboBox();
       motherCombo.addItem("logical");
       motherCombo.addItem("physical");
//       motherCombo.addItem("NULL");
//       motherCombo.addActionListener(this);
    JComboBox eAxisCombo = new JComboBox();
       eAxisCombo.addItem("kXAxis");
       eAxisCombo.addItem("kYAxis");
       eAxisCombo.addItem("kZAxis");
//       eAxisCombo.addItem("kRho");
//       eAxisCombo.addItem("kPhi");

    JComboBox angleCombo = new JComboBox();
       angleCombo.addItem("rad");
       angleCombo.addItem("mrad");
       angleCombo.addItem("deg");

    JComboBox widCombo = new JComboBox();
       widCombo.addItem("mm");
       widCombo.addItem("cm");
       widCombo.addItem("m");
       widCombo.addItem("km");
       widCombo.addItem("micrometer");
       widCombo.addItem("nanoometer");

    JComboBox offsetCombo = new JComboBox();
       offsetCombo.addItem("mm");
       offsetCombo.addItem("cm");
       offsetCombo.addItem("m");
       offsetCombo.addItem("km");
       offsetCombo.addItem("micrometer");
       offsetCombo.addItem("nanoometer");


    TableColumn column[] = new TableColumn[tableIdent.length];
    for (int i=0; i<tableIdent.length; i++){
       column[i] = repTable.getColumn(tableIdent[i]);
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
    JScrollPane scrollpane = new JScrollPane(repTable);
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
       dellRep.setVisible(true);
       dell.setForeground(Color.red);
    }else if(o==create){
       appendRep();
//      createRepDialog.setVisible(true);
//      create.setForeground(Color.red);
    }else if(o==ins){
      insertRepDialog.setVisible(true);
      ins.setForeground(Color.red);
    }else if(o==close){
      setVisible(false);
    }
  }
  public void dellRepCloseAct(){
     dell.setForeground(Color.black);
     dell.repaint();
  }
  public void createRepCloseAct(){
     create.setForeground(Color.black);
     create.repaint();
  }
  public void insertRepCloseAct(){
     ins.setForeground(Color.black);
     ins.repaint();
  }
  public void replicasClear(){
     repTable.setModel(repDataModel =
                             new DefaultTableModel(tableIdent,0));
     repaint();
  } 
  public void appendRep(){
     Object repTmp[] = new Object[10];
      repTmp[0] = "";
      repTmp[1] = "";
      repTmp[2] = "logical";
      repTmp[3] = "";
      repTmp[4] = "kXAxis";
      repTmp[5] = "";
      repTmp[6] = "";
      repTmp[7] = "mm";
      repTmp[8] = "0.0";
      repTmp[9] = "mm";
     repDataModel.addRow(repTmp);
  }
  public void insertRep(){
     Object tmp[] = new Object[10];
      tmp[0] = "";
      tmp[1] = "";
      tmp[2] = "logical";
      tmp[3] = "";
      tmp[4] = "kXAxis";
      tmp[5] = "";
      tmp[6] = "";
      tmp[7] = "mm";
      tmp[8] = "0.0";
      tmp[9] = "mm";
     if(repTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      int repRowCount = repTable.getSelectedRow()+1;
      repDataModel.insertRow(repRowCount,tmp);
     }
  }
  
  public void dellRep(){
     if(repTable.getSelectedRow()==-1){
      frame2 = new JFrame();
      JOptionPane.showMessageDialog(frame2, "Choose the pName","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      repDataModel.removeRow(repTable.getSelectedRow());
      dellRep.setVisible(false);
      dell.setForeground(Color.black);
      repaint();
     }
  }

//C++ source code
  String getCPP(){
    StringBuffer cpp = new StringBuffer("\n// Slicing Translation \n\n");
    if(repDataModel.getRowCount()==0){
      return cpp.toString();
    }else{
     for(int i=0;i<repDataModel.getRowCount();i++){
       String pMother = (String)repDataModel.getValueAt(i,2)+
                                   (String)repDataModel.getValueAt(i,3);

       cpp.append( "G4PVReplica *  physical"+repDataModel.getValueAt(i,0)+
	  "= new G4PVReplica(\"physical"+repDataModel.getValueAt(i,0)+"\",  //name\n" 
	  +"\t\tlogical"+repDataModel.getValueAt(i,1)+",  // its logical\n"  
	  + "\t\t" + pMother+",  // its mother\n"
	  + "\t\t" + repDataModel.getValueAt(i,4)+", // along Axis\n"
	  + "\t\t" + repDataModel.getValueAt(i,5)+", // No of replicas\n"
	  + "\t\t" + repDataModel.getValueAt(i,6)+"*"+repDataModel.getValueAt(i,7)+", //  width\n"
	  + "\t\t" + repDataModel.getValueAt(i,8)+"*"+repDataModel.getValueAt(i,9)+"); // offset\n\n");
    }
   }
   return cpp.toString();
  }
  


}











