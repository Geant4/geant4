// 9.17  C++ source code add coment 
/*
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

class LogicalPanel extends JPanel implements ActionListener {
  private JButton ok, app, ins, del;
  private JComboBox solidSelect;
  private JFrame messageFrame, frame1;
  JTable logTable;
  MaterialScratchTable msTable;
  MaterialCombiTable mcTable;
    //  CreateLogicalDialog createLogDialog;
  DellLogicalDialog dellDialog;
  InsertLogicalDialog insertLogDialog;

  DefaultTableModel logDataModel;
  String tableIdent[] = {"Name","Solid","Material","VisAtb"};
  VolumesFrame volFrame;
  LogicalPanel(VolumesFrame volFrame){
    this.volFrame = volFrame;
    //    createLogDialog = new CreateLogicalDialog(volFrame);
    insertLogDialog = new InsertLogicalDialog(volFrame);
    setLayout(new BorderLayout(2,2));

    logTable = new JTable(logDataModel = new DefaultTableModel(tableIdent,0));
    logTable.setAutoCreateColumnsFromModel(false);
    logTable.setColumnSelectionAllowed(false);
    logTable.setRowSelectionAllowed(true);
    logTable.getTableHeader().setReorderingAllowed(false);
    logTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    TableColumn columns[] = new TableColumn[tableIdent.length];
    for (int i=0; i<tableIdent.length; i++){
      columns[i] = logTable.getColumn(tableIdent[i]);
    }
    DefaultTableCellRenderer logNameCellRenderer = new DefaultTableCellRenderer();
      logNameCellRenderer.setBackground(Color.pink);

    columns[0].setCellEditor(new NameCellEditor());
    columns[0].setCellRenderer(logNameCellRenderer);
    columns[0].setMinWidth(50);
    columns[1].setCellEditor(new SolidCellEditor(volFrame));
    columns[1].setMinWidth(60);
    columns[2].setMinWidth(60);
    columns[3].setMinWidth(50);
    JScrollPane scrollpane = new JScrollPane(logTable);
    scrollpane.setPreferredSize(new Dimension(350, 270));
    scrollpane.setMinimumSize(new Dimension(200, 30));

    JPanel tablePanel = new JPanel();
    tablePanel.setLayout(new BorderLayout(2,2));
     tablePanel.setBorder(new TitledBorder(LineBorder.createBlackLineBorder(),
                                                         " Logical Volume "));
     tablePanel.add("Center", scrollpane);


    JPanel south = new JPanel(new BorderLayout());
     solidSelect = new JComboBox();
     solidSelect.setMaximumRowCount(10);
     for (int i=0; i<CSGItem.solidName.length; i++){
       solidSelect.addItem(new CSGItem(i));
     }
    for (int i=0; i<BREPItem.solidName.length; i++){
      solidSelect.addItem(new BREPItem(i));
    }
     south.add("West",new JLabel("Select a solid "));
     south.add("Center", solidSelect);

    JPanel buttonPanel = new JPanel();
     buttonPanel.setLayout(new GridLayout(1,3));
     buttonPanel.add(app = new JButton("Append"));
     buttonPanel.add(ins = new JButton("Insert"));
     buttonPanel.add(del = new JButton("Delete"));

     south.add("North", buttonPanel);

    JPanel materialPane =new JPanel();
    materialPane.setBorder(new TitledBorder(LineBorder.createBlackLineBorder(),
                                                               " Materials "));
     materialPane.setLayout(new BorderLayout(2,2));
     materialPane.add("North", ok = new JButton("Used Materials"));

     ok.addActionListener(this);
     ok.setForeground(Color.blue);
     ins.addActionListener(this);
     del.addActionListener(this);
     app.addActionListener(this);
     tablePanel.add("North", south);

    add("Center", tablePanel);
    add("South", materialPane);
  }

  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if (o == ok ){ 
       messageFrame = new JFrame();
       ok.setForeground(Color.red);
       int choose = JOptionPane.showConfirmDialog(messageFrame, "Clear and set flags to materials Used in Logical Volum","Used Materials Dialog", JOptionPane.YES_NO_OPTION);
       if(choose == JOptionPane.YES_OPTION){
         volFrame.okCommand();
         ok.setForeground(Color.black); 
       }else{
         ok.setForeground(Color.black);
         return;
       }
    }
    if (o == app){        
        addLogicalRow();
        app.setForeground(Color.red);
        repaint();
        appLogCloseAct();
        return; 
    }
    if (o == ins){
        insertLogDialog.setVisible(true);
        ins.setForeground(Color.red);
        return;      
    }
    if (o == del){
        dellDialog = new DellLogicalDialog(volFrame);
        dellDialog.setVisible(true);
        del.setForeground(Color.red);
       
    }
  }
  public void dellLogCloseAct(){del.setForeground(Color.black);repaint();}
  public void appLogCloseAct(){app.setForeground(Color.black);repaint();}
  public void insertLogCloseAct(){ins.setForeground(Color.black);repaint();}

  public void logClear(){
    logTable.setModel(logDataModel = new DefaultTableModel(tableIdent,0));
    repaint();
  }

  public void addLogicalRow(){
    SolidItem seleItem = (SolidItem)solidSelect.getSelectedItem();
    Object newRow[] = new Object[4];
    newRow[0] = "";
    newRow[2] = "";
    if (seleItem instanceof CSGItem){
      newRow[1] = new CSGItem(seleItem.getSolidType());
    }else{
      newRow[1] = new BREPItem(seleItem.getSolidType());
    }
    newRow[3] = volFrame.visPanel.getSelectedVis();
    if (newRow[3] == null) newRow[3] = "";
    logDataModel.addRow(newRow);
    appLogCloseAct();
    app.setForeground(Color.black);
    repaint();
  }

  public void insertLogicalRow(){
    SolidItem seleItem = (SolidItem)solidSelect.getSelectedItem();
    Object newRow[] = new Object[4];
    newRow[0] = "";
    newRow[2] = "";
    if (seleItem instanceof CSGItem){
      newRow[1] = new CSGItem(seleItem.getSolidType());
    }else{
      newRow[1] = new BREPItem(seleItem.getSolidType());
    }
    newRow[3] = volFrame.visPanel.getSelectedVis();
    if (newRow[3] == null) newRow[3] = "";

     if(logTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the Name","Warning Dialog" , JOptionPane.WARNING_MESSAGE);
     }else{
      int logRowCount = logTable.getSelectedRow()+1;
      logDataModel.insertRow(logRowCount,newRow);
      insertLogDialog.setVisible(false);
      ins.setForeground(Color.black);
      repaint();
     }
  }

  public void dell(){
     if(logTable.getSelectedRow()==-1){
      frame1 = new JFrame();
      JOptionPane.showMessageDialog(frame1, "Choose the Name","Warning Dialog", JOptionPane.WARNING_MESSAGE);
     }else{
      logDataModel.removeRow(logTable.getSelectedRow());
      dellDialog.setVisible(false);
      del.setForeground(Color.black);
      repaint();
     }
  }


/* Data Save Method */
  Vector getSaveData(){
    Vector logicalData = new Vector();
    Component compo[];
    compo = getComponents();
    for (int i = 0; i<compo.length; i++) {
    } 
    return logicalData;
    }
/* Data Load Method */  
  void setLoadData(Vector logicalData){
    removeAll();
    for (int i=0; i<logicalData.size(); i++){
    }
    validate();
  }
 
  String getCPP(){
    StringBuffer cpp = new StringBuffer("\n// Logical  Volumes\n\n");
    Vector data = logDataModel.getDataVector();
    Vector rowData;
    if(data.size()==0){
      return cpp.toString();
    }else{
      for (int i=0; i<data.size(); i++){
        rowData = (Vector)data.elementAt(i);
        String name = (String)rowData.elementAt(0);
        SolidItem solid = (SolidItem)rowData.elementAt(1);
//8-31
         cpp.append(solid.getCPP(name));
         cpp.append("G4LogicalVolume * "+"logical"+name+
                    " = new G4LogicalVolume(solid"+name+", \t //its solid\n"+
                    "\t\t\t\t "+rowData.elementAt(2)+", \t\t //its material\n"+
                    "\t\t\t \"logical"+name+"\" ,\t\t //its name\n"+
                    "\t\t\t 0,0,0);\n\n");
         if (rowData.elementAt(3).equals("")){
           //next data
         }else{
	   cpp.append("logical"+name+"->SetVisAttributes("+
                                        rowData.elementAt(3)+");\n\n");
         }
      }
    }
    return cpp.toString();
  }
}









