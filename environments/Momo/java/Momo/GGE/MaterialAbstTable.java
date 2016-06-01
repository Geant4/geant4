/*
 *
 */
//7.3

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;
import com.sun.java.swing.event.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;

abstract public class MaterialAbstTable extends JTable implements ActionListener {
  protected JButton ins, del;
  protected String tableTitle;
  MaterialFrame parent;
  protected DefaultTableModel model;

  MaterialAbstTable(DefaultTableModel model, MaterialFrame parent){

    super(model);
    this.model = model;
    this.parent = parent;

    setAutoCreateColumnsFromModel(false);
    setColumnSelectionAllowed(false);
    setRowSelectionAllowed(true);
    getTableHeader().setReorderingAllowed(false);
    setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

//Jccombobox items don't corresponds to the GEANT4 units.
//They shall be translated into correct ones in making C++ source codes

    JComboBox stateCombo = new JComboBox();
       stateCombo.addItem("kStateUndefind");
       stateCombo.addItem("kStateSolid");
       stateCombo.addItem("kStateLiquid");
       stateCombo.addItem("kStateGas");
//     stateCombo.addItem("Vaccum");  //deleted 8.13
    JComboBox pressUnitCombo = new JComboBox();
       pressUnitCombo.addItem("pascal");
       pressUnitCombo.addItem("bar");
       pressUnitCombo.addItem("atmosphere");
    JComboBox denl = new JComboBox();
       denl.addItem("mg/cm3");
       denl.addItem("g/cm3");
       denl.addItem("kg/m3");
    JComboBox kl = new JComboBox();
       kl.addItem("kelvin");

//Only used materials have C++ codes
    JComboBox use = new JComboBox();
      use.addItem("    ");
      use.addItem("Used");

    DefaultTableCellRenderer useCellRenderer = new DefaultTableCellRenderer();
      useCellRenderer.setForeground(Color.red);

    DefaultTableCellRenderer nameCellRenderer = new DefaultTableCellRenderer();
      nameCellRenderer.setBackground(Color.pink);
        
    TableColumn column = new TableColumn();
        column = getColumn("Use");
              column.setCellEditor(new DefaultCellEditor(use));
              column.setCellRenderer(useCellRenderer);
              column.setWidth(4);    

        column = getColumn("Name");
              column.setCellEditor(new NameCellEditor());
              column.setCellRenderer(nameCellRenderer);
              column.sizeWidthToFit();
              column.setResizable(true);

        column = getColumn("Density");
              column.setCellEditor( new DoubleCellEditor());
              column.sizeWidthToFit();
              column.setResizable(true);

        column = getColumn("Unit");
              column.setWidth(5);
              column.setCellEditor( new DefaultCellEditor(denl));

        column = getColumn("State");
              column.setCellEditor(new DefaultCellEditor(stateCombo));
              column.setWidth(8);
              column.setResizable(true);

        column = getColumn("Temp");
     	      column.setCellEditor(new DoubleCellEditor());
       	      column.sizeWidthToFit();
              column.setResizable(true);

	column = getColumn(" Unit");
	      column.setCellEditor(new DefaultCellEditor(kl));
              column.setWidth(6);

	column = getColumn("Press");
              column.setCellEditor(new DoubleCellEditor());
              column.sizeWidthToFit();
              column.setResizable(true);

        column = getColumn("Unit ");
              column.setCellEditor(new DefaultCellEditor(pressUnitCombo));
              column.setWidth(6);
  }
  public JPanel makeMaterialPanel(){
    JScrollPane scrollpane = new JScrollPane(this);
    	scrollpane.setPreferredSize(new Dimension(780, 250));
    		scrollpane.setMinimumSize(new Dimension(150, 50));
    JPanel tab = new JPanel();
    	tab.setLayout(new BorderLayout(5,3));

    JLabel mtitle = new JLabel(tableTitle);
                mtitle.setFont(new Font("Serif", Font.BOLD, 12));
                Dimension tsize = mtitle.getPreferredSize();
                mtitle.setSize(tsize);

//this method seems not working as intentioned  8.13
    JPanel buttonPanel = new JPanel();
           buttonPanel.setLayout(new BorderLayout(3,0));
           buttonPanel.add(mtitle, "West");
           ins = new JButton("Create");
            buttonPanel.add(ins, "Center");
            del = new JButton("Delete");
            buttonPanel.add(del, "East");
            tab.add(buttonPanel,"North");
            tab.add(scrollpane, "Center");

            ins.addActionListener(this);
            del.addActionListener(this);
            return tab;

  }
  public void actionPerformed(ActionEvent ae){
    Object o = ae.getSource();
    if (o == ins) {
      create();
      ins.setForeground(Color.red);
      return;
    }
    if (o == del){
      delete();
      del.setForeground(Color.red);
      return;
    }
  }
  public void dellMatCloseAct(){
     del.setForeground(Color.black);
     del.repaint();
  }
  public void createMatCloseAct(){
     ins.setForeground(Color.black);
     ins.repaint();
  }

  abstract void create();
  abstract void delete();

  void stopEditing(){
    if (isEditing()) cellEditor.stopCellEditing();
    clearSelection();
  }
  Object getSelectedName(){
    if (isEditing()) cellEditor.stopCellEditing();
    Vector data = model.getDataVector();
    int row = getSelectedRow();
    if (row < 0) return null;
    return ((Vector)data.elementAt(row)).elementAt(0);
  }
  public void valueChanged(ListSelectionEvent lse){
    super.valueChanged(lse);
    if (lse.getValueIsAdjusting()) return;
    if (lse.getFirstIndex() < 0) return;
    if (lse.getLastIndex() < 0) return;

    System.out.println(lse);
    selected();
  }
  protected abstract void selected();
}
