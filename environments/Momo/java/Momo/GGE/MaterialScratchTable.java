/*
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class MaterialScratchTable extends MaterialAbstTable {
  MaterialScratchItem item;
  JDialog msDialog;
  CreateMaterialScratch cms;
  DellScratchDialog dsd;
  MaterialScratchTable(DefaultTableModel model, MaterialFrame parent){
    super(model, parent);
    tableTitle = "   Scratch ";
    cms = new CreateMaterialScratch(parent);
    dsd = new DellScratchDialog(parent);
    TableColumn column = getColumn("A");
    column.setCellEditor(new IntCellEditor());
    column = getColumn("Z");
    column.setCellEditor(new DoubleCellEditor());
  }

  void create(){ 
     cms.setVisible(true);
     parent.et.setVisible(true);
  }

  void delete(){ dsd.setVisible(true);}
  protected void selected(){
    parent.scratchSelected();
  }
}









