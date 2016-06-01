/*
 *
 */
//7.3


import com.sun.java.swing.*;
//import com.sun.java.swing.event.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;
//import java.awt.event.*;

public class MaterialCombiTable extends MaterialAbstTable {
  MaterialCombiItem item;
  CreateMaterialCombi cmc;
  DellCombiDialog dcd;

  MaterialCombiTable(DefaultTableModel model, MaterialFrame parent) {
    super(model, parent);
    tableTitle = "   Combination ";
    cmc = new CreateMaterialCombi(parent);
    dcd = new DellCombiDialog(parent);
    TableColumn column = getColumn("Elements");
     column.setCellEditor(new RatioCellEditor(parent));
     column.sizeWidthToFit();
     column.setWidth(200);
  }
  void create(){
      cmc.setVisible(true);
      parent.et.setVisible(true);
  }
  void delete(){
      dcd.setVisible(true);
  }
  protected void selected(){
    parent.combiSelected();
  }
}












