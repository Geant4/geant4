/*
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;

class NameCellEditor extends DefaultCellEditor {
  NameCellEditor(){
    super(new JTextField());
//    setClickCountToStart(1);
  }
  public Object getCellEditorValue(){
    String tmp1 = ((String)super.getCellEditorValue()).trim();
    String tmp2 = tmp1.replace(' ', '_');
    return tmp2;
  } 
}
