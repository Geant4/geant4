/*
 *
 */

import com.sun.java.swing.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;

class IntCellEditor extends DefaultCellEditor {
  JFrame frame;
  IntCellEditor(){
    super(new JTextField());
    setClickCountToStart(1);
  }
  public Object getCellEditorValue(){
    Integer value;
    String valueStr = (String)super.getCellEditorValue();
    if (valueStr.equals("")) return "";
    try{
      value = new Integer(valueStr);
    }catch(NumberFormatException e){
    frame = new JFrame();
    JOptionPane.showMessageDialog(frame,"Not a Number", "Warning Dialog",
                                           JOptionPane.WARNING_MESSAGE,null);
     return "";
//      value = new Integer(0);
    }
    return value;
  } 
}



