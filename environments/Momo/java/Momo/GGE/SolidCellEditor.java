import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import com.sun.java.swing.table.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class SolidCellEditor extends Object implements TableCellEditor, ActionListener {
  private Component compo;
  private Vector editorListener;
  private Object value;
  private SolidDialog solidDialog;
  VolumesFrame parent;
  SolidCellEditor(VolumesFrame parent){
    this.parent = parent;
    editorListener = new Vector();
    compo = new JLabel();
  }
  public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int rowIndex, int columnIndex){
    this.value = value;
    System.out.println(isSelected);
    return compo;
  }
  public Object getCellEditorValue(){
    solidDialog.getValues();
    return value;
  } 
  public boolean isCellEditable(EventObject ev){
    return true;
  } 
  public boolean shouldSelectCell(EventObject ev){
    if (value instanceof CSGItem){
      solidDialog = new CSGDialog(this, (CSGItem)value );
    }
    if (value instanceof BREPItem){
      solidDialog = new BREPDialog(this, (BREPItem)value );
    }
    return true;
  } 
  public boolean stopCellEditing(){
    solidDialog.editStop();
    castStop();
    return true;
  } 
  public void cancelCellEditing(){
    solidDialog.editStop();
    castCancel();
  }
  public void csgCastStop(){
    for (int i=0; i<editorListener.size(); i++){
      CellEditorListener cel = (CellEditorListener)editorListener.elementAt(i);
      cel.editingStopped(new ChangeEvent(this));
    }  
  }
  private void castStop(){
    for (int i=0; i<editorListener.size(); i++){
      CellEditorListener cel = (CellEditorListener)editorListener.elementAt(i);
      cel.editingStopped(new ChangeEvent(this));
    }
  }
  private void castCancel(){
    for (int i=0; i<editorListener.size(); i++){
      CellEditorListener cel = (CellEditorListener)editorListener.elementAt(i);
      cel.editingCanceled(new ChangeEvent(this));
    }
  }
  public void addCellEditorListener(CellEditorListener ce){
    editorListener.addElement(ce);
  }
  public void removeCellEditorListener(CellEditorListener ce){
    editorListener.removeElement(ce);
  } 
  public void actionPerformed(ActionEvent ae){
    stopCellEditing();
  }
}






