
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class CreateMaterialCombi extends JDialog implements ActionListener{
   JButton b1, b2, b3;
   JLabel label1, label2;
   private MaterialFrame matFrame;
   private InsertCombiDialog insertCombiDialog;
   public CreateMaterialCombi(MaterialFrame matFrame){
      super(matFrame,"Create Material(Combination)",false);
      this.matFrame = matFrame;
      JPanel createMaterialPanel = new JPanel();
       createMaterialPanel.setLayout(new BorderLayout());
       b1 = new JButton("APPEND");
       b2 = new JButton("INSERT");
       b3 = new JButton("END");
       b1.addActionListener(this);
       b2.addActionListener(this);
       b3.addActionListener(this);
       JPanel labelPanel = new JPanel();
        labelPanel.setLayout(new BorderLayout());
        label1 = new JLabel("  Choose Append or Insert");
        label1.setFont(new Font("Serif",Font.BOLD, 12));
        label2 = new JLabel("  Choose Elements from the ElementsTable");
        label2.setFont(new Font("Serif",Font.BOLD, 12));
        labelPanel.add(label1, BorderLayout.NORTH);
        labelPanel.add(label2, BorderLayout.CENTER);
       createMaterialPanel.add(labelPanel, BorderLayout.NORTH);
       createMaterialPanel.add(b3, BorderLayout.EAST);
       createMaterialPanel.add(b2, BorderLayout.CENTER);
       createMaterialPanel.add(b1, BorderLayout.WEST);
      getContentPane().add(createMaterialPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 100);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        matFrame.appendMC(); 
      }else if(e.getSource() == b2){
        insertCombiDialog = new InsertCombiDialog(matFrame);
        insertCombiDialog.setVisible(true);
      }else if(e.getSource() == b3){
        matFrame.et.setVisible(false);
        matFrame.mcTable.createMatCloseAct();      
        this.setVisible(false);
      } 
   }    
 
}





