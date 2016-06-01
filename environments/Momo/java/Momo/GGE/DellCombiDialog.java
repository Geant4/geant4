
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class DellCombiDialog extends JDialog implements ActionListener{
   JButton b1, b2;
   JLabel label1;
   private MaterialFrame matFrame;
   public DellCombiDialog(MaterialFrame matFrame){
      super(matFrame,"Delete Material (Combination)",false);
      this.matFrame = matFrame;
      JPanel createMaterialPanel = new JPanel();
       createMaterialPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);

        label1 = new JLabel("   Click the MaterialName to delete");
        label1.setFont(new Font("Serif",Font.BOLD, 13));

       createMaterialPanel.add(label1, BorderLayout.NORTH);
       createMaterialPanel.add(b1, BorderLayout.CENTER);
       createMaterialPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createMaterialPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 80);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        matFrame.dellMC();
        this.setVisible(false); 
      }else if(e.getSource() == b2){
        matFrame.mcTable.dellMatCloseAct();
        this.setVisible(false);       
      }
   }
 
}





