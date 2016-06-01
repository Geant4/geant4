
//example
//8.5

import java.awt.*;
import java.awt.event.*;
import com.sun.java.swing.*;
import com.sun.java.swing.event.*;
import java.util.*;

public class InsertCombiDialog extends JDialog implements ActionListener{
   JButton b1, b2;
   JLabel label1, label2;
   private MaterialFrame matFrame;
   public InsertCombiDialog(MaterialFrame matFrame){
      super(matFrame,"Insert Material (Combination)",false);
      this.matFrame = matFrame;
      JPanel createMaterialPanel = new JPanel();
       createMaterialPanel.setLayout(new BorderLayout());
       b1 = new JButton("OK");
       b2 = new JButton("CANCEL");
       b1.addActionListener(this);
       b2.addActionListener(this);
       JPanel labelPanel = new JPanel();
        labelPanel.setLayout(new BorderLayout());
        label1 = new JLabel("   Click the MaterialName after which ");
        label1.setFont(new Font("Serif",Font.BOLD, 12));
        label2 = new JLabel("   you are going to insert");
        label2.setFont(new Font("Serif",Font.BOLD, 13));
        labelPanel.add(label1, BorderLayout.NORTH);
        labelPanel.add(label2, BorderLayout.CENTER);
       createMaterialPanel.add(labelPanel, BorderLayout.NORTH);
       createMaterialPanel.add(b1, BorderLayout.CENTER);
       createMaterialPanel.add(b2, BorderLayout.EAST);
      getContentPane().add(createMaterialPanel, BorderLayout.CENTER);
//      pack();
      setSize(250, 100);

   }
   public void actionPerformed(ActionEvent e){
      if(e.getSource() == b1){
        matFrame.insertMC();
        this.setVisible(false); 
      }else if(e.getSource() == b2){
        this.setVisible(false);       
      }
   }
 
}





