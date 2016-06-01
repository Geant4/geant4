/*
 * Copyright (c) 1997 John Jensen. All rights reserved.
 *
 * This software is FREE FOR COMMERCIAL AND NON-COMMERCIAL USE,
 * provided the following condition is met.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted,
 * provided that any copy or derivative of this software or documentation
 * retaining the name "John Jensen" also retains this condition and the
 * following disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * CopyrightVersion 1.0
 */

package TextEdit;

import java.awt.*;
import java.awt.event.*;
import java.awt.datatransfer.*;

class EditMan extends Object implements ActionListener, ClipboardOwner
{
	// remember stuff

	private TextEdit textEdit;
	private TextCanvas textCanvas;
	private String undoString		= new String("Undo");
	private String redoString		= new String("Redo");
	private String copyString		= new String("Copy");
	private String cutString		= new String("Cut");
	private String pasteString		= new String("Paste");
	private String findString		= new String("Find");
	private String replaceString	= new String("Replace");
	private	MenuItem undoItem;
	private	MenuItem redoItem;
	private	MenuItem copyItem;
	private	MenuItem cutItem;
	private Clipboard clipboard;

	private FindDialog findDialog = null;
	private String findPattern = null;
	private ReplaceDialog replaceDialog = null;
	private String ReplacePattern = null;
	private String ReplaceReplace = null;

	// constructor puts together UI

	public EditMan(TextEdit p)
	{
		textEdit = p;
		textCanvas = textEdit.getCanvas();
		clipboard = textEdit.getToolkit().getSystemClipboard();
	}
	
	public void addItems(Menu menu)
	{
		MenuItem item;
		undoItem = new MenuItem(undoString, new MenuShortcut(KeyEvent.VK_Z,false));
		undoItem.setActionCommand(undoString);
		undoItem.setEnabled(false);
		menu.add(undoItem);
		redoItem = new MenuItem(redoString, new MenuShortcut(KeyEvent.VK_Y,false));
		redoItem.setActionCommand(redoString);
		redoItem.setEnabled(false);
		menu.add(redoItem);
		menu.addSeparator();
		copyItem = new MenuItem(copyString, new MenuShortcut(KeyEvent.VK_C,false));
		copyItem.setActionCommand(copyString);
		copyItem.setEnabled(false);
		menu.add(copyItem);
		cutItem = new MenuItem(cutString, new MenuShortcut(KeyEvent.VK_X,false));
		cutItem.setActionCommand(cutString);
		cutItem.setEnabled(false);
		menu.add(cutItem);
		item = new MenuItem(pasteString, new MenuShortcut(KeyEvent.VK_V,false));
		item.setActionCommand(pasteString);
		menu.add(item);
		menu.addSeparator();
		item = new MenuItem(findString, new MenuShortcut(KeyEvent.VK_F,false));
		item.setActionCommand(findString);
		menu.add(item);
		item = new MenuItem(replaceString, new MenuShortcut(KeyEvent.VK_H,false));
		item.setActionCommand(replaceString);
		menu.add(item);
		menu.addActionListener(this);

		textCanvas.setEditMan(this);
	}
	
	public void updateCopyItems(boolean have_selection)
	{
		copyItem.setEnabled(have_selection);
		cutItem.setEnabled(have_selection);
	}
	
	public void updateUndoItems(boolean have_undo, boolean have_redo)
	{
		undoItem.setEnabled(have_undo);
		redoItem.setEnabled(have_redo);
	}
	
    public void actionPerformed(ActionEvent evt)
	{
        String cmd = evt.getActionCommand();

		if (cmd == null)
			return;

		if (cmd.equals(undoString))
		{
			textCanvas.undo(true);
		}
		if (cmd.equals(redoString))
		{
			textCanvas.undo(false);
		}
		else
        if (cmd.equals(copyString))
		{
			// Implement Copy operation
			String srcData = textCanvas.copy(false);
			if (srcData != null)
			{
                StringSelection contents = new StringSelection(srcData);
                clipboard.setContents(contents, this);
			}
        }
		else
        if (cmd.equals(cutString))
		{
			// Implement Copy operation
			String srcData = textCanvas.copy(true);
			if (srcData != null)
			{
                StringSelection contents = new StringSelection(srcData);
                clipboard.setContents(contents, this);
			}
        }
		else
		if (cmd.equals(pasteString))
		{
			String dstData;
            // Implement Paste operation
            Transferable content = clipboard.getContents(this);
            if (content != null)
			{
                try
				{
                    dstData = (String)content.getTransferData(DataFlavor.stringFlavor);
                }
				catch (Exception e)
				{
                    System.out.println("Could not read clipboard"); 
					return;
                }
				textCanvas.paste(dstData);
           }
        }
		else
		if (cmd.equals(findString))
		{
			if (findDialog == null)
				findDialog = new FindDialog(textEdit,this,findPattern);
			findDialog.show();
		}
		else
		if (cmd.equals(replaceString))
		{
			if (replaceDialog == null)
				replaceDialog = new ReplaceDialog(textEdit,this,ReplacePattern,ReplaceReplace);
			replaceDialog.show();
		}
    }
	
	public void CloseFindDialog(String pat)
	{
		findDialog = null;				
		findPattern = new String(pat);
	}

	public void CloseReplaceDialog(String pat, String rep)
	{
		replaceDialog = null;				
		ReplacePattern = new String(pat);
		ReplaceReplace = new String(rep);
	}

    public void lostOwnership(Clipboard clipboard, Transferable contents)
	{
    }
}

class FindDialog extends Frame implements WindowListener, ActionListener
{
	private Button fbutton,cbutton;
	private EditMan editMan;
	private TextEdit textEdit;
	private TextField pattern;

	public FindDialog(TextEdit te, EditMan em, String pat)
	{
		super("Find");

		setBackground(Color.lightGray);

		textEdit = te;
		editMan = em;

		Panel p1 = new Panel();
		p1.setLayout(new FlowLayout());
		p1.add(new Label("Find:"));


		pattern = new TextField();
		pattern.setColumns(35);

		if (pat != null)
			pattern.setText(pat);

		p1.add(pattern);
		p1.doLayout();
		add("North", p1);

		Panel p2 = new Panel();
		fbutton = new Button("Find Next");
		fbutton.addActionListener(this);
		p2.add(fbutton);
		cbutton = new Button("Cancel");
		cbutton.addActionListener(this);
		p2.add(cbutton);
		add("South",p2);
		setSize(400,110);

		addWindowListener(this);
	}
	
    public void windowDeiconified(WindowEvent event) {}
    public void windowIconified(WindowEvent event) {}
    public void windowActivated(WindowEvent event) {}
    public void windowDeactivated(WindowEvent event) {}
    public void windowOpened(WindowEvent event) {}
    public void windowClosed(WindowEvent event) {}
    public void windowClosing(WindowEvent event)
	{
			editMan.CloseFindDialog(pattern.getText());
			dispose();
    }

	public void actionPerformed(ActionEvent evt)
	{
        if (evt.getSource() == cbutton)
		{
			editMan.CloseFindDialog(pattern.getText());
			dispose();
			return;
		}

        if (evt.getSource() == fbutton)
			if (!textEdit.getCanvas().find(pattern.getText()))
			{
				NotFound nf = new NotFound(textEdit);
				nf.show();
			}
	}

}

class ReplaceDialog extends Frame implements WindowListener, ActionListener
{
	private Button fbutton,rbutton,cbutton;
	private TextField pattern;
	private TextField replace;
	private TextEdit textEdit;
	private EditMan editMan;
	boolean foundOnce = false;

	public ReplaceDialog(TextEdit p, EditMan em, String pat, String rep)
	{
		super("Find and Replace");

		setBackground(Color.lightGray);

		textEdit = p;
		editMan = em;
		Panel p1 = new Panel();
		
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints constraints = new GridBagConstraints();
        p1.setLayout(gridbag);

		Label flabel = new Label("Find:");
        constraints.anchor = GridBagConstraints.NORTHWEST;
		gridbag.setConstraints(flabel, constraints);

		pattern = new TextField();
		pattern.setColumns(35);

		if (pat != null)
			pattern.setText(pat);

		constraints.gridwidth = GridBagConstraints.REMAINDER;
		gridbag.setConstraints(pattern, constraints);

        p1.add(flabel);
        p1.add(pattern);

		Label rlabel = new Label("Replace:");
        constraints.anchor = GridBagConstraints.WEST;
		constraints.gridwidth = 1;
        gridbag.setConstraints(rlabel, constraints);

		replace = new TextField();
		replace.setColumns(35);
		
		if (rep != null)
			replace.setText(rep);

		constraints.gridwidth = GridBagConstraints.REMAINDER;
        gridbag.setConstraints(replace, constraints);

        p1.add(rlabel);
        p1.add(replace);

		add("Center", p1);

		Panel p3 = new Panel();
		fbutton = new Button("Find Next");
		fbutton.addActionListener(this);
		p3.add(fbutton);
		rbutton = new Button("Replace");
		rbutton.addActionListener(this);
		p3.add(rbutton);
		cbutton = new Button("Cancel");
		cbutton.addActionListener(this);
		p3.add(cbutton);
		add("South",p3);
		setSize(400,120);

		addWindowListener(this);
	}
	
    public void windowDeiconified(WindowEvent event) {}
    public void windowIconified(WindowEvent event) {}
    public void windowActivated(WindowEvent event) {}
    public void windowDeactivated(WindowEvent event) {}
    public void windowOpened(WindowEvent event) {}
    public void windowClosed(WindowEvent event) {}
    public void windowClosing(WindowEvent event)
	{
			editMan.CloseReplaceDialog(pattern.getText(),replace.getText());
			dispose();
    }

	public void actionPerformed(ActionEvent evt)
	{
        if (evt.getSource() == cbutton)
		{
			editMan.CloseReplaceDialog(pattern.getText(),replace.getText());
			dispose();
			return;
		}

        if (evt.getSource() == fbutton)
			foundOnce = textEdit.getCanvas().find(pattern.getText());
		else
        if (evt.getSource() == rbutton)
		{
			if (!foundOnce)
			{
				TextCanvas textCanvas = textEdit.getCanvas();
				String selection= textCanvas.copy(false);
				if (selection != null)
					foundOnce = selection.equals(pattern.getText());
			}

			if (foundOnce)
				textEdit.getCanvas().paste(replace.getText());

			foundOnce = textEdit.getCanvas().find(pattern.getText());
		}
		
		if (!foundOnce)
		{
			NotFound nf = new NotFound(textEdit);
			nf.show();
		}
	}

}

class NotFound extends Dialog implements WindowListener, ActionListener
{
	private Button okButton;

	public NotFound(TextEdit te)
	{
		super(te,"",true);

		Panel north = new Panel();
		north.add(new Label("String not found"));
		add("North", north);

		okButton = new Button("OK");
		okButton.addActionListener(this);
		Panel south = new Panel();
		south.add(okButton);
		add("South", south);

		setSize(200,110);

		addWindowListener(this);
	}

    public void windowDeiconified(WindowEvent event) {}
    public void windowIconified(WindowEvent event) {}
    public void windowActivated(WindowEvent event) {}
    public void windowDeactivated(WindowEvent event) {}
    public void windowOpened(WindowEvent event) {}
    public void windowClosed(WindowEvent event) {}
    public void windowClosing(WindowEvent event)
	{
			dispose();
    }

	public void actionPerformed(ActionEvent evt)
	{
        if (evt.getSource() == okButton)
			dispose();
	}

}
