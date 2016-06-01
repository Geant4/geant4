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
import java.util.*;

public class PrintMan extends Thread
{

	private PrintJob pjob;
	private Vector spool;
	private Font font;
	private Frame parent;
	private String jobName;
	private int tabsize;

	private final int PRINT_PROGRESS = 10000;

	public PrintMan(Frame p, String n, Vector lines, Font f, int t)
	{
		int i,max,part,tens;

		parent = p;
		jobName = n;
		font = f;
		tabsize = t;

		max = lines.size();
		spool = new Vector(max);
		tens = max / 10;
		part = 0;

		Progress spoolProgress = null;

		if (max > PRINT_PROGRESS)
		{
			spoolProgress = new Progress(parent,"Saving print copy ... ");
			spoolProgress.show();
		}

		for (i=0;i<max;i++)
		{
			spool.addElement( new String((String)lines.elementAt(i)) );
			if (max > PRINT_PROGRESS)
				part++;
			if (part > tens)
			{
				spoolProgress.update( (int)(100 * i / max) );
				part = 0;
			}
		}

		if (spoolProgress != null)
			spoolProgress.dispose();
	}

	public void run()
	{
		pjob = parent.getToolkit().getPrintJob(parent, jobName, null);

		if (pjob == null)
			return;

		Graphics pg = pjob.getGraphics();

		if (pg != null)
		{
			String nextLine;
			int pageHeight = pjob.getPageDimension().height;
			pg.setFont(font);
			FontMetrics fm = pg.getFontMetrics(font);
			int fontHeight = fm.getHeight();
			int fontDescent = fm.getDescent();
			int curHeight = 0;

			for (int i = 0; i < spool.size(); i++ )
			{
				curHeight += fontHeight;
				if (curHeight > pageHeight)
				{
					pg.dispose();
					pg = pjob.getGraphics();
					if (pg == null)
						return;
					pg.setFont(font);
					curHeight = fontHeight;
				}
				nextLine = detabbed((String)spool.elementAt(i),tabsize);
				pg.drawString(nextLine, 0, curHeight - fontDescent);
			}

			pg.dispose();
		}

		pjob.end();
	}

	private String detabbed(String s, int tabSize)
	{
		if (s.indexOf('\t') < 0)
			return s;

		char c;
		String t = new String("");
		int j=0;
		int tabs;
		int max = s.length();

		for (int i=0; i<max; i++)
		{
			c = s.charAt(i);
			if (c == '\t')
			{
				tabs = tabSize - (j % tabSize);
				j += tabs;
				while (tabs-- > 0) t = t + ' ';
			}
			else
			{
				t = t + c;
				j++;
			}
		}

		return t;
	}

}
