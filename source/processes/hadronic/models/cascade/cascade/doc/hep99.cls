%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{hep99}[1999/02/15 v0.1
 LaTeX2e class file for HEP99]
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse
\newif\if@openbib
\@openbibfalse
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{onecolumn}{\@twocolumnfalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{\@openbibtrue}
\ExecuteOptions{a4paper,10pt,twoside,twocolumn,final}
\ProcessOptions
\input{hep991\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@}
\pretolerance=8000
\tolerance=9500
\hfuzz=0.5pt
\vfuzz=2pt
\hbadness=8000
\vbadness=8000
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{3}
\renewcommand\topfraction{.8}
\setcounter{bottomnumber}{2}
\renewcommand\bottomfraction{.5}
\setcounter{totalnumber}{6}
\renewcommand\textfraction{0}
\renewcommand\floatpagefraction{.8}
\setcounter{dbltopnumber}{3}
\renewcommand\dbltopfraction{.8}
\renewcommand\dblfloatpagefraction{.8}
\def\ps@headings{%
      \def\@oddhead{\small\hfill\thepage}%
      \def\@evenhead{\small\thepage\hfill}%
      \def\@oddfoot{\hfill}
      \def\@evenfoot{\hfill}
      \let\@mkboth\markboth
      \def\sectionmark##1{\markboth {##1}{##1}}%
     }
\def\ps@myheadings{%
      \def\@oddhead{\hfill}%
      \def\@evenhead{\hfill}%
      \def\@oddfoot{\hfill}
      \def\@evenfoot{\hfill}
      \let\@mkboth\markboth
      \def\sectionmark##1{\markboth {##1}{##1}}%
     }

\def\title#1{\gdef\@title{{\LARGE\bf #1\par}}}
\def\@title{}

\def\maketitle{\thispagestyle{myheadings}\twocolumn[\@maketitle]}
\def\@maketitle{\begin{centred}
   \begin{minipage}{30pc}
   \begin{centred}
   \@title
   \vspace*{24\p@}%
   \@author
   \vspace*{11\p@}%
   \@address
   \end{centred}
   \vspace*{11\p@}%
   \@abstract
   \end{minipage}
   \end{centred}
   \vspace{22\p@ plus 22\p@ minus 6\p@}}

\def\author#1{\gdef\@author{{\Large\rm #1\par}}}
\def\@author{}

\def\address#1{\gdef\@address{{\small\rm #1\par}}}
\def\@address{}
\let\affil=\address

\def\abstract#1{\gdef\@abstract{\begin{centred}{\bf\abstractname}
    \end{centred}{\normalsize\rm#1\par}}}
\def\@abstract{}

\setcounter{secnumdepth}{3}
\newcounter {part}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand\thepart          {\Roman{part}}
\renewcommand\thesection       {\arabic{section}}
\renewcommand\thesubsection    {\thesection.\arabic{subsection}}
\renewcommand\thesubsubsection {\thesubsection.\arabic{subsubsection}}
\renewcommand\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\renewcommand\thesubparagraph  {\theparagraph.\arabic{subparagraph}}

\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\reset@font\normalsize\bfseries}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\reset@font\normalsize\itshape}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {3ex\@plus 1ex \@minus .2ex}%
                                     {-1.5ex \@plus -.2ex}%
                                     {\reset@font\small\bfseries}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\reset@font\small\bfseries}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\reset@font\small\bfseries}}

\def\@sect#1#2#3#4#5#6[#7]#8{\ifnum #2>\c@secnumdepth
     \let\@svsec\@empty\else
     \refstepcounter{#1}\edef\@svsec{\csname the#1\endcsname. }\fi
     \@tempskipa #5\relax
      \ifdim \@tempskipa>\z@
        \begingroup #6\relax
          \noindent{\hskip #3\relax\@svsec}{\interlinepenalty \@M #8\par}%
        \endgroup
       \csname #1mark\endcsname{#7}\addcontentsline
         {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                      \protect\numberline{\csname the#1\endcsname. }\fi
                    #7}\else
        \def\@svsechd{#6\hskip #3\relax  %% \relax added 2 May 90
                   \@svsec #8\csname #1mark\endcsname
                      {#7}\addcontentsline
                           {toc}{#1}{\ifnum #2>\c@secnumdepth \else
                             \protect\numberline{\csname the#1\endcsname. }\fi
                       #7}}\fi
     \@xsect{#5}}
%
\def\@ssect#1#2#3#4#5{\@tempskipa #3\relax
   \ifdim \@tempskipa>\z@
     \begingroup #4\noindent{\hskip #1}{\interlinepenalty \@M #5\par}\endgroup
   \else \def\@svsechd{#4\hskip #1\relax #5}\fi
    \@xsect{#3}}

\setlength\leftmargini  {2em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}

\setlength\leftmarginv  {.5em}
\setlength\leftmarginvi {.5em}

\setlength\leftmargin    {\leftmargini}
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\arabic{enumi}}
\renewcommand\theenumii{\alph{enumii}}
\renewcommand\theenumiii{\roman{enumiii}}
\renewcommand\theenumiv{\Alph{enumiv}}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{\normalfont\bfseries --}
\newcommand\labelitemiii{$\m@th\ast$}
\newcommand\labelitemiv{$\m@th\cdot$}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}

\newenvironment{verse}
               {\let\\=\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item[]}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item[]}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item[]}
               {\endlist}
\newenvironment{titlepage}
    {%
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \if@compatibility
        \setcounter{page}{0}
      \else
        \setcounter{page}{1}%
      \fi}%
    {\if@restonecol\twocolumn \else \newpage \fi
      \setcounter{page}{1}%
    }
\newcommand\appendix{\par
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  \renewcommand\thesection{\Alph{section}}}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand\theequation{\arabic{equation}}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\small{\bf #1.} #2}%
  \ifdim \wd\@tempboxa >\hsize
    \small{\bf #1.} #2\par
  \else
    \hbox to\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\newcommand{\cal}{\protect\pcal}
\newcommand{\pcal}{\@fontswitch{\relax}{\mathcal}}
\newcommand{\mit}{\protect\pmit}
\newcommand{\pmit}{\@fontswitch{\relax}{\mathnormal}}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg {2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\newcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{\uppercase{\contentsname}}{\uppercase{\contentsname}}}%
    \@starttoc{toc}%
    }
\newcommand\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{\@secpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \setlength\@tempdima{3em}%
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hbox to\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}
      \fi
    \endgroup
  \fi}
\newcommand\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty{\@secpenalty}%
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hbox to\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\newcommand\l@subsection   {\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand\l@paragraph    {\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand\l@subparagraph {\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename
      \@mkboth{\uppercase{\listfigurename}}%
              {\uppercase{\listfigurename}}}%
    \@starttoc{lof}%
    }
\newcommand\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename
      \@mkboth{\uppercase{\listtablename}}{\uppercase{\listtablename}}}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure

\newenvironment{centred}{\begincentred}{\endcentred}
\def\begincentred{\parindent=0pt\parskip=0pt
   \topsep=0pt\partopsep=0pt\parsep=0pt\itemsep=0pt\trivlist
   \centering\item[]}
\def\endcentred{\endtrivlist}

\newdimen\bibindent
\bibindent=1.5em
\newenvironment{thebibliography}[1]
     {\section*{\refname
        \@mkboth{\uppercase{\refname}}{\uppercase{\refname}}}%
      \list{\@biblabel{\arabic{enumiv}}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \if@openbib
              \advance\leftmargin\bibindent
              \itemindent -\bibindent
              \listparindent \itemindent
              \parsep \z@
            \fi
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\arabic{enumiv}}}%
      \if@openbib
        \renewcommand\newblock{\par}
      \else
        \renewcommand\newblock{\hskip .11em \@plus.33em \@minus.07em}%
      \fi
      \raggedright\topsep=0\p@\parsep=0\p@\partopsep=0\p@\itemsep=0\p@
      \clubpenalty4000\widowpenalty4000%
      \sfcode`\.=\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{}
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\uppercase{\indexname}}%
                        {\uppercase{\indexname}}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem  {\par\hangindent 40\p@}
\newcommand\subitem   {\par\hangindent 40\p@ \hspace*{20\p@}}
\newcommand\subsubitem{\par\hangindent 40\p@ \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}

\def\footnoterule{}%
\setcounter{footnote}{1}
\long\def\@makefntext#1{\parindent 1em\noindent 
 \makebox[1em][l]{\footnotesize\rm$\m@th{\fnsymbol{footnote}}$}%
 \footnotesize\rm #1}
\def\@makefnmark{\hbox{${\fnsymbol{footnote}}\m@th$}}
\def\@thefnmark{\fnsymbol{footnote}}
\def\footnote{\protect\pfootnote}
\def\pfootnote{\@ifnextchar[{\@xfootnote}{\stepcounter{\@mpfn}%
       \begingroup\let\protect\noexpand
       \xdef\@thefnmark{\thempfn}\endgroup
     \@footnotemark\@footnotetext}}
\def\@xfootnote[#1]{\setcounter{footnote}{#1}%
   \addtocounter{footnote}{-1}\footnote}
\def\@fnsymbol#1{\ifcase#1\or \dagger\or \ddagger\or \S\or
   \|\or \P\or ^{+}\or ^{\tsty *}\or ^{\tsty\#}\or 
   ^{\tsty \bullet}\or ^{\tsty\circ}\or ^{\tsty\diamond}\or 
   ^{\tsty\star}\or ^\times\or ^+_+\or 
   ^{\tsty*}_{\tsty*}\or ^\#_\#\or
   ^{\tsty\bullet}_{\tsty\bullet}\or
   ^{\tsty\circ}_{\tsty\circ}\or ^{\tsty\diamond}_{\tsty\diamond}\or 
   ^{\tsty\star}_{\tsty\star}\else\@ctrerr\fi\relax}

\newcommand\ftnote{\protect\pftnote}
\newcommand\pftnote[1]{\setcounter{footnote}{#1}%
   \addtocounter{footnote}{-1}\footnote}
\newcommand{\fntext}[1]{\def\@makefnmark{\hbox{#1}}\@footnotetext}


\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\refname{References}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\appendixname{Appendix}
\newcommand\abstractname{Abstract}
\newcommand\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{headings}
\pagenumbering{arabic}      % Arabic page numbers
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi

% Redefinition of \@citex to reduce the space between
% reference numbers in e.g. [3, 4] from full space to thin space
%
\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout{\string\citation{#2}}\fi
  \let\@citea\@empty
  \@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{,\penalty\@m\thinspace}%
     \def\@tempa##1##2\@nil{\edef\@citeb{\if##1\space##2\else##1##2\fi}}%
     \expandafter\@tempa\@citeb\@nil
     \@ifundefined{b@\@citeb}{{\reset@font\bf ?}\@warning
       {Citation `\@citeb' on page \thepage \space undefined}}%
     \hbox{\csname b@\@citeb\endcsname}}}{#1}}
%
\def\@cite#1#2{[{#1\if@tempswa ,\,#2\fi}]}
%
\def\thereferences{\list{}{\topsep=0\p@\parsep=0\p@
 \partopsep=0\p@\itemsep=0\p@\labelsep=0\p@\itemindent=-18\p@
\labelwidth=0\p@\leftmargin=18\p@
}\small\rm\raggedright
\def\newblock{\ }
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax
}
%
\let\endthereferences=\endlist
%
% Macro to used for references in the Harvard system
%
\newenvironment{harvard}{\list{}{\topsep=0\p@\parsep=0\p@
\partopsep=0\p@\itemsep=0\p@\labelsep=0\p@\itemindent=-18\p@
\labelwidth=0\p@\leftmargin=18\p@
}\small\rm\raggedright
\def\newblock{\ }
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax}{\endlist}
%
\def\refs{\begin{harvard}}
\def\endrefs{\end{harvard}}
%
\newcommand{\CQG}{{\em Class. Quantum Grav.} }
\newcommand{\CTM}{{\em Combust. Theory Modelling\/} }
\newcommand{\DSE}{{\em Distrib. Syst. Engng\/} }
\newcommand{\EJP}{{\em Eur. J. Phys.} } 
\newcommand{\HPP}{{\em High Perform. Polym.} }              % added 4/5/93
\newcommand{\IP}{{\em Inverse Problems\/} }
\newcommand{\JHM}{{\em J. Hard Mater.} }                    % added 4/5/93
\newcommand{\JO}{{\em J. Opt.} }
\newcommand{\JOA}{{\em J. Opt. A: Pure Appl. Opt.} }
\newcommand{\JOB}{{\em J. Opt. B: Quantum Semiclass. Opt.} }
\newcommand{\JPA}{{\em J. Phys. A: Math. Gen.} }
\newcommand{\JPB}{{\em J. Phys. B: At. Mol. Phys.} }      %1968-87
\newcommand{\jpb}{{\em J. Phys. B: At. Mol. Opt. Phys.} } %1988 and onwards
\newcommand{\JPC}{{\em J. Phys. C: Solid State Phys.} }   %1968--1988
\newcommand{\JPCM}{{\em J. Phys.: Condens. Matter\/} }    %1989 and onwards
\newcommand{\JPD}{{\em J. Phys. D: Appl. Phys.} }
\newcommand{\JPE}{{\em J. Phys. E: Sci. Instrum.} }
\newcommand{\JPF}{{\em J. Phys. F: Met. Phys.} }
\newcommand{\JPG}{{\em J. Phys. G: Nucl. Phys.} }         %1975--1988
\newcommand{\jpg}{{\em J. Phys. G: Nucl. Part. Phys.} }   %1989 and onwards
\newcommand{\MSMSE}{{\em Modelling Simulation Mater. Sci. Eng.} }
\newcommand{\MST}{{\em Meas. Sci. Technol.} }                 %1990 and onwards
\newcommand{\NET}{{\em Network: Comput. Neural Syst.} }
\newcommand{\NJP}{{\em New J. Phys.} }
\newcommand{\NL}{{\em Nonlinearity\/} }
\newcommand{\NT}{{\em Nanotechnology} }
\newcommand{\PAO}{{\em Pure Appl. Opt.} }
\newcommand{\PM}{{\em Physiol. Meas.} }                        % added 4/5/93
\newcommand{\PMB}{{\em Phys. Med. Biol.} }
\newcommand{\PPCF}{{\em Plasma Phys. Control. Fusion\/} }      % added 4/5/93
\newcommand{\PSST}{{\em Plasma Sources Sci. Technol.} }
\newcommand{\PUS}{{\em Public Understand. Sci.} }
\newcommand{\QO}{{\em Quantum Opt.} }
\newcommand{\QSO}{{\em Quantum Semiclass. Opt.} }
\newcommand{\RPP}{{\em Rep. Prog. Phys.} }
\newcommand{\SLC}{{\em Sov. Lightwave Commun.} }               % added 4/5/93 
\newcommand{\SMS}{{\em Smart Mater. Struct.} }      
\newcommand{\SST}{{\em Semicond. Sci. Technol.} }
\newcommand{\SUST}{{\em Supercond. Sci. Technol.} }
\newcommand{\WRM}{{\em Waves Random Media\/} }
%
% Other commonly quoted journals
%
\newcommand{\AC}{{\em Acta Crystallogr.} }
\newcommand{\AM}{{\em Acta Metall.} }
\newcommand{\AP}{{\em Ann. Phys., Lpz.} }
\newcommand{\APNY}{{\em Ann. Phys., NY\/} }
\newcommand{\APP}{{\em Ann. Phys., Paris\/} }
\newcommand{\CJP}{{\em Can. J. Phys.} }
\newcommand{\JAP}{{\em J. Appl. Phys.} }
\newcommand{\JCP}{{\em J. Chem. Phys.} }
\newcommand{\JJAP}{{\em Japan. J. Appl. Phys.} }
\newcommand{\JP}{{\em J. Physique\/} }
\newcommand{\JPhCh}{{\em J. Phys. Chem.} }
\newcommand{\JMMM}{{\em J. Magn. Magn. Mater.} }
\newcommand{\JMP}{{\em J. Math. Phys.} }
\newcommand{\JOSA}{{\em J. Opt. Soc. Am.} }
\newcommand{\JPSJ}{{\em J. Phys. Soc. Japan\/} }
\newcommand{\JQSRT}{{\em J. Quant. Spectrosc. Radiat. Transfer\/} }
\newcommand{\NC}{{\em Nuovo Cimento\/} }
\newcommand{\NIM}{{\em Nucl. Instrum. Methods\/} }
\newcommand{\NP}{{\em Nucl. Phys.} }
\newcommand{\PL}{{\em Phys. Lett.} }
\newcommand{\PR}{{\em Phys. Rev.} }
\newcommand{\PRL}{{\em Phys. Rev. Lett.} }
\newcommand{\PRS}{{\em Proc. R. Soc.} }
\newcommand{\PS}{{\em Phys. Scr.} }
\newcommand{\PSS}{{\em Phys. Status Solidi\/} }
\newcommand{\PTRS}{{\em Phil. Trans. R. Soc.} }
\newcommand{\RMP}{{\em Rev. Mod. Phys.} }
\newcommand{\RSI}{{\em Rev. Sci. Instrum.} }
\newcommand{\SSC}{{\em Solid State Commun.} }
\newcommand{\ZP}{{\em Z. Phys.} }
%%%%%%%%%%%%%%%%%%%%%
% Tables rules      %
%%%%%%%%%%%%%%%%%%%%%

\newcommand{\boldarrayrulewidth}{1\p@} 
% Width of bold rule in tabular environment.

\def\bhline{\noalign{\ifnum0=`}\fi\hrule \@height  
\boldarrayrulewidth \futurelet \@tempa\@xhline}

\def\@xhline{\ifx\@tempa\hline\vskip \doublerulesep\fi
      \ifnum0=`{\fi}}

\newcommand{\ms}{\noalign{\vspace{3\p@ plus2\p@ minus1\p@}}}
%
% Rules for tables with extra space around
%
\newcommand{\br}{\ms\bhline\ms}
\newcommand{\mr}{\ms\hline\ms}
%
\def\0{\hbox{\phantom{0}}}%
\def\m{\hbox{\llap{$-$}}}%
 
\def\endcolumn{\parfillskip=0pt\par\newpage
   \noindent\parfillskip=0pt plus 1fil}

\endinput
%% End of file `hep99.cls'.

