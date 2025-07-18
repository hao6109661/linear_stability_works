%% muthesis Class for LaTeX2e.
%% Last update 12th August 2013: amended evensidemargin for duplex

%%
%% Time-stamp: <July 23, 1997 11:33:00 BST graham@cs.man.ac.uk>
%%
%% This class file has been cobbled together for use of postgraduate
%% students in the Department of Computer Science at the University of
%% Manchester. Any comments or problems about its use should be made
%% via the local newsgroup man.cs.tex.
%%
%% The \school command can be used to define as a school other than
%% Mathematics, but no support is provided for users outside the
%% School of Mathematics.
%%
%% History
%%    Acquired for uk.ac.oxford.prg by Stephen Page <sdpage> on 16 Apr
%%    86 from Howard Trickey <trickey@arpa.su-aimvax> , then hacked by
%%    Michael Fisher into Manchester University style.
%%
%%    Hacked from muthesis.sty by Rhod Davies (rmd@uk.ac.man.cs.r5) to
%%    match the latest ridiculous rules from the exams office (March
%%    1990) as given below.
%%
%%    Further modified by David Carlisle (DPC) September 1993.
%%
%%    Further modified by David Carlisle (DPC) September 1994.
%%      Add [PhD] option (default) and [MSc] option.
%%      Make into LaTeX2e Class file.
%%
%%    Further modified by Graham Gough at various dates
%%    Added [MPhil option.

%%    Modified William Dyke (WBD) July 1997 so that declaration is
%%    correct for use in departments other than Computer Science

%%    Modified by Chris Paul March 2006 to reflect new University 
%%    Schools and Faculty setup

%%    Modified by Chris Paul August 2013 to support duplex printed margins

%% Note that this makes no attempt to deal with multi-volume thesis -
%% if you do one that long you deserve everything you get.
%% 1) Page numbering starts on the title page and is arabic (i.e. title
%%       page is no.1), all pages are numbered (declaration, figures
%%       etc included), but the title page itself shall not show the
%%       page number!
%%
%%
%% 2) The preface will consist of the following sections in the
%%       following order!
%%       i  )    Title Page
%%       ii )    Contents page (actually we shall be fancy and have a
%%               table of contents, list of tables and list of figures.
%%               This is the correct place to put these extras)
%%       iii)    A short Abstract
%%       iv)     Declaration of honesty
%%       v)      Statement of qualifications and research
%%       vi)     Anything else you want, such as acknowledgements
%%
%%
%% 3) The ordering of the items on the title page has been changed
%%       slightly from the old muthesis style - it is now
%%               Statement - A thesis ..... Faculty of ...
%%               Year of Presentation
%%               Candidate's name
%%               Candidate's department
%%
%%
%% 4) Certain members of this university never having heard of
%%       typesetting, (i.e. the ones who will be the first against the
%%       wall when the revolution comes), the whole thing (apart from
%%       footnotes and quotations must be doublespaced. By a strang
%%       coincidence, a set of the thesis rules fell through a time
%%       warp from the future and landed on my desk. These note that
%%       single spacing is allowed in typeset thesis, with a footnote
%%       that those who were in favour of double spacing were the first
%%       against the wall when the revolution came.
%%               If everyone who submits a thesis goes and argues with
%%       them about this they may eventually get the hint, but this is
%%       unlikely. To get this we set baselinestretch to 1.5 which
%%       looks about right, and they have never noticed the difference.
%%       By all means set this to something sensible while producing
%%       draft copies, but remember DOUBLE SPACE before you hand in.
%%
%%
%% 5) The binding (left hand) margin must be no less than 40mm, other
%%       margins must be no less than 15mm, with page numbers falling
%%       within these margins. The binding margin is thus rather large,
%%       even when you take the binding (normally about 10mm) into account.
%%       Consequently you cannot reduce the other margins without
%%       making the page look unbalanced. Consequently we have a 40mm
%%       binding margin and 15 mm margin on all other sides. NOTE, this
%%       does not allow for any marginal notes. I recommend that you
%%       avoid using these, since a) there is not much space for them
%%       and they will format badly and b) they will violate the margin
%%       rules.
%%
%%
%% The rest of these notes will not be found in the written submission
%% rules, but come from converstions with the exams office.
%%
%%
%% 6)    They have not particular rules about the format of the
%%       bibliography, so long as you are consistent in your naming
%%       scheme. I would recommend the alpha style.
%%
%%
%% 7)    They do not seem to have any rules about appendices.
%%
%%
%% 8)    As before, any extra things you want to put in -
%%       acknowledgements etc, should come after the mandatory preface
%%       sections and before the first chapter.
%%
%% This style can usefully be used in conjunction with the doublespace
%% style which introduces a single spaced environment.
%%
%% It is not compatible with the drafthead style.


\ProvidesClass{muthesis}[2013/08/12 MU Thesis Class]

\DeclareOption{PhD}{\def\degreetitle{Doctor of Philosophy}}
\DeclareOption{MSc}{\def\degreetitle{Master of Science}}
\DeclareOption{MPhil}{\def\degreetitle{Master of Philosophy}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\def\@msc{Master of Science}

\ExecuteOptions{PhD}

\ProcessOptions


\LoadClass[a4paper]{report}

%% Manchester University PhD thesis style -- modifications to the report style
%%
\long\def\comment#1{}
\comment{

  Example of use:
    \documentclass[12pt,MSc]{muthesis}
    \begin{document}
    \title{How to Write Theses\\
            With Two Line Titles}
    \author{John Henry Candidate}
    \principaladviser{John Parker}
    \firstreader{John Green}
    \secondreader{John BigBooty\\(Another Department)}

    \beforeabstract
    \prefacesection{Abstract}
        This thesis answers all the questions you were afraid to ask ...
    \afterabstract
    \prefacesection{Acknowledgements}
        I would like to thank...
    \afterpreface

    \chapter{Introduction}
         ...
    \chapter{Conclusions}
         ...
    \appendix
    \chapter{A Long Proof}
         ...
%%%%    \addcontentsline{toc}{chapter}{Bibliography} %%% NO!
%%%% This \addcontentsline never worked anyway (it made the page
%%%% number in the table of contents one too small). It is not needed
%%%% at all with the sept 93 version. DPC.
    \bibliographystyle{alpha}
    \bibliography{mybib}


Documentation:
    It sets the margins and interline spacing and
    disallows page breaks at hyphens.

    The \beforeabstract command creates the title page, a copyright
    page (default), and a signature page (optional), the tables of
    contents, tables and figures. The signature page is not mentioned
    in the rules for submission, and by implication are not allowed.

    The user then includes the preface using the
    \prefacesection{Abstract} command.

    Next should come an \afterabstract command to insert the declaration.

    Then the user should put preface section(s), using the
    \prefacesection{section title} command. The first of these should
    be the statement of education and research, others may be
    acknowledgements etc.

    The \afterpreface command then  sets things up to start
    the main body (continuing the sequence of page numbering, but
    putting page headings etc.).

    Having one of \beforepreface or \afterpreface without the other
    may cause problems, so udon't do it.

    The following commands can control what goes in the front matter
    material:

        \title{thesis title}
        \author{author's name}
        \school{author's School}
        \faculty{author's Faculty}
        \principaladviser{the principal advisor's name}
           (or \principaladvisor, if you prefer advisor spelled with o)
        \firstreader{the first reader's name}
        \secondreader{the second reader's name}
        \submitdate{month year in which submitted to GPO}
                - date LaTeX'd if omitted
        \copyrightyear{year degree conferred (next year if submitted in Dec.)}
                - year LaTeX'd (or next year, in December) if omitted
        \copyrighttrue or \copyrightfalse
                - produce or don't produce a copyright page (true by default)
        \figurespagetrue or \figurespagefalse
                - produce or don't produce a List of Figures page
                  (true by default)
        \tablespagetrue or \tablespagefalse
                - produce or don't produce a List of Tables page
                  (true by default)
        \sigpagtrue or \sigpagfalse
                - produce or don't produce a signature page
                  (again false by default, as Manchester theses don't
                   require this. However, it's left in so that you can
                   your own `examiners page' if you like!)

        If you set all these to false, then perhaps you don't want any
        preface at all, so perhaps you hould leave out both \beforepreface
        and \afterpreface -).

This style uses interline spacing that is 1.5 times normal (about double
spaced), except in the figure and table environments where normal spacing
is used.

That can be changed by doing:
    \renewcommand{\baselinestretch}{1.6}
(or whatever you want instead of 1.6)
Unfortunately, that command probably won't take effect unless it comes
before the \begin{document} or after something like \small \normalsize.

The example above shows the 12pt substyle being used.  This seems to give
acceptable looking results, but it may be omitted to get smaller print.

}
%% We need >15mm  margins except on the binding edge, where it is >40mm
\oddsidemargin 15truemm \evensidemargin -8truemm
\marginparwidth 40pt \marginparsep 10pt
\topmargin -11.5truemm \headsep 7truemm
\textheight 245truemm \textwidth 154truemm

%% Disallow page breaks at hyphens (this will give some underfull vbox's,
%% so an alternative is to use \brokenpenalty=100 and manually search
%% for and fix such page breaks)
\brokenpenalty=10000

%% Use 1.5 times the normal baseline-to-baseline skip (approx double space)
\def\stretch{1.5}
\renewcommand{\baselinestretch}{\stretch}

%% Redefine the macro used for floats (including figures and tables)
%% so that single spacing is used.
\let\latex@xfloat=\@xfloat
\def\@xfloat #1[#2]{
        \latex@xfloat #1[#2]
        \def\baselinestretch{1}
        \@normalsize \normalsize
}

%% Redefine the macro used for footnotes to use single spacing
\long\def\@footnotetext#1{\insert\footins{\def\baselinestretch{1}\footnotesize
    \interlinepenalty\interfootnotelinepenalty
    \splittopskip\footnotesep
    \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
    \hsize\columnwidth \@parboxrestore
   \edef\@currentlabel{\csname p@footnote\endcsname\@thefnmark}\@makefntext
    {\rule{\z@}{\footnotesep}\ignorespaces
      #1\strut}}}

%% \author, \title are defined in report; here are the rest of the
%% front matter defining macros
%% Note: If you really want, you can define your own macro for, say
%% \supervisor , just by copying the definitions of \principaladvisor
%%
\def\school#1{\gdef\@school{#1}}
\def\faculty#1{\gdef\@faculty{#1}}
\def\principaladviser#1{\gdef\@principaladviser{#1}}
\def\advis@r{Advisor} % default spelling---we are English!
%% Obviously, we could add a \supervisor macro
\def\principaladvisor#1{\gdef\@principaladviser{#1}\gdef\advis@r{Advisor}}
\def\firstreader#1{\gdef\@firstreader{#1}}
\def\secondreader#1{\gdef\@secondreader{#1}}
\def\submitdate#1{\gdef\@submitdate{#1}}
\def\copyrightyear#1{\gdef\@copyrightyear{#1}} % \author, \title in report
\def\@title{}
\def\@author{}
\def\@school{}
\def\@faculty{}
\def\@principaladviser{}\def\@firstreader{}\def\@secondreader{}
\def\@submitdate{\number\the\year}
%% Stanford says that December theses should have the next year as the
%% copyright year
\ifnum\month=12
    \@tempcnta=\year \advance\@tempcnta by 1
    \edef\@copyrightyear{\number\the\@tempcnta}
\else
    \def\@copyrightyear{\number\the\year}
\fi
\newif\ifcopyright \newif\iffigurespage \newif\iftablespage \newif\ifsigpage
\copyrighttrue \figurespagetrue \tablespagetrue \sigpagefalse


\ifx\degreetitle\@msc
    \def\@thesis{dissertation}
    \else
    \def\@thesis{thesis}
\fi

%% If you're doing an MSc, I'm sure you can guess which parts to change
%% below; if not then you shouldn't be doing one.
\def\titlep{%
        \thispagestyle{empty}%
        \null\vskip1in%
        \begin{center}
                \Huge\uppercase\expandafter{\@title}
        \end{center}
        \vfill
        \begin{center}
           \sc A first-year report submitted to the University of Manchester\\
           for the Doctor of Philosophy\\
           in the Faculty of \expandafter{\@faculty}
        \end{center}
        \vfill
        \begin{center}
                \@submitdate
        \end{center}
        \vfill
        \begin{center}
                {\bf\@author}\\
                \ifx\@empty\@school
                \else
                Department of \expandafter{\@school}
                \fi
        \end{center}\vskip.5in\newpage}

\def\copyrightpage{%
\ifx\degreetitle\@msc
\prefacesection{Intellectual Property Statement}
\else
\prefacesection{Copyright Statement}
\fi

\begin{list}{\labelitemi}{\leftmargin=1em}
\item[{\bf i.}] The author of this \expandafter{\@thesis} (including any
appendices and/or schedules to this \expandafter{\@thesis}) owns certain
copyright or related rights in it (the ``Copyright'') and s/he has given
The University of Manchester certain rights to use such Copyright,
including for administrative purposes.
\item[{\bf ii.}] Copies of this \expandafter{\@thesis}, either in full or in
extracts and whether in hard or electronic copy, may be made {\bf only} in
accordance with the Copyright, Designs and Patents Act 1988 (as amended) and
regulations issued under it or, where appropriate, in accordance with
licensing agreements which the University has
\ifx\degreetitle\@msc
entered into.
\else
from time to time.
\fi
This page must form part of any such copies made.
\item[{\bf iii.}] The ownership of certain Copyright, patents, designs, trade
marks and other intellectual property (the ``Intellectual Property'') and any
reproductions of copyright works in the \expandafter{\@thesis}, for example
graphs and tables  (``Reproductions''), which may be described in this
\expandafter{\@thesis}, may not be owned by the author and may be owned by
third parties. Such Intellectual Property and Reproductions cannot and must
not be made available for use without the prior written permission of the
owner(s) of the relevant Intellectual Property and/or Reproductions.
\item[{\bf iv.}] Further information on the conditions under which disclosure,
publication and commercialisation of this \expandafter{\@thesis}, the
Copyright and any Intellectual Property and/or Reproductions described in
it may take place is available in the University IP Policy (see
http://documents.manchester.ac.uk/DocuInfo.aspx?DocID=487), in any relevant
\ifx\degreetitle\@msc
Dissertation
\else
Thesis
\fi
 restriction declarations deposited in the University Library, The University
Library's regulations (see
http://www.manchester.ac.uk/library/ab\-outus/regul\-ations) and in The 
University's 
\ifx\degreetitle\@msc
Guidance on Presentation of Dissertations.
\else
Policy on Presentation of Theses.
\fi
\end{list}
\vfill\newpage}

\long\def\signature#1{%
\begin{center}
\begin{minipage}{4in}
\parindent=0pt
I certify that I have read this thesis and that in my opinion
it is fully adequate, in scope and in quality, as a dissertation for the degree
of \degreetitle.
\par
\vspace{.5in}
\hbox to 4in{\hfil\shortstack{\vrule width 3in height 0.4pt\\#1}}
\end{minipage}
\end{center}}

\def\signaturepage{%
        \signature{\@principaladviser\\(Principal \advis@r)}
        \vfill
        \signature\@firstreader
        \vfill
        \signature\@secondreader
        \vfill
%% Again, the real format will have to be added here%
        \begin{center}
        \begin{minipage}{4in}
        Approved for the University Committee on Graduate Studies:\par
        \vspace{.5in}
        \hbox to 4in{\hfil\shortstack{\vrule width 3in height 0.4pt\\
                        Dean of Graduate Studies \& Research}}
        \end{minipage}
        \end{center}}

\def\declaration{No portion of the work referred to in the
\expandafter{\@thesis} has been submitted in support of an application for
another degree or qualification of this or any other university or other
institute of learning.}
\def\declpage{%
\prefacesection{Declaration}
        \begin{center}
        \begin{minipage}{4in}
            \parindent=0pt
            \declaration 
            \par
        \end{minipage}%
        \end{center}%
        \vfill}

\def\beforeabstract{
        \pagenumbering{arabic}
        \pagestyle{plain}
        \titlep
        \ifsigpage\signaturepage\fi
        \tableofcontents
	\vfill
	\centerline{Word count \wordcount}
        \newpage
        \iftablespage
                \addvspace{10pt}
                \listoftables
                \newpage
        \fi
        \iffigurespage
                \addvspace{10pt}
                \listoffigures
                \newpage
        \fi
        \addcontentsline{toc}{chapter}{Abstract}
        \ifx\degreetitle\@msc
             \chapter*{Abstract}
             \else
             \begin{center}\LARGE\bf The University of Manchester\end{center}
             \noindent{\bf\@author}\\
             \noindent{\bf\degreetitle}\\
             \noindent{\bf\@title}\\
             \noindent{\bf\today}\\[0.1\baselineskip]
             \def\baselinestretch{1.0}\normalsize
        \fi
        }

\def\afterabstract{
%% \newpage inserted by DPC 15/9/93
%% make the declaration come on the page after the abstract.
        \def\baselinestretch{\stretch}\normalsize
        \newpage
        \declpage
        \ifcopyright\copyrightpage\fi % moved to conform with 95 regs
        }

%% The declaration page above is the standard
%% declaration, but if you want the other sort ( where you admit
%% that parts of your thesis have been used elsewhere) then you'll
%% have to hack it yourself -- just replace the text by what you'd like
%%
\def\prefacesection#1{%
        \chapter*{#1}
        \addcontentsline{toc}{chapter}{#1}}

%% The \addvspace business below is to avoid the extra space that gets put
%% in the tables of figures and tables between chapters (it makes no sense
%% now that the numbering is continuous)
\def\afterpreface{\newpage
        \pagestyle{headings}}

\pagestyle{headings}

%% This produces the chapter title at the top of each page. It can be
%% changed in the obvious way, e.g., \pagestyle{empty} or \pagestyle{plain}
%%

%%% DPC 15/9/93 Make the bibliography generate an entry in the table
%%% of contents.

\def\thebibliography#1{\chapter*{\bibname\@mkboth
  {\uppercase{\bibname}}{\uppercase{\bibname}}}%
  \addcontentsline{toc}{chapter}{\bibname}%
  \list
  {\@biblabel{\arabic{enumiv}}}{\settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth
    \advance\leftmargin\labelsep
    \usecounter{enumiv}%
    \let\p@enumiv\@empty
    \def\theenumiv{\arabic{enumiv}}}%
    \def\newblock{\hskip .11em plus.33em minus.07em}%
    \sloppy\clubpenalty4000\widowpenalty4000
    \sfcode`\.=\@m}

\def\endthebibliography{%
  \def\@noitemerr{\@warning{Empty `thebibliography' environment}}%
  \endlist}
