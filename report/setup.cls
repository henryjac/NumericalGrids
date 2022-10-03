\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{setup}[Class for homework]

\LoadClass[titlepage]{article}

\RequirePackage{scalerel,stackengine}

\RequirePackage{pdfpages}
\RequirePackage{parskip}
\setlength{\parindent}{0pt}

\RequirePackage{minted}
\RequirePackage{xcolor}
\RequirePackage{subfigure}
\RequirePackage{svg}
\RequirePackage{float}
\RequirePackage{graphicx}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage[margin=2cm]{geometry}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
\RequirePackage{amssymb}
\RequirePackage{amsthm}
\RequirePackage{mathtools}

% Commands
\newcommand{\pth}[1]{\left(#1\right)}
\newcommand{\brk}[1]{\left[#1\right]}
\newcommand{\brc}[1]{\left\{#1\right\}}
\newcommand{\eps}{\varepsilon}
\newcommand{\mtr}[1]{\begin{bmatrix}#1\end{bmatrix}}
\newcommand{\Tx}{\Tilde{x}}
\newcommand{\case}[1]{\begin{cases}#1\end{cases}}
\newcommand{\R}{\mathbb{R}}
\def\ddfrac#1#2{\displaystyle\frac{\displaystyle #1}{\displaystyle #2}}

\newcommand{\pdv}[3][]{\frac{\partial^{#1}#2}{\partial #3^{#1}}}
\newcommand{\pddv}[3]{\frac{\partial^2#1}{\partial #2\partial #3}}
