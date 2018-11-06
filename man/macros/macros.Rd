% chemical formulas
\newcommand{\CO2}{\ifelse{latex}{\eqn{\mathrm{CO_{2}}}}{\ifelse{html}{\out{CO<sub>2</sub>}}{CO2}}}
\newcommand{\NH3}{\ifelse{latex}{\eqn{\mathrm{NH_{3}}}}{\ifelse{html}{\out{NH<sub>3</sub>}}{NH3}}}
\newcommand{\H2S}{\ifelse{latex}{\eqn{\mathrm{H_{2}S}}}{\ifelse{html}{\out{H<sub>2</sub>S}}{H2S}}}
\newcommand{\H2O}{\ifelse{latex}{\eqn{\mathrm{H_{2}O}}}{\ifelse{html}{\out{H<sub>2</sub>O}}{H2O}}}
\newcommand{\O2}{\ifelse{latex}{\eqn{\mathrm{O_{2}}}}{\ifelse{html}{\out{O<sub>2</sub>}}{O2}}}
\newcommand{\Hplus}{\ifelse{latex}{\eqn{\mathrm{H^{+}}}}{\ifelse{html}{\out{H<sup>+</sup>}}{H+}}}
\newcommand{\eminus}{\ifelse{latex}{\eqn{e^{-}}}{\ifelse{html}{\out{<i>e</i><sup>-</sup>}}{e-}}}
\newcommand{\Mgplus2}{\ifelse{latex}{\eqn{\mathrm{Mg^{+2}}}}{\ifelse{html}{\out{Mg<sup>+2</sup>}}{Mg+2}}}
\newcommand{\H3PO4}{\ifelse{latex}{\eqn{\mathrm{H_{3}PO_{4}}}}{\ifelse{html}{\out{H<sub>3</sub>PO<sub>4</sub>}}{H3PO4}}}
\newcommand{\Fe2O3}{\ifelse{latex}{\eqn{\mathrm{Fe_{2}O_{3}}}}{\ifelse{html}{\out{Fe<sub>2</sub>O<sub>3</sub>}}{Fe2O3}}}
\newcommand{\CH4}{\ifelse{latex}{\eqn{\mathrm{CH_{4}}}}{\ifelse{html}{\out{CH<sub>4</sub>}}{CH4}}}
\newcommand{\NH4plus}{\ifelse{latex}{\eqn{\mathrm{NH_{4}^{+}}}}{\ifelse{html}{\out{NH<sub>4</sub><sup>+</sup>}}{NH4+}}}
\newcommand{\H2}{\ifelse{latex}{\eqn{\mathrm{H_{2}}}}{\ifelse{html}{\out{H<sub>2</sub>}}{H2}}}

% subscripted numbers
\newcommand{\s2}{\ifelse{latex}{\eqn{_2}}{\ifelse{html}{\out{<sub>2</sub>}}{2}}}
\newcommand{\s3}{\ifelse{latex}{\eqn{_3}}{\ifelse{html}{\out{<sub>3</sub>}}{3}}}
\newcommand{\s4}{\ifelse{latex}{\eqn{_4}}{\ifelse{html}{\out{<sub>4</sub>}}{4}}}
\newcommand{\s}{\ifelse{latex}{\eqn{_{#1}}}{\ifelse{html}{\out{<sub>#1</sub>}}{#1}}}

% superscript numbers
\newcommand{\S}{\ifelse{latex}{\eqn{^{#1}}}{\ifelse{html}{\out{<sup>#1</sup>}}{^#1}}}

% other common variables
\newcommand{\T}{\ifelse{latex}{\eqn{T}}{\ifelse{html}{\out{<I>T</I>}}{T}}}
\newcommand{\P}{\ifelse{latex}{\eqn{P}}{\ifelse{html}{\out{<I>P</I>}}{P}}}
% italicize any letter
\newcommand{\I}{\ifelse{latex}{\eqn{#1}}{\ifelse{html}{\out{<I>#1</I>}}{#1}}}
\newcommand{\logK}{\ifelse{latex}{\eqn{\log K}}{\ifelse{html}{\out{log<i>K</i>}}{log K}}}
\newcommand{\logfO2}{\ifelse{latex}{\eqn{\log f_{\mathrm{O_{2}}}}}{\ifelse{html}{\out{log<i>f</i><sub>O<sub>2</sub></sub>}}{logfO2}}}
\newcommand{\logfH2}{\ifelse{latex}{\eqn{\log f_{\mathrm{H_{2}}}}}{\ifelse{html}{\out{log<i>f</i><sub>H<sub>2</sub></sub>}}{logfH2}}}
\newcommand{\logaH2}{\ifelse{latex}{\eqn{\log a_{\mathrm{H_{2}}}}}{\ifelse{html}{\out{log<i>a</i><sub>H<sub>2</sub></sub>}}{logaH2}}}
\newcommand{\logaH2O}{\ifelse{latex}{\eqn{\log a_{\mathrm{H_{2}O}}}}{\ifelse{html}{\out{log<i>a</i><sub>H<sub>2</sub>O</sub>}}{logaH2O}}}
\newcommand{\loga}{\ifelse{latex}{\eqn{\log a}}{\ifelse{html}{\out{log<i>a</i>}}{log a}}}
\newcommand{\ZC}{\ifelse{latex}{\eqn{Z_\mathrm{C}}}{\ifelse{html}{\out{<I>Z</I><sub>C</sub>}}{ZC}}}
% use \nH2O{&#772;} to call this macro (the html code can't be defined in the macro,
% which interprets '#' followed by a number as a placeholder for an argument)
\newcommand{\nH2O}{\ifelse{latex}{\eqn{\bar{n}_{\mathrm{H_{2}O}}}}{\ifelse{html}{\out{<i>n</i>#1<sub>H<sub>2</sub>O</sub>}}{nH2O}}}
\newcommand{\Cp}{\ifelse{latex}{\eqn{C_P}}{\ifelse{html}{\out{<I>C<sub>P</sub></I>}}{Cp}}}

% other symbols
\newcommand{\degC}{\ifelse{latex}{\eqn{^{\circ}}C}{\ifelse{html}{\out{&deg;}C}{°C}}}
\newcommand{\le}{\ifelse{latex}{\eqn{\le}}{\ifelse{html}{\out{&le;}}{≤}}}
\newcommand{\ge}{\ifelse{latex}{\eqn{\ge}}{\ifelse{html}{\out{&ge;}}{≥}}}
\newcommand{\Psat}{\ifelse{latex}{\eqn{P_{\mathrm{SAT}}}}{\ifelse{html}{\out{<i>P</i><sub>SAT</sub>}}{Psat}}}
\newcommand{\Delta}{\ifelse{latex}{\eqn{\Delta}}{\ifelse{html}{\out{&Delta;}}{Δ}}}
\newcommand{\aacute}{\ifelse{latex}{\out{\'{a}}}{\ifelse{html}{\out{&aacute;}}{á}}}
\newcommand{\eacute}{\ifelse{latex}{\out{\'{e}}}{\ifelse{html}{\out{&eacute;}}{é}}}
