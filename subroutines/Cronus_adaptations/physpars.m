% This function was published within Cronus v2.1 by Marrero et al. 2016
% and addpted for WeCode by Richard Ott, 2021
%
%
% pp=physpars()
%
% Sets up basic physics parameters for all nuclides. Does not require
% inputs.
%
function pp=physpars()

% scaling_model = 'st';
global scaling_model


% decay constants
pp.lambda36Cl=2.3028e-6;
pp.lambda10Be=4.998e-7;
pp.lambda14C=1.21e-4;
pp.lambda3He=0.0;
pp.lambda21Ne=0.0;
pp.lambda26Al=9.83e-7;

%parameters for chlorine IDMS
pp.Wn=3.54527e+1;
pp.Wo=3.4963e+1;
pp.Sn=3.12e+0;
pp.So=2.93e+2;
pp.An=pp.Sn/(pp.Sn+1);
pp.Ao=pp.So/(pp.So+1);

%
% The table below corresponds to the table on the TABLE page in CHLOE.  The data is from Fabryka-Martin 1988, and was updated by Vasily Alfimov using the
% "Atlas of Neutron Resonances" by S.F. Mughabghab, Elsevier Science, 2006.The rows correspond to the following elements (top to bottom):
%  O, H, C, Na, Mg, Al, Si, P, K, Ca, Ti, Mn, Fe, Cl, B, Sm, Gd, U, Th, Cr, & Li.  Li & Cr are from the original Fabryka-Martin table and were not
%  updated.
% The columns correspond to the following categories: col 1= Ai=Atomic weight of element (g/mol),
%  col2 = average log decrement of energy per neutron collison with element i [-], col 3= neutron scattering cross-section of element i[i],
%  col 4= thermal neutron absorption cross-section of element i [-], col 5= dilute resonance integral for element i [-],
%  col 6= Mass stopping power of element i for alpha particles of a given energy [-],
%  col 7= neutron yield of element i per ppm U in radioequilibrium [-], col 8= Neutron yield of element i per ppm Th in radioequilibrium,
%  col 9 = Km, 602/atomic weight of sample, used to convert to at/g from
%  ppm. col 10= Stoichiometry ratio of oxide (oxygen to element i) [-], col 11=atomic number [-],
% col 12=Average capture probability relative to oxygen P(Z) (ref: von
% Egidy & Hartman 1982, value for C from back calculation from Heisinger
% 2002 compound factor calculations).
%
pp.table=[...
16          0.120004104		3.761		0.00019	0.0002693   539     0.23	0.079	0		0   8   1.00;...% O
1.01		1               20.49		0.3326	0           1542	0       0		0		0.5 1   0.00;...% H
12.01		0.157760474		4.74		0.0035	0.0018      573     0.45	0.18	13.691	2   6   0.36;...% C
22.99		0.084543589		3.038		0.517	0.311		454     14.5	6.8		19.42	0.5 11  1.00;...% Na
24.31		0.08009077		3.414		0.0666	0.038		463     5.8     2.6		14.94	1   12  0.93;...% Mg
26.98		0.072337427		1.413		0.231	0.17		449     5.1     2.6		11.812	1.5 13  0.76;...% Al
28.09		0.069559975		2.044		0.171	0.082		455     0.69	0.335	10.01	2   14  0.84;...% Si
30.97		0.063210393		3.134		0.165	0.079		444     0       0		8.48	2.5 15  1.04;...% P
39.1		0.050295528		2.04		2.1		1           432     0.45	0.305	12.78	0.5 19  1.54;...% K
40.08		0.049086179		2.93		0.43	0.233		436     0       0		10.73	1   20  1.90;...% Ca
47.87		0.041208508		4.09		6.41	3.1         367     0       0		7.53	2   22  2.66;...% Ti
54.94		0.035968204		2.06		13.36	13.4		351     0       0		8.486	1   25  2.73;...% Mn
55.85		0.035390922		11.35		2.56	1.36		353     0.19	0.205	7.54	1.5 26  3.28;...% Fe
35.45		0.055371497		15.8		33.14	13.83		420     0       0		16.98	0   17  1.32;...% Cl
10.81		0.174236264		4.27		767		343         537     62.3	19.2	55.68	0   5   0.25;...% B
150.36      0.013242694		38  		9640	1400		0       0       0		4.004	0   62  4.4;...% Sm
157.25      0.012664667		172     	41560	390         0       0       0		3.828	0   64  5.8;...% Gd
238.03      0.008378873		9.08		2.68	277         0       0       0		2.529	0   92  4.7;...% U
232.04      0.008594581		13.55		7.34	83.3		0       0       0		2.594	0   90  3.0;...% Th
52.0        0.038           3.38        3.1     1.6         0.0     0.0     0.0     11.578  0.0 24  2.98;...% Cr
6.9         0.264           0.95        70.5    0.0         548     21.1    9.6     86.731  0.0 3   0.18]; % Li

pp.ag=6.02e23;             % Avagadro's number.
pp.Sigmaetha=0.0548;
pp.Sigmatha=0.060241;
pp.Sigmasca=0.3773;
pp.Aa=14.5;
%
% Published values were used for these parameters.
%
pp.PsTi0=3.98002E-22;
%pp.PsTi0=3.98002E-22*3.5; Value from Fink (2000)
pp.PsFe0=1.72421e-22; %value from Stone (2005) abstract
pp.YstarCa=2.0e-24;
pp.YstarK=3.85e-24;
pp.Psimu0=175.0;

pp.Natoms3 = 2.006e22;
pp.Natoms10 = 2.006e22;
pp.Natoms14 = 2.006e22;
pp.Natoms26 = 1.003e22;

%k_negpartial must be multiplied
%by fstar(below)
pp.k_negpartial14 = 0.704 * 0.1828;
pp.delk_negpartial14 = 0.704 * 0.1828 * 0.0011;
pp.k_negpartial36K=0.802;
pp.k_negpartial36Ca=0.8486;
pp.k_negpartial10=(0.704 * 0.1828)./1.106;
pp.sigmak_neg10 = (0.704 * 0.1828 * 0.0003)./1.106;
pp.k_negpartial26=0.296 * 0.6559;
pp.sigmak_neg26 = 0.296 * 0.6559 * 0.002;

pp.sigma190_14 = 0.45e-27;
pp.delsigma190_14 = 0.25e-27;
pp.sigma190_10 = (0.094e-27)./1.106;
pp.sigmasigma190_10 = (0.013e-27)./1.106;
pp.sigma190_26 = 1.41e-27;
pp.sigmasigma190_26 = 0.17e-27;
pp.sigma190_36Ca=1.40e-27;
pp.sigmasigma190_36Ca=0.30e-27;
%
% More parameters.
%
pp.pEtha=0.56;
pp.La=3.9;
pp.Xa=0.1346;
pp.Phimuf0=700000;
pp.Ys=0.44;
pp.Yf=1.0;

%This value needed for muons and is from Nat's geomag scaling model values
pp.SPhiInf = 416.492; % Changed 12/13/11 to reflect updated SPhi values from Usoskin et al 2011

%
% Calibrated production rates.  See calibratexx.out.
%
%Production rates for all nuclides/all scaling schemes. Muon parameters calibrated for Al, Be, and
%Cl.
% References:
% CRONUScalc 2.0 refers to the original version of this code as
% published in 2016. The individual paper reference is given afterwards by
% abbreviation as indicated below. Any deviations from this are noted.
% CRONUScalc 2.1 refers to production rates changed/recalibrated after updates in 2020.

% Borchers: Borchers, B., Marrero, S., Balco, G., Caffee, M.,
% Goehring, B., Lifton, N., Nishiizumi, K., Phillips, F., Schaefer, J.,
% and Stone, J., 2016, Geological calibration of spallation production rates
% in the CRONUS-Earth project: Quaternary Geochronology, v. 31, p. 188–198,
% doi:10.1016/j.quageo.2015.01.009.

% Marrero: S. Marrero, F.M. Phillips, M. Caffee, J. Gosse. (2016b).
% CRONUS-Earth cosmogenic chlorine-36 calibration. Quaternary Geochronology,
% 31, 199-219, doi:10.1016/j.quageo.2015.10.002.

% Phillips: F. M. Phillips, D. C. Argento, G. Balco, M. W. Caffee, J. Clem,
% T. J. Dunai, R. Finkel, B. Goehring, J. C. Gosse, A. M. Hudson, A. J. T. Jull,
% M. A. Kelly, M. Kurz, D. Lal, N. Lifton, S. M. Marrero, K. Nishiizumi,
% R. C. Reedy, J. Schaefer, J. O. H. Stone, T. Swanson, M. G. Zreda. (2016a).
% The CRONUS-Earth Project: A Synthesis. Quaternary Geochronology, 31, 119-154,
% doi:10.1016/j.quageo.2015.09.006.

% Moore: A. K. Moore and D. E. Granger. (2019). Calibration of the
% production rate of cosmogenic 36Cl from Fe. Quaternary Geochronology, 51,
% 87-98, doi: 10.1016/j.quageo.2019.02.002.

switch scaling_model
case 'de'

	%DE scaling

	pp.PsAl=26.3;%CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=4.3;%CRONUScalc 2.0 - Phillips
	pp.PsBe=3.69;%CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.58;%CRONUScalc 2.0 - Phillips
    %pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
	pp.PsCa0=55.639824; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 55.596487)
	pp.sigmaPsCa0=12.4;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 13.4)
	pp.PsK0=127.652890;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 127.329329)
	pp.sigmaPsK0=28;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 31)
	pp.Pf0=671.187286;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 722.4699)
	pp.sigmapf0=250;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 181)
	pp.PsC=12.49;%CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=0.0;%CRONUScalc 2.0 - Phillips
	pp.PsHe=123;%CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 17; %CRONUScalc 2.0 - Phillips
    pp.PsNe=15.31; % CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.78; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.84e-3; %CRONUScalc 2.0 - Phillips
	pp.fstar26=10.9e-3; %CRONUScalc 2.0 - Phillips
	pp.fstar36K=0.05291409;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.05232411)
	pp.fstar36Ca=0.0134791;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.01367434)
	pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.sigma036K=9.966357e-30;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 10.102183e-30)
	pp.sigma036Ca=8.260059e-30;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.279102e-30)
	pp.sigma010=0.251e-30; %CRONUScalc 2.0 - Phillips
	pp.sigma026=4.21e-30; %CRONUScalc 2.0 - Phillips

%pp.StAl = 27; % for Lal/Stone age approximation in Attn Length Calc
%pp.StBe = 4; % for Lal/Stone age approximation in Attn Length Calc
%pp.StClCa = 54; % for Lal/Stone age approximation in Attn Length Calc
%pp.StC = 12.4; % for Lal/Stone age approximation in Attn Length Calc
%pp.StC = 12.4; % for Lal/Stone age approximation in Attn Length Calc
%pp.StHe = 120; % for Lal/Stone age approximation in Attn Length Calc

case 'du'
% DU scaling

	pp.PsAl=27.3; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=4.4;% CRONUScalc 2.0 – Phillips
	pp.PsBe=3.70; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.58; % CRONUScalc 2.0 – Phillips
    %pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
	pp.PsCa0=54.981815; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 54.983145)
	pp.sigmaPsCa0=12.2; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 12.9)
	pp.PsK0=128.290537; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 128.006482)
	pp.sigmaPsK0=28; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 30)
	pp.Pf0=669.535792; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 719.6554)
	pp.sigmapf0=252; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 182)
	pp.PsC=12.4; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=1.2;% CRONUScalc 2.0 – Phillips
	pp.PsHe=123; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 14; % CRONUScalc 2.0 – Phillips
    pp.PsNe=15.35; % CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.78; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.85e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar26=10.9e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar36K=0.05418697; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.05359793)
	pp.fstar36Ca=0.01376199; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.0139712)
	pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.sigma036K=9.799104e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 9.934741e-30)
	pp.sigma036Ca=8.308403e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.330115e-30)
	pp.sigma010=0.251e-30; % CRONUScalc 2.0 – Phillips
	pp.sigma026=4.21e-30; % CRONUScalc 2.0 – Phillips

case 'li'
% LI scaling
	pp.PsAl=28.7; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=4.6; % CRONUScalc 2.0 – Phillips
	pp.PsBe=4.06; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.57;% CRONUScalc 2.0 – Phillips
    %pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
	pp.PsCa0=60.133328; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 60.140146)
	pp.sigmaPsCa0=11.9; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.2)
	pp.PsK0=142.085283; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 141.677893)
	pp.sigmaPsK0=28; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 19)
	pp.Pf0=723.732503; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 778.6495)
	pp.sigmapf0=263; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 183)
	pp.PsC=13.4; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=0.4;% CRONUScalc 2.0 – Phillips
	pp.PsHe=131; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 19; % CRONUScalc 2.0 – Phillips
    pp.PsNe=16.86; % CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.82; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.86e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar26=11.1e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar36K=0.05349911; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.05288851)
	pp.fstar36Ca=0.01337862; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.01357641)
	pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.sigma036K=9.879147e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 10.017604e-30)
	pp.sigma036Ca=8.242888e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.262502e-30)
	pp.sigma010=0.251e-30; % CRONUScalc 2.0 – Phillips
	pp.sigma026=4.19e-30; % CRONUScalc 2.0 – Phillips


case 'lm'
% LM scaling
	pp.PsAl=27.9; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=2.7; % CRONUScalc 2.0 – Phillips
	pp.PsBe=4.00; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.32; % CRONUScalc 2.0 – Phillips
    %pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
	pp.PsCa0=51.686259;   % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 51.664818)
	pp.sigmaPsCa0=3.3;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 4.9)
	pp.PsK0=150.004199;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 150.715117)
	pp.sigmaPsK0=10;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 14)
	pp.Pf0=645.291055;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 691.0873)
    pp.sigmapf0=232;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 186)
	pp.PsC=12.2; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=0.0; % CRONUScalc 2.0 – Phillips
	pp.PsHe=117; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 13; % CRONUScalc 2.0 – Phillips
    pp.PsNe=16.59; % CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.01; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.87e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar26=11.0e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar36K=0.06211791;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.06146167)
	pp.fstar36Ca=0.01396003;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.01416637)
	pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.sigma036K=8.719879e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.861757e-30)
	pp.sigma036Ca=8.340333e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.36167e-30)
	pp.sigma010=0.251e-30; % CRONUScalc 2.0 – Phillips
	pp.sigma026=4.22e-30; % CRONUScalc 2.0 – Phillips

case 'st'
% ST scaling
	pp.PsAl=27.9; % CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=2.8;% CRONUScalc 2.0 – Phillips
	pp.PsBe=4.01;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.33;% CRONUScalc 2.0 – Phillips
    %pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
	pp.PsCa0=52.166027;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 52.160727)
	pp.sigmaPsCa0=3.7; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 5.2)
	pp.PsK0=150.075197; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 149.790172)
	pp.sigmaPsK0=11; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 15)
	pp.Pf0=649.405271; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 695.5976)
	pp.sigmapf0=230; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 185)
	pp.PsC=12.2;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=0.0; % CRONUScalc 2.0 – Phillips
	pp.PsHe=118;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 18; % CRONUScalc 2.0 – Phillips
    pp.PsNe=16.63; % CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.04; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.87e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar26=10.9e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar36K=0.05897637; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.05830165)
	pp.fstar36Ca=0.0133582; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.0135710)
	pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.sigma036K=9.103510e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 9.253363e-30)
	pp.sigma036Ca=8.2331016e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.254806e-30)
	pp.sigma010=0.251e-30; % CRONUScalc 2.0 – Phillips
	pp.sigma026=4.22e-30; % CRONUScalc 2.0 – Phillips


case 'sf'
% SF scaling
	pp.PsAl=28.6;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=3.3;% CRONUScalc 2.0 – Phillips
	pp.PsBe=4.09;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.35;% CRONUScalc 2.0 – Phillips
    %pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
	pp.PsCa0=55.853169;  % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 56.24897)
	pp.sigmaPsCa0=4.5; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 4.6)
	pp.PsK0=153.493296; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 153.296853)
	pp.sigmaPsK0=12.4; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 12)
	pp.Pf0=695.858691; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 742.7535)
	pp.sigmapf0=196; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 179)
	pp.PsC=12.7;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=0.0;% CRONUScalc 2.0 – Phillips
	pp.PsHe=119;% CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 19; % CRONUScalc 2.0 – Phillips
    pp.PsNe=16.96;% CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.12; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.89e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar26=11.7e-3; % CRONUScalc 2.0 – Phillips
	pp.fstar36K=0.0579802; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.05736245)
	pp.fstar36Ca=0.01331762; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.01358239)
	pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.sigma036K=9.262416e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 9.398043e-30)
	pp.sigma036Ca=8.228633e-30; % CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.261573e-30)
	pp.sigma010=0.252e-30; % CRONUScalc 2.0 – Phillips
	pp.sigma026=4.10e-30; % CRONUScalc 2.0 – Phillips


case 'sa'
% SA scaling
	pp.PsAl=28.5; %CRONUScalc 2.0 - Borchers
	pp.sigmaPsAl=3.1; %CRONUScalc 2.0 - Phillips
	pp.PsBe=3.92; %CRONUScalc 2.0 - Borchers
	pp.sigmaPsBe=0.31; %CRONUScalc 2.0 - Phillips
	%pp.PsFe0=??; %CRONUScalc 2.1 - Moore
    %pp.sigmaPsFe0=??; %CRONUScalc 2.1 - Moore
    pp.PsCa0=55.560068; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 56.0)
    pp.sigmaPsCa0=4.2; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 4.1)
    pp.PsK0=155.654793;%CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 155)
    pp.sigmaPsK0=11.8; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 11)
    pp.Pf0=714.313247; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 759)
    pp.sigmapf0=191; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 180)
	pp.PsC=12.8; %CRONUScalc 2.0 - Borchers
	pp.sigmaPsC=0.0; %CRONUScalc 2.0 - Phillips
	pp.PsHe=115; %CRONUScalc 2.0 - Borchers
	pp.sigmaPsHe = 19; %CRONUScalc 2.0 - Phillips
	pp.PsNe=16.26; % CRONUScalc 2.0 - From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =1.96; % CRONUScalc 2.0 - Value from Balco's paper, error propagated.
	pp.fstar14=0.137; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth
	pp.fstar10=1.89e-3; %CRONUScalc 2.0 - Phillips
	pp.fstar26=12.1e-3; %CRONUScalc 2.0 - Phillips
    pp.fstar36Ca=0.01328490; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.0135)
    pp.fstar36K=0.05833281; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 0.0577)
	pp.sigma036Ca=8.220769e-30; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 8.25e-30)
    pp.sigma036K=9.217989e-30; %CRONUScalc 2.1; CRONUScalc 2.0 - Marrero (Prev Value: 9.36e-30)
    pp.sigma010=0.252e-30; %CRONUScalc 2.0 - Phillips
	pp.sigma026=4.03e-30; %CRONUScalc 2.0 - Phillips
    pp.sigma014=8.79321E-30; % CRONUScalc 2.0 - Not calibrated as part of CRONUS-Earth

end
