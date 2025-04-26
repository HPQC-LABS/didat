clear('all')
%% Fundamental Constants

bohrRadii2angstroms=0.52917721092; % 5.2917721092(17)�10?11 from https://en.wikipedia.org/wiki/Bohr_radius
hartrees2wavenumbers=219474.6313705; % 219?474.631?3705(15) cm?1 from https://en.wikipedia.org/wiki/Hartree
eV2wavenumbers=8065.54468111324; % from http://www.highpressurescience.com/onlinetools/conversion.html
%%      << THEORETICAL LONG-RANGE POTENTIALS FOR b-STATE >>
%% 1.   Define the r array (independent variable) for theoretical long-range potential, and Set r_e, D_e, Delta E (all for b-state)
r=1:0.01:500;

%deltaE=0.33534652830;                                      % (7,7) Li2 from Sansonetti 1995 PRA 52 pg 2682
deltaE=0.096;                                               % From Sansonetti 1995

re=2.589769594404D+00;
re=2.589769594404D+00;

reA=3.10792898; %re for A-state

De=1.218218000000D+04;                                     % (7,7) Li2 from Le Roy 2009 JCP 131 ar 204309 (same as in 2011 Le Roy 109 pg 435 since this state was fixed in that analysis)
De=1.218218000000D+04;                                     % (6,6) Li2 from Dattani 2011 JMS 268 pg 199 , Supplementary Material Table III.A

V_lim=0;
%%     << DAMPING >>
firstIonizationPotential_hydrogen=13.5984;            % eV , for ^1 H, ground state. from http://physics.nist.gov/cuu/Archive/2002RMP.pdf 2005 Mohr, Rev. Mod. Phys. 77
firstIonizationPotential_Li_2S=43487.15940/8065.54445; % cm-1/8065.54445 = eV , experimental for 7Li, from http://pra.aps.org/pdf/PRA/v75/i5/e052503 2007 Bushaw, PRA, 75, article: 05203, Isotope shift between 7Li and 6Li is also given in abstract, but is small.
firstIonizationPotential_Li_3P=firstIonizationPotential_Li_2S-30925.517/5.765593821614D+03 ; % eV , subtracted D1(cog)/eV2wavenumber for 7Li from firstIonizationPotential_Li_S according to Bob LeRoy's suggestion in his e-mail to Nike on 7 September 2013 18:35

rhoA=(firstIonizationPotential_Li_2S/firstIonizationPotential_hydrogen)^(2/3);
rhoB=(firstIonizationPotential_Li_3P/firstIonizationPotential_hydrogen)^(2/3);
rhoAB=2*rhoA*rhoB/(rhoA+rhoB);
rhoAB=0.5; % should be  0.031366494940143 for 3b state, but this is just flat!
%% 2.   Set C3,C6,C8,C9,C10,C11

% C3

%C3=3.5777*10^5;                                                         % 7 Li with relativistic correction Tang 2010 PRA 81 pg 042521
C3=dispersionCoefficientsFromAtomicUnits2spectroscopic(3.3314e-2,0,3,1); % 2007 Zhang
C3p_triplet_u=C3/2;                                                      % by symmetry 
C3_Sigma_triplet_u=-C3;
C3_Pi_singlet_u=-C3p_triplet_u;

% C6

%C6=1.0005417*10^7;                                                                  % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
C6=dispersionCoefficientsFromAtomicUnits2spectroscopic(3.8236e4,0,6,1);              % 2007 Zhang

%C6p_triplet_u=6.7815899*10^6;                                                       % C6p_singlet_g=0;  % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
C6p_triplet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(2.0282e4,0,6,1);   % C6p_singlet_g=0;  % 2007 Zhang

C6_Sigma_triplet_u=C6;
C6_Pi_singlet_u=C6p_triplet_u;

% C8

%C8=3.6995320*10^8;                                                                  % C8 sigma % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
C8=dispersionCoefficientsFromAtomicUnits2spectroscopic(2.4870e7,0,8,1);              % C8 sigma % 2007 Zhang

%C8p_triplet_u=6.5543451*10^7;                                                       % C8p_triplet_u=0;  % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
C8p_triplet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(3.7222e5,0,8,1);   % C8p_triplet_u=0;  % 2007 Zhang

C8_Sigma_triplet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(2.3183e7,0,8,1);
C8_Pi_singlet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(7.8976e5,0,8,1);

% C9

% C9=1.633994329773709e+08;              % C9 sigma          % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
% C9p_triplet_u=3.694341201013722e+07;   % C9p_triplet_u=0;  % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
% 
% C9_Sigma_triplet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(2.2880e5,0,9,1);
% C9_Pi_singlet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(5.173e4,0,9,1);
% 
% % C10
% 
% C10=1.137377*10^(10);                  % C10=0;            % infinite mass Li with no relativistic correction Zhang 2007 PRA 75 pg 042509
% C10p_triplet_u=3.4707458*10^6;         % C10p_triplet_u=0; % infinite mass Li with no relativistic correction Zhang 2007 PRA 75 pg 042509
% 
% C10_Sigma_triplet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(3.0096e7,0,10,1);
% C10_Pi_singlet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(9.1839e5,0,10,1);
% 
% 
% % C11
% 
% C11=-1.185909325855541e+10;            % C11 sigma         % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
% C11p_triplet_u=5.303594489323599e+09;  % C11p_triplet_u=0; % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
% 
% C11_Sigma_triplet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(-5.930e7,0,11,1);
% C11_Pi_singlet_u=dispersionCoefficientsFromAtomicUnits2spectroscopic(2.652e7,0,11,1);

%% 2_u state , no diagonalization
Cm_power=[           3          6             8             ];
Cm_Sigma_singlet_u=[C3         C6            C8             ];
Cm_Sigma_triplet_u=[C3_Sigma_triplet_u C6_Sigma_triplet_u C8_Sigma_triplet_u];
Cm_Pi_singlet_u=[C3_Pi_singlet_u C6_Pi_singlet_u C8_Pi_singlet_u];
Cm_Pi_triplet_u=[C3p_triplet_u C6p_triplet_u C8p_triplet_u];

V=V_lim-zeros(length(Cm_power),length(r));
V_undamped=V;

for ii=1:length(Cm_power)
    for jj=1:ii;
        Dm=(1-exp(-((3.30*rhoAB*r+0.423*(rhoAB*r).^2)/Cm_power(jj)))).^(Cm_power(jj)-1); % has dimension of r
        V(ii,:)=V(ii,:)-(Cm_Pi_triplet_u(jj).*Dm./r.^Cm_power(jj));
        V_undamped(ii,:)=V_undamped(ii,:)-(Cm_Pi_triplet_u(jj)./r.^Cm_power(jj));
    end
end

%V=V+deltaE;

%% 1_u state , 3x3 diagonalization

M=zeros(3,3,length(r));
M(2,2,:)=deltaE;
M(3,3,:)=deltaE;

for ii=1:length(Cm_power)
Dm=(1-exp(-((3.30*rhoAB*r+0.423*(rhoAB*r).^2)/Cm_power(ii)))).^(Cm_power(ii)-1); % has dimension of r
M(1,1,:)=squeeze(M(1,1,:)).' -(1/3)*          ( (    Cm_Sigma_triplet_u(ii)*Dm + Cm_Pi_singlet_u(ii)*Dm + Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
M(1,2,:)=squeeze(M(1,2,:)).' -(1/(3*sqrt(2)))*( ( -2*Cm_Sigma_triplet_u(ii)*Dm + Cm_Pi_singlet_u(ii)*Dm + Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
M(1,3,:)=squeeze(M(1,3,:)).' -(1/sqrt(6))*    ( (                              - Cm_Pi_singlet_u(ii)*Dm + Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;

M(2,1,:)=M(1,2,:);
M(2,2,:)=squeeze(M(2,2,:)).' -(1/6)*          ( (  4*Cm_Sigma_triplet_u(ii)*Dm + Cm_Pi_singlet_u(ii)*Dm + Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
M(2,3,:)=squeeze(M(2,3,:)).' -(1/(2*sqrt(3)))*( (                              - Cm_Pi_singlet_u(ii)*Dm + Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;

M(3,1,:)=M(1,3,:);
M(3,2,:)=M(2,3,:);
M(3,3,:)=squeeze(M(3,3,:)).' -(1/2)*          ((                                 Cm_Pi_singlet_u(ii)*Dm + Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
end

for ii=1:length(r)
    eigenvalues=eig(M(:,:,ii));
    eig_1(ii)=(eigenvalues(1));  eig_2(ii)=eigenvalues(2); eig_3(ii)=eigenvalues(3); 
    V_1_b(ii)=min([eig_1(ii) eig_2(ii)]);
    V_1_2a(ii)=max([eig_1(ii) eig_2(ii)]);
end

Mtrue3x3_b_B_2a=[eig_1' eig_2' eig_3'];                                               

plot(r,Mtrue3x3_b_B_2a)
%% 0^+_u state , 2x2 diagonalization
M=zeros(2,2,length(r));
M(2,2,:)=deltaE;

for ii=1:length(Cm_power)
Dm=(1-exp(-((3.30*rhoAB*r+0.423*(rhoAB*r).^2)/Cm_power(ii)))).^(Cm_power(ii)-1); % has dimension of r
M(1,1,:)=squeeze(M(1,1,:)).' -(1/3)*      ( (  Cm_Sigma_singlet_u(ii)*Dm + 2*Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
M(1,2,:)=squeeze(M(1,2,:)).' +(sqrt(2)/3)*( (  Cm_Sigma_singlet_u(ii)*Dm -   Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;

M(2,1,:)=M(1,2,:);
M(2,2,:)=squeeze(M(2,2,:)).' -(1/3)*      ( ( 2*Cm_Sigma_singlet_u(ii)*Dm +  Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
end

for ii=1:length(r)
    eigenvalues=eig(M(:,:,ii));
    eig_1(ii)=(eigenvalues(1));  eig_2(ii)=eigenvalues(2); 
    V_0plus_A(ii)=min([eig_1(ii) eig_2(ii)]);
    V_0plus_b(ii)=max([eig_1(ii) eig_2(ii)]);
end

Mtrue2x2_A_b=[eig_1' eig_2'];                                               %  max(max(Mtrue-MtrueNoRetardation))=0.110466420650482. plot(r,abs(Mtrue(:,1)-MtrueNoRetardation)) peaks at like 2.5 Angstroms and 0.11 cm-1, then monotonically decreases with r. Shouldn't it nicrease with r ? I changed e-4 to e-1 and now the general trend is still a decrease, but it's oscillatory.

%% 0^-_u state , 2x2 diagonalization

M=zeros(2,2,length(r));
M(2,2,:)=deltaE;

for ii=1:length(Cm_power)
Dm=(1-exp(-((3.30*rhoAB*r+0.423*(rhoAB*r).^2)/Cm_power(ii)))).^(Cm_power(ii)-1); % has dimension of r
M(1,1,:)=squeeze(M(1,1,:)).' -(1/3)*      ( (  Cm_Sigma_triplet_u(ii)*Dm + 2*Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
M(1,2,:)=squeeze(M(1,2,:)).' +(sqrt(2)/3)*( (  Cm_Sigma_triplet_u(ii)*Dm -   Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;

M(2,1,:)=M(1,2,:);
M(2,2,:)=squeeze(M(2,2,:)).' -(1/3)*      ( ( 2*Cm_Sigma_triplet_u(ii)*Dm +  Cm_Pi_triplet_u(ii)*Dm)./r.^Cm_power(ii) ) ;
end

for ii=1:length(r)
    eigenvalues=eig(M(:,:,ii));
    eig_1(ii)=(eigenvalues(1));  eig_2(ii)=eigenvalues(2); 
    V_0minus_b(ii)=min([eig_1(ii) eig_2(ii)]);
    V_0minus_2a(ii)=max([eig_1(ii) eig_2(ii)]);
end

Mtrue2x2_b_2a=[eig_1' eig_2'];                                               %  max(max(Mtrue-MtrueNoRetardation))=0.110466420650482. plot(r,abs(Mtrue(:,1)-MtrueNoRetardation)) peaks at like 2.5 Angstroms and 0.11 cm-1, then monotonically decreases with r. Shouldn't it nicrease with r ? I changed e-4 to e-1 and now the general trend is still a decrease, but it's oscillatory.



%%      << ab initio GLOBAL POTENTIALS >>
%% 12.  Set r array and V array for (6,6) Li2 potential from Schmidt-Mink 1985 Chemical Physics. 
% N16round__66Li2bstatePotential_from2009JCP=...
% [3.25 50899; 
% 3.50  38091; 
% 4.00  21577; 
% 4.50  14192; 
% 5.00  12782; 
% 5.50  15031; 
% 6.00  19374; 
% 6.50  24786; 
% 7.00  30597; 
% 7.50  36361; 
% 8.00  41784; 
% 9.00  50902; 
% 10.5  59529; 
% 12.0  63554; 
% 13.5  65352; 
% 15.0  66244; 
% 17.0  66877; 
% 20.0  67356; 
% 30.0  67859;
% 
% r_66Li2bstate_from2009JCP=N16round__66Li2bstatePotential_from2009JCP(:,1);
% V_66Li2bstate_from2009JCP=N16round__66Li2bstatePotential_from2009JCP(:,2);
%% 12.  Set r array and V array for (6,6) Li2 potential from Schmidt-Mink 1985 Chemical Physics. 
% p6q3=load('N08p6q3Rf3_20_uAx.10.txt');
% r_p6q3=p6q3(:,1);
% V_p6q3=p6q3(:,2);

abInitio_3C=load('s3_1pig.txt');
r_abInitio_3C=abInitio_3C(:,1)*bohrRadii2angstroms;
V_abInitio_3C=abInitio_3C(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

C_08p6q8r07_00=load('3C_08p6q8r07_00.10.txt');
r_3C_08p6q8r07_00=C_08p6q8r07_00(:,1);
V_3C_08p6q8r07_00=C_08p6q8r07_00(:,2);

abInitio_3B=load('s3_1piu.txt');
r_abInitio_3B=abInitio_3B(:,1)*bohrRadii2angstroms;
V_abInitio_3B=abInitio_3B(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

abInitio_3A=load('s3_1siu+.txt');
r_abInitio_3A=abInitio_3A(:,1)*bohrRadii2angstroms;
V_abInitio_3A=abInitio_3A(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

abInitio_3d=load('s3_3pig.txt');
r_abInitio_3d=abInitio_3d(:,1)*bohrRadii2angstroms;
V_abInitio_3d=abInitio_3d(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

C_3d_04p6q9r07_20=load('3d_04p6q9r07_20.10.txt');
r_3d_04p6q9r07_20=C_3d_04p6q9r07_20(:,1);
V_3d_04p6q9r07_20=C_3d_04p6q9r07_20(:,2);

abInitio_3b=load('s3_3piu.txt');
r_abInitio_3b=abInitio_3b(:,1)*bohrRadii2angstroms;
V_abInitio_3b=abInitio_3b(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

b_p6q7r05_90=load('3b_p6q7r05_90.10.txt');
r_3b_p6q7r05_90=b_p6q7r05_90(:,1);
V_3b_p6q7r05_90=b_p6q7r05_90(:,2);

B_p6q6r06_20=load('3B_p6q6r06_20.10.txt');
r_3B_p6q6r06_20=B_p6q6r06_20(:,1);
V_3B_p6q6r06_20=B_p6q6r06_20(:,2);

abInitio_3c=load('s3_3sig+.txt');
r_abInitio_3c=abInitio_3c(:,1)*bohrRadii2angstroms;
V_abInitio_3c=abInitio_3c(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

abInitio_6X=load('s6_1sig+.txt');
r_abInitio_6X=abInitio_6X(:,1)*bohrRadii2angstroms;
V_abInitio_6X=abInitio_6X(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC

abInitio_6a=load('s6_3siu+.txt');
r_abInitio_6a=abInitio_6a(:,1)*bohrRadii2angstroms;
V_abInitio_6a=abInitio_6a(:,2)*hartrees2wavenumbers-(-14.94710546*hartrees2wavenumbers+3.83389013*eV2wavenumbers); % 14903.2967364 = 6Li(D1cog), 14.94710546 a.u. = S + S dissociation limit in Table 2 of 2014 Musial JCTC
%% 16.  Vibratioanl energies from LeRoy 2009 JCP from scienide2:ndattani/diat/Li2/A-X/2013-A-X/final_files_from_2009-JCP-vol131-ar204309/66first/16round.7
v0=-11958.58920;
v(1)=-11589.29680;

%%      << FIG 1 >>
%% 15.  Plot (6,6) Li2 potential from Semczuk 2013 PRA
close('all')
figure1=figure(1);

set(gca,'Position',[0.101190476190476 0.105359342915811 0.885 0.888357289527721])
set(gcf,'Color','w')

minX=1;maxX=20;
splineMeshSize_for_r=0.00001; % I just chose it to be of the same precision as r_e defined in cell 1. Splines can only be necessary for short-range, since long-range seems to have points very close together.
splineMesh_for_r=minX:splineMeshSize_for_r:maxX;

V_ab_initio_spline_6X=spline(r_abInitio_6X,V_abInitio_6X,splineMesh_for_r);
V_ab_initio_spline_6a=spline(r_abInitio_6a,V_abInitio_6a,splineMesh_for_r);
V_ab_initio_spline_3A=spline(r_abInitio_3A,V_abInitio_3A,splineMesh_for_r);
V_ab_initio_spline_3B=spline(r_abInitio_3B,V_abInitio_3B,splineMesh_for_r);
V_ab_initio_spline_3C=spline(r_abInitio_3C,V_abInitio_3C,splineMesh_for_r);
V_ab_initio_spline_3b=spline(r_abInitio_3b,V_abInitio_3b,splineMesh_for_r);
V_ab_initio_spline_3c=spline(r_abInitio_3c,V_abInitio_3c,splineMesh_for_r);
V_ab_initio_spline_3d=spline(r_abInitio_3d,V_abInitio_3d,splineMesh_for_r);

minY=-6000;maxY=300;
axis([minX,maxX,minY,maxY])
line([minX maxX], [0 0],'Color','k','LineWidth',3)

hold('on')

% [Vmin_ab_initio, indexOfVmin_ab_initio]=min(V_ab_initio_spline);
% plotPointsSymmetricallyOnPotential(vibrationalEnergies(1:2:end),V_ab_initio_spline,splineMesh_for_r,indexOfVmin_ab_initio);

axis([minX,maxX,minY,maxY])
plot(r_abInitio_6X,V_abInitio_6X,'Color',[0,255,100]./255,'Linewidth',3);
plot(r_abInitio_3A,V_abInitio_3A,'Color','r','Linewidth',3);
plot(r_abInitio_3B,V_abInitio_3B,'Color','k','Linewidth',3);
plot(r_abInitio_3C,V_abInitio_3C,'Color','c','Linewidth',3);

figure2=figure(2);

set(gca,'Position',[0.101190476190476 0.105359342915811 0.885 0.888357289527721])
set(gcf,'Color','w')

minX=1;maxX=20;
minY=-10000;maxY=300;
axis([minX,maxX,minY,maxY])
line([minX maxX], [0 0],'Color','k','LineWidth',3)

splineMeshSize_for_r=0.00001; % I just chose it to be of the same precision as r_e defined in cell 1. Splines can only be necessary for short-range, since long-range seems to have points very close together.
splineMesh_for_r=minX:splineMeshSize_for_r:maxX;

hold('on')

plot(r_abInitio_6a,V_abInitio_6a,'Color','b','Linewidth',3);
plot(r_abInitio_3b,V_abInitio_3b,'Color','m','Linewidth',3);
plot(r_abInitio_3c,V_abInitio_3c,'Color','g','Linewidth',3);
plot(r_abInitio_3d,V_abInitio_3d,'Color',[100,100,100]./255,'Linewidth',3);

% [Vmin_ab_initio, indexOfVmin_ab_initio]=min(V_ab_initio_spline);
% plotPointsSymmetricallyOnPotential(vibrationalEnergies(1:2:end),V_ab_initio_spline,splineMesh_for_r,indexOfVmin_ab_initio);

axis([minX,maxX,minY,maxY])
plot(r_abInitio_6a,V_abInitio_6a,'Color','b','Linewidth',3);
plot(r_abInitio_3b,V_abInitio_3b,'Color','m','Linewidth',3);
plot(r_abInitio_3c,V_abInitio_3c,'Color','g','Linewidth',3);
plot(r_abInitio_3d,V_abInitio_3d,'Color',[100,100,100]./255,'Linewidth',3);

%%
close('all');figure3=figure(3);hold('on')
set(gca,'Position',[0.101190476190476 0.114161849710983 0.892953008436171 0.879554782732549])
set(gcf,'Color','w')

minX=2;maxX=15;
minY=-6000;maxY=300;
axis([minX,maxX,minY,maxY])
line([minX maxX], [0 0],'Color','k','LineWidth',3)
re=3.98;

plot(r_abInitio_3b,V_abInitio_3b,'Color',[0,255,100]./255,'Linewidth',10);
plot(r_3b_p6q7r05_90,V_3b_p6q7r05_90,'k','LineWidth',3)

v0=-5654.932449854602055; 
v( 1)=-5493.765135128664951;
v( 2)=-5385.710984514686970;
v( 3)=-5241.269953247034209;
v( 4)=-5083.976206177168933;
v( 5)=-4914.85;
v( 6)=-4737.691247106215087;
v( 7)=-4554.829936403927604;
v( 8)=-4368.240502184650722;
v( 9)=-4179.292096277415112;
v(10)=-3988.756231213321371;
v(11)=-3796.972235927646125;
v(12)=-3604.218750080568952;
v(13)=-3411.033149015913750;
v(14)=-3218.297919388399805;
v(15)=-3027.108781656958399;
v(16)=-2838.477232599816034;
v(17)=-2652.930601841349471;
v(18)=-2470.302596380747673;
v(19)=-2290.105676375571875;
v(20)=-2112.256988503596288;
v(21)=-1937.405398146945536;
v(22)=-1766.631273238948864;
v(23)=-1600.818435795562496;
v(24)=-1440.124051819757312;
v(25)=-1284.340696145073408;
v(26)=-1134.583144744511360;
v(27)= -996.778888028432000;
v(28)= -914.777113202062848;
v(29)= -882.576799450836608;
v(30)= -845.026425674076544;
v(31)= -800.327833982061696;
v(32)= -753.924510206392448;
v(33)= -705.900692148498688;
v(34)= -654.935744351130112;
v(35)= -602.077023609842304;
v(36)= -548.234951970257472;
v(37)= -493.706925401166272;
v(38)= -438.926074728071424;
v(39)= -384.474432243322304;
v(40)= -330.903787772794240;
v(41)= -278.814572462308608;
v(42)= -228.913912424278848;
v(43)= -182.018345069401152;
v(44)= -139.060176467093648;
v(45)= -101.066168262798576;
v(46)=  -69.047731355521066;
v(47)=  -43.749645705548336;
v(48)=  -25.299916111188852;
v(49)=  -13.022995406988458;
v(50)=   -5.669747141689355;
v(51)=   -1.855271158639904;
v(52)=   -0.324750667836605;
v(53)=   -0.002393731491646;
vibrationalEnergies=[v0 v];

[Vmin_ab_initio_3b, indexOfVmin_ab_initio_3b]=min(V_ab_initio_spline_3b);
%plotPointsSymmetricallyOnPotential(vibrationalEnergies(1:2:end),V_ab_initio_spline_3b,splineMesh_for_r,indexOfVmin_ab_initio_3b,200);

vibrationalLevelsMeasured=[0:53];
for vibrationalLevel=vibrationalLevelsMeasured+1
line([interp1(V_3b_p6q7r05_90(r_3b_p6q7r05_90<re),r_3b_p6q7r05_90(r_3b_p6q7r05_90<re),vibrationalEnergies(vibrationalLevel)) interp1(V_3b_p6q7r05_90(r_3b_p6q7r05_90>re),r_3b_p6q7r05_90(r_3b_p6q7r05_90>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','b','LineWidth',3)
r_innerTurningPoint(vibrationalLevel)=interp1(V_3b_p6q7r05_90(r_3b_p6q7r05_90<re),r_3b_p6q7r05_90(r_3b_p6q7r05_90<re),vibrationalEnergies(vibrationalLevel));
scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
r_outerTurningPoint(vibrationalLevel)=interp1(V_3b_p6q7r05_90(r_3b_p6q7r05_90>re),r_3b_p6q7r05_90(r_3b_p6q7r05_90>re),vibrationalEnergies(vibrationalLevel));
scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
end

% vibrationalLevelsMeasured=[0:53];
% for vibrationalLevel=vibrationalLevelsMeasured+1
% line([interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel)) interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','r','LineWidth',3)
% r_innerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel));
% scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% r_outerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel));
% scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% end

plot(r(700:end),V(end,700:end),'r','LineWidth',6)

annotation(figure3,'textbox',[0.620 0.760 0.327 0.087],'Color','r','String','$V(r) = - \frac{D_3(r)C_3}{r^3} - \frac{D_6(r)C_6}{r^6} - \frac{D_8C_8(r)}{r^8}$','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure3,'textbox',[0.421 0.609 0.107 0.0636],'String','MLR$_{5.9}^{6,7}(17)$','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure3,'textbox',[0.293 0.223 0.107 0.0636],'String','Li$_2\left(3b,3^3\Pi_u\right)$','LineStyle','none','Interpreter','latex','FontSize',36,'FontName','Helvetica','FitBoxToText','off');
annotation(figure3,'textbox',[0.438 0.532 0.107 0.0636],'Color',[0,255,100]./255,'String','Original','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure3,'line',[0.366764275256223 0.412884333821376],[0.63775968992248 0.637209302325581],'LineWidth',3);
scatter(5.89,-2750,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
scatter(6.5,-2750,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
annotation(figure3,'arrow',[0.620790629575403 0.588579795021962],[0.850162790697674 0.924031007751938],'LineWidth',3);
annotation(figure3,'arrow',[0.922401171303075 0.953147877013177],[0.768992248062015 0.689922480620155],'LineWidth',3);

xlabel('Internuclear distance \AA','Interpreter','Latex','FontSize',36)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$V(r)$ cm$^{-1}$','Interpreter','Latex','FontSize',36)
box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',20);

% % %%%% INSET %%%%% % % 

insetAxesHandle=axes('Position',[0.63177159590044 0.264178120123497 0.343779050175466 0.46140327522534]);hold('on')

minX=0;maxX=0.001;
%minY=C3p_triplet_u-100;maxY=C3p_triplet_u+1000;
minY=0;maxY=100000;
axis([minX,maxX,minY,maxY])

plot(1./r.^3,-r.^3.*V(end,:),'r','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
plot(1./r.^3,-r.^3.*V_undamped(end,:),'c','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
plot(1./r_3b_p6q7r05_90.^3,-r_3b_p6q7r05_90.^3.*V_3b_p6q7r05_90,'k','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
%plot(1./splineMesh_for_r.^3,-splineMesh_for_r.^3.*V_ab_initio_spline_3b,'b','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

for vibrationalLevel=0:53 % the spline goes bad after v=56, r=11.16 Angstroms. The ab initio points are at r=9.5, 10.1, 10.6, 11.2, then 28
%line([1./r_outerTurningPoint(vibrationalLevel+1).^3 1./r_outerTurningPoint(vibrationalLevel+1).^3 ], [minY maxY])
r_innerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_innerTurningPoint(vibrationalLevel+1).^3,-r_innerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
r_outerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_outerTurningPoint(vibrationalLevel+1).^3,-r_outerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
end

%plot(1./r_abInitio_3b.^3,-r_abInitio_3b.^3.*V_abInitio_3b,'Color',[0,255,100]./255,'LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

xlabel('$1/r^3$ \AA$^{-3}$','Interpreter','Latex','FontSize',26)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$-r^3V(r)$ \AA$^{3}$cm$^{-1}$','Interpreter','Latex','FontSize',26)

box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);
%% B-state (1u)
close('all');figure4=figure(4);hold('on')
set(gca,'Position',[0.101190476190476 0.114161849710983 0.892953008436171 0.879554782732549])
set(gcf,'Color','w')

minX=2;maxX=15;
minY=-6000;maxY=300;
axis([minX,maxX,minY,maxY])
line([minX maxX], [0 0],'Color','k','LineWidth',3)
re=3.165;

plot(r_abInitio_3B,V_abInitio_3B,'Color',[0,255,100]./255,'Linewidth',10);
plot(r_3B_p6q6r06_20,V_3B_p6q6r06_20,'k','LineWidth',3)

v0=-5259.702267161964301;  
v( 1)=-5047.140849881975555; 
v( 2)=-4838.569786159541763; 
v( 3)=-4634.049062932902416; 
v( 4)=-4433.598339738490722; 
v( 5)=-4237.211653070495231; 
v( 6)=-4044.873993594876993; 
v( 7)=-3856.578296560157924; 
v( 8)=-3672.341252509541391; 
v( 9)=-3492.216542011438378; 
v(10)=-3316.304488645995661; 
v(11)=-3144.757464404133316; 
v(12)=-2977.780421310530528; 
v(13)=-2815.625477755599604; 
v(14)=-2658.578476551557742; 
v(15)=-2506.933989605654915; 
v(16)=-2360.954118872159142; 
v(17)=-2220.807831851416722; 
v(18)=-2086.495694733651712; 
v(19)=-1957.782947061481216; 
v(20)=-1834.183256277384192; 
v(21)=-1715.028175716195840; 
v(22)=-1599.606072525047552; 
v(23)=-1487.301883932875520; 
v(24)=-1377.677567197110272; 
v(25)=-1270.486719610525440; 
v(26)=-1165.651954085971712; 
v(27)=-1063.232673372752512; 
v(28)= -963.396914960888704; 
v(29)= -866.400745208404096; 
v(30)= -772.574265319793280; 
v(31)= -682.311966587884160; 
v(32)= -596.064857021040256; 
v(33)= -514.331528103759552; 
v(34)= -437.644846655112448; 
v(35)= -366.550363916987072; 
v(36)= -301.572455022926592; 
v(37)= -243.165988378659008; 
v(38)= -191.656817726308320; 
v(39)= -147.183982888330304; 
v(40)= -109.665218902316736; 
v(41)=  -78.804168171930848; 
v(42)=  -54.136957611555160; 
v(43)=  -35.092203999862157; 
v(44)=  -21.035204396617212; 
v(45)=  -11.285153823763400; 
v(46)=   -5.113465586640047; 
v(47)=   -1.737975355206322; 
v(48)=   -0.324712343973808;
v(49)=   -0.007461527144271;
vibrationalEnergies=[v0 v];

[Vmin_ab_initio_3B, indexOfVmin_ab_initio_3B]=min(V_ab_initio_spline_3B);
%plotPointsSymmetricallyOnPotential(vibrationalEnergies(1:2:end),V_ab_initio_spline_3b,splineMesh_for_r,indexOfVmin_ab_initio_3b,200);

vibrationalLevelsMeasured=[0:49];
for vibrationalLevel=vibrationalLevelsMeasured+1
line([interp1(V_3B_p6q6r06_20(r_3B_p6q6r06_20<re),r_3B_p6q6r06_20(r_3B_p6q6r06_20<re),vibrationalEnergies(vibrationalLevel)) interp1(V_3B_p6q6r06_20(r_3B_p6q6r06_20>re),r_3B_p6q6r06_20(r_3B_p6q6r06_20>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','b','LineWidth',3)
r_innerTurningPoint(vibrationalLevel)=interp1(V_3B_p6q6r06_20(r_3B_p6q6r06_20<re),r_3B_p6q6r06_20(r_3B_p6q6r06_20<re),vibrationalEnergies(vibrationalLevel));
scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
r_outerTurningPoint(vibrationalLevel)=interp1(V_3B_p6q6r06_20(r_3B_p6q6r06_20>re),r_3B_p6q6r06_20(r_3B_p6q6r06_20>re),vibrationalEnergies(vibrationalLevel));
scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
end

% vibrationalLevelsMeasured=[0:53];
% for vibrationalLevel=vibrationalLevelsMeasured+1
% line([interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel)) interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','r','LineWidth',3)
% r_innerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel));
% scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% r_outerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel));
% scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% end

plot(r(700:end),V(end,700:end),'r','LineWidth',6)

annotation(figure4,'textbox',[0.620 0.760 0.327 0.087],'Color','r','String','$V(r) = - \frac{D_3(r)C_3}{r^3} - \frac{D_6(r)C_6}{r^6} - \frac{D_8C_8(r)}{r^8}$','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure4,'textbox',[0.421 0.609 0.107 0.0636],'String','MLR$_{6.2}^{6,6}(4)$','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure4,'textbox',[0.293 0.223 0.107 0.0636],'String','Li$_2\left(3B,3^1\Pi_u\right)$','LineStyle','none','Interpreter','latex','FontSize',36,'FontName','Helvetica','FitBoxToText','off');
annotation(figure4,'textbox',[0.438 0.532 0.107 0.0636],'Color',[0,255,100]./255,'String','Original','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure4,'line',[0.366764275256223 0.412884333821376],[0.63775968992248 0.637209302325581],'LineWidth',3);
scatter(5.89,-2750,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
scatter(6.5,-2750,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
annotation(figure4,'arrow',[0.620790629575403 0.588579795021962],[0.850162790697674 0.924031007751938],'LineWidth',3);
annotation(figure4,'arrow',[0.922401171303075 0.953147877013177],[0.768992248062015 0.689922480620155],'LineWidth',3);

xlabel('Internuclear distance \AA','Interpreter','Latex','FontSize',36)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$V(r)$ cm$^{-1}$','Interpreter','Latex','FontSize',36)
box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',20);

% % %%%% INSET %%%%% % % 

insetAxesHandle=axes('Position',[0.63177159590044 0.264178120123497 0.343779050175466 0.46140327522534]);hold('on')

minX=0;maxX=0.001;
%minY=C3p_triplet_u-100;maxY=C3p_triplet_u+1000;
minY=0;maxY=100000;
axis([minX,maxX,minY,maxY])

plot(1./r.^3,-r.^3.*V(end,:),'r','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
%plot(1./r.^3,-r.^3.*V_undamped(end,:),'c','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
plot(1./r_3B_p6q6r06_20.^3,-r_3B_p6q6r06_20.^3.*V_3B_p6q6r06_20,'k','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
%plot(1./splineMesh_for_r.^3,-splineMesh_for_r.^3.*V_ab_initio_spline_3b,'b','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

for vibrationalLevel=0:49 % the spline goes bad after v=56, r=11.16 Angstroms. The ab initio points are at r=9.5, 10.1, 10.6, 11.2, then 28
%line([1./r_outerTurningPoint(vibrationalLevel+1).^3 1./r_outerTurningPoint(vibrationalLevel+1).^3 ], [minY maxY])
r_innerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3B(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_innerTurningPoint(vibrationalLevel+1).^3,-r_innerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
r_outerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3B(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_outerTurningPoint(vibrationalLevel+1).^3,-r_outerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
end

%plot(1./r_abInitio_3b.^3,-r_abInitio_3b.^3.*V_abInitio_3b,'Color',[0,255,100]./255,'LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

xlabel('$1/r^3$ \AA$^{-3}$','Interpreter','Latex','FontSize',26)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$-r^3V(r)$ \AA$^{3}$cm$^{-1}$','Interpreter','Latex','FontSize',26)

box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',16);

%% C-state (1g)
close('all');figure5=figure(5);hold('on')
set(gca,'Position',[0.101190476190476 0.114161849710983 0.892953008436171 0.879554782732549])
set(gcf,'Color','w')

minX=2;maxX=15;
minY=-4200;maxY=300;
axis([minX,maxX,minY,maxY])
line([minX maxX], [0 0],'Color','k','LineWidth',3)
re=3.137;

plot(r_abInitio_3C,V_abInitio_3C,'Color',[0,255,100]./255,'Linewidth',10);
plot(r_3C_08p6q8r07_00,V_3C_08p6q8r07_00,'k','LineWidth',3)

v0=-3957.890301282239761;
v( 1)=-3745.119157025943878;
v( 2)=-3524.020966461889202;
v( 3)=-3324.020966461889202;
v( 4)=-3116.570131136164491;
v( 5)=-2911.754405554478126;
v( 6)=-2710.006065926872907;
v( 7)=-2511.771813604795170;
v( 8)=-2317.519675049029047;
v( 9)=-2127.744091638082304;
v(10)=-1942.968271627786496;
v(11)=-1763.742033056433920;
v(12)=-1590.632282885245952;
v(13)=-1424.202382217587712;
v(14)=-1264.977470546079744;
v(15)=-1113.398784719593088;
v(16)= -969.785127262634752;
v(17)= -834.340424003393024;
v(18)= -707.248873597836544;
v(19)= -588.846535677036288;
v(20)= -479.738057606813440;
v(21)= -380.571382655633792;
v(22)= -291.279706608250048;
v(23)= -210.843993944130240;
v(24)= -139.563348452164896;
v(25)=  -75.962574346390944;
v(26)=  -59.183547111123872;
v(27)=  -44.917842694088616;
v(28)=  -32.408391003076368;
v(29)=  -13.147887047208048;
vibrationalEnergies=[v0 v];

[Vmin_ab_initio_3C, indexOfVmin_ab_initio_3C]=min(V_ab_initio_spline_3C);
%plotPointsSymmetricallyOnPotential(vibrationalEnergies(1:2:end),V_ab_initio_spline_3b,splineMesh_for_r,indexOfVmin_ab_initio_3b,200);

vibrationalLevelsMeasured=[0:30];
for vibrationalLevel=vibrationalLevelsMeasured+1
line([interp1(V_3C_08p6q8r07_00(r_3C_08p6q8r07_00<re),r_3C_08p6q8r07_00(r_3C_08p6q8r07_00<re),vibrationalEnergies(vibrationalLevel)) interp1(V_3C_08p6q8r07_00(r_3C_08p6q8r07_00>re),r_3C_08p6q8r07_00(r_3C_08p6q8r07_00>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','b','LineWidth',3)
r_innerTurningPoint(vibrationalLevel)=interp1(V_3C_08p6q8r07_00(r_3C_08p6q8r07_00<re),r_3C_08p6q8r07_00(r_3C_08p6q8r07_00<re),vibrationalEnergies(vibrationalLevel));
scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
r_outerTurningPoint(vibrationalLevel)=interp1(V_3C_08p6q8r07_00(r_3C_08p6q8r07_00>re),r_3C_08p6q8r07_00(r_3C_08p6q8r07_00>re),vibrationalEnergies(vibrationalLevel));
scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
end

% vibrationalLevelsMeasured=[0:53];
% for vibrationalLevel=vibrationalLevelsMeasured+1
% line([interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel)) interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','r','LineWidth',3)
% r_innerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel));
% scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% r_outerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel));
% scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% end

plot(r(700:end),V(end,700:end),'r','LineWidth',6)

annotation(figure4,'textbox',[0.822781844802344 0.843815028901732 0.112796486090776 0.070924855491331],'Color',[1 0 0],'String','$V(r) = - u(r)$','LineStyle','none','Interpreter','latex','FontSize',24,'FitBoxToText','off');
annotation(figure4,'textbox',[0.393162518301611 0.606109826589595 0.107 0.0636],'String','MLR$_{7.0}^{6,8}(8)$','LineStyle','none','Interpreter','latex','FontSize',24,'FitBoxToText','off');
annotation(figure4,'textbox',[0.293 0.223 0.107 0.0636],'String','Li$_2\left(3C,3^1\Pi_g\right)$','LineStyle','none','Interpreter','latex','FontSize',36,'FontName','Helvetica','FitBoxToText','off');
annotation(figure4,'textbox',[0.4 0.532 0.107 0.0636],'Color',[0,255,100]./255,'String','Original','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure4,'line',[0.326500732064422 0.385797950219619],[0.638728323699422 0.63728323699422],'LineWidth',3);
scatter(5.39,-1900,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
scatter(6,-1900,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
annotation(figure4,'arrow',[0.94289897510981 0.943631039531479],[0.855491329479769 0.812138728323699],'LineWidth',3);
annotation(figure4,'arrow',[0.823572474377745 0.569546120058565],[0.893063583815029 0.893063583815029],'LineWidth',3);

xlabel('Internuclear distance \AA','Interpreter','Latex','FontSize',36)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$V(r)$ cm$^{-1}$','Interpreter','Latex','FontSize',36)
box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',20);

% % %%%% INSET - LE ROY SPACE %%%%% % % 

insetAxesHandle=axes('Position',[0.624450951683748 0.58643245538361 0.34 0.234376793171301]);hold('on')
 
minX=0;maxX=0.001;
%minY=C3p_triplet_u-100;maxY=C3p_triplet_u+1000;
minY=0;maxY=100000;
axis([minX,maxX,minY,maxY])

plot(1./r.^3,-r.^3.*V(end,:),'r','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
%plot(1./r.^3,-r.^3.*V_undamped(end,:),'c','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
plot(1./r_3C_08p6q8r07_00.^3,-r_3C_08p6q8r07_00.^3.*V_3C_08p6q8r07_00,'k','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
%plot(1./splineMesh_for_r.^3,-splineMesh_for_r.^3.*V_ab_initio_spline_3b,'b','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

for vibrationalLevel=0:53 % the spline goes bad after v=56, r=11.16 Angstroms. The ab initio points are at r=9.5, 10.1, 10.6, 11.2, then 28
%line([1./r_outerTurningPoint(vibrationalLevel+1).^3 1./r_outerTurningPoint(vibrationalLevel+1).^3 ], [minY maxY])
r_innerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3C(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_innerTurningPoint(vibrationalLevel+1).^3,-r_innerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
r_outerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3C(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_outerTurningPoint(vibrationalLevel+1).^3,-r_outerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
end

%plot(1./r_abInitio_3b.^3,-r_abInitio_3b.^3.*V_abInitio_3b,'Color',[0,255,100]./255,'LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

xlabel('$1/r^3$ \AA$^{-3}$','Interpreter','Latex','FontSize',20)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$-r^3V(r)$ \AA$^{3}$cm$^{-1}$','Interpreter','Latex','FontSize',20)

box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',14);

% % %%%% INSET %%%%% % % 

insetAxesHandle=axes('Position',[0.624450951683748 0.243946906250664 0.34 0.234376793171301]);hold('on')

minX=7;maxX=10;
minY=-100;maxY=-65;
axis([minX,maxX,minY,maxY])

%   plot(r_abInitio_3C,V_abInitio_3C,'Color',[0,255,100]./255,'LineWidth',10,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
%plot(splineMesh_for_r,V_ab_initio_spline_3C,'b','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
plot(r_3C_08p6q8r07_00,V_3C_08p6q8r07_00,'k','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
scatter(r_abInitio_3C,V_abInitio_3C,150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k') % contains C3, C6, C8, C9, C10, C11

for vibrationalLevel=25:25 % the spline goes bad after v=56, r=11.16 Angstroms. The ab initio points are at r=9.5, 10.1, 10.6, 11.2, then 28
text(r_outerTurningPoint(vibrationalLevel+1)-0.5,vibrationalEnergies(vibrationalLevel+1)+4,strcat('$v=',num2str(vibrationalLevel),'$'),'Interpreter','latex','Color','b','FontSize',20)
line([r_innerTurningPoint(vibrationalLevel+1) r_outerTurningPoint(vibrationalLevel+1) ], [vibrationalEnergies(vibrationalLevel+1) vibrationalEnergies(vibrationalLevel+1)],'Color','b','LineWidth',3)
r_innerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3C(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel+1));
scatter(r_innerTurningPoint(vibrationalLevel+1),vibrationalEnergies(vibrationalLevel+1),'MarkerFaceColor','b','MarkerEdgeColor','b');
r_outerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3C(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel+1));
scatter(r_outerTurningPoint(vibrationalLevel+1),vibrationalEnergies(vibrationalLevel+1),'MarkerFaceColor','b','MarkerEdgeColor','b');
end

xlabel('$r$ \AA','FontSize',20,'Interpreter','latex');
ylabel('$V(r)$ cm$^{-1}$','FontSize',20,'Interpreter','latex');

box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',14);

%% 3d state (2g)

close('all');figure6=figure(6);hold('on')
set(gca,'Position',[0.101190476190476 0.114161849710983 0.892953008436171 0.879554782732549])
set(gcf,'Color','w')

minX=2;maxX=15;
minY=-6200;maxY=300;
axis([minX,maxX,minY,maxY])
line([minX maxX], [0 0],'Color','k','LineWidth',3)
re=3.18;

plot(r_abInitio_3d,V_abInitio_3d,'Color',[0,255,100]./255,'Linewidth',10);
plot(r_3d_04p6q9r07_20,V_3d_04p6q9r07_20,'k','LineWidth',3)

v0=-5909.082870709067720; 
v( 1)=-5641.796305197752190; 
v( 2)=-5377.524117379149175; 
v( 3)=-5116.452670555453551; 
v( 4)=-4858.751050544754435; 
v( 5)=-4604.574956177269087; 
v( 6)=-4354.070515379085919; 
v( 7)=-4107.378161105932122; 
v( 8)=-3864.636683065670240; 
v( 9)=-3625.987562434165739; 
v(10)=-3391.579698179508796; 
v(11)=-3161.574647280548106; 
v(12)=-2936.152530982379176; 
v(13)=-2715.518811966155681; 
v(14)=-2499.912233637356167; 
v(15)=-2289.614350938369626; 
v(16)=-2084.961304408457984; 
v(17)=-1886.358855954255872; 
v(18)=-1694.302335629324288; 
v(19)=-1509.404298133581312; 
v(20)=-1332.434943932483072; 
v(21)=-1164.385230868400384; 
v(22)=-1006.574505456192768; 
v(23)= -860.858922231548800; 
v(24)= -730.125234005433344; 
v(25)= -620.025713485010048; 
v(26)= -556.5838; %from regular lev (new_stolyarov/analytic_no_stolyarov.9)
v(27)= -536.279739461353984; 
v(28)= -502.412950961499776; 
v(29)= -464.979122336951168; 
v(30)= -425.955105861283072; 
v(31)= -385.425289908528128; 
v(32)= -344.186189032638336; 
v(33)= -302.862201947904064; 
v(34)= -262.005598583399328; 
v(35)= -222.176157230817248; 
v(36)= -183.954413800329920; 
v(37)= -114.815825508788192; 
v(38)= -114.815825508788192; 
v(39)= -94.815825508788192 ;
v(40)=  -69.005984215226432; 
v(41)=  -39.005984215884736; 
v(42)=  -23.171796637275224; 
v(43)=  -12.113609603112186; 
v(44)=   -5.236141959367891; 
v(45)=   -1.618761491885656; 
vibrationalEnergies=[v0 v];

[Vmin_ab_initio_3d, indexOfVmin_ab_initio_3d]=min(V_ab_initio_spline_3d);
%plotPointsSymmetricallyOnPotential(vibrationalEnergies(1:2:end),V_ab_initio_spline_3b,splineMesh_for_r,indexOfVmin_ab_initio_3b,200);

vibrationalLevelsMeasured=[0:44];
for vibrationalLevel=vibrationalLevelsMeasured+1
line([interp1(V_3d_04p6q9r07_20(r_3d_04p6q9r07_20<re),r_3d_04p6q9r07_20(r_3d_04p6q9r07_20<re),vibrationalEnergies(vibrationalLevel)) interp1(V_3d_04p6q9r07_20(r_3d_04p6q9r07_20>re),r_3d_04p6q9r07_20(r_3d_04p6q9r07_20>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','b','LineWidth',3)
r_innerTurningPoint(vibrationalLevel)=interp1(V_3d_04p6q9r07_20(r_3d_04p6q9r07_20<re),r_3d_04p6q9r07_20(r_3d_04p6q9r07_20<re),vibrationalEnergies(vibrationalLevel));
scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
r_outerTurningPoint(vibrationalLevel)=interp1(V_3d_04p6q9r07_20(r_3d_04p6q9r07_20>re),r_3d_04p6q9r07_20(r_3d_04p6q9r07_20>re),vibrationalEnergies(vibrationalLevel));
scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','b','MarkerEdgeColor','b');
end

% vibrationalLevelsMeasured=[0:53];
% for vibrationalLevel=vibrationalLevelsMeasured+1
% line([interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel)) interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel))], [vibrationalEnergies(vibrationalLevel) vibrationalEnergies(vibrationalLevel)],'Color','r','LineWidth',3)
% r_innerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel));
% scatter(r_innerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% r_outerTurningPoint(vibrationalLevel)=interp1(V_ab_initio_spline_3b(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel));
% scatter(r_outerTurningPoint(vibrationalLevel),vibrationalEnergies(vibrationalLevel),'MarkerFaceColor','r','MarkerEdgeColor','r');
% end

plot(r(500:end),V(end,500:end),'r','LineWidth',6)

annotation(figure5,'textbox',[0.822781844802344 0.843815028901732 0.112796486090776 0.070924855491331],'Color',[1 0 0],'String','$V(r) = - u(r)$','LineStyle','none','Interpreter','latex','FontSize',24,'FitBoxToText','off');
annotation(figure5,'textbox',[0.393162518301611 0.606109826589595 0.107 0.0636],'String','MLR$_{7.2}^{6,9}(4)$','LineStyle','none','Interpreter','latex','FontSize',24,'FitBoxToText','off');
annotation(figure5,'textbox',[0.293 0.223 0.107 0.0636],'String','Li$_2\left(3d,3^3\Pi_g\right)$','LineStyle','none','Interpreter','latex','FontSize',36,'FontName','Helvetica','FitBoxToText','off');
annotation(figure5,'textbox',[0.4 0.532 0.107 0.0636],'Color',[0,255,100]./255,'String','Original','LineStyle','none','Interpreter','latex','FontSize',24);
annotation(figure5,'line',[0.326500732064422 0.385797950219619],[0.638728323699422 0.63728323699422],'LineWidth',3);
scatter(5.395,-2900,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
scatter(6.005,-2900,200,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k')
annotation(figure5,'arrow',[0.94289897510981 0.943631039531479],[0.855491329479769 0.812138728323699],'LineWidth',3);
annotation(figure5,'arrow',[0.823572474377745 0.569546120058565],[0.893063583815029 0.893063583815029],'LineWidth',3);

xlabel('Internuclear distance \AA','Interpreter','Latex','FontSize',36)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$V(r)$ cm$^{-1}$','Interpreter','Latex','FontSize',36)
box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',20);

% % %%%% INSET - LE ROY SPACE %%%%% % % 

insetAxesHandle=axes('Position',[0.624450951683748 0.58643245538361 0.34 0.234376793171301]);hold('on')
 
minX=0;maxX=0.001;
%minY=C3p_triplet_u-100;maxY=C3p_triplet_u+1000;
minY=0;maxY=100000;
axis([minX,maxX,minY,maxY])

plot(1./r.^3,-r.^3.*V(end,:),'r','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
%plot(1./r.^3,-r.^3.*V_undamped(end,:),'c','LineWidth',5,'LineStyle','-') % contains C3, C6, C8
plot(1./r_3d_04p6q9r07_20.^3,-r_3d_04p6q9r07_20.^3.*V_3d_04p6q9r07_20,'k','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
%plot(1./splineMesh_for_r.^3,-splineMesh_for_r.^3.*V_ab_initio_spline_3b,'b','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

for vibrationalLevel=0:44 % the spline goes bad after v=56, r=11.16 Angstroms. The ab initio points are at r=9.5, 10.1, 10.6, 11.2, then 28
%line([1./r_outerTurningPoint(vibrationalLevel+1).^3 1./r_outerTurningPoint(vibrationalLevel+1).^3 ], [minY maxY])
r_innerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3d(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_innerTurningPoint(vibrationalLevel+1).^3,-r_innerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
r_outerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3d(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel+1));
scatter(1./r_outerTurningPoint(vibrationalLevel+1).^3,-r_outerTurningPoint(vibrationalLevel+1).^3.*vibrationalEnergies(vibrationalLevel+1),150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k');
end

%plot(1./r_abInitio_3b.^3,-r_abInitio_3b.^3.*V_abInitio_3b,'Color',[0,255,100]./255,'LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11

xlabel('$1/r^3$ \AA$^{-3}$','Interpreter','Latex','FontSize',20)                                                                                    % (r) taken out as per request by Kirk Madison in email on 25/8/2012
ylabel('$-r^3V(r)$ \AA$^{3}$cm$^{-1}$','Interpreter','Latex','FontSize',20)

box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',14);

% % %%%% INSET %%%%% % % 

insetAxesHandle=axes('Position',[0.624450951683748 0.243946906250664 0.34 0.234376793171301]);hold('on')

minX=6;maxX=8.2;
minY=-630;maxY=-510;
axis([minX,maxX,minY,maxY])

%   plot(r_abInitio_3C,V_abInitio_3C,'Color',[0,255,100]./255,'LineWidth',10,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
%plot(splineMesh_for_r,V_ab_initio_spline_3C,'b','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
plot(r_3d_04p6q9r07_20,V_3d_04p6q9r07_20,'k','LineWidth',3,'LineStyle','-') % contains C3, C6, C8, C9, C10, C11
scatter(r_abInitio_3d,V_abInitio_3d,150,'MarkerFaceColor',[0,255,100]./255,'MarkerEdgeColor','k') % contains C3, C6, C8, C9, C10, C11

% FOR SOME REASON THE BELOW GIVES INDEX EXCEEDS ARRAY DIMENSIONS IF I CHANGE TO MLR WHILE USING SAME RE
for vibrationalLevel=25:27 % the spline goes bad after v=56, r=11.16 Angstroms. The ab initio points are at r=9.5, 10.1, 10.6, 11.2, then 28
text(r_outerTurningPoint(vibrationalLevel+1)-0.255,vibrationalEnergies(vibrationalLevel+1)+6,strcat('$v=',num2str(vibrationalLevel),'$'),'Interpreter','latex','Color','b','FontSize',18)
line([r_innerTurningPoint(vibrationalLevel+1) r_outerTurningPoint(vibrationalLevel+1) ], [vibrationalEnergies(vibrationalLevel+1) vibrationalEnergies(vibrationalLevel+1)],'Color','b','LineWidth',3)
r_innerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3d(splineMesh_for_r<re),splineMesh_for_r(splineMesh_for_r<re),vibrationalEnergies(vibrationalLevel+1));
scatter(r_innerTurningPoint(vibrationalLevel+1),vibrationalEnergies(vibrationalLevel+1),'MarkerFaceColor','b','MarkerEdgeColor','b');
r_outerTurningPoint(vibrationalLevel+1)=interp1(V_ab_initio_spline_3d(splineMesh_for_r>re),splineMesh_for_r(splineMesh_for_r>re),vibrationalEnergies(vibrationalLevel+1));
scatter(r_outerTurningPoint(vibrationalLevel+1),vibrationalEnergies(vibrationalLevel+1),'MarkerFaceColor','b','MarkerEdgeColor','b');
line([r_innerTurningPoint(vibrationalLevel+1) r_outerTurningPoint(vibrationalLevel+1) ], [vibrationalEnergies(vibrationalLevel+1) vibrationalEnergies(vibrationalLevel+1)],'Color','b','LineWidth',3)
end

xlabel('$r$ \AA','FontSize',20,'Interpreter','latex');
ylabel('$V(r)$ cm$^{-1}$','FontSize',20,'Interpreter','latex');

box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',2,'FontSize',14);

function [convertedValue,convertedUncertainty]=dispersionCoefficientsFromAtomicUnits2spectroscopic(originalValue,originalUncertainty,order,atomicUnits2SpectroscopicUnitsORspectroscopicUnits2atomicUnits)

% atomicUnits2SpectroscopicUnitsORspectroscopicUnits2atomicUnits = 1 to get spectroscopic units, and -1 to get atomic units

BohrRadius=5.2917721092e-1; % Angstronms , from  2010 CODATA , http://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0

RinfinityTIMESTwo=219474.6313605; % this is just what i had in my excel spreadsheet .. couldn't find it in these units online, should be recalculated some time

RinfinityTIMESTwo=219474.6313705; % wavenumbers , http://en.wikipedia.org/wiki/Hartree says this number came from 2010 CODATA, but the CODATA recommended value seems to be in Joules, and in case we don't trust wikipedia's conversion:

RinfinityTIMESTwo=219474.631333725; % 4.35974434(19)e-18 Joules /((6.62606957e-34 Js)*(29979245800 cm/s)) , 4.35974434(19)e-18 = raw value from 2010 CODATA http://physics.nist.gov/cgi-bin/cuu/Value?hr

conversionFactor=RinfinityTIMESTwo*(BohrRadius^order); 

convertedValue=originalValue*conversionFactor^atomicUnits2SpectroscopicUnitsORspectroscopicUnits2atomicUnits;     

convertedUncertainty=originalUncertainty*conversionFactor^atomicUnits2SpectroscopicUnitsORspectroscopicUnits2atomicUnits;     

end
