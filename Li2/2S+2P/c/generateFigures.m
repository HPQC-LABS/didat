clear('all')
%% Fundamental Constants

bohrRadii2angstroms=0.52917721092; % 5.2917721092(17)×10?11 from https://en.wikipedia.org/wiki/Bohr_radius
hartrees2wavenumbers=219474.6313705; % 219?474.631?3705(15) cm?1 from https://en.wikipedia.org/wiki/Hartree
eV2wavenumbers=8065.54468111324; % from http://www.highpressurescience.com/onlinetools/conversion.html

%%      << THEORETICAL LONG-RANGE POTENTIALS FOR c-STATE >>
%% 1.   Define the r array (independent variable) for theoretical long-range potential, and Set r_e, D_e, Delta E (all for b-state)
r=1:0.01:500;

deltaE=0.3353246;                                                                           % (7,7) Li2 from Sansonetti 1995 PRA 52 pg 2682
deltaE=0.3353246;                                                                           % (6,6) Li2 from Sansonetti 1995 PRA 52 pg 2682

re=3.065436D+00;
re=3.065436D+00;

De=7.0934926000000D+03;                                                                     % (7,7) Li2 from Le Roy 2009 JCP 131 ar 204309 (same as in 2011 Le Roy 109 pg 435 since this state was fixed in that analysis)
De=7.0934926000000D+03;                                                                     % (6,6) Li2 from Dattani 2011 JMS 268 pg 199 , Supplementary Material Table III.A

V_lim=0;
%%     << DAMPING >>
firstIonizationPotential_hydrogen=13.5984;                                                  % eV , for ^1 H, ground state. from http://physics.nist.gov/cuu/Archive/2002RMP.pdf 2005 Mohr, Rev. Mod. Phys. 77
firstIonizationPotential_Li_S=43487.15940/8065.54445;                                       % cm-1/8065.54445 = eV , experimental for 7Li, from http://pra.aps.org/pdf/PRA/v75/i5/e052503 2007 Bushaw, PRA, 75, article: 05203, Isotope shift between 7Li and 6Li is also given in abstract, but is small.
firstIonizationPotential_Li_P=firstIonizationPotential_Li_S-14903.2967364042/8065.54445;    % eV , subtracted D1(cog)/eV2wavenumber for 7Li from firstIonizationPotential_Li_S according to Bob LeRoy's suggestion in his e-mail to Nike on 7 September 2013 18:35

rhoA=(firstIonizationPotential_Li_S/firstIonizationPotential_hydrogen)^(2/3);
rhoB=(firstIonizationPotential_Li_P/firstIonizationPotential_hydrogen)^(2/3);
rhoAB=2*rhoA*rhoB/(rhoA+rhoB);
rhoAB=0.5; % rounded from 0.4647 since there are approximations in the above: S-state ionization energy is for 7Li, P-state ionization energy is crudely approximated (I believe)

%% 2.   Set C3,C6,C8,C9,C10,C11

% C3

%C3=dispersionCoefficientsFromAtomicUnits2spectroscopic(11.0007,0,3,1);                      % 7 Li with relativistic correction Tang 2010 PRA 81 pg 042521
C3=3.576828e5;                       % 6 Li with relativistic correction Tang (from Jim Mitroy's e-mail to Nike on  4 June 2013 03:31)
C3_Sigma_triplet_g=C3;
C3_Pi_singlet_g=-C3/2;
C3_Pi_triplet_g=C3/2;

% C6

C6=1.00059e7;                       % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
C6=1.00059e7;                       % 6 Li with no relativistic correction Tang 2009 PRA 79 pg 06712

C6_Pi_triplet_g=6.78183e6;            % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
%C6_Pi_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(1407.20,0,6,1);           % 6 Li with no relativistic correction Tang 2009 PRA 79 pg 06712. If we want the true value that will get back the atomic units value when converted back, we need to keep all digits. Rounding is just for reporting, so that the reader doesn't think there's more significant digits (if using significance arithmetic, ie if the uncertainty/error is not known) or so that the presentation looks neater when the uncertainty/error is in fact known (no use agonizing about digit 7 if we're statistically concerned about digit 6

C6_Sigma_triplet_g=C6;
C6_Pi_singlet_g=C6_Pi_triplet_g;

% C8

%C8_Sigma_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(274128,0,8,1);                        % C8 sigma % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
C8_Sigma_triplet_g=3.69965e8;                        % C8 sigma % 6 Li with no relativistic correction Tang 2009 PRA 79 pg 06712

C8_Pi_triplet_g=6.55441e8;             % C8p_triplet_u=0;  % 7 Li with no relativistic correction Tang 2009 PRA 79 pg 06712
%C8_Pi_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(103053,0,8,1);            % C8p_triplet_u=0;  % 6 Li with no relativistic correction Tang 2009 PRA 79 pg 06712

C8_Pi_singlet_g=1.39076e8;          % 7Li
C8_Pi_singlet_g=1.39076e8;          % 6Li

% C9

C9=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,9,1);                      % C9 sigma          % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
C9_Pi_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,9,1);            % C9p_triplet_u=0;  % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
C9_Sigma_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,9,1);
C9_Pi_singlet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,9,1);

% C10

C10=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,10,1);                    % C10=0;            % infinite mass Li with no relativistic correction Zhang 2007 PRA 75 pg 042509
C10_Pi_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,10,1);         % C10p_triplet_u=0; % infinite mass Li with no relativistic correction Zhang 2007 PRA 75 pg 042509

C10_Sigma_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,10,1);
C10_Pi_singlet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,10,1);

% C11

C11=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,11,1);                    % C11 sigma         % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502
C11_Pi_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,11,1);          % C11p_triplet_u=0; % infinite mass Li with no relativistic correction Tang 2011 PRA 84 pg 052502

C11_Sigma_triplet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,11,1);
C11_Pi_singlet_g=dispersionCoefficientsFromAtomicUnits2spectroscopic(0,0,11,1);

%% 2_u state , no diagonalization
Cm_power=[           3          6             8             9             10             11           ];
%Cm_Sigma_singlet_g=[C3         C6            C8            C9            C10            C11           ];
Cm_Sigma_triplet_g=[C3_Sigma_triplet_g C6_Sigma_triplet_g+C3_Sigma_triplet_g^2/(4*De) C8_Sigma_triplet_g C9_Sigma_triplet_g C10_Sigma_triplet_g C11_Sigma_triplet_g];
Cm_Pi_singlet_g=[C3_Pi_singlet_g C6_Pi_singlet_g C8_Pi_singlet_g C9_Pi_singlet_g C10_Pi_singlet_g C11_Pi_singlet_g];
Cm_Pi_triplet_g=[C3_Pi_triplet_g C6_Pi_triplet_g C8_Pi_triplet_g C9_Pi_triplet_g C10_Pi_triplet_g C11_Pi_triplet_g];

V=V_lim-zeros(length(Cm_power)+1,length(r));

V=V+deltaE;

%% 1_u state , 3x3 diagonalization

M=zeros(3,3,length(r));
M(2,2,:)=deltaE;
M(3,3,:)=deltaE;

for ii=1:length(Cm_power)
Dm=(1-exp(-((3.30*rhoAB*r+0.423*(rhoAB*r).^2)/Cm_power(ii)))).^(Cm_power(ii)-1); % has dimension of r
Dm=Dm./Dm; % make equal to 1
M(1,1,:)=squeeze(M(1,1,:)).' -(1/3)*          ( (    Cm_Sigma_triplet_g(ii)*Dm + Cm_Pi_singlet_g(ii)*Dm + Cm_Pi_triplet_g(ii)*Dm)./r.^Cm_power(ii) ) ;
M(1,2,:)=squeeze(M(1,2,:)).' -(1/(3*sqrt(2)))*( ( -2*Cm_Sigma_triplet_g(ii)*Dm + Cm_Pi_singlet_g(ii)*Dm + Cm_Pi_triplet_g(ii)*Dm)./r.^Cm_power(ii) ) ;
M(1,3,:)=squeeze(M(1,3,:)).' -(1/sqrt(6))*    ( (                              - Cm_Pi_singlet_g(ii)*Dm + Cm_Pi_triplet_g(ii)*Dm)./r.^Cm_power(ii) ) ;

M(2,1,:)=M(1,2,:);
M(2,2,:)=squeeze(M(2,2,:)).' -(1/6)*          ( (  4*Cm_Sigma_triplet_g(ii)*Dm + Cm_Pi_singlet_g(ii)*Dm + Cm_Pi_triplet_g(ii)*Dm)./r.^Cm_power(ii) ) ;
M(2,3,:)=squeeze(M(2,3,:)).' -(1/(2*sqrt(3)))*( (                              - Cm_Pi_singlet_g(ii)*Dm + Cm_Pi_triplet_g(ii)*Dm)./r.^Cm_power(ii) ) ;

M(3,1,:)=M(1,3,:);
M(3,2,:)=M(2,3,:);
M(3,3,:)=squeeze(M(3,3,:)).' -(1/2)*          ((                                 Cm_Pi_singlet_g(ii)*Dm + Cm_Pi_triplet_g(ii)*Dm)./r.^Cm_power(ii) ) ;
end

for ii=1:length(r)
    eigenvalues=eig(M(:,:,ii));
    V_1_c(ii)=min(eigenvalues);
    V_2(ii)=max(eigenvalues);
end

%% 15.  Plot (6,6) Li2 potential from Semczuk 2013 PRA
close('all')
figure1=figure(1);
hold('on')
plot(r,V_1_c);
plot(r,V_2)
set(gca,'Position',[0.101190476190476 0.105359342915811 0.885 0.885])
set(gcf,'Color','w')
axis([50 200 -1e0 1e0])