% ==============================================================================================================
%   Created by Renan Vieira on 04/19.
%   Copyright (c) 2017 Renan Vieira. All rights reserved.
%   Code might be used as long as the author is properly credited.
%   This script is an attempt to implement the solutions for vertical wells
%   in multilayer reservoirs.
% ==============================================================================================================

clear; clc; close all     %clearing the console and the memory, and closing all graphical windows===============
%%
% setting the unit conversion constants - DONT CHANGE!!!
global alphap alphat
alphap=19.03;
alphat=0.0003484;
% setting the time vector size (must be an even number, to include falloff)
dim=142/2;
% setting the first timestep
t0=9.6d-6;
deltat=10.0^(0.1);

%%
% setting the input parameters
% initializing the number of layers
nlayers=2;
% setting oil and water viscosities
mio=2.3*ones(1,nlayers);
miw=0.52;%*ones(1,nlayers);
% setting the injection flow-rate
qinj=100.0;
% setting the wellbore radius
rw=0.108;

%%
% initializing the layer properties vectors
kj=zeros(1,nlayers);
hj=kj;%zeros(1,nlayers);
phij=kj;%zeros(1,nlayers);
cr=8e-5*ones(1,nlayers);
co=1.138e-4*ones(1,nlayers);
cw=4.04e-5*ones(1,nlayers);
ct=kj;
rskinj=rw*ones(1,nlayers);
kskinj=kj;%zeros(1,nlayers);

% setting layer properties
hj(1)=15.0;
hj(2)=10.0;
% hj(3)=7.0;
kskinj(1)=200.0;
kskinj(2)=200.0;
% kskinj(3)=300.0;
for j=1:nlayers
    % setting layer permeability
    if j == 1
        kj(j)=600.0;
    end
    
    if j == 2
        kj(j)=600.0;
    end
    
%     if j == 3
%         kj(j)=800.0;
%     end
    
    % setting layer thickness
%     hj(j)=20.0;
    % setting layer skin zone permeability
%     kskinj(j)=kj(j)*0.25;
    % setting layer skin zone radius
    if (j == 1)
        rskinj(j)=rw+0.0;
    end
    
    if (j == 2)
        rskinj(j)=rw+0.0;
    end
    
%     if (j == 3)
%         rskinj(j)=rw+0.2;
%     end
    % setting layer porosity
    phij(j)=0.25;
end

%%
% setting the relative permeability data
% initializing the initial water saturation
swi=0.25;
% initializing the residual oil saturation
sor=0.28;

% water relative permeability data
krw(1)=0.000000000; krw(2)=0.004730000; krw(3)=0.009458000;
krw(4)=0.014188000; krw(5)=0.018916000; krw(6)=0.023646000;
krw(7)=0.028374000; krw(8)=0.033110000; krw(9)=0.037873000;
krw(10)=0.042757000; krw(11)=0.047613000; krw(12)=0.052476000;
krw(13)=0.058981000; krw(14)=0.066197000; krw(15)=0.073700000;
krw(16)=0.082294000; krw(17)=0.095461000; krw(18)=0.110594000;
krw(19)=0.127831000; krw(20)=0.148597000; krw(21)=0.173000000;

% oil relative permeability data
kro(1)=0.540300000; kro(2)=0.500145000; kro(3)=0.459990000;
kro(4)=0.419835000; kro(5)=0.379685000; kro(6)=0.339530000;
kro(7)=0.299375000; kro(8)=0.261046000; kro(9)=0.223727000;
kro(10)=0.188473000; kro(11)=0.154791000; kro(12)=0.122027000;
kro(13)=0.095422000; kro(14)=0.071638000; kro(15)=0.051815000;
kro(16)=0.034833000; kro(17)=0.024335000; kro(18)=0.015993000;
kro(19)=0.009207000; kro(20)=0.003933000; kro(21)=0.000000000;

% transposing the relative permeability vectors
krw=krw';
kro=kro';

%%
% computing the endpoint mobilities
lohat=kro(1)./mio;
lwhat=krw(length(kro))./miw;
% computing the fractional flow and total mobility datanlay
sw=zeros(length(kro),nlayers);
lambdat=sw;
dfw=sw;
for j=1:nlayers
    [sw(:,j),lambdat(:,j),dfw(:,j),ct(j)]=fill_data(length(kro),krw,kro,swi,sor,mio(j),miw,co(j),cw(j),cr(j));
end

% computing layer hidraulic diffusivities
etaj=kj.*kro(1)./phij./mio./ct;
% computing layer skin factors
Sj=(kj./kskinj-1).*log(rskinj./rw);

%%
% computing pressure data
[t,tp,pwf,deltapo,deltapl,qj] = compute_pwf(nlayers,dim,t0,deltat,rw,qinj,kj,hj,phij,etaj,ct,kskinj,rskinj,Sj,mio,miw,lohat,lwhat,sw,dfw,lambdat);
% computing pressure derivative data
[dpwf]=compute_derivative(t,pwf,tp);
% % checking the results with the line-source solution
% 
% a=alphap*qinj/2/kj(1)/hj(1)/lohat;
% po=a.*log(t)+a*log(4*alphat*kj*lohat/exp(0.5522)/phij/ct/rw/rw);
%%
% plotting the results during injection
figure(1)
if (t(end)==tp)
    g1=loglog(t(1:dim),pwf(1:dim), t(1:dim),dpwf(1:dim));
else
    g1=loglog(t(1:dim/2),pwf(1:dim/2), t(1:dim/2),dpwf(1:dim/2));
end
g1(1).LineWidth=1.5;
g1(2).LineWidth=1.5;
a=round(log10(t(1))); a=max([10^(a) 1e-4]); 
b=floor(log10(tp))+1; b=10^b;
c=min(min((pwf)),min(abs(dpwf))); c=10^log10(floor(c));
d=floor(log10(max(pwf)))+1; d=10^d;
axis([a,b,c,d])
title('Pressure and Pressure Derivative during Injection')
xlabel('t (h)')
ylabel('P, deltaP (kgf/cm²)')
grid on
% plotting results during falloff
% figure(2)
% g2=loglog(t(dim/2+1:dim)-tp,deltapo(dim/2+1:dim), t(dim/2+1:dim)-tp,dpwf(dim/2+1:dim));
% g2(1).LineWidth=1.5;
% g2(2).LineWidth=1.5;
% axis([a,b,c,d])
% title('Pressure and Pressure Derivative during Falloff')
% xlabel('t (h)')
% ylabel('P, deltaP (kgf/cm²)')

% res=[t(11:dim/2) pwf(11:dim/2) dpwf(11:dim/2) deltapo(11:dim/2) deltapl(11:dim/2)];
% res=[res; t(dim/2+11:dim) pwf(dim/2+11:dim) dpwf(dim/2+11:dim) deltapo(dim/2+11:dim) deltapl(dim/2+11:dim)];

%%
% clearing the console
clc 
% clearing some auxiliary variables
clear t0 deltat dim co cr cw ii j compsw swi sor alphap alphat a b c d g1 g2
