
function R=rho_wright(S,Th,P);

%
% function R=rho_dwright(S,Th,P);

% EOS of seawater following D. G. Wright (1997, JAOT);

% Calculate in Sigma (R, in kg/m^3, related to in situ density rho) 

% using salnity (S, in psu) and potential temperature 

% (Th, in degree C) and pressure (P, in dbar, 

% approximately meters in water depth).

%
% All elements can be scalars, vectors, or matrices,

% but should be the same size.



if (nargin<3)
disp('Insufficient input arguments')

end



if (size(S)~=size(Th))
   
disp('T and S sizes do not match');
   
elseif (size(Th)~=size(P)) 
   
disp('size of P does not match that of T and S')

end  



% a0 a1 a2: values in DGW table1 *1e3


a0=7.133718e-1;
a1=2.724670e-4;
a2=-1.646582e-4;


% b0--b5: values in DGW table1 *1e-5


b0=5.613770e3;
b1=3.600337e1;
b2=-3.727194e-1;
b3=1.660557e-3;
b4=6.844158;
b5=-8.389457e-2;


%c0--c5: values in DGW table1 *1e-2 


c0=1.6098930e3;
c1=8.427815;
c2=-6.931554e-2;
c3=3.869318e-4;
c4=-1.664201;
c5=-2.765195e-2;


T=Th;
T2=T.*T;
T3=T2.*T;
ST=S.*T;

Alp0=T*a1+S*a2+a0;
P1=T*b1+T2*b2+T3*b3+S*b4+ST*b5+b0;

Ram=T*c1+T2*c2+T3*c3+S*c4+ST*c5+c0;

P0=P/10; 
% convert dbar to bar

P01=P0+P1;
dum=Alp0.*P01;
dum=dum+Ram;
R=P01./dum;   
% in g/cm^3


R=R-1;

R=R*1000+0.025;     
% in kg/m^3



