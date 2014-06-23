function ypred=f_OI_gaspari_cohn(x1,x2,Ls,x1obs,x2obs,yobs,sigobs)

%% Specify observation error matrix
nobs=length(yobs);
R=eye(nobs)*sigobs^2; 

%% Calculate H matrix using nearest neighbour
H=zeros(nobs,length(x1));
for k=1:nobs
    z=x1obs(k)-x1+1i*(x2obs(k)-x2);
    [~,I]=min(abs(z));
    H(k,I)=1;
end

%% Calculate intragrid distances
z=x1(:)+1i*x2(:);
Z=repmat(z,1,length(z)); 
Dist=abs(Z-Z.');

%% Calculate B matrix using Gaspari-Cohn correlation function
B=0*Dist;
I1=find(Dist<Ls);           
I2=find(Dist>=Ls&Dist<2*Ls); 
z=Dist(I1)/Ls; B(I1)=-z.^5/4 +z.^4/2+z.^3*5/8-z.^2*5/3    +1;
z=Dist(I2)/Ls; B(I2)= z.^5/12-z.^4/2+z.^3*5/8+z.^2*5/3-5*z+4-2/3./z;

%% Calculate predicted y values from B, H, R and yobs
temp=(H*B*H'+R)\yobs(:);
ypred=B*H'*temp;

