function [tc,sc]=TS_convect1_fast(tt,ss,z1,dz1)

% input: tt -- potential temperature (deg C)
%        ss -- salinity (psu)
%        z1 -- depth (m)
%       dz1 -- level thickness (m)

tc2=tt(:); sc2=ss(:); z2=z1(:); dz2=dz1(:);

convect=1;

while (convect==1)
    %check stability
    stable=1;
    

    % FAST METHOD
    nk=length(dz2);
    pp=0.5*(z2(1:nk-1)+z2(2:nk));
    R1=rho_wright(sc2(1:nk-1), tc2(1:nk-1),pp);
    R2=rho_wright(sc2(2:nk), tc2(2:nk),pp);
    for k=1:length(dz2)-1
        k1=k+1;
        if R1(k)>R2(k)
            stable=0;
            break;
        end
    end

    % ORIGINAL METHOD
    %     for k=1:length(dz2)-1
    %         k1=k+1;
    %         pp=(z2(k)+z2(k1))/2;
    %         %zz=dz2(k)+dz2(k1)*(k1-k);
    %         r1=rho_wright(sc2(k),tc2(k),pp);
    %         r2=rho_wright(sc2(k1),tc2(k1),pp);
    %         if (r1> r2)
    %             stable=0;
    %             break;
    %         end
    %     end


    %end checking

    if (stable==1)        %stable
        convect=0;

    elseif (stable==0)    %unstable, then convective

        % note: the td is changed to use td(index)=x rather
        % than td=[td x] for speedup purposes only

        %      td=[];sd=[];zd=[];dzd=[];
        td=zeros(1,length(z2)); sd=td; zd=td; dzd=td;

        k=1;
        index=0;
        
        while k<=length(z2);
            index=index+1;
            k1=min(k+1,length(z2));
            pp=(z2(k)+z2(k1))/2;
            zz=dz2(k)+dz2(k1)*(k1-k);
            S=sc2([k k1]);T=tc2([k k1]);P=[pp; pp];
            r=rho_wright(S(:), T(:), P(:));
%             r1=rho_wright(sc2(k),tc2(k),pp);
%             r2=rho_wright(sc2(k1),tc2(k1),pp);

             if (r(1)>r(2))
                tdum=(tc2(k)*dz2(k)+tc2(k1)*dz2(k1))/zz;
                sdum=(sc2(k)*dz2(k)+sc2(k1)*dz2(k1))/zz;
                %           td=[td tdum];sd=[sd sdum];
                %           zd=[zd pp];dzd=[dzd zz];
                td(index)=tdum;  zd(index)=pp;
                sd(index)=sdum; dzd(index)=zz;
                k=k+2;
            else
                %           td=[td tc2(k)];sd=[sd sc2(k)];
                %           zd=[zd z2(k)]; dzd=[dzd dz2(k)];
                td(index)=tc2(k);  zd(index)=z2(k);
                sd(index)=sc2(k); dzd(index)=dz2(k);
                k=k+1;
            end
        end
        %       tc2=td; sc2=sd;z2=zd;dz2=dzd;
        tc2=td(1:index); sc2=sd(1:index);z2=zd(1:index);dz2=dzd(1:index);
    end  %end if stable

end  % end convection

dz22=dz2;
k=1;l=0;
while k<=length(dz2) && l<length(z1)
    l=l+1;
    tc(l)=tc2(k);
    sc(l)=sc2(k);
    dz22(k)=dz22(k)-dz1(l);
    %[k,l,dz22(k)]
    %disp([l k dz22(k)])
    if dz22(k)<=1e-2
        k=k+1;
    end
end

%subplot(121);plot(tt,-z1,'r',tc,-z1,'b');
%subplot(122);plot(ss,-z1,'r',sc,-z1,'b');

%subplot(121);plot(tt,-z1,'r',tc2,-z2,'b');
%subplot(122);plot(ss,-z1,'r',sc2,-z2,'b');

