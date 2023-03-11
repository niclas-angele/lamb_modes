function [L,M]=dispersionreal(w,maxi,h)
%%%% give the real wavenumbers of a elastic plate of width 2h at frequency
%%%% w. maxi gives the maximum window length of search (usualy 5*w). Return
%%%% L the list of symmetric wavenumbers and M the list of antisymmetric
%%%% wavenumbers. 
    %dilation to come back to a plate of width h=1
    H=h; 
    w=w*h;
    maxi=maxi*h;
    ct=0.5;
    cl=0.9;
    h=1;
    p=@(k,w) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k,w) sqrt(w.^2./(ct^2)-k.^2);
    %antisymetric relation 
    a=@(k,w) (q(k,w).^2-k.^2).^2.*cos(q(k,w)*h)+(4*k.^2.*p(k,w)).*q(k,w).*cos(p(k,w)*h).*sin(q(k,w)*h)./sin(p(k,w)*h);
    %symmetric relation 
    s=@(k,w) (q(k,w).^2-k.^2).^2.*cos(p(k,w)*h)+(4*k.^2.*q(k,w)).*p(k,w).*cos(q(k,w)*h).*sin(p(k,w)*h)./sin(q(k,w)*h);
    %number of search points 
    point=min(max(floor(maxi*10),150),500); 
    k1=linspace(-maxi,maxi,point); 
    E=[]; %symmetric list 
    F=[]; %antisymmetric list 
    %rewritting of s and a for the min search 
    sym=@(x) log(abs(s(x,w)));
    ant=@(x) log(abs(a(x,w)));
    for i=1:point
            opts = optimset('Display','off'); 
            t=fminsearch(sym,k1(i),opts);
            u=fminsearch(ant,k1(i),opts);
            E=[E,t];
            F=[F,u];
    end
    %remove all the abberant points 
    L=triE(E);
    M=triE(F);
    %sort using the right going classification 
    L=ordre(L);
    M=ordre(M);
    %remove the left going using the condition \partial_k \omega >0 
    ind=[];
    for i=1:length(L)
        if imag(L(i))==0
            symb=@(x) log(abs(s(x,w+0.01)));
            symc=@(x) log(abs(s(x,w-0.01)));
            testb=fminsearch(symb,L(i),opts);
            testc=fminsearch(symc,L(i),opts);
            if testb<testc 
                ind=[ind,i];
            end
        end
    end
    L(ind)=[];
    ind=[];
    for i=1:length(M)
        if imag(M(i))==0
            antb=@(x) log(abs(a(x,w+0.01)));
            antc=@(x) log(abs(a(x,w-0.01)));
            testb=fminsearch(antb,M(i),opts);
            testc=fminsearch(antc,M(i),opts);
            if testb<testc 
                ind=[ind,i];
            end
        end
    end
    M(ind)=[];
    %come back to the plate of width h
    L=L/H; 
    M=M/H; 
    %remove the vanishing values 
    if abs(L(length(L)))<10^(-8)
        L(length(L))=[]; 
    end
    if abs(M(length(M)))<10^(-8)
        M(length(M))=[]; 
    end
end

function L=triE(E)
%%%% remove duplicates and lonely points in a list E
    l=length(E);
    LL=[]; 
    count=[]; %count of doublons 
    for i=1:l
        test=0;
        for j=1:length(LL)
            if abs(LL(j)-E(i))<10^(-2) %E(i) is a doublon ?
                count(j)=count(j)+1; %if yes, we count it 
                test=1;
            end
        end
        if test==0  %if not doublon, we add it to the list 
            LL=[LL,E(i)];
            count=[count,1]; %and initialize its count to 1
        end
    end
    L=[]; % nouvelle liste pour enlever les points abérants
    for i=1:length(LL)
        if count(i)>10 %we keep the point if it appears 10 times 
            L=[L,LL(i)];
        end
    end
end

function L=ordre(E)
%%% sort values of E using the right going classification 
    for i=1:length(E)
        if imag(E(i))<10^(-8)
            E(i)=real(E(i));
        end
    end
    M=zeros(length(E),2);
    for i=1:length(E)
        M(i,1)=real(E(i));
        M(i,2)=imag(E(i));
    end
    M=sortrows(M,[2,-1]);
    for i=1:length(E)-1
        if abs(M(i,1)+M(i+1,1))<10^(-6)
            if real(M(i,1))<0
                M([i,i+1],:)=M([i+1,i],:);
            end
        end
    end
    L=zeros(1,length(E));
    for i=1:length(E)
        L(i)=M(i,1)+1i*M(i,2);
    end
end
