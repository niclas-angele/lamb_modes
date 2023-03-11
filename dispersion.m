function [L,M]=dispersion(w,max)
%%%% give dispersion curves in an elastic plate associated to the frequency
%%%% w. L contain symmetric wavenumbers and M antisymmetric ones. max is
%%%% the length of the search windox (usualy 5*w). 
    ct=0.5;
    cl=0.9;
    h=1;
    p=@(k,w) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k,w) sqrt(w.^2./(ct^2)-k.^2);
    a=@(k,w) (q(k,w).^2-k.^2).^2.*sin(p(k,w)*h).*cos(q(k,w)*h)./(4*k.^2.*p(k,w))+q(k,w).*cos(p(k,w)*h).*sin(q(k,w)*h);
    s=@(k,w) (q(k,w).^2-k.^2).^2.*sin(q(k,w)*h).*cos(p(k,w)*h)./(4*k.^2.*q(k,w))+p(k,w).*cos(q(k,w)*h).*sin(p(k,w)*h);
    %search points
    point=min(floor(max*10),400);
    k1=linspace(-max,max,point); %real part of k
    k2=linspace(-0.3,2*max,floor(point/2)); %imaginary part of k 
    E=[]; 
    F=[]; 
    %new relation to find minimum values 
    sym=@(x) log(abs(s(x(1)+1i*x(2),w)));
    ant=@(x) log(abs(a(x(1)+1i*x(2),w)));
    for i=1:point
        for j=1:(point/2)
            opts = optimset('Display','off'); 
            t=fminsearch(sym,[k1(i),k2(j)],opts);
            u=fminsearch(ant,[k1(i),k2(j)],opts);
            if max>t(1) && t(1)>-max
                if t(2)>-0.0001
                    E=[E,t(1)+1i*t(2)];
                end 
            end
            if max>u(1) && u(1)>-max
                if u(2)>-0.0001 
                      F=[F,u(1)+1i*u(2)];
                end 
            end
        end
    end
    %we delete doublons points 
    L=triE(E);
    M=triE(F);
    %we sort according to the right going classification 
    L=ordre(L);
    M=ordre(M);
    %we remove left going modes (\partial_k \omega<0)
    ind=[];
    for i=1:length(L)
        if imag(L(i))==0
            symb=@(x) log(abs(s(x,w+0.0001)));
            symc=@(x) log(abs(s(x,w-0.0001)));
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
end

function L=triE(E)
%%%% fonction qui trie une liste donnée en enlevant les points qui sont
%%%% tout seuls et les doublons pour les points au même endroits 
    l=length(E);
    LL=[]; %liste intermédiaire où on enlève les doublons 
    count=[]; %compte du nombre de doublons qu'il y avait 
    for i=1:l
        test=0;
        for j=1:length(LL)
            if abs(LL(j)-E(i))<10^(-5) %test pour dire si E(i) est le doublon de qqun
                count(j)=count(j)+1; %si oui on le compte
                test=1; %et on dit que c'était le doublon de qqun
            end
        end
        if test==0  %si c'est pas un doublon, on l'ajoute à notre liste
            LL=[LL,E(i)];
            count=[count,1]; %et il est présent une fois 
        end
    end
    L=[]; % nouvelle liste pour enlever les points abérants
    for i=1:length(LL)
        if count(i)>2 % on veut que le point apparaisse au moins 3 fois pour le garder 
            L=[L,LL(i)];
        end
    end
end

function L=ordre(E)
%%% fonction qui range les valeurs par ordre crois imaginaire et dec réel 
    for i=1:length(E)
        if imag(E(i))<10^(-8) %&& -10^(-4)<imag(E(i))
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
