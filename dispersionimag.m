function [L,M]=dispersionimag(w,maxi,h)
%%%% similar to disperionn, but extract only the imaginary part of k  
    H=h; 
    w=w*h;
    maxi=maxi*h;
    ct=0.5;
    cl=0.9;
    h=1;
    p=@(k,w) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k,w) sqrt(w.^2./(ct^2)-k.^2);
    a=@(k,w) (q(k,w).^2-k.^2).^2.*sin(p(k,w)*h).*cos(q(k,w)*h)+4*k.^2.*p(k,w).*q(k,w).*cos(p(k,w)*h).*sin(q(k,w)*h);
    s=@(k,w) (q(k,w).^2-k.^2).^2.*sin(q(k,w)*h).*cos(p(k,w)*h)+4*k.^2.*q(k,w).*p(k,w).*cos(q(k,w)*h).*sin(p(k,w)*h);
    point=min(max(floor(maxi*10),100),300);  
    k1=linspace(0,maxi,point);
    E=[]; 
    F=[]; 
    sym=@(x) log(abs(s(1i*x,w)));
    ant=@(x) log(abs(a(1i*x,w)));
    for i=1:point
            opts = optimset('Display','off'); 
            t=fminsearch(sym,k1(i),opts);
            u=fminsearch(ant,k1(i),opts);
            E=[E,1i*t];
            F=[F,1i*u];
    end 
    L=triE(E);
    M=triE(F);
    L=ordre(L);
    M=ordre(M);
    ind=[];
    for i=1:length(L)
        if imag(L(i))==0
            symb=@(x) log(abs(s(x,w+0.0000001)));
            symc=@(x) log(abs(s(x,w-0.0000001)));
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
    if length(L)>0 && abs(L(1))<10^(-5)
        L(1)=[]; 
    end
    if length(M)>0 && abs(M(1))<10^(-5)
        M(1)=[]; 
    end
    L=L/H; 
    M=M/H; 
end

function L=triE(E) 
    l=length(E);
    LL=[]; 
    count=[]; 
    for i=1:l
        test=0;
        for j=1:length(LL)
            if abs(LL(j)-E(i))<10^(-3) 
                count(j)=count(j)+1;
                test=1; 
            end
        end
        if test==0  
            LL=[LL,E(i)];
            count=[count,1]; 
        end
    end
    L=LL; 
end

function L=ordre(E)
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
