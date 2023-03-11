function [L]=dispersionn(w)
%%%% function which gives the values for a given frequency w of 
%%%% the k of the Lamb wave dispersion curve (sym and antisym). We set a 
%%%% max search window for k, such that real(k)<max and imag(k)<5max. 
%%%% To search we place ourselves on a fine grid and we start from each 
%%%% point of the grid to look for the zeros of the dispersion relation.
    ct=0.5;
    cl=0.9;
    h=1;
    p=@(k,w) sqrt(w.^2./(cl^2)-k.^2);
    q=@(k,w) sqrt(w.^2./(ct^2)-k.^2);
    s=@(k,w) (q(k,w).^2-k.^2).^2.*sin(q(k,w)*h).*cos(p(k,w)*h)./(4*k.^2.*q(k,w))+p(k,w).*cos(q(k,w)*h).*sin(p(k,w)*h);
    point=floor(50); %number of search starting points min
    k1=linspace(0,10,point); %real part of k
    k2=linspace(-0.2,10,point); %imaginary part of k
    E=[]; %list of symmetric points
    sym=@(x) log(abs(s(x(1)+1i*x(2),w)));
    for i=1:point
        for j=1:point
            opts = optimset('Display','off'); 
            t=fminsearch(sym,[k1(i),k2(j)],opts);
            %we check that the mins found are in our search window then 
            % if so we add them to the list
            if 10>t(1) && t(1)>0
                if t(2)>-0.0001 %&& 5*max>t(2)
                    E=[E,t(1)+1i*t(2)];
                end 
            end
        end
    end
    %we remove the aberrant points and the duplicates
    L=triE(E);
    %we sort according to the classification of right going waves
    L=ordre(L);
end

function L=triE(E)
%%%% function that sorts a given list by removing the points that are 
%%%% alone and duplicates for points in the same places
    l=length(E);
    LL=[]; %intermediate list where duplicates are removed
    count=[]; %count of the number of duplicates there were
    for i=1:l
        test=0;
        for j=1:length(LL)
            if abs(LL(j)-E(i))<10^(-5) %test to say if E(i) is the duplicate of someone
                count(j)=count(j)+1; %if so, we count it
                test=1; %otherwise, it is the duplicate of someone
            end
        end
        if test==0  %if it is not a duplicate, we add it to the list
            LL=[LL,E(i)];
            count=[count,1]; %and it is counted for 1
        end
    end
    L=[]; % new list to remove outliers
    for i=1:length(LL)
        if count(i)>2 % we want the point to appear at least 3 times to keep it
            L=[L,LL(i)];
        end
    end
end

function L=ordre(E)
%%% function that arranges values in ascending imaginary order and descending real order
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


