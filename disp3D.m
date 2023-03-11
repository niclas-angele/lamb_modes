 %%%% To do the general 3D plot: 

 %%%% interesting frequencies and zoom around large variations
W=sort([linspace(0.2,4.5,20), linspace(1.29,1.33,5),linspace(1.4,1.58,10),linspace(3.1,3.2,7),linspace(4.15,4.25,7)]); 
%storage of 30 lamb modes 
M=zeros(length(W),30); 

%Rayleigh Lamb solutions for each frequency in W
for i=1:length(W)
    w=W(i); 
    [S]=dispersionn(w);
    M(i,1:length(S))=S; 
end
L=M;

%%%% extraction of each mode. s: reel, i: inhomogeneous, e: envanescent 
%%%S0
figure
S=[0;M(:,1)];
plot3(imag(S),real(S),[0,W],'Color','#7E2F8E')
axis([0,9,0,9,0,4.5])
L(:,1)=[];
grid on
hold on 

%S1
S=[0.86;M(10:49,2)];
plot3(imag(S),real(S),[1.33,W(10:49)],'Color','#7E2F8E')
L(10:49,1:27)=M(10:49,3:29);

%I1
S=[1.1+1i*2.1;L(1:9,1);0.86];
plot3(imag(S),real(S),[0,W(1:9),1.33],'Color','#EDB120')

%S1b
S=[0.86;L(10:12,1);0]; 
plot3(imag(S),real(S),[1.33,W(10:12),1.41],'Color','#7E2F8E')

%E1
S=[0;L(13:21,1);0]; 
plot3(imag(S),real(S),[1.41,W(13:21),1.57],'Color','#77AC30')

%S2
S=[0;L(22:49,1)];
plot3(imag(S),real(S),[1.57,W(22:49)],'Color','#7E2F8E')
L(:,1)=[];

%S3
S=[0;L(32:49,1)];
plot3(imag(S),real(S),[3.14,W(32:49)],'Color','#7E2F8E')
L(32:49,1:27)=L(32:49,2:28);

%E2b
S=[1i*1.67;L(29:31,1);0];
plot3(imag(S),real(S),[3.09,W(29:31),3.14],'Color','#77AC30')
L(29:31,1:27)=L(29:31,2:28);

%E2
S=[1i*1.67;L(29:46,1);0];
plot3(imag(S),real(S),[3.09,W(29:46),4.24],'Color','#77AC30')

%S4
S=[0;L(47:49,1)];
plot3(imag(S),real(S),[4.24,W(47:49)],'Color','#7E2F8E')

%I2
S=[1.55+1i*5.35;L(1:28,1);0.5+2*1i;1i*1.67];
plot3(imag(S),real(S),[0,W(1:28),3,3.09],'Color','#EDB120')
L(:,1)=[]; 

%E3b
S=[5.13*1i;L(41:49,1)]; 
plot3(imag(S),real(S),[4.13,W(41:49)],'Color','#77AC30')
L(41:49,1:26)=L(41:49,2:27); 

%E3
S=[5.13*1i;L(41:49,1)]; 
plot3(imag(S),real(S),[4.13,W(41:49)],'Color','#77AC30')

%I3
S=[1.77+1i*8.55;L(1:40,1);1i*5.13];
plot3(imag(S),real(S),[0,W(1:40),4.13],'Color','#EDB120')
L(:,1)=[]; 

%I4
S=[L(39:49,1)];
plot3(imag(S),real(S),[W(39:49)],'Color','#EDB120')

xlabel('imag(k_n)h')
ylabel('real(k_n)h')
zlabel('\omega h')

%%%%% critical points
%L and T
plot3(0,0.86,1.33,'.m','Markersize',30)
plot3(0,0,1.41,'.b','Markersize',30)
plot3(0,0,1.57,'.r','Markersize',30)
plot3(0,0,3.14,'.r','Markersize',30)
plot3(0,0,4.24,'.b','Markersize',30)
plot3(0,0.86,1.33,'.r','Markersize',30)
plot3(0,0,1.41,'.r','Markersize',30)
plot3(0,0,1.57,'.r','Markersize',30)
plot3(0,0,3.14,'.r','Markersize',30)
plot3(0,0,4.24,'.r','Markersize',30)
%ZGV
plot3(1.67,0,3.09,'.r','Markersize',30)
plot3(5.13,0,4.13,'.r','Markersize',30)
 view(500,15)


%%%%% zoom around L point
figure 
grid on 
%new frequency set 
W=linspace(1.3, 1.5,50);
%storage of lamb modes 
Im=W*0; 
Re=W*0; 
Re2=W*0; 
%real Rayleigh Lamb solutions for each frequency in W
for i=1:length(W)
    w=W(i); 
    [ks1,ka1]=dispersionimag(w,5*w,1);
    if length(ks1)>0 
        Im(i)=imag(ks1(1)); 
    end
end
%imaginary Rayleigh Lamb solutions for each frequency in W
for i=1:length(W)
    w=W(i);
    [ks,ka]=dispersionreal(w,5*w,1);
    Re(i)=ks(length(ks));
    if length(ks)>1
        Re2(i)=ks(length(ks)-1); 
    end
end

%%%%extraction of modes 
Re(1:7)=[]; 
Re2(1:7)=[];
Re2(22:43)=Re(22:43); 
Re(22:43)=[];
Re=[-0.86,Re,0]; 
Re2=[0.86,Re2]; 

%%%plots 
plot3(0*Re,-Re,[1.326,W(8:28),1.412],'b')
hold on 
plot3(0*Re,Re,[1.326,W(8:28),1.412],'b')
plot3(Re2*0,Re2,[1.326,W(8:50)],'b')
plot3(Re2*0,-Re2,[1.326,W(8:50)],'b')
plot3(Im(28:50),Im(28:50)*0, [1.412,W(29:50)],'r')
plot3(-Im(28:50),Im(28:50)*0, [1.412,W(29:50)],'r')
plot3(0,0,1.412,'.m','Markersize',30)

%%% ZGV zoom
figure 
grid on 
W=linspace(1.3,1.35,40);
M=zeros(length(W),30); 
for i=1:length(W)
    w=W(i); 
    [S]=dispersionn(w);
    M(i,1:length(S))=S; 
end
L=M;

S1=[0.858;M(22:40,2)];
I1=[M(1:21,2);0.858]; 
S2=[0.858;M(22:40,3)]; 
W1=[1.326,W(22:40)];
W2=[W(1:21),1.326];


plot3(imag(S1),real(S1),W1,'r')
hold on 
grid on
plot3(imag(S2),real(S2),W1,'r')
plot3(imag(I1),real(I1),W2,'b')
plot3(-imag(I1),real(I1),W2,'b')
plot3(imag(S2),-real(S2),W1,'r')
plot3(imag(I1),-real(I1),W2,'b')
plot3(-imag(I1),-real(I1),W2,'b')
plot3(imag(S1),-real(S1),W1,'r')
plot3(0,0.858,1.326,'.m','Markersize',30)
plot3(0,-0.858,1.326,'.m','Markersize',30)

