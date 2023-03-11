function [S,A]=dispersionh(w,max,h)
%%% same function as dispersion but return the dispersion relation for a
%%% general waveguide of width 2h
    [S,A]=dispersion(w*h,max*h); 
    S=S/h; 
    A=A/h;
end