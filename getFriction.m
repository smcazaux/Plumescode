% friction calc.
function [friction] = getFriction(tau, r, L, A_d, throat, c, dx)
%friction = tau.*pi.*2.*r;%.*L./sqrt(A_d(throat)); % the "L./sqrt(A_d(throat)" is to non-dimensionalize 
friction = tau.*c.*dx;
end