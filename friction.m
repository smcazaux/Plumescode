% friction calc.
function friction = getFriction(tau, r, L, A)
friction = tau.*pi.*2.*r.*L./sqrt(A(throat)); % the "L./sqrt(A(throat)" is to non-dimensionalize 
end