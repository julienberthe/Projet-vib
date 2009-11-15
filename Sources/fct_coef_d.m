%%Dérivée de la fonction définissant les coefficients de la superposition modale
function fct=fct_coef_d(ed,t)
for i=1:size(ed.A,2)
    fct(i)=-ed.A(i)*ed.W(i)*cos(ed.W(i)*t)+ed.B(i)*ed.W(i)*sin(ed.W(i)*t);
end