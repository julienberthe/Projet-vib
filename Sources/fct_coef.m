%%Fonction définissant les coefficients de la superposition modale
function toto=fct_coef(ed,t)
for iter=1:size(ed.A)
    toto(iter)=ed.A(iter)*cos(ed.W(iter)*t)+ed.B(iter)*sin(ed.W(iter)*t)+ed.C(iter);
end