%détermination du déplacement "statique"
function toto=Statique_EF(matrice,donnee);

 toto.Us=matrice.K_ef^-1*matrice.F;
%toto.Us=matrice.Ks\matrice.Fs';
toto.Ug=toto.Us';
toto.U(1:size(donnee.Elem,2)+1)=0;

for j=2:size(donnee.Elem,2)
    toto.U(j:j+1)=toto.Ug(j-1:j);
end
end