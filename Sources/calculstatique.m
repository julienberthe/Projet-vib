function toto=Statique_EF(matrice);

toto.Us=((matrice.Ks)^-1)*matrice.Fs;

toto.U(1)=0;

for j=2:size(donnee.Elem,2)
    toto.U(j:j+1)=toto.U(j:j+1)+toto.Us(j-1:j);
end
end