function ModePropre=CalculModePropre(matrice,donnee)

[VecteurPropre,ValeurPropre]=eig(matrice.M^-1*matrice.K_ef);
%  keyboard
for i=1:size(ValeurPropre,1)
	vp(i)=ValeurPropre(i,i);
end

[VP,index]=sort(vp,'ascend');

ModePropre.n=size(vp,2);
ModePropre.Matrice=VecteurPropre(:,index);
	
for i=1:size(ValeurPropre,1)
	ModePropre.Vecteur{i}	=ModePropre.Matrice(:,i);
	ModePropre.Valeur(i)	=VP(i);
end


%%%Normalisation des vecteurs propres
%%Caclul de l'intégrale int(rho*S*phi^2)

%for i=1:size(ValeurPropre,1)
%    int(i)=0;
%    for j=1:donnee.nelem-1
%        int(i)=int(i)+donnee.dx/2*donnee.mat.rho*donnee.mat.S*(ModePropre.Matrice(j,i)^2+ModePropre.Matrice(j+1,i)^2);
%    end
%end
%Normalisation des vecteurs propres
%for i=1:size(ValeurPropre,1)
%    ModePropre.Matrice(:,i)=ModePropre.Matrice(:,i)/int(i);
%end
end