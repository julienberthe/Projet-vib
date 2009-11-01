function toto=Construction_EF(donnee);


%calcul de la matrice de rigidite
	K_elem=[1 -1;-1 1];			%matrice elementaire
	toto.K_ef=zeros(size(donnee.Elem,2)+1);		%initialisation de la matrice de rigidite globale

	%assemblage
	for j=1:size(donnee.Elem,2) %ici on fait une boucle sur tous les elements
		toto.K_ef(j:j+1,j:j+1)=toto.K_ef(j:j+1,j:j+1)+	(donnee.Elem{j}.young*donnee.Elem{j}.S/donnee.Elem{j}.dx)*K_elem;
	end




%calcul de la matrice de masse
	M_elem=[1/2 0;0 1/2];			%matrice elementaire
	toto.M=zeros(size(donnee.Elem,2)+1);		%initialisation de la matrice de masse globale

	%assemblage
	for j=1:size(donnee.Elem,2) %ici on fait une boucle sur tous les elements
		toto.M(j:j+1,j:j+1)=toto.M(j:j+1,j:j+1)+(donnee.Elem{j}.S*donnee.Elem{j}.dx*donnee.Elem{j}.rho)*M_elem;
	end

		
%calcul de la matrice d'amortissement
	toto.A=donnee.mat.alpha*toto.K_ef+ donnee.mat.beta*toto.M;
	
	
end