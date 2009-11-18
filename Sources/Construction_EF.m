function toto=Construction_EF(donnee);


%calcul de la matrice de rigidite
	K_elem=[1 -1;-1 1];			%matrice elementaire
	toto.K_ef=zeros(size(donnee.Elem,2)+1);		%initialisation de la matrice de rigidite globale

	%assemblage
	for j=1:size(donnee.Elem,2) %ici on fait une boucle sur tous les elements
        %donnee.Elem{j}.young*donnee.Elem{j}.S/donnee.Elem{j}.dx
		toto.K_ef(j:j+1,j:j+1)=toto.K_ef(j:j+1,j:j+1)+	(donnee.Elem{j}.young*donnee.Elem{j}.S/donnee.Elem{j}.dx)*K_elem;
	end




%calcul de la matrice de masse
	M_elem=[1/3 1/6;1/6 1/3];			%matrice elementaire
	toto.M=zeros(size(donnee.Elem,2)+1);		%initialisation de la matrice de masse globale

	%assemblage
	for j=1:size(donnee.Elem,2) %ici on fait une boucle sur tous les elements
        
		toto.M(j:j+1,j:j+1)=toto.M(j:j+1,j:j+1)+(donnee.Elem{j}.S*donnee.Elem{j}.dx*donnee.Elem{j}.rho)*M_elem;
	end

		
%calcul de la matrice d'amortissement
	toto.A=donnee.mat.alpha*toto.K_ef+ donnee.mat.beta*toto.M;
    
    %construction du vecteur F pour le cas statique
    toto.F=zeros(donnee.nelem,1);
    toto.F(donnee.nelem)=donnee.mat.S;
 
	
end