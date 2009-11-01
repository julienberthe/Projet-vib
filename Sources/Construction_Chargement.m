function toto=Construction_Chargement(chargement,donnee)

%initialisation
toto=chargement;
F=zeros(donnee.nelem,donnee.npas+1);

%construction du chargement
switch chargement.type

	case 'echelon en bout de poutre'
		%on cherche les points les plus pres des abscisses adimensionnees
		pas1=find((abs(donnee.t-donnee.T*chargement.parametre{2}(1))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(1)))));
		pas2=find((abs(donnee.t-donnee.T*chargement.parametre{2}(2))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(2)))));
		F(donnee.nelem,pas1:pas2)=chargement.parametre{1};
		
	case 'echelon en bout de poutre doux'
		%on cherche les points les plus pres des abscisses adimensionnees
		pas1=find((abs(donnee.t-donnee.T*chargement.parametre{2}(1))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(1)))));
		pas2=find((abs(donnee.t-donnee.T*chargement.parametre{2}(2))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(2)))));
		if pas2<pas1+5
			disp('je ne peux pas construire un bel effort rampe-echelon...')
			F(donnee.nelem,pas1:pas2)=chargement.parametre{1};
		else
			F(donnee.nelem,pas1:pas1+5)=(0:5)*chargement.parametre{1}/4;
			F(donnee.nelem,pas1+5:pas2)=chargement.parametre{1};
		end
	case 'creneau en bout de poutre doux'
		%on cherche les points les plus pres des abscisses adimensionnees
		pas1=find((abs(donnee.t-donnee.T*chargement.parametre{2}(1))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(1)))));
		pas2=find((abs(donnee.t-donnee.T*chargement.parametre{2}(2))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(2)))));
		if pas2<pas1+5
			disp('je ne peux pas construire un bel effort rampe-echelon...')
			F(donnee.nelem,pas1:pas2)=chargement.parametre{1};
		elseif pas2+5>donnee.npas
			disp('je ne peux pas construire un bel effort rampe-echelon...')
			F(donnee.nelem,pas1:pas2)=chargement.parametre{1};
		else
			F(donnee.nelem,pas1:pas1+5)=(0:5)*chargement.parametre{1}/5;
			F(donnee.nelem,pas1+5:pas2)=chargement.parametre{1};
			F(donnee.nelem,pas2:pas2+5)=(5:-1:0)*chargement.parametre{1}/5;
		end
	case 'harmonique'
		F(donnee.nelem,:)=chargement.parametre{1}*sin(2*pi*chargement.parametre{2}*donnee.t);
	
end
toto.F=F;
end