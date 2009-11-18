%%Fichier de résolution EF

function [ModePropre,U,Eps]=EF(chargement,matrice,donnee,option)


switch option.resolution
    case 'modale'
        disp('Résolution sur la base des modes propres')
        disp('A	Calcul des modes et valeurs propres');
        ModePropre=CalculModePropre(matrice,donnee);
        
    

        disp('B	Détermination du mode statique')    
        SolutionStatique=Statique_EF(matrice,donnee);
        % figure
        %plot(donnee.x,SolutionStatique.U)
    
        disp('C	Resolution du probleme EF sur la base des modes propres');
        [U,Eps]=Resolution_EFSM(chargement,donnee,ModePropre,SolutionStatique,option);
    
    
    case 'directe'        
        disp('Resolution directe');
        [U,Eps]=Resolution_EFD(chargement,donnee,matrice,option);
        ModePropre=0;
end
