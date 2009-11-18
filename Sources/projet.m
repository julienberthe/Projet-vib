%--------------------------------------------------------------------------
%                programme ecrit par Bibou et Loulou
%
%                 "petit" projet numerique du TACS
%                      EF classiques et malins
%
%==========================================================================


%Definition du probleme:

% on considere un poutre elastique de longueur L, de section S et de module
% de Young E
close all;
clear; clc;
tps1=tic;
disp('I	Construction du probleme mecanique EF');
	%parametres materiaux
	donnee.mat.L=1;			%longueur:
	donnee.mat.S=10^(-4);		%section:
	donnee.mat.E=220*10^9;		%module de young
	donnee.mat.rho=70000;		%masse volumique
	donnee.mat.alpha=0;		%coefficient d'amortissement C=alpha*K+beta*M
	donnee.mat.beta=0;
	%parametres utiles par la suite
	donnee.changementbase='non';	
	donnee.statique='non';
	
	%Parametres de la methode EF
	donnee.nelem = 160;	%nombre d'elements
	donnee.npas  = 500;	%nombre de pas de temps

	%caracteristique du probleme
	donnee.T= 0.1;	%temps d'etude
    
    %mise en donnee
	for i=1:donnee.nelem
		donnee.Elem{i}.xinit=(i-1)*donnee.mat.L/donnee.nelem;
		donnee.Elem{i}.xfinal=i*donnee.mat.L/donnee.nelem;
		donnee.Elem{i}.dx=donnee.mat.L/donnee.nelem;
		donnee.Elem{i}.S=donnee.mat.S;
		donnee.Elem{i}.young=donnee.mat.E;
		donnee.Elem{i}.rho=donnee.mat.rho;
	end

	donnee.dt =donnee.T/donnee.npas; 
	donnee.t  =[0:donnee.dt:donnee.T];
	donnee.dx =donnee.mat.L/donnee.nelem;
	donnee.x  =[0:donnee.dx:donnee.mat.L];
    donnee.f = donnee.mat.E*donnee.mat.S;

	%construction des differentes matrice du probleme EF
disp('II	Construction des matrices');
    matrice=Construction_EF(donnee);
    
    %Prise en comtpe des CL (methode de substitution)
disp('IIb   Conditionnement des matrices par subsitution')
    [matriceS,donneeS]=Substitution(matrice,donnee);
    
    
disp('III	Construction du chargement');
	% type de chargement accessible: echelon en bout de poutre, creneau en bout de poutre, harmonique... pour plus d'info voir Construction_Chargement.m
	chargement.type = 'echelon en bout de poutre';
	chargement.parametre{1}=10000;				%amplitude
	chargement.parametre{2}=[0.2 0.5];		%quand ?
    %chargement.type = 'harmonique';
    %chargement.parametre{1}=10;		%amplitude
    %chargement.parametre{2}=60;		%frequence
	chargement=Construction_Chargement(chargement,donnee);

	%affichage du chargement, pour voir les options d'affichage voir Affichage.m
	option.type='en fonction du temps';
	option.titre='effort en bout de poutre en fonction du temps';
    Affichage(chargement.F(donnee.nelem,:),donnee,option)
	
    
    disp('========================================');
    disp('==============Résolution EF=============');
    disp('========================================');
    option.schema='euler AR';  %%champs possibles 'euler AR', 'euler AV' (non implémenté) ou 'newmark' (non implémenté)
    option.resolution='modale'; %%champs possibles 'modale' ou 'directe'
    %nombre de modes propres à pendre en compte
    option.Nb_ef=5;
    [ModePropre,U,Eps]=EF(chargement,matriceS,donneeS,option);
    
    
    %%%%%Affichages divers et variés
    option.type='4 Modes Propres';
    option.Mode=[1 2 3 4];
    option.titre=sprintf('Mode propre numero %d',option.Mode);

    Affichage(ModePropre,donneeS,option);
    %Affichage(ModePropree,donneee,option);
        
    option.type='en fonction du temps';
	option.titre='déplacement extrémité de la poutre en fct du tps';
	Affichage(U.U(donnee.nelem+1,:),donneeS,option)
    option.titre='vitesse extrémité de la poutre en fct du tps';
    %Affichage(U.V(donnee.nelem+1,:),donnee,option)
    
    
    option.titre='accélération extrémité de la poutre en fct du tps';
    %Affichage(U.A(donnee.nelem+1,:),donnee,option)
    
    option.type='animation en fonction du temps';
	option.titre='Déplacement des points en fct du temps';
    
    option.save='non';    %%champs possibles 'non', 'film' ou 'images'
    option.dossier='echelon';
    %option.nbm=ModePropre.Nb_ef;
	%Affichage(U,donnee,option)

    option.type='3D';
    option.titre='Déplacement de la poutre en fonction du temps et de labscisse';
    Affichage(U.U,donnee,option);
    option.titre='Vitesse de la poutre en fonction du temps et de labscisse';
    Affichage(U.V,donnee,option);
    
   tps2=toc(tps1);
   disp('Temps de calcul: ')
   disp(tps2)