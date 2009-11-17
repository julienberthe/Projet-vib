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
	donnee.nelem = 200;	%nombre d'elements
	donnee.npas  = 100;	%nombre de pas de temps

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

	%construction des differentes amtrice du probleme EF
disp('II	Construction de la matrice');
    matrice=Construction_EF(donnee);
	
disp('III	Construction du chargement');
	% type de chargement accessible: echelon en bout de poutre, creneau en bout de poutre, harmonique... pour plus d'info voir Construction_Chargement.m
	chargement.type = 'echelon en bout de poutre';
	chargement.parametre{1}=10000;				%amplitude
	chargement.parametre{2}=[0.2 0.6];		%quand ?
    %chargement.type = 'harmonique';
    %chargement.parametre{1}=10;		%amplitude
	%chargement.parametre{2}=60;		%frequence
	chargement=Construction_Chargement(chargement,donnee);

	%affichage du chargement, pour voir les options d'affichage voir Affichage.m
	option.type='en fonction du temps';
	option.titre='effort en bout de poutre en fonction du temps';
Affichage(chargement.F(donnee.nelem,:),donnee,option)
	
    
disp('IV	Calcul des modes et valeurs propres');
	ModePropre=CalculModePropre(matrice,donnee);
    matricee=matrice;
    
    matricee.K_ef(:,1)=[];
    matricee.K_ef(1,:)=[];
    matricee.M(:,1)=[];
    matricee.M(1,:)=[];
    donneee=donnee;
    donneee.x(1)=[];
    
    ModePropree=CalculModePropre(matricee,donnee);
	option.type='4 Modes Propres';
	option.Mode=[1 2 3 4];
	option.titre=sprintf('Mode propre numero %d',option.Mode);
    ModePropre=ModePropree;
	%Affichage(ModePropre,donnee,option);
    %Affichage(ModePropree,donneee,option);
    

    disp('IVb Détermination du mode statique')    
SolutionStatique=Statique_EF(matrice,donnee);
   % figure
    %plot(donnee.x,SolutionStatique.U)
    
disp('V	Resolution du probleme EF sur la base des modes propres');
%nombre de modes propres à pendre en compte
ModePropre.Nb_ef=150;
option.type='euler AR';
[U,Eps]=Resolution_EF(chargement,donnee,ModePropre,SolutionStatique,option);
    
    option.type='en fonction du temps';
	option.titre='déplacement extrémité de la poutre en fct du tps';
	Affichage(U.U(donnee.nelem+1,:),donnee,option)
    option.titre='vitesse extrémité de la poutre en fct du tps';
    
    %Affichage(U.V(donnee.nelem+1,:),donnee,option)
    option.titre='accélération extrémité de la poutre en fct du tps';
    %Affichage(U.A(donnee.nelem+1,:),donnee,option)
    
    option.type='animation en fonction du temps';
	option.titre='Déplacement des points en fct du temps';
	Affichage(U,donnee,option)
