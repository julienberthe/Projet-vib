function [toto1,toto2]=Resolution_EF(chargement,donnee,ModePropre,SolutionStatique,option)
%%initialisation
nbmode=ModePropre.Nb_ef; %nombre de modes propres � prendre en compte
gg=zeros(nbmode,donnee.npas+1);
ggd=zeros(nbmode,donnee.npas+1);
ggdd=zeros(nbmode,donnee.npas+1);
uu=zeros(donnee.nelem+1,donnee.npas+1);
uud=zeros(donnee.nelem+1,donnee.npas+1);
uudd=zeros(donnee.nelem+1,donnee.npas+1);


%%Nouvelle matrice des vecteurs propres (tenant compte d'un nombre r�duit
%%de mode propre pris en compte)
mat_vp=ModePropre.Matrice(:,1:nbmode);

%%Traitement du chargement � l'extremit
%Initialisation
F=chargement.F(donnee.nelem,:);

UFs=zeros(donnee.nelem+1,donnee.npas+1);
UFsd=zeros(donnee.nelem+1,donnee.npas+1);
UFsdd=zeros(donnee.nelem+1,donnee.npas+1);

Fd(1)=0;
Fdd(1)=0;
for i=1:donnee.npas
    %Calcul de la d�riv�e en temps du chargement
    Fd(i+1)=(F(i+1)-F(i))/donnee.dt;
    %Calcul de la d�riv�e seconde en temps du chargement
    Fdd(i+1)=(Fd(i+1)-Fd(i))/donnee.dt;
end

%Calcul du terme de prise en charge de la charge statique
for i=1:donnee.npas
    UFs(:,i)=F(i).*SolutionStatique.U';
    UFsd(:,i)=Fd(i).*SolutionStatique.U';
    UFsdd(:,i)=Fdd(i).*SolutionStatique.U';
end


%%r�solution de l'�quation diff�rentielle
%%dai/dtdt+wi�*dai/dt=rho*S*dsig/dtdt*int(phi*us)dx


%calcul de l'integrale int(phi*us)dx � l'aide de la m�thode des trap�zes
for i=1:nbmode
    int=-donnee.dx/2*(mat_vp(1,i)*SolutionStatique.U(1)+mat_vp(donnee.nelem,i)*SolutionStatique.U(donnee.nelem))+donnee.dx*dot(mat_vp(:,i),SolutionStatique.U);
end

switch option.type
    case 'euler AR'
    %%R�solution de l'�quation diff�rentielle en temps par usage d'un sch�ma
    %%Euler arri�re (implicite
    %Conditions initiales
    gg(:,1)=zeros(nbmode,1);
    gg(:,2)=zeros(nbmode,1);
    ggd(:,1)=zeros(nbmode,1);
    ggd(:,2)=zeros(nbmode,1);
    ggdd(:,1)=zeros(nbmode,1);
        for i=1:(donnee.npas-1)
            for j=1:nbmode
                %Calcul des coefficients (deplacement)
                gg(j,i+2)=(2*gg(j,i+1)-gg(j,i)+F(i+2)-2*F(i+1)-F(i))/(1+ModePropre.Valeur(j)*donnee.dt^2);
                %Calcul des d�riv�es (vitesse)
                ggd(j,i+2)=(gg(j,i+2)-gg(j,i+1))/donnee.dt;
                %Calcul des d�riv�es seconde (acc�l�ration)
                ggdd(j,i+2)=(ggd(j,i+2)-ggd(j,i+1))/donnee.dt;
            end
        end

end



%%calcul du deplacement, de la vitesse et de l'acc�l�ration
for i=1:donnee.npas+1
   uu(:,i)=mat_vp*gg(:,i);
   uud(:,i)=mat_vp*ggd(:,i);
   uudd(:,i)=mat_vp*ggdd(:,i);
end

%initialisation
U=zeros(donnee.nelem+1,donnee.npas+1);
V=zeros(donnee.nelem+1,donnee.npas+1);
A=zeros(donnee.nelem+1,donnee.npas+1);

U=uu+UFs;
V=uud+UFsd;
A=uudd+UFsdd;
%sauvegarde des r�sultats
toto1.U=U;
toto1.V=V;
toto1.A=A;


toto2=uu/donnee.mat.L;