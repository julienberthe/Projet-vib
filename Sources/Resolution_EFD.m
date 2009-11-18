function [toto1,toto2]=Resolution_EFD(chargement,donnee,matrice,option)


%%Résolution directe matricielle + schéma itératif en temps
switch option.schema

    case 'euler AR'  %schéma implicite

        %%création des matrices
        %initialisation
        MM=zeros(2*donnee.nelem);
        KK=zeros(2*donnee.nelem);
        F=zeros(2*donnee.nelem,donnee.npas+1);
        q=zeros(2*donnee.nelem,donnee.npas+1);
        %création
        MM(1:donnee.nelem,1:donnee.nelem)=eye(donnee.nelem,donnee.nelem);
        MM(donnee.nelem+1:2*donnee.nelem,donnee.nelem+1:2*donnee.nelem)=matrice.M;
        KK(1:donnee.nelem,donnee.nelem+1:2*donnee.nelem)=-eye(donnee.nelem,donnee.nelem);
        KK(donnee.nelem+1:2*donnee.nelem,1:donnee.nelem)=matrice.K_ef;
        F(donnee.nelem+1:2*donnee.nelem,:)=chargement.F;
        %Calcul préliminaire
        C=(MM^-1*KK)/donnee.dt;

        %%établissement du schéma itératif
        for i=1:donnee.npas
            q(:,i+1)=q(:,i)-C*q(:,i)+MM^-1*F(:,i);
        end
        U=q(1:donnee.nelem,:);
        V=q(donnee.nelem+1:2*donnee.nelem,:);
        A=0;




    case 'euler AV'  %schéma explicite
        %non implémenté

    case 'newmark'
        %non implémenté
end
toto1.U=zeros(1:donnee.nelem+1,donnee.npas+1);
toto1.V=zeros(1:donnee.nelem+1,donnee.npas+1);
%tot1.A=zeros(1:donnee.nelem+1,donnee.npas);
            %sauvegarde des résultats
    toto1.U(2:donnee.nelem+1,:)=U;
    toto1.V(2:donnee.nelem+1,:)=V;
    toto1.A=A;

    toto2=toto1.U/donnee.mat.L;