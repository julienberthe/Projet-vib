function [toto1,toto2]=Resolution_EFD(chargement,donnee,matrice,option)


%%R�solution directe matricielle + sch�ma it�ratif en temps
switch option.schema

    case 'euler AR'  %sch�ma implicite

        %%cr�ation des matrices
        %initialisation
        MM=zeros(2*donnee.nelem);
        KK=zeros(2*donnee.nelem);
        F=zeros(2*donnee.nelem,donnee.npas+1);
        q=zeros(2*donnee.nelem,donnee.npas+1);
        %cr�ation
        MM(1:donnee.nelem,1:donnee.nelem)=eye(donnee.nelem,donnee.nelem);
        MM(donnee.nelem+1:2*donnee.nelem,donnee.nelem+1:2*donnee.nelem)=matrice.M;
        KK(1:donnee.nelem,donnee.nelem+1:2*donnee.nelem)=-eye(donnee.nelem,donnee.nelem);
        KK(donnee.nelem+1:2*donnee.nelem,1:donnee.nelem)=matrice.K_ef;
        F(donnee.nelem+1:2*donnee.nelem,:)=chargement.F;
        %Calcul pr�liminaire
        C=(eye(2*donnee.nelem,2*donnee.nelem)+(MM^-1*KK)*donnee.dt)^-1;

        %%�tablissement du sch�ma it�ratif
        for i=1:donnee.npas
            q(:,i+1)=C*(q(:,i)+donnee.dt*MM^-1*F(:,i+1));
        end
        U=q(1:donnee.nelem,:);
        V=q(donnee.nelem+1:2*donnee.nelem,:);
        A=0;




    case 'euler AV'  %sch�ma explicite
        %non impl�ment�

    case 'newmark'
        %non impl�ment�
end
toto1.U=zeros(donnee.nelem+1,donnee.npas+1);
toto1.V=zeros(donnee.nelem+1,donnee.npas+1);
%tot1.A=zeros(1:donnee.nelem+1,donnee.npas);
            %sauvegarde des r�sultats
    toto1.U(2:donnee.nelem+1,:)=U;
    toto1.V(2:donnee.nelem+1,:)=V;
    toto1.A=A;

    toto2=toto1.U/donnee.mat.L;