function [toto1,toto2]=Resolution_EFSM(chargement,donnee,ModePropre,SolutionStatique,option)

        %%Résolution sur la base modale + mode statique

        %%initialisation
        nbmode=option.Nb_ef; %nombre de modes propres à prendre en compte
        gg=zeros(nbmode,donnee.npas+1);
        ggd=zeros(nbmode,donnee.npas+1);
        ggdd=zeros(nbmode,donnee.npas+1);
        uu=zeros(donnee.nelem+1,donnee.npas+1);
        uud=zeros(donnee.nelem+1,donnee.npas+1);
        uudd=zeros(donnee.nelem+1,donnee.npas+1);


        %%Nouvelle matrice des vecteurs propres (tenant compte d'un nombre réduit
        %%de mode propre pris en compte)
        mat_vp=zeros(donnee.nelem+1,nbmode);
        mat_vp(2:donnee.nelem+1,:)=ModePropre.Matrice(:,1:nbmode);
        %mat_vp=ModePropre.Matrice(:,1:nbmode);


        %%Traitement du chargement à l'extremit
        %Initialisation
        F=chargement.F(donnee.nelem,:);

        UFs=zeros(donnee.nelem+1,donnee.npas+1);
        UFsd=zeros(donnee.nelem+1,donnee.npas+1);
        UFsdd=zeros(donnee.nelem+1,donnee.npas+1);

        Fd(1)=0;
        Fdd(1)=0;
        for i=1:donnee.npas
            %Calcul de la dérivée en temps du chargement
            Fd(i+1)=(F(i+1)-F(i))/donnee.dt;
            %Calcul de la dérivée seconde en temps du chargement
            Fdd(i+1)=(Fd(i+1)-Fd(i))/donnee.dt;
        end

        %Calcul du terme de prise en charge de la charge statique
        for i=1:(donnee.npas+1)
            UFs(:,i)=F(i).*SolutionStatique.U';
            UFsd(:,i)=Fd(i).*SolutionStatique.U';
            UFsdd(:,i)=Fdd(i).*SolutionStatique.U';
        end


        %%résolution de l'équation différentielle
        %%dai/dtdt+wi²*dai/dt=rho*S*dsig/dtdt*int(phi*us)dx/int(rho*S*phi*phi)dx


        %calcul de l'integrale int(phi*us)dx à l'aide de la méthode des trapèzes
        for i=1:nbmode
            int(i)=0;
            for j=1:donnee.nelem
                int(i)=int(i)+donnee.dx/2*donnee.mat.rho*donnee.mat.S*(mat_vp(j,i)*SolutionStatique.U(j)+mat_vp(j+1,i)*SolutionStatique.U(j+1));
            end
        end

        %calcul de l'integrale int(rho*S*phi*phi)dx à l'aide de la méthode des trapèzes
        for i=1:nbmode
            intt(i)=0;
            for j=1:donnee.nelem
                intt(i)=intt(i)+donnee.dx/2*donnee.mat.rho*donnee.mat.S*(mat_vp(j,i)^2+mat_vp(j+1,i)^2);
            end
        end
        %sintt
        switch option.schema
            case 'euler AR'
            %%Résolution de l'équation différentielle en temps par usage d'un schéma
            %%Euler arrière (implicite)
            %Conditions initiales
            gg(:,1)=zeros(nbmode,1);
            gg(:,2)=zeros(nbmode,1);
            ggd(:,1)=zeros(nbmode,1);
                for i=1:(donnee.npas-1)
                    for j=1:nbmode
                        %Calcul des coefficients (deplacement)
                        gg(j,i+2)=(2*gg(j,i+1)-gg(j,i)+(F(i+2)-2*F(i+1)+F(i))*int(j)/intt(j))/(1+ModePropre.Valeur(j)*donnee.dt^2);
                    end
                end



                for i=1:(donnee.npas)
                   %Calcul des dérivées (vitesse)
                   ggd(:,i+1)=(gg(:,i+1)-gg(:,i))/donnee.dt;
                   %Calcul des dérivées seconde (accélération)
                   ggdd(:,i+1)=(ggd(:,i+1)-ggd(:,i))/donnee.dt;
                end
            case 'euler AV'
                %%Résolution de l'équation différentielle en temps par usage d'un schéma
            %%Euler avant (explicite)
            %Conditions initiales
            gg(:,1)=zeros(nbmode,1);
            gg(:,2)=zeros(nbmode,1);
            gg(:,3)=zeros(nbmode,1);
            ggd(:,1)=zeros(nbmode,1);
                for i=1:(donnee.npas-2)
                    for j=1:nbmode
                        %Calcul des coefficients (deplacement)
                        gg(j,i+3)=gg(j,i+2)*(2-donnee.dt^2*ModePropre.Valeur(j))-gg(j,i+1)+(F(i+2)-2*F(i+1)+F(i))*int(j)/intt(j);
                    end
                end



                for i=1:(donnee.npas)
                   %Calcul des dérivées (vitesse)
                   ggd(:,i+1)=(gg(:,i+1)-gg(:,i))/donnee.dt;
                   %Calcul des dérivées seconde (accélération)
                   ggdd(:,i+1)=(ggd(:,i+1)-ggd(:,i))/donnee.dt;
                end

        end


        %%calcul du deplacement, de la vitesse et de l'accélération
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

    
        %sauvegarde des résultats
        toto1.U=U;
        toto1.V=V;
        toto1.A=A;

    toto2=toto1.U/donnee.mat.L;