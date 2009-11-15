function [toto1,toto2]=Resolution_EF(chargement,donnee,ModePropre,matrice,mat_vp)
%%initialisation
nbmode=ModePropre.Nb_ef; %nombre de modes propres à prendre en compte
gg=zeros(nbmode,donnee.npas+1);
ggd=zeros(nbmode,donnee.npas+1);
uu=zeros(donnee.nelem+1,donnee.npas+1);


%%redimensionnement matrice des vecteurs propres
    mat_vp=ModePropre.Matrice(:,1:nbmode);
%mat_vp

%%Détermination des coefficients de l'équation différentielle
    MM=mat_vp'*matrice.M*mat_vp;
    KK=mat_vp'*matrice.K_ef*mat_vp;
    F=zeros(donnee.nelem+1,donnee.npas+1);
    F(2:donnee.nelem+1,:)=chargement.F(:,:);

    for i=1:donnee.npas+1
    FF(:,i)=mat_vp'*F(:,i);
    end
    ed.WW=inv(MM)*KK;
   
    for i=1:nbmode
        ed.W(i)=sqrt(ed.WW(i,i));        
    end
    
%%on cherche u au noeud de la forme ui=Ai*cos(wi*t)+Bi*sin(wi*t)+Ci
%%on distingue pour cela 3 phases de calcul
%%%%%%  si 0<t<T1 F(t)=0  donc u=0
%%%%%%  si T1<t<T2 F(t)=F=cste  donc u varie
%%%%%%  si T2<t  F(t)=0 donc u varie de nouveau


%on cherche les points les plus pres des abscisses adimensionnees
		pas1=find((abs(donnee.t-donnee.T*chargement.parametre{2}(1))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(1)))));
		pas2=find((abs(donnee.t-donnee.T*chargement.parametre{2}(2))) == min((abs(donnee.t-donnee.T*chargement.parametre{2}(2)))));

%%%%%T1<t<T2
        for i=1:nbmode
            ed.C(i)=FF(i,pas1)/ed.W(i);
            ed.A(i)=-ed.C(i)*cos(ed.W(i)*pas1);
            ed.B(i)=-ed.C(i)*sin(ed.W(i)*pas1);
        end
        
        for j=pas1:pas2
            gg(:,j)=fct_coef(ed,j);
            ggd(:,j)=fct_coef_d(ed,j);
        end
%%%%%%T2<t
        for i=1:nbmode
            ed.C(i)=0;
            ed.A(i)=(sin(ed.W(i)*pas2)*ggd(i,pas2)-ed.W(i)*cos(ed.W(i)*pas2)*gg(i,pas2))/ed.W(i);
            ed.B(i)=(cos(ed.W(i)*pas2)*ggd(i,pas2)+ed.W(i)*sin(ed.W(i)*pas2)*gg(i,pas2))/ed.W(i);
        end
        for j=pas2:donnee.T
            gg(:,j)=fct_coef(ed,j);
            ggd(:,j)=fct_coef_d(ed,j);
        end
       
%%calcul du deplacement
for i=1:donnee.npas
   uu(:,i)=mat_vp*gg(:,i);
end
toto1=uu;
toto2=uu/donnee.mat.L;