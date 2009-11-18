function Affichage(Champs,donnee,option)

fenetre=figure;

switch option.type
    case 'animation en fonction du temps'
       minx=donnee.x(1);
       maxx=donnee.x(donnee.nelem+1);
       chpsmin=min(Champs.U(:));
       chpsmax=max(Champs.U(:));
       %création du repertoire de sauvegarde
       if strcmp(option.save,'film')|strcmp(option.save,'images')
                named=([option.dossier,'_el',num2str(donnee.nelem),'_npas',num2str(donnee.npas),'_mode',num2str(option.nbm)]);
                mkdir(named);
       end
       
        for i=1:donnee.npas
            
            plot(donnee.x,Champs.U(:,i));
            axis([minx maxx chpsmin chpsmax]);
            title(option.titre);
            xlabel(sprintf('x'));
            ylabel('amplitude');
            F(i)=getframe;
            %sauvegarde si souhaité
            if strcmp(option.save,'images')
                namef=([named,'/','dep_',num2str(i)]);
                saveas(gcf, namef , 'jpg');
            end
        end
        %affichage de l'animation 
        if strcmp(option.save,'film')
            namev=([named,'/anim']);
            namev
            movie(fenetre,F)
            movie2avi(F,namev,'compression','Indeo3','fps',15)
        end
        
        
	case 'en fonction du temps'
		if size(Champs,2)==donnee.npas+1
			plot(donnee.t,Champs);
		else
			plot(donnee.t(1:donnee.npas),Champs);
		end
%  		axis([0 donnee.mat.L -0.0001 0.0001]);
		title(option.titre);
		xlabel(sprintf('t'));
		ylabel('amplitude');
        
	case 'pasapas'
		for i=1:donnee.npas+1
			if size(Champs,1)==donnee.nelem+1
			plot(donnee.mat.L*(0:1/donnee.nelem:1),Champs(:,i));
			else
			plot(donnee.mat.L*(1/donnee.nelem:1/donnee.nelem:1),Champs(:,i));
			end
%  			axis([0 donnee.mat.L -0.001 0.001]);
			title(option.titre);
			xlabel(sprintf('x (au temps t=%ds)',(i-1)*donnee.T/donnee.npas));
			ylabel('amplitude');
			pause(1/10);
		end
	case '3D'
		if size(Champs,1)==donnee.nelem+1
			%champs exprime au noeuds, interpolation lineaire entre chaque noeud
			surf(donnee.t,donnee.x,Champs);
		else
			%champs exprime par elements (constant par element)
			x_dbl(1:2:2*donnee.nelem-1)	=donnee.x(1:donnee.nelem);
			x_dbl(2:2:2*donnee.nelem)	=donnee.x(2:donnee.nelem+1);
			Champs_dbl(1:2:2*donnee.nelem-1,:)	=Champs;
			Champs_dbl(2:2:2*donnee.nelem,:)	=Champs;
			surf(donnee.t,x_dbl,Champs_dbl);
%  			size(x_dbl);
%  			size(donnee.t);
%  			size(Champs_dbl);
		end
		title(option.titre);
		xlabel('t');
		ylabel('x');
		zlabel('amplitude');
		shading interp
	case 'Un Mode Propre'
		plot(donnee.x,Champs.Vecteur{option.Mode});
		xlabel('x');
		ylabel('Amplitude')
		title(option.titre);
	case '4 Modes Propres'
		title(option.titre)
		subplot 221
			plot(donnee.x,Champs.Vecteur{option.Mode(1)});
			xlabel('x');
			ylabel('Amplitude')
			title(sprintf('Mode propre numero %d',option.Mode(1)))
		subplot 222
			plot(donnee.x,Champs.Vecteur{option.Mode(2)});
			xlabel('x');
			ylabel('Amplitude')
			title(sprintf('Mode propre numero %d',option.Mode(2)))
		subplot 223
			plot(donnee.x,Champs.Vecteur{option.Mode(3)});
			xlabel('x');
			ylabel('Amplitude')
			title(sprintf('Mode propre numero %d',option.Mode(3)))
		subplot 224
			plot(donnee.x,Champs.Vecteur{option.Mode(4)});
			xlabel('x');
			ylabel('Amplitude')
			title(sprintf('Mode propre numero %d',option.Mode(4)));
	case 'Tous les Modes Propres'
		title(option.titre);
			N=Champs.n;
			for i=1:floor(N/4)
				option.type='4 Modes Propres';
				option.Mode=(i-1)*4+[1 2 3 4];
				option.titre=sprintf('Mode propre numero %d',option.Mode);
				Affichage(Champs,donnee,option);
			end
		figure
		title(option.titre)
		if N>floor(N/4)*4;
			subplot 221
				plot(donnee.x,Champs.Vecteur{floor(N/4)*4+1});
				xlabel('x');
				ylabel('Amplitude')
				title(sprintf('Mode propre numero %d',floor(N/4)*4+1));
			if N>floor(N/4)*4+1;
				subplot 222
					plot(donnee.x,Champs.Vecteur{floor(N/4)*4+2});
					xlabel('x');
					ylabel('Amplitude')
					title(sprintf('Mode propre numero %d',floor(N/4)*4+2));
				if N>floor(N/4)*4+2;
					subplot 223
						plot(donnee.x,Champs.Vecteur{floor(N/4)*4+3});
						xlabel('x');
						ylabel('Amplitude')
						title(sprintf('Mode propre numero %d',floor(N/4)*4+3));
				end
			end
		end
	end
	
	
end