% close all;
% clear all;

fid = fopen('data_mgd2DKSE','r');
    output = textscan(fid,'%*s %f %*s %f %*s %f');
fclose(fid);

col1 = output{1}; col2 = output{2}; col3 = output{3};

M = col1(1); N = col2(1); dt = col3(1);
Tfinal = col1(2);
kappa = col1(3); gamma = col2(3); s = col3(3);
sizeLvec = col1(4); L_1 = col2(4); stepL = col3(4);
a = col1(5);

L(1:sizeLvec)=0.0;
for i = 1:sizeLvec
    L(i) = L_1 + round(i-1)*stepL;
end

string_M = int2str(M);
string_N = int2str(N);

if M < 10
	string_M = ['000', string_M];
elseif M < 100
	string_M = ['00', string_M];
elseif M < 1000
	string_M = ['0', string_M];
elseif M < 10000
    string_M = string_M;
end

if N < 10
	string_N = ['000', string_N];
elseif N < 100
	string_N = ['00', string_N];
elseif N < 1000
	string_N = ['0', string_N];
elseif N < 10000
    string_N = string_N;
end

for j1 = 1:sizeLvec
    
                  
Lint = floor(L(j1)); Lfrac = L(j1) - Lint;
aint = floor(a); afrac = a - aint;
kappaint = floor(abs(kappa)); kappafrac = abs(kappa) - kappaint;

Lint = num2str(Lint); aint = num2str(aint); kappaint = num2str(kappaint);

  Lfrac1 = floor(Lfrac*10);
  Lfrac2 = floor(Lfrac*100)-Lfrac1*10;
  Lfrac3 = floor(Lfrac*1000)-Lfrac1*100 - Lfrac2*10;
  
  afrac1 = floor(afrac*10);
  afrac2 = floor(afrac*100)-afrac1*10;
  afrac3 = floor(afrac*1000)-afrac1*100 - afrac2*10;
  
  kappafrac1 = floor(kappafrac*10);
  kappafrac2 = floor(kappafrac*100)-kappafrac1*10;
  kappafrac3 = floor(kappafrac*1000)-kappafrac1*100 - kappafrac2*10;

 Lfrac = [num2str(Lfrac1),num2str(Lfrac2),num2str(Lfrac3)];  

 afrac = [num2str(afrac1),num2str(afrac2),num2str(afrac3)];  

 kappafrac = [num2str(kappafrac1),num2str(kappafrac2),num2str(kappafrac3)];  

         % Load files:
         
         x_k1 = load(['x_k1_M=', string_M, '_N=', string_N, '.txt']);      
         y_k2 = load(['y_k2_M=', string_M, '_N=', string_N, '.txt']);  
         L2norminp = load(['L2norm_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_kappa=',kappaint,'-',kappafrac,'.txt']);        
         Linfnorminp = load(['Linfnorm_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_kappa=',kappaint,'-',kappafrac,'.txt']);
         profiles = load(['profiles_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_kappa=',kappaint,'-',kappafrac,'.txt']);     
         
         x = x_k1(:,1); y = y_k2(:,1); k1 = x_k1(:,2); k2 = y_k2(:,2);
         uinitial = profiles(1:2*M,:); u = profiles(2*M+1:4*M,:); Fu_final = profiles(4*M+1:6*M,:);
         tdata = L2norminp(:,1); L2norm = L2norminp(:,2); Linfnorm = Linfnorminp(:,2);

 
         L1 = L(j1); L2 = L1^a ;
         
        % Plot of solution initial condition:        
        figure; 
        subplot(2,3,1)
        surfc(y*L2/L(end)^a,x*L1/L(end),uinitial);
        shading interp
        xlabel('$y$','Interpreter','LaTex','Fontsize',24); xlim([0 L2]);
        ylabel('$x$','Interpreter','LaTex','Fontsize',24); ylim([0 L1]);
        zlabel('$u_0$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'ZLabel'),'Rotation',0);
        
        % Plot of solution u:        
        %figure; 
        subplot(2,3,2)
        surfc(y*L2/L(end)^a,x*L1/L(end),u);
        shading interp
        xlabel('$y$','Interpreter','LaTex','Fontsize',24); xlim([0 L2]);
        ylabel('$x$','Interpreter','LaTex','Fontsize',24); ylim([0 L1]);
        zlabel('$u$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'ZLabel'),'Rotation',0);
        
        % Plot final spectrum Fu_final
        
        %figure; 
        subplot(2,3,3)
        surf(k2,k1,log10(Fu_final)); colorbar
        shading interp
          
        % Plot of L2norm squared:
        %figure; 
        subplot(2,3,4)
        plot(tdata,L2norm,'-k');
        xlabel('t','Interpreter','LaTex','Fontsize',24);
        ylabel('$||u||_{2}^2$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'YLabel'),'Rotation',0);
      
        % Plot of Linfnorm:
        %figure; 
        subplot(2,3,5)
        plot(tdata,Linfnorm,'-k');
        xlabel('t','Interpreter','LaTex','Fontsize',24);
        ylabel('$||u||_{\infty}$','Interpreter','LaTex','Fontsize',24);
        set(get(gca,'YLabel'),'Rotation',0);
     
               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculate norms of unstable modes
        clear unstableabssq
        unstable = load(['unstable_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_kappa=',kappaint,'-',kappafrac,'.txt']);  
        control = load(['control_L=',Lint,'-',Lfrac,'_a=',aint,'-',afrac,'_kappa=',kappaint,'-',kappafrac,'.txt']);  
        
       
        unstableL2 = sqrt(sum(unstable.^2,2));
        
        figure(20001)
        hold on
        plot3 = plot(tdata,L2norm,'-k','LineWidth',2);
        set(plot3(1),'DisplayName','$L_0^2$-norm of $\eta^*$',...
    'Color',[0 0 0]);
        plot4 = plot(tdata(1:40:end-1,1),unstableL2(1:40:end),'--k','LineWidth',2);
        set(plot4(1),...
    'DisplayName','$L_0^2$-norm of $\tilde{\eta}_{\textrm{u}}^*$',...
    'MarkerSize',0.5,...
    'Color',[0 0 0]);
        xlabel('$t$','Interpreter','LaTex','Fontsize',24);
        ylabel('$L_0^2$-norm','Interpreter','LaTex','Fontsize',24);
        
     
        unstable = unstable(:,3:end);
        unstableabssq = unstable(:,1:2:end).^2 + unstable(:,2:2:end).^2;
        countunstable = size(unstableabssq,2); unstableHssq=0;
        unstableL2redo = sqrt(2*sum(unstableabssq,2));
        for jkl = 1:countunstable
            unstableHssq = unstableHssq + 2*(jkl*pi/L2)^(2*s)*unstableabssq(:,jkl);
        end
        unstableHs = sqrt(unstableHssq);
        controlL2 = sqrt(sum(control.^2,2));
        
        
        figure(20003)
        hold on
        plot2 = plot(tdata(1:end-1,1),unstableHs,'-k','LineWidth',2);
        set(plot2(1),...
    'DisplayName','$H_0^s$-norm of $\tilde{\eta}_{\textrm{u}}^*$',...
    'Color',[0 0 0]);
        plot3 = plot(tdata(1:40:end-1,1),controlL2(1:40:end,1),'--b','LineWidth',2);
        set(plot3(1),'DisplayName','$L_0^2$-norm of $\tilde{\zeta}^*$','LineStyle','--',...
    'Color',[0 0 1]);
        
        
        figure(20005)
        surf(y*L2/L(end)^a,x*L1/L(end),u);
        shading interp
        set(gca,'FontSize',16,'FontWeight','bold');
        xlabel('$y$','Interpreter','LaTex','Fontsize',24); xlim([0 L2]);
        ylabel('$x$','Interpreter','LaTex','Fontsize',24); ylim([0 L1]);
        zlabel('$\eta^*(T)$','Interpreter','LaTex','Fontsize',24);
        %set(get(gca,'ZLabel'),'Rotation',0);
        set(gca,'YDir','reverse');
        view([-38,74]);
        
        
        % Calculate cost functional value !!! timestep for plot is 0.025
        Cost = 0.005*sum(unstableHs(1:end-1).^2)/2 + unstableHs(end).^2/2 + gamma*0.005*sum(controlL2(1:end-1,1).^2)/2;
       
        
        
   
end
