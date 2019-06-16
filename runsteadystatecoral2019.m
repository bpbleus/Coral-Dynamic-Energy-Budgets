
% DYNAMIC EQUATIONS OF CORAL MODEL + CALCIFICATION + PARTIAL OA EFFECTS

% NOTE THAT OA EFFECTS ARE ONLY ON CALCIFICATION
% OA EFFECTS ON CO2 ENRICHMENT (FOR PHOTOSYNTHESIS) AND EFFECT
% OF H+ TOXICITY (IMPACT OF PH ON GENERAL PHYSIOLOGY) ARE EXCLUDED

% NOTE THAT ENVIRONMENTAL RESOURCES (DIN, DIC,IRRADIANCE) ARE CONSTANT
% THIS AIDS UNDERSTANDING THE DYNAMICS OF HOST/SYMBIONT BUT IT IS 'WRONG'
% (the assumption implied is that some entity (Neptune?) performs miracles 
% by immediately replacing a nutrient molecule whenever 1 was taken up) 

% Notation closely follows paper; 
% 'X' are concentration measures, i.e. mass per system volume
% 1st letter refers to symbol; 2nd and beyond are mostly subscripts

% ALL MOLAR QUANTITIES ARE IN muMOL; TIME UNIT IS DAY
close all


% ENVIRONMENTAL PARAMETERS AND CONDITIONS
% to delete Jlf [(4:10:84)*10e6];%[(0.01:0.02:0.89)*10e7];
Xc      =   2000; 	% ambient inorganic C density (min of range)
XcU     =   Xc; 	% ambient inorganic C density (max of range)
interv  =   25;      % number of env. conditions in range
Xn      =   0.1;	% ambient inorganic N density (min of range)
XnU     =   4;	% ambient inorganic N density (max of range)
X1      =   0;      % ambient prey density   (min of range)
X1U     =   0;      % ambient prey density (max of range)
Jlf     =   2e6;    % flux of photons; 43.2e6: 500 mumol q/m2.s
JlfU    =   6e7;    % flux of photons; 43.2e6: 500 mumol q/m2.s
%hrc     =   0.05;   % not used, hrc is solved.
h       =   1;      % dilution rate of solubles 
frJlfD  =   1;      % fraction of Jlf potentially causing damage
                    
% vector of environmental parameters and conditions

envp = [Xc Xn X1 Jlf XcU XnU X1U JlfU h interv frJlfD];

% DEFAULT PARAMETERS
% 1. YIELDS & ELEMENTAL COEFFICIENTS
nhgly    =  2.67;   % H:C ratio in glycerol
nogly    =  1.0;    % O:C ratio in in glycerol
nngly    =  0.0;    % N:C ratio in in glycerol
nhx1    =   1.8;    % H:C ratio in prey
nox1    =   0.4;    % O:C ratio in prey
nnx1    =   0.2;    % N:C ratio in prey
nhp1    =   1.8;    % H:C ratio in feces from prey
nop1    =   0.4;    % O:C ratio in feces from prey
nnp1    =   0.05;   % N:C ratio in feces from prey
nhes    =   1.8;    % H:C ratio in symbiont reserves
noes    =   0.4;    % O:C ratio in symbiont reserves
nnes    =   0.15;   % N:C ratio in symbiont reserves
nhvs    =   1.8;    % H:C ratio in symbiont structure
novs    =   0.4;    % O:C ratio in symbiont structure
nnvs    =   0.15;   % N:C ratio in symbiont structure
nheh    =   1.8;    % H:C ratio in host reserves
noeh    =   0.4;    % O:C ratio in host reserves
nneh    =   0.2;    % N:C ratio in host reserves
nhvh    =   1.8;    % H:C ratio in host structure
novh    =   0.4;    % O:C ratio in host structure
nnvh    =   0.2;    % N:C ratio in host structure

yg3pl   =   0.05;    % yield of g3p from photons, tentative value
yesg3p  =   0.8;    % yield of symbiont reserves from g3p
yvses   =   0.8;    % yield of symbiont structure from reserves
yehx1   =   0.4;    % yield of host reserves from prey
yp1x1   =   0.5;    % yield of feces from prey
yehgly  =   0.8;    % yield of host reserves from glycerol
yvheh   =   0.8;    % yield of host structure from host reserves
ymeh    =   0;      % yield of mucus from host reserves; 
                    % 0 switches off mucus production
yrheh   =   1;      % yield of reproductive matter from host reserves

% vector of standard stoichiometric parameters 

stoichp = [nhgly   nogly  nngly   nhx1   nox1    nnx1   nhp1  nop1 ...
           nnp1   nhes    noes    nnes    nhvs   novs    nnvs ...    
           nheh   noeh    nneh    nhvh   novh   nnvh ...    
        yg3pl   yesg3p  yvses  yehx1  yp1x1   yehgly yvheh ymeh   yrheh];
    
% DEFAULT PARAMETERS
% 2. DEB PARAMETERS (note: fluxes with large value do not affect dynamics)     
jg3pasm =   6;      % max specific g3p synthesis rate
jIp     =   50e6;   % photoinhibition scaling parameter (production phase)
jIb     =	50e6; 	% photoinhibition scaling parameter (binding phase)
switchIp=   0;      % on/off switch photoinhibition (production phase)
switchIb=   0;      % on/off switch photoinhibition (binding phase)
 % if 0 < switch parameter < 1, interpretation is fraction of Jlf 
 % potentially causing photoinhibition
jIDs 	=	1e10;      % photodamage scaling parameter symbiont
jIDh  	=	1e10;      % photodamage scaling parameter host
Jlfsnec =   1e10;      % no effect irradiation symbiont
Jlfhnec =   1e10;      % no effect irradiation host
jesasm  =   1e8;    % max specific g3p assimilation rate      
jesds   =   0.2;    % specific maintenance rate in symbiont
jxm     =   1;      % max spec feeding rate of host
jcfhm   =   6e5;      % max spec inorganic C uptake rate of host
jnfhm   =   0.1;    % max spec ammonium uptake rate of host
jehdh   =   0.2;    % spec maintance rate of host
jeha3hm =   1e8;    % max spec glycerol assimilation rate of host
kes     =   1.2;    % storage turn over rate in symbiont
rhoas   =   1;      % probability that reserve SU binds G3P
Xkc     =   400;  	% saturation constant for carbonate uptake; 
Xkn     =   1;      % saturation constant for ammonia uptake;
Xk1     =   1;      % saturation constant for catching prey    
keh     =   0.1;    % storage turn over rate in host
kaph    =   1;      % fraction of commitment flux to maturation/reproduction
ds      =   1/(5e-6);% phi_s: symbiont structure/symbiont cross sectional area 

% vector of DEB and effect parameters (except OA)
debp =[jg3pasm jIp jIb switchIp switchIb jesasm jesds jxm ...
        jcfhm jnfhm jehdh jeha3hm kes rhoas Xkc Xkn Xk1 keh kaph ds ...
        jIDs jIDh Jlfsnec Jlfhnec];
    
% CALCIFICATION and OA PARAMETERS

ycaas =     0.00;   % yield of CaCO3 due to assimilation symbiont
ycaah =     0.01;   % yield of CaCO3 due to assimilation host
ycagh =     0.01;   % yield of CaCO3 due to growth host
ycadh =     0.03;   % yield of CaCO3 due to dissipation (maintenance) host
omegaenv =  4;      % aragonite saturation state
komega =    0.1;    % OA effect scaling parameter, muM^-1
jcamin =    0;      % Dissolution rate (default 0)
% add H+ toxicity

% vector of calcification and OA parameters

calcp = [ycaas ycaah ycagh ycadh omegaenv komega jcamin];

% START VALUES                    
mes0    =   0.5;    % initial guess symbiont reserve density
hrc0    =   0.01;   % initial guess specific growth rate
Xvs     =   1;      % implied, arbitrary  
meh0    =   0.5;  	% initial guess host reserve density
Xvh0    =   1;      % initial guesshost structure density

% vector of integration parameters

initp = [mes0 hrc0 Xvs meh0 Xvh0];

%--------------------- CALCULATIONS-------------------------------------
%--------------------- START HERE --------------------------------------

[Y,fluxes,remark] =steadystatecoral(stoichp,debp,calcp,initp,envp);

%----- PREPARE RESULTS FOR PLOTTING ----------------------
% Not all calculations are used for present plots;
% I keep them in so it is easy to add more results to the plots
% Change to 'subplot' for easier 'play'; removed it to get plots to publish

mes=squeeze(Y(1,:,:));hrc=squeeze(Y(2,:,:));meh=squeeze(Y(3,:,:));       
Xvh=squeeze(Y(4,:,:));Cin=squeeze(Y(5,:,:));Nin=squeeze(Y(6,:,:));

biomasss=Xvs.*mes+Xvs;
biomassh=Xvh.*meh+Xvh;
ratiostructures=Xvs./Xvh;

jg3pas=squeeze(fluxes(1,:,:));jesas=squeeze(fluxes(2,:,:));
jeha3h=squeeze(fluxes(3,:,:));jeha12h=squeeze(fluxes(4,:,:));
timeCs=squeeze(fluxes(5,:,:));timeLs=squeeze(fluxes(6,:,:));
timeNs=squeeze(fluxes(7,:,:));timeGLYs=squeeze(fluxes(8,:,:));
timeNh=squeeze(fluxes(9,:,:));timeGLYh=squeeze(fluxes(10,:,:));
jvsgs=squeeze(fluxes(11,:,:));jvhgh=squeeze(fluxes(12,:,:));
jglyf=squeeze(fluxes(13,:,:));joh=squeeze(fluxes(14,:,:));
joresps=squeeze(fluxes(15,:,:)); jophotos=squeeze(fluxes(16,:,:));
% calculate capacities and limitation coefficients
capphoto=jg3pas/jg3pasm;
limcoeffphoto=log10(timeLs./timeCs);
capas=jesas/jesasm;
limcoeffas=1./(timeGLYs.*jesas);
capa3h=jeha3h/jeha3hm;
limcoeffa3h=1./(timeGLYh.*jeha3h);

%calculate yields
yglyfc=jglyf.*Xvh./(jg3pas.*Xvs);
yehc=jeha3h.*Xvh./(jg3pas.*Xvs);
yesc=jesas./jg3pas;

autoCinhost=jeha3h./(jeha12h+jeha3h);

symbiomass=Xvs.*mes+Xvs; hostbiomass=Xvh.*meh+Xvh;

domainN=(Xn:(XnU-Xn)/interv:XnU);     
domainJlf=(Jlf:(JlfU-Jlf)/interv:JlfU)/(24*3600);
x=domainN; y=domainJlf;
%[x,y]=meshgrid(domainN,domainJlf);

pos=[0.75 0.12 0.05 0.8];
figure(1)
	h1=imagesc(x,y,hrc');  
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 
    h2=colorbar;
    text(1.4*max(x),0.5*max(y),'d^{-1}','Rotation',270,...
            'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
	text(0.5*max(x),0.95*max(y),'Specific growth rate',...
        'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
        'Color',[1 1 1])
 figure(2) 
  	%title('Symbiont to host ratio (\{propto} surface area)')
	h1=imagesc(x,y,ratiostructures');
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 
    h2=colorbar;
    text(1.4*max(x),0.5*max(y),'mol C/ mol C (structural biomass)','Rotation',270,...
            'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
	text(0.5*max(x),0.95*max(y),'Symbiont density (\propto surface area)',...
        'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
        'Color',[1 1 1])

figure(3) 
	%title('Symbiont density (biomass)')
    h1=imagesc(x,y,(symbiomass./(symbiomass+hostbiomass))');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 
    h2=colorbar;
    text(1.4*max(x),0.5*max(y),'mol symbiont C/ mol total C','Rotation',270,...
            'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
    text(0.5*max(x),0.95*max(y),'Symbiont biomass density',...
        'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
        'Color',[1 1 1]) 
figure(4) 
	%title('Dark Respiration (\propto surface area)')
    h1=imagesc(x,y,(0.5*(joresps.*Xvs./Xvh+joh))');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 
    h2=colorbar;
    text(1.4*max(x),0.5*max(y),'mol O_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
            'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
    text(0.5*max(x),0.95*max(y),'O_{2} consumption rate (\propto surface area)',...
        'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
        'Color',[1 1 1])       
figure(5) 
	%title('Net respiration rate (\propto surface area)')
    h1=imagesc(x,y,(0.5*((joresps-jophotos).*Xvs./Xvh+joh))');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 
    h2=colorbar;
    text(1.4*max(x),0.5*max(y),'mol O_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
            'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
    text(0.5*max(x),0.95*max(y),'Net respiration rate (\propto surface area)',...
        'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
        'Color',[1 1 1])          
figure(6) 
	%title('Photosynthesis (\propto surface area)')
    h1=imagesc(x,y,jg3pas.*Xvs./Xvh');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 
	h2=colorbar;
    text(1.3*max(x),0.5*max(y),'mol CO_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
            'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')    
    text(0.5*max(x),0.95*max(y),'Photosynthesis rate (\propto surface area)',...
        'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
        'Color',[1 1 1])   
    
 figure(61) 
	%title('Photosynthesis (\propto surface area)')
    h1=imagesc(x,y,limcoeffas');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
               'LineWidth',1.5,'FontSize',24) 

%     set(h2,'YTick',yt)
%     text(1.3*max(x),0.5*max(y),'mol CO_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
%             'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')    
%     text(0.5*max(x),0.95*max(y),'Photosynthesis rate (\propto surface area)',...
%         'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
%         'Color',[1 1 1])   
 
% figure(62) 
% 	%title('Photosynthesis (\propto surface area)')
%     h1=imagesc(x,y,limcoeffa3h');   
% 	axis xy
% 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% 	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
%         'FontSize',36,'FontWeight','demi');
%     set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
%                'LineWidth',1.5,'FontSize',24) 
% 	h2=colorbar;
%     asp = [1 10 1];
%     yt = (2:2:8);
%     set(h2,'Position',pos,'FontSize',24)
% %     set(h2,'YTick',yt)
% %     text(1.3*max(x),0.5*max(y),'mol CO_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
% %             'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')    
% %     text(0.5*max(x),0.95*max(y),'Photosynthesis rate (\propto surface area)',...
% %         'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
%         'Color',[1 1 1])   
    
    

% figure(7) 
% 	%title('Fraction of photosynthate assimilated by host')
%     h1=imagesc(x,y,yehc');   
% 	axis xy
% 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% 	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
%         'FontSize',36,'FontWeight','demi');
%     set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
%                'LineWidth',1.5,'FontSize',24) 
% 	h2=colorbar;
%     asp = [1 10 1];
%     set(h2,'Position',pos,'PlotBoxAspectRatio',asp,'FontSize',24)
%     %set(h2,'YTick',yt,'FontSize',24)
%     text(1.3*max(x),0.5*max(y),'Fraction of photosynthate assimilated by host','Rotation',270,...
%             'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')    
%     text(0.5*max(x),0.95*max(y),'Autotrophic C to host',...
%         'FontSize',24,'FontWeight','demi','HorizontalAlignment','center',...
%         'Color',[1 1 1])              

    
    
%         figure(4) 
% 	%title('Dark Respiration')
%     h1=imagesc(x,y,(0.5*(joresps.*Xvs+joh.*Xvh)./(symbiomass+hostbiomass))');   
% 	axis xy
% 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% 	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
%         'FontSize',36,'FontWeight','demi');
%     set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
%                'LineWidth',1.5,'FontSize',24) 
% 	h2=colorbar;
%     asp = [1 10 1];
%     set(h2,'Position',pos,'PlotBoxAspectRatio',asp,'YTick',yt,'FontSize',24)
%     %set(h2,'YTick',yt,'FontSize',24)
%     text(6,500,'Dark respiration rate (mol O_{2}] mol biomass C^{-1} day^{-1})','Rotation',270,...
%             'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
%         
% figure(5) 
% 	%title('Light Respiration')
%     h1=imagesc(x,y,(0.5*((joresps-jophotos).*Xvs+joh.*Xvh)./(symbiomass+hostbiomass))');   
% 	axis xy
% 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% 	ylabel('Irradiation, \mumol q. m^{-2}. s^{-1}',...
%         'FontSize',36,'FontWeight','demi');
%     set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 200 400 600],...
%                'LineWidth',1.5,'FontSize',24) 
% 	h2=colorbar;
%     asp = [1 10 1];
%     set(h2,'Position',pos,'PlotBoxAspectRatio',asp,'YTick',yt,'FontSize',24)
%     %set(h2,'YTick',yt,'FontSize',24)
%     text(6,500,'Respiration rate (mol O_{2}] mol biomass C^{-1} day^{-1})','Rotation',270,...
%             'FontSize',36,'FontWeight','demi','HorizontalAlignment','center')
        
% figure(9)
%   %title('')
%         h1=mesh(x,y,(0.5*joh.*Xvh)');
%         h2=zlabel('Host respiration in \mumol O_{2}]/ day');

% figure(10)
%   %title('')
%         h1=mesh(x,y,(0.5*symresp.*Xvh)');
%         h2=zlabel('Symbiont respiration in \mumol O_{2}]/ day');
%         h3=xlabel('DIN in \muM N');
%         h4=ylabel('Irradiation in \mumol q. m^{-2}. s^{-1}');
%         %axis([0 4 0 1000 -0.2 .6])
%         %set(gca,'ztick',[-.2 0 .2 0.4 0.6])
%      	%axis ij
%         
%         %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%          %   'LineWidth',1.5,'FontSize',24)       
%         set([h2 h3 h4],'FontSize',24,'FontWeight','demi')
%         set([h1],'LineWidth',1)
%         %annotation('textbox','String',{'B'},'LineStyle','none',...
%          %   'Position',[0.12 0.7001 0.07031 0.07537],...
%           %  'FitHeightToText','on','FontSize',24);
% figure(11)
%   title('Coral respiration per total biomass (w/o photo)')
%         h1=mesh(x,y,(0.5*(symresp.*Xvs+hostresp.*Xvh)./(symbiomass+hostbiomass))');
%         h2=zlabel('mol O_{2}]/ mol C day');
%         h3=xlabel('DIN in \muM N');
%         h4=ylabel('Irradiation in \mumol q. m^{-2}. s^{-1}');
%         %axis([0 4 0 1000 -0.2 .6])
%         %set(gca,'ztick',[-.2 0 .2 0.4 0.6])
%      	%axis ij
%         
%         %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%          %   'LineWidth',1.5,'FontSize',24)       
%         set([h2 h3 h4],'FontSize',24,'FontWeight','demi')
%         set([h1],'LineWidth',1)
%         %annotation('textbox','String',{'B'},'LineStyle','none',...
%          %   'Position',[0.12 0.7001 0.07031 0.07537],...
%           %  'FitHeightToText','on','FontSize',24);
% figure(12)
%   title('Coral respiration per total biomass (w photosynthesis)')
%         h1=mesh(x,y,(0.5*((symresp-symphotoo2).*Xvs+hostresp.*Xvh)./(symbiomass+hostbiomass))');
%         h2=zlabel('mol O_{2}]/ mol C day');
%         h3=xlabel('DIN in \muM N');
%         h4=ylabel('Irradiation in \mumol q. m^{-2}. s^{-1}');
%         %axis([0 4 0 1000 -0.2 .6])
%         %set(gca,'ztick',[-.2 0 .2 0.4 0.6])
%      	%axis ij
%         
%         %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%          %   'LineWidth',1.5,'FontSize',24)       
%         set([h2 h3 h4],'FontSize',24,'FontWeight','demi')
%         set([h1],'LineWidth',1)
%         %annotation('textbox','String',{'B'},'LineStyle','none',...
%          %   'Position',[0.12 0.7001 0.07031 0.07537],...
%           %  'FitHeightToText','on','FontSize',24);

        


