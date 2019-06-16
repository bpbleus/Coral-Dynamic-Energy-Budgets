
% STEADY STATES CORAL MODEL + CALCIFICATION + PHOTO AND PARTIAL OA EFFECT
% DAMAGE MODEL IS MULTIPLICATIVE LINEAR RZ, LINEAR RD MODEL IN
% T.Klanjscek etal./JournalofTheoreticalBiology404(2016)361?374
% THIS MODEL IS REFFERED TO AS 'TinTox' in the code

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

% ALL MOLAR QUANTITIES ARE IN muMOL; TIME UNIT IS DAY (except TinTox
% module)
close all


% ENVIRONMENTAL PARAMETERS AND CONDITIONS
% to delete Jlf [(4:10:84)*10e6];%[(0.01:0.02:0.89)*10e7];
Xc      =   2000; 	% ambient inorganic C density (min of range)
XcU     =   Xc; 	% ambient inorganic C density (max of range)
interv  =   50;      % number of env. conditions in range
                   % IF THIS IS A SMALL NUMBER: FAST, BUT CRUDE PLOTS 
Xn      =   0.2;	% ambient inorganic N density (min of range)
XnU     =   4;	% ambient inorganic N density (max of range)
X1      =   0;      % ambient prey density   (min of range)
X1U     =   0;      % ambient prey density (max of range)
Jlf     =   3e6;    % flux of photons; 43.2e6: 500 mumol q/m2.s
JlfU    =   9e7;    % flux of photons; 43.2e6: 500 mumol q/m2.s
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
jIp     =   75;   % photoinhibition scaling parameter (production phase)
jIb     =	75; 	% photoinhibition scaling parameter (binding phase)
switchIp=   0;      % on/off switch photoinhibition (production phase)
switchIb=   0;      % on/off switch photoinhibition (binding phase)
 % if 0 < switch parameter < 1, interpretation is fraction of Jlf 
 % potentially causing photoinhibition
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

% vector of DEB and effect parameters (except OA and photodamage)
debp =[jg3pasm jIp jIb switchIp switchIb jesasm jesds jxm ...
        jcfhm jnfhm jehdh jeha3hm kes rhoas Xkc Xkn Xk1 keh kaph ds];

% photodamage parameters
% DIFFERENT UNITS: time scale is minutes, ROS in picomol (/mumol C) 
% damage in nanomol (/mumol C)
damageswitch =  1;      % if 1, then TinTox on, else off
kZs     =   36000/1440;	% ROS degradation rate, min^-1
yZD     =   1e-6*0.2;	% yield of ROS from damage,   
yDZ     =   1e-6*1;     % yield of damage from ROS
rhojlf  =   1e9*1e-4;	% probability quantum creates ROS, pmol ROS/mol quanta
rhoD    =   1e6*1e-5;	% probability quantum creates damage, nmol damage/ mol quanta
vZs     =   72000/1440;	% ROS neutralization rate coefficient, min^-1
vDs     =   0.5;        % Damage repair rate coefficient, min^-1
KmZs    =   1e9*1e-8;	% Half saturation constant of ROS neutralization 
                     	% induction, mumol ROS/ mol C
KmDs    =   1e6*1e-4;	% Half saturation constant of damage repair induction
                        % mmol damage/ mol C
yeD     =   1e-3*(1-yvses); 	% mol reserve C needed to repeair 1 mmol damage
damp=[damageswitch kZs yZD yDZ rhojlf rhoD vZs vDs KmZs  KmDs yeD];    
    
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
mes0    =   0.1688;    % initial guess symbiont reserve density
hrc0    =   0.0045;   % initial guess specific growth rate
Xvs     =   4;      % implied, arbitrary  
meh0    =   2.0545;  	% initial guess host reserve density
Xvh0    =   4.2;      % initial guesshost structure density
mZs0    =   1;   % initial guess ROS density (pmol ROS/mumol C)
mDs0    =   1;   % initial guess damage density (nmol D/ mumol C)


% first solve steady state damage and ROS densities (only depend on Jlf)
JJlf=(Jlf:(JlfU-Jlf)/interv:JlfU); 
nn=length(JJlf);
YTinTox=NaN*ones(2,nn);
syms SYMmDs SYMmZs;
options=optimset('Display','Off','TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',5e4);
if damageswitch==1 % find steady state damage and ROD densities at lowest jlf
    jlf=JJlf(1)/ds/1440;	
	dmZsSYM=(1+yZD*SYMmDs)*rhojlf.*jlf-SYMmZs*...
                (kZs+vZs*SYMmZs/(KmZs+SYMmZs));
	dmDsSYM=rhoD*jlf+yDZ*kZs*SYMmZs-vDs*SYMmDs^2/(KmDs+SYMmDs);
   	f = matlabFunction([dmZsSYM;dmDsSYM],'vars',{[ SYMmZs SYMmDs]});
    [Yi,fval,exitflag,output] = fsolve(f,[mZs0,mDs0],options);
    mZs = Yi(1);mDs=Yi(2);YTinTox(:,1)=Yi;
   	if exitflag<1  % will print problematic solutions
     	exitflag
    	mZs
       	mDs
    end
else
    YTinTox(:,:)=0;
end
if damageswitch==1 % find damage and ROS levels at subsequent jlf levels
    for ii=2:nn
        Jlfii=JJlf(ii);
        jlf=Jlfii/ds/1440;
        dmZsSYM=(1+yZD*SYMmDs)*rhojlf.*jlf-SYMmZs*...
                    (kZs+vZs*SYMmZs/(KmZs+SYMmZs));
        dmDsSYM=rhoD*jlf+yDZ*kZs*SYMmZs-vDs*SYMmDs^2/(KmDs+SYMmDs);
        f = matlabFunction([dmZsSYM;dmDsSYM],'vars',{[ SYMmZs SYMmDs]});
        [Yi,fval,exitflag,output] = fsolve(f,[mZs,mDs],options);
        mZs = Yi(1);mDs=Yi(2);YTinTox(:,ii)=Yi;
                if exitflag<1  % will print problematic solutions
                        exitflag
                        mZs
                        mDs
                end
    end
end
    

% vector of intitial values

initp = [mes0 hrc0 Xvs meh0 Xvh0];

%-------------- CALCULATIONS STEADY STATE VALUES CORAL MODEL-----------
%--------------------- START HERE --------------------------------------

[Y,fluxes,remark] =steadystatecoralphoto(stoichp,debp,...
        damp,calcp,initp,envp,YTinTox);

%----- PREPARE RESULTS FOR PLOTTING ----------------------
% Not all calculations are used for present plots;
% I keep them in so it is easy to add more results to the plots
% Change to 'subplot' for easier 'play'; removed it to get plots to publish

mes=squeeze(Y(1,:,:));hrc=squeeze(Y(2,:,:));meh=squeeze(Y(3,:,:));       
Xvh=squeeze(Y(4,:,:));Cin=squeeze(Y(5,:,:));Nin=squeeze(Y(6,:,:));
% mes(1,1)
% meh(1,1)
% Xvh(1,1)
mZs=YTinTox(1,:);mDs=YTinTox(2,:);

symbiomass=Xvs.*mes+Xvs; Ns=Xvs.*mes*nnes+Xvs*nnvs;
hostbiomass=Xvh.*meh+Xvh; Nh =Xvh.*meh*nneh+Xvh*nneh;
ratiostructures=Xvs./Xvh;

jg3pas=squeeze(fluxes(1,:,:));jesas=squeeze(fluxes(2,:,:));
jeha3h=squeeze(fluxes(3,:,:));jeha12h=squeeze(fluxes(4,:,:));
timeCs=squeeze(fluxes(5,:,:));timeLs=squeeze(fluxes(6,:,:));
timeNs=squeeze(fluxes(7,:,:));timeGLYs=squeeze(fluxes(8,:,:));
timeNh=squeeze(fluxes(9,:,:));timeGLYh=squeeze(fluxes(10,:,:));
jvsgs=squeeze(fluxes(11,:,:));jvhgh=squeeze(fluxes(12,:,:));
jglyf=squeeze(fluxes(13,:,:));jo2s_resp =squeeze(fluxes(14,:,:));
jo2s_photo=squeeze(fluxes(15,:,:)); jo2h_auto=squeeze(fluxes(16,:,:));
jo2h_hetero =squeeze(fluxes(17,:,:));
% calculate capacities and limitation coefficients
capphoto=jg3pas/jg3pasm;
limcoeffphoto=log10(timeLs./timeCs);
capas=jesas/jesasm;
limcoeffas=log10(timeNs./timeGLYs);
capa3h=jeha3h/jeha3hm;
limcoeffa3h=log10(timeNh./timeGLYh);

%calculate yields
yglyfc=jglyf.*Xvh./(jg3pas.*Xvs);
yehc=jeha3h.*Xvh./(jg3pas.*Xvs);
yesc=jesas./jg3pas;

autoCinhost=jeha3h./(jeha12h+jeha3h);

jesds=jesds+1440*yeD*vDs*mDs.^2./(KmDs+mDs);
   
domainN=(Xn:(XnU-Xn)/interv:XnU);     
domainJlf=(Jlf:(JlfU-Jlf)/interv:JlfU)/(24*3600);
x=domainN; y=domainJlf;
%[x,y]=meshgrid(domainN,domainJlf);

pos=[0.75 0.12 0.05 0.8];

figure(2) 
  	%title('Symbiont to host ratio (\{propto} surface area)')
	h1=imagesc(x,y,ratiostructures');
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiance, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[ 1 2 3 4],'ytick',[200 400 600 800 1000],...
               'LineWidth',1.5,'FontSize',24,'FontWeight','demi') 
	h2=colorbar;
colormap(flipud(jet))
    text(1.33*max(x),0.5*max(y),'mol C/ mol C (structural biomass)','Rotation',270,...   
        'FontSize',32,'FontWeight','demi','HorizontalAlignment','center')    
    text(0.5*max(x),0.95*max(y),'Symbiont density ',...
        'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
        'Color',[1 1 1])  
	text(0.5*max(x),0.9*max(y),'\propto surface area',...
        'FontSize',28,'FontWeight','bold','HorizontalAlignment','center',...
        'Color',[1 1 1]) 
% 
% % figure(3) 
% % 	%title('Symbiont density (biomass)')
% %     h1=imagesc(x,y,(symbiomass./(symbiomass+hostbiomass))');   
% % 	axis xy
% % 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% % 	ylabel('Irradiance, \mumol q. m^{-2}. s^{-1}',...
% %         'FontSize',36,'FontWeight','demi');
% %     set(gca,'xtick',[ 1 2 3 4],'ytick',[200 400 600 800 1000],...
% %                'LineWidth',1.5,'FontSize',24,'FontWeight','demi') 
% % 	h2=colorbar;
% %     asp = [1 10 1];
% %  	caxis([0.095 0.165]);
% %     yt = (0.1:0.02:0.18);
% %     %caxis([0.21 0.36]);
% %     %yt = (0.22:0.03:0.36);
% %     set(h2,'Position',pos,'PlotBoxAspectRatio',asp,'FontSize',24)
% %     set(h2,'YTick',yt)
% %     text(1.33*max(x),0.5*max(y),'mol symbiont C/ mol total C','Rotation',270,...
% %     	'FontSize',32,'FontWeight','demi','HorizontalAlignment','center')    
% %     text(0.5*max(x),0.95*max(y),'Symbiont biomass density',...
% %         'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
% %         'Color',[1 1 1 ])  
% % 
% figure(4) 
% 	%title('Dark Respiration (\propto surface area)')
%     h1=imagesc(x,y,(jo2s_resp.*Xvs./Xvh+jo2h_hetero)');   
% 	axis xy
% 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% 	ylabel('Irradiance, \mumol q. m^{-2}. s^{-1}',...
%         'FontSize',36,'FontWeight','demi');
%     set(gca,'xtick',[ 1 2 3 4],'ytick',[200 400 600 800 1000],...
%                'LineWidth',1.5,'FontSize',24,'FontWeight','demi') 
% 	h2=colorbar;
%     asp = [1 10 1];
%     caxis([0.2 0.65]);
%     %caxis([0.26 0.51]);
%     yt = (0.3:0.05:0.5);
%     set(h2,'Position',pos,'PlotBoxAspectRatio',asp,'FontSize',24)
%    %	set(h2,'YTick',yt)
%     text(1.33*max(x),0.5*max(y),'mol O_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
%         'FontSize',32,'FontWeight','demi','HorizontalAlignment','center')    
%     text(0.5*max(x),0.95*max(y),'"Dark" respiration rate',...
%         'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
%         'Color',[1 1 1])  
% 	text(0.5*max(x),0.9*max(y),'\propto surface area',...
%         'FontSize',28,'FontWeight','bold','HorizontalAlignment','center',...
%         'Color',[1 1 1]) 
% figure(5) 
% 	%title('Net respiration rate (\propto surface area)')
%     h1=imagesc(x,y,(-1*((jo2s_resp-jo2s_photo).*Xvs./Xvh+jo2h_auto+jo2h_hetero))');   
% 	axis xy
% 	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
% 	ylabel('Irradiance, \mumol q. m^{-2}. s^{-1}',...
%         'FontSize',36,'FontWeight','demi');
%     set(gca,'xtick',[ 1 2 3 4],'ytick',[200 400 600 800 1000],...
%                'LineWidth',1.5,'FontSize',24,'FontWeight','demi') 
% 	h2=colorbar;
%     text(1.33*max(x),0.5*max(y),'mol O_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
%         'FontSize',32,'FontWeight','demi','HorizontalAlignment','center')    
%     text(0.5*max(x),0.95*max(y),'Net O_{2} production',...
%         'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
%         'Color',[0 0 0])  
% 	text(0.5*max(x),0.9*max(y),'\propto surface area',...
%         'FontSize',28,'FontWeight','bold','HorizontalAlignment','center',...
%         'Color',[0 0 0]) 
figure(6) 
	%title('Photosynthesis (\propto surface area)')
    h1=imagesc(x,y,(jo2s_photo.*Xvs./Xvh)');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiance, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[ 1 2 3 4],'ytick',[200 400 600 800 1000],...
               'LineWidth',1.5,'FontSize',24,'FontWeight','demi') 
	h2=colorbar;
    text(1.33*max(x),0.5*max(y),'mol O_{2} mol host structural C^{-1} d^{-1}','Rotation',270,...
            'FontSize',32,'FontWeight','demi','HorizontalAlignment','center')    
    text(0.5*max(x),0.95*max(y),'Photosynthesis rate',...
        'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
        'Color',[0 0 0])  
	text(0.5*max(x),0.9*max(y),'\propto surface area',...
        'FontSize',28,'FontWeight','bold','HorizontalAlignment','center',...
        'Color',[0 0 0]) 


figure(66) 
	%title('Photosynthesis (\propto surface area)')
    h1=imagesc(x,y,jo2s_photo');   
	axis xy
	xlabel('DIN, \muM N','FontSize',36,'FontWeight','demi');
	ylabel('Irradiance, \mumol q. m^{-2}. s^{-1}',...
        'FontSize',36,'FontWeight','demi');
    set(gca,'xtick',[ 1 2 3 4],'ytick',[200 400 600 800 1000],...
               'LineWidth',1.5,'FontSize',24,'FontWeight','demi') 
	h2=colorbar;
 
    text(1.33*max(x),0.5*max(y),'mol O_{2} mol symbiont structural C^{-1} d^{-1}','Rotation',270,...
            'FontSize',32,'FontWeight','demi','HorizontalAlignment','center')    
    text(0.5*max(x),0.95*max(y),'Photosynthesis rate',...
        'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
        'Color',[1 1 1])  
	text(0.5*max(x),0.9*max(y),'per unit symbiont structural biomass',...
        'FontSize',28,'FontWeight','bold','HorizontalAlignment','center',...
        'Color',[1 1 1]) 

    


