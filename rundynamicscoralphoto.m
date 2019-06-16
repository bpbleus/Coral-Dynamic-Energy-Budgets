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

startSS =   1;      % if 1 initial values are steady state values

Xc      =   2000; 	% initial ambient inorganic C density
XcU     =   Xc; 	% upper limit when starting at SS
Xc2      =   Xc; 
Xn      =   0.2;	% initial ambient inorganic N density
XnU     =   Xn*1.00001;	% ambient inorganic N density (max of range)
Xn2      =   Xn;	% Xn after shift (dynamic, no start at SS)
X1      =   0;      %0.1; 	% initial prey density 
X1U     =   0; 
X12      =   X1;
qms100  =   8.64e6;
Jlf     =   3e6;    % flux of photons; 43.2e6: 500 mumol q/m2.s
JlfU    =   1.001*Jlf;%0.0001*qms100;    % flux of photons; 43.2e6: 500 mumol q/m2.s
Jlf2     =   1.001*Jlf;%0*qms100;    % flux of photons; 43.2e6: 500 mumol q/m2.s
hrc     =   0.01;   % change later for SS harvesting rate = mortality and/or dilution rate
hrc2     =   hrc;   % set later harvesting rate = mortality and/or dilution rate

interv  =   2;      % number of env. conditions in range

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
nnes    =   0.2;   % N:C ratio in symbiont reserves
nhvs    =   1.8;    % H:C ratio in symbiont structure
novs    =   0.4;    % O:C ratio in symbiont structure
nnvs    =   0.15;   % N:C ratio in symbiont structure
nheh    =   1.8;    % H:C ratio in host reserves
noeh    =   0.4;    % O:C ratio in host reserves
nneh    =   0.15;    % N:C ratio in host reserves
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
yeD     =   1e-3*0.1; 	% mol reserve C needed to repeair 1 mmol damage

% Solve steady state ROS and damage levels; assumptions: growth dilution =
% 0 and pseudosteady state is quickly reached

options=optimset('Display','Off','TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',5e4);
if damageswitch==1
    syms SYMmDs SYMmZs;
    jlf=Jlf/ds/1440;
    dmZsSYM=(1+yZD*SYMmDs)*rhojlf.*jlf-SYMmZs*(kZs+vZs*SYMmZs/(KmZs+SYMmZs));
    dmDsSYM=rhoD*jlf+yDZ*kZs*SYMmZs-vDs*SYMmDs^2/(KmDs+SYMmDs);
    f = matlabFunction([dmZsSYM;dmDsSYM],'vars',{[ SYMmZs SYMmDs]});
    [Yi,fval,exitflag,output] = fsolve(f,[1e-0,1e-0],options);
    mZs = Yi(1);
    mDs=Yi(2);
    if exitflag<1  % will print problematic solutions
                exitflag
                mZs
                mDs
    end
else
    mZs = 0;mDs=0;
end

% vector of photodamage parameters

damp=[damageswitch kZs yZD yDZ rhojlf rhoD vZs vDs KmZs  KmDs yeD mZs mDs];
    
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

% INTEGRATION PARAMETERS                    
if startSS == 1	   
    h       =   1;      % dilution rate of solubles 
    frJlfD  =   1;      % fraction of Jlf potentially causing damage
    mes0    =   0.2;    % initial guess symbiont reserve density
    hrc0    =   0.042;   % initial guess specific growth rate
    Xvs     =   4;      % implied, arbitrary  
    meh0    =   2.91;  	% initial guess host reserve density
    Xvh0    =   3;      % initial guesshost structure density
    mZs0    =   1;   % initial guess ROS density (pmol ROS/mumol C)
    mDs0    =   1;   % initial guess damage density (nmol D/ mumol C
    YTinTox =   ones(2,interv+1);
    initp = [mes0 hrc0 Xvs meh0 Xvh0];
    envpSS = [Xc Xn X1 Jlf XcU XnU X1U JlfU h interv frJlfD];
    [Y,fluxes,remark] =steadystatecoralphoto(stoichp,debp,...
        damp,calcp,initp,envpSS,YTinTox);
    mes0    =   Y(1,1,1);   % initial symbiont reserve density
    Xvs0    =   Xvs	% initial symbiont structure density
    hrc     =   Y(2,1,1);
    meh0    =   Y(3,1,1);      % initial host reserve density
    Xvh0    =   Y(4,1,1); 	% initial host structure density
    envp(5) =   hrc;
    hrc2    =   hrc;%Y(2,end,end);
    t0      =   0;      % start time
    tend    =   50;	% end time of integration  
    t02     =   tend;
    tend2   =   20+tend; % end time shift up
else
    mes0    =   0.18;   % initial symbiont reserve density
    Xvs0    =   4;	% initial symbiont structure density  
    meh0    =   2;      % initial host reserve density
    Xvh0    =   8; 	% initial host structure density
    t0      =   0;      % start time
    tend    =   50;
    t02      =   tend;      % start time
    tend2    =   t02+tend;	% end time of integration
end

% vector of environmental parameters and conditions

envp = [Xc Xn X1 Jlf hrc];
% vector of integration parameters

intp = [mes0 Xvs0 meh0 Xvh0 t0 tend];

%--------------------- CALCULATIONS-------------------------------------
%--------------------- START HERE --------------------------------------

[T,Y,fluxes] = dynamicscoralphoto(stoichp,debp,damp,calcp,intp,envp);

mes=Y(:,1); Xvs=Y(:,2); meh=Y(:,3); Xvh=Y(:,4); Cin=Y(:,5); Nin=Y(:,6);

% 'UNCOMMENT' THE NEXT 30+ LINES TO ACTIVATE SHIFT-UP/ DOWN OF 
% ENVIRONMENTAL AND/OR CALCIFICATION/OA PARAMETERS 

ycaas2 =     ycaas;   % yield of CaCO3 due to assimilation symbiont
ycaah2 =     ycaah;   % yield of CaCO3 due to assimilation host
ycagh2 =     ycagh;   % yield of CaCO3 due to growth host
ycadh2 =     ycadh;   % yield of CaCO3 due to dissipation (maintenance) host
omegaenv2 =  omegaenv;      % aragonite saturation state
komega2 =    komega;    % OA effect scaling parameter, muM^-1
jcamin2 =    jcamin;      % Dissolution rate (default 0)
% add H+ toxicity
calcp2 = [ycaas2 ycaah2 ycagh2 ycadh2 omegaenv2 komega2 jcamin2];

mes02    =   mes(end);  % initial symbiont reserve density
Xvs02    =   Xvs(end);	% initial symbiont structure density  
meh02    =   meh(end); 	% initial host reserve density
Xvh02    =   Xvh(end); 	% initial host structure density

intp2    = [mes02 Xvs02 meh02 Xvh02 t02 tend2];
envp2 =[Xc2 Xn2 X12 Jlf2 hrc2];

if damageswitch==1
    syms SYMmDs SYMmZs;
    jlf=Jlf/ds/1440;
    dmZsSYM=(1+yZD*SYMmDs)*rhojlf.*jlf-SYMmZs*(kZs+vZs*SYMmZs/(KmZs+SYMmZs));
    dmDsSYM=rhoD*jlf+yDZ*kZs*SYMmZs-vDs*SYMmDs^2/(KmDs+SYMmDs);
    f = matlabFunction([dmZsSYM;dmDsSYM],'vars',{[ SYMmZs SYMmDs]});
    [Yi,fval,exitflag,output] = fsolve(f,[1e-0,1e-0],options);
    mZs = Yi(1)
    mDs=Yi(2)

    if exitflag<1  % will print problematic solutions
                exitflag
                mZs
                mDs
    end
else
    mZs = 0;mDs=0;
end
damp2=[damageswitch kZs yZD yDZ rhojlf rhoD vZs vDs KmZs  KmDs yeD mZs mDs];

[T2,Y2,fluxes2] = dynamicscoralphoto(stoichp,debp,damp2,calcp2,intp2,envp2);

T=[T;T2]; mes=[mes;Y2(:,1)];Xvs=[Xvs;Y2(:,2)]; meh=[meh;Y2(:,3)];
    Xvh=[Xvh;Y2(:,4)];Cin=[Cin;Y2(:,5)];Nin=[Nin;Y2(:,6)];
Xvh(1:50:end);
fluxes=[fluxes;fluxes2];

symbiomass=Xvs.*mes+Xvs; hostbiomass=Xvh.*meh+Xvh;

jcalight    = fluxes(:,1); % calcification rate (per Xv) 
jcadark     = fluxes(:,2); % dark calcification rate (per Xv)
jo2s_resp   = fluxes(:,3);  % symbiont O2 consumption
jo2s_photo  = fluxes(:,4);  % symbiont O2 production
jo2h_auto	= fluxes(:,5); % host O2 consumption due to photosynthate assimilation 
jo2h_hetero	= fluxes(:,6); % host O2 consumption except due to photosynthate ass
jg3pas  	= fluxes(:,7);
jvsgs  	= fluxes(:,8);
jvhgh  	= fluxes(:,9);
wildcard  	= fluxes(:,10);

 figure(1)
%[h1,h2,h3]=plotyy(T,yoesgs*jvsgs/yvses,T,yoehgh*jvhgh/yvheh);

    h1=plot(T,meh,'k','LineWidth',2);
    hold on
    h1=plot(T,mes);
    set(h1,'LineStyle','--','LineWidth',2,'Color','k')
    h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
    h3=ylabel('Reserve densities','FontSize',36,'FontWeight','demi');
    %axis([0 4 0 1000 0 0.2])	
	set(gca,'FontSize',36,'LineWidth',2)
        
figure(2)
    h1=plot(T,Xvs./Xvh,'k','LineWidth',3);
    h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
    h3=ylabel('Symbiont density (Xs/Xv)','FontSize',36,'FontWeight','demi');

	set(gca,'FontSize',36,'LineWidth',2)
    hold on

% figure(33)
%     h1=plot(T,symbiomass./(symbiomass+hostbiomass),'LineWidth',2);
%    	h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
%     h3=ylabel('Symbiont biomass density (mol C/ mol C)',...
%         'FontSize',36,'FontWeight','demi');
%     %axis([0 4 0 1000 0 0.2])
% 	set(gca,'FontSize',36,'LineWidth',2)

% figure(14)
%    % h1=plot(T,jo2s_photo.*Xvs./Xvh+jo2h_auto,'k','LineWidth',3);
%   %  hold on
%   
%   
%   neto2=-1*((jo2s_resp-jo2s_photo).*Xvs./Xvh+jo2h_auto+jo2h_hetero);
%     h1=plot(T,neto2,'r','LineWidth',3);
% % 
%           text(0.5*max(T),4.25,'Net O_{2} production rate',...
%         'FontSize',32,'FontWeight','bold','HorizontalAlignment','center',...
%         'Color',[0 0 0])
%     text(0.5*max(T),4,'\propto surface area',...
%         'FontSize',28,'FontWeight','bold','HorizontalAlignment','center',...
%         'Color',[0 0 0])
%     hold on
% %     set(h1,'LineStyle','--','LineWidth',3,'Color','k')
% %   	h1=plot(T,jo2s_resp.*Xvs./Xvh+jo2h_hetero));
% %     set(h1,'LineStyle',':','LineWidth',3,'Color','k')
%     h2=xlabel('Time, d','FontSize',36,'FontWeight','demi');
%     h3=ylabel('mol O_{2} mol C^{-1} d^{-1}','FontSize',36,'FontWeight','demi');
% axis([0 2 0 4.5])	
% set(gca,'FontSize',36,'LineWidth',2,'xtick',[0 0.5 1 1.5 2],'ytick',[1 2 3 4]) 
% figure(15)
% h1=plot(T,-1*((jo2s_resp-jo2s_photo).*Xvs./Xvh+jo2h_auto+jo2h_hetero)...
%             ,'k','LineWidth',3);
%     h2=xlabel('Time, d','FontSize',36,'FontWeight','demi');
%     h3=ylabel('O_{2} evolution, mol O_{2} mol C^{-1} d^{-1}','FontSize',36,'FontWeight','demi');
%    % axis([0 4 0 1000 0 0.2])	
% 	set(gca,'FontSize',36,'LineWidth',2)
 figure(16)   
 plot(T,jvsgs,T,jvhgh)
% figure(5)
%     h1=plot(T,jg3pas.*Xvs./Xvh,'LineWidth',2);
%    	h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
%     h3=ylabel('Symbiont biomass density (mol C/ mol C)',...
%         'FontSize',36,'FontWeight','demi');
%     %axis([0 4 0 1000 0 0.2])
% 	set(gca,'FontSize',36,'LineWidth',2)  
%  figure(6)
%     h1=plot(T,jo2s_resp,'LineWidth',2);
%    	h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
%     h3=ylabel('Symbiont biomass density (mol C/ mol C)',...
%         'FontSize',36,'FontWeight','demi');
%     %axis([0 4 0 1000 0 0.2])
% 	set(gca,'FontSize',36,'LineWidth',2)     
%     
%  figure(7)
%     h1=plot(T,wildcard,T,jg3pas);%'LineWidth',2);
%    	h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
% %     h3=ylabel('Symbiont biomass density (mol C/ mol C)',...
% %         'FontSize',36,'FontWeight','demi');
% %     %axis([0 4 0 1000 0 0.2])
% % 	set(gca,'FontSize',36,'LineWidth',2)    
% %        
% % figure(8)
% %     h1=plot(T,jo2h_auto,'LineWidth',2);
% %    	h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
% %     h3=ylabel('Symbiont biomass density (mol C/ mol C)',...
% %         'FontSize',36,'FontWeight','demi');
% %     %axis([0 4 0 1000 0 0.2])
% % 	set(gca,'FontSize',36,'LineWidth',2) 
% figure(66)
%     h1=plot(T,Xvh,T,Xvs);
%    	h2=xlabel('Time (d)','FontSize',36,'FontWeight','demi');
%     h3=ylabel('Xv','FontSize',36,'FontWeight','demi');
%     %axis([0 4 0 1000 0 0.2])
% 	set(gca,'FontSize',36,'LineWidth',2)
% % 
% 
% 
% 



% figure(1)
%     h1=plot(T,jcalight,T,jcadark);  
%     h2=xlabel('Time (d)');
%     h3=ylabel('Specific calcification rate');%(mol CaCO_{3}. mol biomass C^{-1}. day^{-1})
%     %axis([0 4 0 1000 0 0.2])
%     %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%     %        'LineWidth',1.5,'FontSize',24)       
% 	set([h2 h3],'FontSize',24,'FontWeight','demi')
% 	set([h1],'LineWidth',1)
% 
% figure(2)
%     h1=plot(T,Cin,T,Nin,'--');  
%     h2=xlabel('Time (d)');
%     h3=ylabel('Internal C and N pools');
%     title('CONTROL: POOLS MUST REMAIN CLOSE TO 0')
%     %axis([0 4 0 1000 0 0.2])
%     %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%     %        'LineWidth',1.5,'FontSize',24)       
% 	set([h2 h3],'FontSize',24,'FontWeight','demi')
% 	set([h1],'LineWidth',1)
% 
% figure(3)
%     h1=plot(T,mes,T,meh);  
%     h2=xlabel('Time (d)');
%     h3=ylabel('Reserve densities');
%     %axis([0 4 0 1000 0 0.2])
%     %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%     %        'LineWidth',1.5,'FontSize',24)       
% 	set([h2 h3],'FontSize',24,'FontWeight','demi')
% 	set([h1],'LineWidth',1)
%         
% figure(4)
%     h1=plot(T,Xvs./Xvh);
%    	h2=xlabel('Time (d)');
%     h3=ylabel('Symbiont density (Xs/Xv)');
%     %axis([0 4 0 1000 0 0.2])
%     %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%     %        'LineWidth',1.5,'FontSize',24)       
% 	set([h2 h3],'FontSize',24,'FontWeight','demi')
% 	set([h1],'LineWidth',1)  
% 
% figure(5)
%     h1=plot(T,symbiomass./(symbiomass+hostbiomass));
%    	h2=xlabel('Time (d)');
%     h3=ylabel('Symbiont biomass density (mol C/ mol C)');
%     %axis([0 4 0 1000 0 0.2])
%     %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%     %        'LineWidth',1.5,'FontSize',24)       
% 	set([h2 h3],'FontSize',24,'FontWeight','demi')
% 	set([h1],'LineWidth',1)  
%     
% figure(6)
%     h1=plot(T,Xvh,T,Xvs);
%    	h2=xlabel('Time (d)');
%     h3=ylabel('Xv');
%     %axis([0 4 0 1000 0 0.2])
%     %set(gca,'xtick',[0 1 2 3 4 ],'ytick',[0 250 500 750 1000],...
%     %        'LineWidth',1.5,'FontSize',24)       
% 	set([h2 h3],'FontSize',24,'FontWeight','demi')
% 	set([h1],'LineWidth',1)  
 

    