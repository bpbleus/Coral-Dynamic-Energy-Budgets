function [Y,fluxes,remark] = steadystatecoral(stoichp,debp,calcp,initp,envp)
% calculates ultimate densities of C and N containing compounds in an 
% exponentially growing culture given the density of 1 (density of 
% symbiont structure, which sets both C and N)
% Similar to the calculation steady states in partial recycling reactor; 
% IN:   parameter vector p and start values pp 
% OUT:  matrix Y: mes, hrc, meh, Xvh, Cin, Nin, Xgly, Xp1
%                 (Cin and Nin are dummy state vars with value approx 0)
%       Array fluxes: see around line 190-200
%       matrix remark: list of problems encountered in minimazation
% Erik Muller, 2/08, last edited 4/7/08
% Comments added 1/21/09


% NOTE ON NOTATION
% Notation closely follows paper; 
% 'X' are concentration measures, i.e. mass per system volume
% 1st letter refers to symbol; 2nd and beyond are mostly subscripts

% ALL MOLAR QUANTITIES ARE IN muMOL; TIME UNIT IS DAY

% YIELDS & ELEMENTAL COEFFICIENTS
nhgly = stoichp(1);nogly = stoichp(2);nngly = stoichp(3); nhx1 = stoichp(4);
nox1 = stoichp(5); nnx1 = stoichp(6); nhp1 = stoichp(7);nop1 = stoichp(8);
nnp1 = stoichp(9); nhes = stoichp(10); noes = stoichp(11);nnes = stoichp(12);
nhvs = stoichp(13);novs = stoichp(14);nnvs = stoichp(15);nheh = stoichp(16);
noeh = stoichp(17);nneh = stoichp(18);nhvh = stoichp(19);novh = stoichp(20);
nnvh = stoichp(21);yg3pl = stoichp(22);yesg3p = stoichp(23);
yvses = stoichp(24);yehx1 = stoichp(25);yp1x1 = stoichp(26);
yehgly = stoichp(27);yvheh = stoichp(28);ymeh = stoichp(29);yrheh = stoichp(30);

% DEB & EFFECT PARAMETERS (except OA)
jg3pasm = debp(1);jesasm = debp(6);jesds = debp(7);jxm = debp(8);
jcfhm = debp(9);jnfhm = debp(10);jehdh = debp(11);jeha3hm = debp(12);
kes = debp(13);rhoas = debp(14);Xkc = debp(15); Xkn = debp(16);
Xk1 = debp(17); keh = debp(18);kaph = debp(19);ds = debp(20);
jIDs = debp(21);jIDh = debp(22);Jlfsnec = debp(23); Jlfhnec = debp(24);

jIp = debp(2);jIb = debp(3);switchIp = debp(4); switchIb = debp(5);

% CALCIFICATION PARAMETERS
ycaas = calcp(1);ycaah = calcp(2);ycagh = calcp(3);ycadh = calcp(4);
omegaenv = calcp(5);komega = calcp(6);jcamin = calcp(7);

% ENVIRONMENTAL PARAMETERS 
Xc=envp(1);Xn=envp(2);X1=envp(3);Jlf=envp(4); 
XcU=envp(5);XnU=envp(6);X1U=envp(7);JlfU=envp(8);h=envp(9);interv=envp(10);
frJlfD = envp(11);

% DERIVED STOICHIOMETRIC PARAMETERS (y*c and yesc are variable)
ycesgs  =   1-yvses;
ynesgs  =   nnes-yvses*nnvs;
ycesds  =   1;
ynesds  =   nnes;
ycx1    =   1-yehx1-yp1x1;
ynx1    =   nnx1-yehx1*nneh-yp1x1*nnp1;
ycgly   =   1-yehgly;
yngly   =   yehgly*nneh;
ycehgh  =   1-yvheh;
ynehgh  =   nneh-yvheh*nnvh;
ycehdh  =   1-ymeh;
ynehdh  =   nneh;
yox1    =   1-0.5*nox1+0.25*nhx1-0.75*nnx1-yehx1*...
            (0.25*nheh-0.75*nneh-0.5*noeh+1)-...
            yp1x1*(0.25*nhp1-0.75*nnp1-0.5*nop1+1);
yogly   =   1.17-yehgly*(0.25*nheh-0.75*nneh-0.5*noeh+1);
yoehgh  =   0.5*(-noeh+2+0.5*nheh-1.5*nneh-yvheh*(2+novh+0.5*nhvh-1.5*nnvh));
yoehdh  =   1+0.25*nheh-0.75*nneh-0.5*noeh;
yoesgs  =   0.5*(-noes+2+0.5*nhes-1.5*nnes-yvses*(2+novs+0.5*nhvs-1.5*nnvs));
yoesds  =   1+0.25*nhes-0.75*nnes-0.5*noes;

% start values                    
mes0=initp(1);hrc0=initp(2);Xvs=initp(3);meh0=initp(4);Xvh0=initp(5);

% inorganic C & N in host. SHOULD NOT ACCUMULATE!! These are 'implicit' state variables
% whose values remain close to 0
Cin0=0.00001;Nin0=0.00001;
jcex=1e6;jnex=1e6;

XXn=(Xn:(XnU-Xn)/interv:XnU);     
n=length(XXn);
JJlf=(Jlf:(JlfU-Jlf)/interv:JlfU); 
nn=length(JJlf);
Y=NaN*ones(6,n,nn);
fluxes=NaN*ones(16,n,nn);
options=optimset('Display','Off','TolX',1e-6,'TolFun',1e-6,'MaxFunEvals',5e4);
remark=[];

% SOLVE STEADY STATE VALUES
y=[mes0,hrc0,meh0,Xvh0,Cin0,Nin0];
for ii=1:nn
   Jlf=JJlf(ii);
 	for i=1:n    
        Xn=XXn(i); 
        if i==1         % this loop ensures 'best' initial guesses
            [Yi,fval,exitflag,output]=fsolve(@coral, y,options);
            y=Yi;y2=y;
            if exitflag<1   % will print problematic ambient conditions
                exitflag
                Xn
                Jlf
                remark=[remark;exitflag Xn Jlf];
            end 
        else
            [Yi,fval,exitflag,output]=fsolve(@coral, y2,options);
            y2=Yi;
            if exitflag<1  % will print problematic ambient conditions
                exitflag
                Xn
                Jlf
                remark=[remark;exitflag Xn Jlf];
            end 
        end       
        Y(:,i,ii)=Yi;  
        fluxes(:,i,ii)=[jg3pas;jesas;jeha3h;jeha12h;B1;C1;B2;C2;B3;C3;...
                jvsgs;jvhgh;jglyf;joh;joresps; jophotos];
        if Yi(2)<0
            Y(:,i,ii)=NaN;
            fluxes(:,i,ii)=NaN;
        end        
    end
end



%-----------MODEL CODE STARTS HERE ---------------------------
        
function FF=coral(y)
  
mes=y(1);hrc=y(2);meh=y(3);Xvh=y(4);Cin=y(5);Nin=y(6);

% inorganic fluxes from the environment
    jlf=Jlf/ds;
    jlfD=frJlfD*jlf;
    jcfh=jcfhm*Xc/(Xkc+Xc);
    jnfh=jnfhm*Xn./(Xkn+Xn);
% Growth rates
    jesdsD=jesds*(1+(Jlfsnec/ds<jlfD).*(jlfD-Jlfsnec/ds)/jIDs);
    jehdhD=jehdh*(1+(Jlfhnec/ds<jlfD).*(jlfD-Jlfhnec/ds)/jIDh);
    jvsgs=(kes*mes-jesdsD)/(mes+1/yvses);
    jvhgh=(kaph*keh*meh-jehdhD)/(kaph*meh+1/yvheh);
% Food assimilation and feeding rates of host
    jx1=jxm*X1/(X1+Xk1);
    jeha12h = jx1*yehx1;     
% internal C and N fluxes 
    jcasstar=(1-yesg3p)*kes*mes/yesg3p;
    jnas=nnes*kes*mes;
    jccs=ycesgs*jvsgs/yvses+ycesds*jesdsD;
    jncs=ynesgs*jvsgs/yvses+ynesds*jesdsD;
    jcah=(ycx1/yehx1-ycgly/yehgly)*jeha12h+ycgly*keh*meh/yehgly;
    jnah=(ynx1/yehx1+yngly/yehgly)*jeha12h-yngly*keh*meh/yehgly;
    jna12h=ynx1*jeha12h/yehx1;
    jna3h=jna12h-jnah;
    jcch=ycehgh*jvhgh/yvheh+ycehdh*jehdhD;
    jnch=ynehgh*jvhgh/yvheh+ynehdh*jehdhD;  
% calculate nitrogen arrival fluxes at assimilation SUs
% Assumption: symbiont needs more N than it generates via catabolism
 	jstarcpluss=(jcfh+jcch+jcah)*Xvh/Xvs+jccs+jcasstar;
	N2sym=(jnfh+jnch+jna12h-jna3h);
	N2sym=(N2sym>0)*N2sym*Xvh;
	N2host=(N2sym+(jncs-jnas)*Xvs);
	N2host=(N2host>0)*N2host;
	jstarnpluss=N2sym/Xvs+jncs;
 	jstarnplush=(jnfh+jnch+jna12h);                    
% glyceraldehyde-3-phosphate production rate as function of light and CO2
	inhibp = (1+switchIp*jlf/jIp);
	inhibb = (1+switchIb*jlf/jIb);
	A1=1/jg3pasm ;
	B1=1/jstarcpluss;
	C1=1/(yg3pl*jlf);
	D1=1/(yg3pl*jlf+jstarcpluss);
	jg3pas=1/(A1*inhibp+(B1+C1-D1)*inhibb);
% Reserve production rate symbiont
	A2=1/jesasm;
	B2=nnes/(yesg3p*jstarnpluss); 
	C2=1/(yesg3p*rhoas*jg3pas);
	D2=nnes/(yesg3p*jstarnpluss+nnes*yesg3p*rhoas*jg3pas);
	jesas=1/(A2+B2+C2-D2);
% glycerol assimilation flux host 
	jglyas=jg3pas-kes*mes/yesg3p;
	A3=1/jeha3hm;
	C3=1/(yehgly*jglyas*Xvs/Xvh); 
	B3=yngly/(yehgly*jstarnplush);
	D3=yngly/(yehgly*yngly*jglyas*Xvs/Xvh+yehgly*jstarnplush);
	jeha3h=1/(A3+B3+C3-D3);
% Net assimilation host  
	jehah = jeha12h+jeha3h;
% spec glycerol excretion rate (into medium)
	jglyah=jglyas*Xvs/Xvh;
	jglyf=jglyah-jeha3h/yehgly;
% specific respiration rates
    yesc=jesas./jg3pas;
    yoc=1.17+0.5*yesc*(0.5*nhes-1.5*nnes-0.33-noes);
    joh=yox1*jeha12h/yehx1+yogly*jeha3h/yehgly+yoehgh*jvhgh/yvheh...
            +yoehdh*jehdhD;
	joresps=(yoesgs*jvsgs/yvses+yoesds*jesdsD);
    jophotos=yoc.*jesas./yesc;    

% ODEs      
	% symbiont reserve density
	dmes=(jesas-kes*mes);
	% symbiont structure
	dXvs=(jvsgs-hrc)*Xvs;
	% host reserve density
	dmeh=(jehah-keh*meh);
	% host structure
	dXvh=(jvhgh-hrc)*Xvh;
	% internal C and N concentrations; these are 'hidden' or dummy state variables
    % that need to be calculated in order to get around the problem 
    % of solving implicit equations
    % Their values need to remain close to 0
	dCin=(jcfh+jcch+jcah)*Xvh+(jcasstar+jccs-jg3pas)*Xvs-jcex*Cin;
	dNin=(jnfh+jnch+jna12h-jna3h)*Xvh+N2host-N2sym-jnex*Nin;  
FF=[dXvh;dXvs;dmes;dmeh;dCin;dNin];

end
end