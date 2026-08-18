% =============================================================================
% 3D Fiber-Based Cell Migration Model
% =============================================================================
% Final stage of the modeling pipeline.
%
% This script uses cell-mechanical and adhesion parameters obtained from the
% upstream signaling/cytoskeletal model to simulate single-cell migration
% through a 3D fibrous extracellular matrix (ECM).
%
% Main inputs include:
%   - koff / kon      : adhesion bond dissociation and association rates
%   - polarity        : probability controlling persistence of migration
%   - kcell           : effective cellular spring constant
%   - Vpseudo         : pseudopod extension velocity
%   - Fmax            : maximum contractile force
%
% The model explicitly represents:
%   1. ECM fiber organization and crosslinks
%   2. Integrin-ligand binding
%   3. Pseudopod extension and retraction
%   4. Actomyosin-driven contraction
%   5. 3D cell displacement and trajectory generation
%
% Model outputs include 3D trajectories, invasiveness, average velocity,
% persistence length, and phase-dependent migration behavior.
% =============================================================================

% Fiber based cell migration code - Single cell persistent random walk

clear all
%delete(gcp('nocreate'))
 c = parcluster;
% c.NumWorkers = 20;
% saveProfile(c);
%parpool(c)
%c

% =============================================================================
% 1. SELECT SIMULATION CONDITION
% =============================================================================
% 'para' selects one of the predefined ECM / cell-state parameter sets below.
% Each condition contains a specific combination of collagen concentration,
% stiffness, adhesion kinetics, polarity, cell stiffness, and protrusion rate.
para=4;%
C_gel_1=[1 2 4 1 2 4];
Stiffness_1=[30 150 600 30 150 600];
AI_1=0;
%% 
name=["AI0Black1mg_","AI0Black2mg","AI04Black4mg","AI0White1mg", "AI0White2mg","AI0White4mg"];

degredation=0;
im=0;

num_points=50;

countnumber = 0;
contractilitytime = 0;


gadd_fiber=[];
%
Vinst=0;
counter=0;

%Black --> para=1
%White --> para=2



% =============================================================================
% 2. CELL-MECHANICAL INPUTS
% =============================================================================
% These parameter arrays define the cell state used by the migration model.
% In the complete pipeline, these quantities are informed by the upstream
% signaling/cytoskeletal model and are then passed into the 3D migration model.
%
% Fmax      : maximum contractile force
% tsearch   : characteristic pseudopod search/extension time
% koff/kon  : adhesion bond dissociation/association rates
% polarity  : probability controlling directional persistence
% kcell     : effective cell spring constant
% Vpseudo   : pseudopod extension speed
%
Fmax_1=[2.08E-09	2.09E-09	2.10E-09	2.19E-09	2.19E-09	2.21E-09];
tsearch_1=[307.6696705	314.1086752	316.7259194	664.7687962	654.6108944	647.1731002];
koff_1=[0.593734929	1.066852099	0.857855499	1.143291349	0.976954049	0.360345891];
kon_1=[0.000951007	0.001011239	0.001126061	0.002128564	0.002171587	0.002270644];
polarity_1=[0.515404157	0.513305829	0.507505906	0.276669635	0.299357317	0.316814797];
kcell_1=[7.30018E-10	7.53873E-10	8.19183E-10	5.46182E-09	5.49301E-09	5.55463E-09];
Vpseudo_1=[0.098080321	0.100883534	0.106202486	0.183502008	0.18426506	0.187550201];





Bmax_1=140 ;
Bhalf =40;%....................Number of bound ligands at half max. force
Bmin=20;%...............Min bound ligands needed for pseudopod outgrowth f




Integrin_1=327528.750;











 datapoints=1; %increase
 Avg_fold=linspace(-1,1,5);



cellinvasive=zeros(10,10);





filename = strcat(name(para),date);
rng(datenum(date));%......................................For repeatability
tic

% =============================================================================
% 3. MODEL DEFINITIONS AND SIMULATION TIME
% =============================================================================
% Define simulation duration, time step, cell geometry, ECM properties, and
% matrices used to store trajectories and migration metrics.
%
%Model Definitions
cellAx = zeros(1,3,1);
gpd = zeros(1,3,1);
pV2 = zeros(1,3);
gPV = zeros(1,3,1);
gFiberDir = zeros(1,3,1);

runTime =12*3600;%...........................Total runtime for model (sec)v
dt =1;%.......................................Size of each time step (sec)
elmtSize = 0.01;%............Size of binding site element (um) - max is 0.3
samples =10;%..................................Sample size to get averages
dataPnts = 1;%.............................Number of parameter data points



%Initiate matrices for MSD calculation
time = ((0:1+runTime/dt-1)');%...................Time at each cell position
% tMax =11;%.............Maximum time for MSD fitting
% tMsd = dt*((1:tMax)')/3600;%.....................Time steps for MSD fitting



% =============================================================================
% 4. CONSTANT CELL AND ECM ATTRIBUTES
% =============================================================================
% These values define geometric and physical properties that remain fixed
% during a simulation, including cell size, collagen geometry, viscosity,
% bond stiffness, and temperature.
%
%Constant Cell Attributes
ve=[1,0,0];
prevve=[1,0,0];
Lsearch =0.3;%...............Length of pseudopod tip in contact with fiber
Dcell = 15;%...........................................Cell's diameter (um)
Vcell = (4/3)*pi*(Dcell/2)^3;%.........................Cell's volume (um^3)



%Bmax =130;%.............Min bound ligands needed for pseudopod contraction 
retractionspeed=0.05; % (um/sec)
%Constant ECM Attributes
Ltropo = 0.3;%.................................Length of tropocollagen (um)
Dtropo = 0.0015;%............................Diameter of tropocollagen (um)
Lfiber = 75;%.............................................Fiber length (um)
otAngle = 0;%...................Average angle between ECM fibers and x-axis
vRef = [1,0,0];%.........................Fiber orientation referance vector
eta = 1e-11;%............................Extracellular viscosity (N*s/um^2)
kI = 0.25e-9;%........................................Bond stiffness (N/um)
%kon = 1.4e-3;%...........................Bond association rate  (um^2/s^-1)
kB = 1.38e-17;%.........................Boltzmann's contant (um^2*kg/s^2*K)
T = 311;%...................................................Temperature (K)

%Output Parameter Preallocations
avgVel = zeros(dataPnts,samples);%.......................Average cell speed
LpR = zeros(dataPnts,samples);%.................Persistence length by <R^2>
D = zeros(dataPnts,samples);%...................Random motility coefficient
alp = zeros(dataPnts,samples);%.......................................Alpha
R1 = zeros(dataPnts,samples);%.....................Velocity Goodness of fit
R2 = zeros(dataPnts,samples);%..........................MSD Goodness of fit
R3 = zeros(dataPnts,samples);%..........................LpR Goodness of fit
MSD=zeros(9,samples); %...........................Mean Square Displacement

%=========================================================================%
% =============================================================================
% 5. ASSIGN CONDITION-SPECIFIC PARAMETERS
% =============================================================================
% Select the cell and ECM parameters associated with the current condition.
% This step connects the upstream predicted cell state to the migration model.
%
%Parametric Sweep
%=========================================================================%
for par=1:datapoints
    Fmax = Fmax_1(para);
    %p=[2,0.0058,tsearch_1(para),0,6];%2.5,0.0058 6,0.0132
    koff =koff_1(para);%.............................Bond dissociation rate .1-100 (s^-1) 
    dRecp =Integrin_1;%..............................Integrin receptors per cell 
    dR = dRecp/saellipsoid(15/2);%...Integrin receptor density (receptors/um^2)
    kcell = kcell_1(para);%.........................................Cell spring constant
    polarize=polarity_1(para);%................................................Polarize cell
    Vpseudo =Vpseudo_1(para);%...................Velocity of extending pseudopod (um/sec)
    tC = Lsearch/Vpseudo;%............................Binding site contact time
    Bmax = Bmax_1;%....Min bound ligands needed for pseudopod contraction, added to loop 12/12/2024
    kon= kon_1(para);%.......................Bond association rate  (um^2/s^-1), added to loop 12/12/2024


    C_gel=C_gel_1(para);%..................................Collegen concentration
    fiberDens = 0.0021*C_gel+ 0.0006;%.......................................Fiber densit
    avgSearchTime = tsearch_1(para);%.......................Pseudopod extension time
    AI =AI_1;%.............................................Fiber Alignment %0.22
    Fiber_D = sqrt((C_gel*6.022e20*(Ltropo+0.067)*Dtropo^2)/...
    (1e12*fiberDens*805000*0.7*Lfiber* 0.9));%.......Average fiber diameter
    RGD =6;%.......................................RGD Peptides/monomer
    dL = (RGD*.7*.9)/((Ltropo+0.067)*Dtropo);%.............RGD ligand density
   
   
  
    sigma = sqrt(-1642*reallog(AI));%.Angle deviation from reference vector   
    if isinf(sigma)
        sigma = 1e10;
    end
    %sigma=48
    crslnkPerL = ((4.439*sigma)/(sigma+4.55))/150;%**..........Crosslinks per fiber %high denisty:150 lowdensity:160
    %crslnkPerL=0.0296;%0.0296
    crslnkPerF = Lfiber*crslnkPerL;%Crosslink number per Fiber length
    crslnkDens = fiberDens*crslnkPerF%..................Crosslink density
    
    %Equation from Lin et. al. 2015
    if crslnkPerF < 3.5
        gel_stiffness = 1039.9*crslnkPerF-1992.9
    else
        gel_stiffness = 5247.9*crslnkPerF-16274
    end%.................................................Gel stiffness (Pa)  

    gel_stiffness=Stiffness_1(para);

    kecm = (1/crslnkDens^(1/3))*...
        gel_stiffness/1e12%.....................ECM spring constant (N/um)
  


   

%=========================================================================%
% =============================================================================
% 6. STOCHASTIC SAMPLE SWEEP
% =============================================================================
% Repeat the same biological condition for multiple stochastic realizations.
% Each sample experiences independently generated local fibers, adhesion
% events, pseudopod choices, and migration trajectories.
%
%Sample Sweep - repeats simulation for number of given samples
%=========================================================================%
    for samp=1:samples

        ma=0;
        countlen=1;

        cellPosClu = [0     0     0;
                      12.31 7.5   7.5;
                      12.31 -7.5  0;
                      12.31 7.5  -7.5;
                      25   -7.5   -7.5;
                      25   7.5  0;
                      25   -7.5  7.5;
                      38 7.5   7.5;
                      38 -7.5  0;
                      38 7.5  -7.5];
        cellPosClu=zeros(samp,3);

        %Generates cell's initial position/polarity
       
%             cellPos(1,:,cell) = oriPos*(2*rand(1,3)-1);%Random init
%             cellPos(1,:,cell) = [cell*15 rand/2 rand/2];%Linear init
        %Cluster init

        disp([samp par])

        % if crslnkPerF < 3.5
        %     gel_stiffness = 1039.9*crslnkPerF-1992.9
        % else
        %     gel_stiffness = 5247.9*crslnkPerF-16274
        % end%.................................................Gel stiffness (Pa)  

         gel_stiffness=Stiffness_1(para);


        kecm = (1/crslnkDens^(1/3))*...
            gel_stiffness/1e12;%.....................ECM spring constant (N/um)
        
     
        
        %if samp > samples/2
           % polarize = 0.5;
        %else 
            %polarize = 0.2;
        %end
        
        %=================================================================%
        % =============================================================================
% 7. INITIALIZE CELL POSITION AND MIGRATION PHASE
% =============================================================================
% cellPos stores the 3D cell trajectory.
% phaseRT records the migration phase at each time point:
%       1 = retracting
%       2 = pseudopod outgrowth
%       3 = contraction / forward movement
%
%Position and Phase Matrices
        %=================================================================%
        cellPos = zeros(1,3);%.........................Cell position matrix
        phaseRT = zeros;%....................Pseudopod phase runtime matrix
        pseudoVect = zeros(1,3);%.................Pseudopod directoin vector
        
        %=================================================================%
        %Tracking variables
        %=================================================================%
        Lpseudo = 0;%............................Tracks length of pseudopod
        newAngle = 1;%.....................Tracks when cell finds new fiber
        retracting = 1;%............................Tracks retracting phase
        outgrowth = 0;%..............................Tracks outgrowth phase
        contracting = 0;%..........................Tracks contracting phase
        t = 0;%........................Tracks pseudopod extension frequency
        pos = 1;%....................Tracks cell position at each time step
        
        %=================================================================%
        %Temporary Variables
        %=================================================================%
        fiberDir = 0;%.........................Unit vector of guiding fiber
        searchTime = 0;%...........................Pseudopod extension time
        numElmts = 0;%...........................Discrete elements of fiber
        Bsites = 0;%..............................Distributed binding sites
        numSearch = 0;%.............Number of elements searched per step
        bonds = 0;%.........................................Filopodia bonds
        Bpseudo = 0;%...................Number bonds along entire pseudopod
        dMove = 0;%......................Distance cell moves in contraction
        localFibers = 0;%....................Number of fibers touching cell
        cell = 0;%.............................................Axes of cell
        dR = 0;%...........................................Receptor density
        polAngle = 0;........................................Polarity angle
        dl = 0;
        a = 0;
        b = 0;
        
        %========================================fig=========================%
        %Generates Cell's Initial Position (+/- 0.5um from origin)
        %=================================================================%
        %cellPos(pos,:) = (2*rand(1,3)-1);%...................Initial position
        cellPos(1,:) = cellPosClu(samp,:);
        %Calculates Fiber's Angle from Reference Vector
        if sigma == 1e10
            fiberAngle = (2*rand-1)*360;%...........Randomly aligned fibers
        else
            fiberAngle = (180*(randi(2)-1))+...
                normrnd(otAngle,sigma);%.....................Aligned fibers
        end
        
        %Sets Initial Dinstance to Fiber Inersection
        Ltemp  = exprnd(1/crslnkPerL);
        if Ltemp <= Lsearch
            Ltemp = Lsearch;%............Limited to length of pseudopod tip
        end
        
        %Sets Initial Polarity
        vPrev = 2*rand(1,3)-1;  %cell polarity
        vPrev = vPrev/norm(vPrev);
        
        %=================================================================%
        % =============================================================================
% 8. TIME-RESOLVED 3D MIGRATION SIMULATION
% =============================================================================
% At each time step the cell senses nearby fibers, chooses a direction,
% forms integrin-mediated adhesions, extends or retracts a pseudopod, and
% contracts when sufficient adhesion strength is reached.
%
%Calculates path of cell
        %=================================================================%
        for i=dt:dt:runTime
            
            % figure(4)
            % hold on
            % plot3([0,10*ve(1,1)],[0,10*ve(1,2)],[0,10*ve(1,3)])
            % 
         

            
            % -------------------------------------------------------------------------
% LOCAL ECM SENSING AND FIBER SELECTION
% -------------------------------------------------------------------------
% Generate the local collagen fibers available to the cell. Fiber number,
% orientation, spacing, and crosslink distance depend on ECM properties.
%
%Generates new fiber direction
            if newAngle 

                %Generates local fibers in contact with cell
                localFibers = poissrnd(Vcell*fiberDens);
                numFibers = localFibers;
                if localFibers < 1 || outgrowth
                    localFibers = 1;
                end
                if numFibers < 1
                    numFibers = 1;
                end
                
                %Reset fiber, angle, and length matrices
                fiber = zeros(1,3);%.................Fiber direction matrix
                fb = zeros(1,3);
                angle = zeros(1,1);%.....................Fiber angle matrix
                acAng = zeros(1,1);%...............Fiber acute angle matrix
                numSearch = round(Lsearch/elmtSize);%.....number of Elements
                L = zeros(localFibers,1);%.........Distance to intersection
                cs = zeros;
                sn = zeros;
                

                for j=1:localFibers
                       L(j) = exprnd(1/(crslnkDens)^(1/3));
                        if L(j) <= Lsearch
                            L(j)  = Lsearch;
                        end
           

                     

                 
                    if sigma == 1e10
                        newDir = 2*rand(1,3)-1;%........New fiber direction
                        fiber(j,:) = newDir(1,:)/norm(newDir);%..New fibers

                         %Calculates aligned fiber angle
                    else
                        % fiberAngle =generatedangles()*(randi([0, 1]) * 2 - 1);
                        % if fiberAngle<-90 
                        %     fiberAngle=fiberAngle+180;
                        % 
                        % end
                        % 
                        % if fiberAngle>90 
                        %      fiberAngle=fiberAngle-180;
                        % 
                        % 
                        % end
                        % 
                        % 
                        % v1 = L(j)*cosd(fiberAngle)* vRef;%'x' component of new vector
                        % vRand = 2*rand(1,3)-1;%Random vector to get new vector rotation
                        % xProd = (cross(vRef,vRand));%Gets new vector's rotation about ref vector
                        % v2 =(L(j)*sind(fiberAngle)*xProd)/norm(xProd);%'y' component of new vector
                        % newDir = (v1+v2);%............New fiber direction
                        % fiber(j,:) = newDir(1,:)/norm(newDir);%Unit vector for new fibers


                        fiberAngle = normrnd(otAngle,sigma);
                        v1 = (L(j)*cosd(fiberAngle)*vRef)/norm(vRef);
                        vRand = 2*rand(1,3)-1;
                        xProd = (cross(vRef,vRand));
                        v2 =(L(j)*sind(fiberAngle)*xProd)/norm(xProd);
                        newDir = (v1 + v2);%............New fiber direction
                        fiber(j,:) = newDir(1,:)/norm(newDir);%..New fibers
                    end        

                

                  

                    
                    %Get fiber vectors in direction of previous direction
                    fb(j,:) = fiber(j,:);
                    an = acosd(dot(vPrev,fb(j,:))/(norm(vPrev)*norm(fb(j,:))));
                    if an > 90
                        fb(j,:) =-fb(j,:);
                    end
                end %local fibers loop

                
                

                
                %Elongation Vector
                ve = sum(fb,1)+vPrev/norm(vPrev);
                
                

              
                
                % ve and prev all cell vectors
                %Cell follows angle between ve and prev direction
                veAngle = acosd(dot(vPrev,ve)/(norm(vPrev)*norm(ve)));
                if veAngle > 90
                    ve =-ve;
                end


                %Get ratio of a0 to b0
                fbp = vertcat(vPrev,fb);%...........All fibers cell touches
                for j=1:size(fbp,1)
                    cs(j,1) = dot(ve,fbp(j,:))/(norm(ve)*norm(fbp(j,:)));              
                    sn(j,1) = sind(acosd(dot(ve,fbp(j,:))/(norm(ve)*...
                        norm(fbp(j,:))))); 
                end
                a0 = sum(cs);   
                b0 = sum(sn);
                
                
                

                bondsTotal = localFibers*RGD*(0.7*(Dcell/2)/Ltropo)...
                    *(0.9*Fiber_D/Dtropo)*pi/2;%......Bonds at trailing end
                keq = (bondsTotal*kI*kecm)/(bondsTotal*kI+kecm);

                %Define axes of the cell
                b = ((3*Vcell*(kcell+keq))/(4*pi*((a0/b0)*keq+kcell)))^(1/3);
                c = b;
                a = (3*Vcell)/(4*pi*b^2);
                cell = [a,b,c];%........Cell axes 
               
                
                %Calculate cell aahicity
                A_cell = saellipsoid(cell);%..........Cell surface area 
                s = (pi^(1/3)*(6*Vcell)^(2/3))/A_cell;%.Cell sphericity
                dR = dRecp/A_cell;%...........Membrane integrin density
                polAngle = 180*s^2;......................Polarity angle

                %Fiber angles with respect to previous direction
                for j=1:localFibers
                    angle(j,1) = acosd(dot(fiber(j,:),ve)/...
                        (norm(fiber(j,:))*norm(ve)));  
                end
                
                %Gets acute angle for all fibers
                for j=1:size(angle,1)
                    if angle(j) > 90
                        acAng(j,1) = 180-angle(j);
                    else
                        acAng(j,1) = angle(j);
                    end
                end

                %acAng: angle between cell and fibers
               
                % -------------------------------------------------------------------------
% POLARITY-DEPENDENT DIRECTION SELECTION
% -------------------------------------------------------------------------
% Compare available fiber directions with the cell's current orientation.
% The polarity parameter controls how strongly the cell preserves its
% previous migration direction rather than reversing along a fiber.
%
%Determines which fiber to follow based on cell's polarity 
                
                
                 cell_pol = 90*normrnd(0,s);
                 %polanglesave(i,samp)=90*normrnd(0,s);
                



                
                
                ang = abs(acAng-cell_pol);   
               
                
                min_idx = ang == min(ang);
                l=0;
                if sum(min_idx)>1 
                   for im=1:length(min_idx)
                       if min_idx(im)==1  && l==0
                           l=l+1;
                           min_idx(im)=1;
                       else
                          min_idx(im)=0;
                       end
                   end
                    
                end

                min_idx;

                



               
                
                %Extension vector
                fiberDir = fiber(min_idx,:);

                
 
                %Limits reversal of polarity if new pseudopod
                   if retracting && angle(min_idx) > 90
                            if rand <= polarize
                                fiberDir = -1*fiberDir;
                            end

                   elseif retracting && angle(min_idx) <90
                            if rand >= polarize
                                fiberDir= -1*fiberDir;
                            end

                        %Determines which way cell moves at crosslink
                   elseif outgrowth && angle(min_idx) > 90 
                            if rand <= 0.99
                                fiberDir= -1*fiberDir;
                            end
                   elseif outgrowth && angle(min_idx) < 90 
                            if rand >= 0.9
                                fiberDir= -1*fiberDir;
                            end
                   end

                

                
  

                
                %Gets distance to next crosslink
                crossL = L(min_idx);%............Distance to next crosslink
                %g_crossL(i,samp)=crossL;
                numElmts = floor(crossL/elmtSize);%......Elements to search 
                
                
                % -------------------------------------------------------------------------
% INTEGRIN-LIGAND ADHESION FORMATION
% -------------------------------------------------------------------------
% Estimate the number of available adhesion sites along the selected fiber.
% Bond formation depends on receptor/ligand density and the kinetic rates
% kon and koff.
%
%Distributes # of adhesion sites in each element - assumes 
                %cell touches half of fiber surface
                Fdm = normrnd(Fiber_D,0.02);%......Randomize fiber diameter
                if (Fdm <= 0)
                    Fdm = 0.01;
                end
                dmax = min([dL,dR]);%.Greater of ligand or receptor density
                Pb = ((kon*dmax)/(kon*dmax+koff));
                    %(1-exp(-(kon*dmax+koff)*tC))%......Binding probability 
                Bsites = poissrnd(Pb*RGD*0.7*0.9*(elmtSize/(Ltropo+.067))*(Fdm/Dtropo)*...
                    (pi/2),numElmts,1);%..........Distributed binding sites
                Bpseudo = Bpseudo + sum(Bsites(1:numSearch-1));
                
                %Sets time until new pseudopod extends
                if (retracting)
                    searchTime = exprnd(avgSearchTime);%.....Extension time        
                    if searchTime > 2*avgSearchTime
                        searchTime = 2*avgSearchTime;%..........Extension limit
                    end
                end
                
                %Reset values
                newAngle = 0;
                retracting = 0;
                outgrowth = 0;
            end    %newangle loop

            
            %Calculates bonds formed for each time step
            if not(contracting)
                dl = 0;  
                while(numSearch < numElmts && dl < dt*Vpseudo)              
                    numSearch = numSearch+1;
                    dl = dl+elmtSize;
                    bonds =sum(Bsites((numSearch-(round(Lsearch/elmtSize)-1)):numSearch));
                    Bpseudo = Bpseudo+Bsites(numSearch);
                  
                end
            end
            
            %=============================================================%
            % =========================================================================
% 9. MIGRATION STATE TRANSITIONS
% =========================================================================
% The cell switches between three migration states according to the number
% of bonds formed at the pseudopod.
%
% RETRACTING PHASE - pseudopod will retract if weak bonds are 
            %formed, or if new pseudopod starts extending
            %=============================================================%
            if bonds < Bmin || (bonds < Bmax && t > searchTime)
                retracting = 1;
                phaseRT(pos,1) = 1;
                
                %Resets values
                Lpseudo = 0;.........................Reset pseudopod length
                t = 0;%................................Reset time searching
                bonds = 0;%.................Reset bonds at tip of pseudopod
                Bpseudo = 0;%...................Reset bonds along pseudopod
                pseudoVect = zeros(1,3);%.........Reset pseudopod direction
                newAngle = 1;%.............................Select new fiber
                cellPos(pos+1,:) = cellPos(pos,:);%.....Store cell position
                
            %=============================================================%
            % -------------------------------------------------------------------------
% OUTGROWTH PHASE
% -------------------------------------------------------------------------
% Intermediate adhesion strength supports actin-driven pseudopod extension,
% but is not yet sufficient to generate productive contraction.
%
%OUTGROWTH PHASE - pseudopod will continue growing as long as 
            %enough stable bonds form to allow for actin polymerization
            %=============================================================%
            elseif bonds >= Bmin && bonds < Bmax


     

                outgrowth = 1;
                phaseRT(pos,1) =2;

                
                %Tracks distance and time pseudopod extends along 
                %current fiber
                pseudoVect = pseudoVect+dl*fiberDir;
                Lpseudo = norm(pseudoVect);%............Length of pseudopod
                t = t+dt;%...................Increment time spent searching 
                
                %Get new fiber if cell reaches cross-link
                if numSearch >= numElmts
                    newAngle = 1;
                end
                cellPos(pos+1,:) = cellPos(pos,:);%.....Store cell position
                
            %=============================================================%
            % -------------------------------------------------------------------------
% CONTRACTING PHASE
% -------------------------------------------------------------------------
% Once enough bonds are formed, adhesions can resist actomyosin contractile
% force. The pseudopod then shortens and pulls the cell body forward.
%
%CONTRACTING PHASE - pseudopod will contract when enough bonds 
            %are formed to withstand the acto-myosin contractile force
            %=============================================================%
            
            elseif bonds >= Bmax 
                phaseRT(pos,1) = 3;

                if rand>=0.99
                        newAngle=1; %give random
                end

                %t=t+dt;
                
                contractilitytime = contractilitytime+1;
               

                

                
                % -------------------------------------------------------------------------
% FORCE BALANCE AND CELL DISPLACEMENT
% -------------------------------------------------------------------------
% Calculate adhesion force, bond friction, viscous drag, instantaneous
% migration velocity, and the resulting cell displacement.
%
%Calculates final pseudopod direction, magnitude, 
                %and contractile force
                if not(contracting)
                    %countlen=countlen+1;
                   
                    
                    contracting = 1;%.................Set contracting phase
                    Ltouching = mean(cell)/4;%Length of fibers touching cell
                    bondsTotal = localFibers*RGD*(0.7*Ltouching/Ltropo)...
                        *(0.9*Fiber_D/Dtropo)*pi/2;%Bonds at trailing end
                    
                    pseudoVect = pseudoVect+dl*fiberDir;
                    Lpseudo = norm(pseudoVect);%..Final length of pseudopod
                    
                    if Lpseudo > 1
                        countnumber = countnumber + 1;%for number of protrusions
                        gLpseudo(countnumber,samp)= norm(pseudoVect);

                    end
                    
                    
                    

                    
                    
                    
                 
                    F0 = (Fmax*Bpseudo)/(Bhalf+Bpseudo);%....Adhesion force

                    kT_eff = (F0*kecm^2*Lpseudo^3)/...
                        (F0+kecm*Lpseudo)^2;%...Effective equivalent of kBT
                    fb = bondsTotal*((kecm*kI)/((kecm+kI)*koff))...
                        *exp(-(kT_eff)/(bondsTotal*kB*T));%Bond friction 
                                   
                    if a < b
                        beta = b/a;%.........................Inverse sphericity
                        Kprime = ((4/3)*(beta^2-1))/(((2*beta^2-1)/...
                            sqrt(beta^2-1))*reallog(beta+sqrt(beta^2-1))-...
                            beta);%.....................Drag adjustment factor
                    else
                        beta = a/b;
                        Kprime = ((4/3)*(beta^2-1))/((beta*(beta^2-2)/...
                            sqrt(beta^2-1))*atan(sqrt(beta^2-1))+...
                            beta);%.....................Drag adjustment factor
                    end
                    
                    fv = 6*pi*eta*a*Kprime;%...............Viscous friction
                    if ~isreal(fv)
                        disp(fv)
                    end
                        
                    Vinst = (F0*(F0+kecm*Lpseudo)-F0^2*...
                       reallog(F0+kecm*Lpseudo))/((fb+fv)*kecm*Lpseudo)-...
                       (F0*(F0+kecm*0)-F0^2*reallog(F0+kecm*0))/...
                       ((fb+fv)*kecm*Lpseudo);%......Insantaneous velocity
                    dMove = Vinst*dt;%.................Distance moved in dt   
                    
                end

                

                
                dMove = min([Lpseudo,dMove]);%.........Eliminates overshoot
                
                Lpseudo = Lpseudo-dMove;%..Decrement Lspeudo each time step                
                %Stores new cell position and direction
                cellPos(pos+1,:) = dMove*(pseudoVect./norm(pseudoVect))+...
                    cellPos(pos,:);
                vPrev = cellPos(pos+1,:)-cellPos(pos,:);

                
               

                
                
                
                %Resets values after cell is done contracting
                if Lpseudo <= 0
                    
               
                    searchTime = exprnd(avgSearchTime);
                    if searchTime > 2*avgSearchTime
                        searchTime = 2*avgSearchTime;
                    end
                    Lpseudo = 0;
                    t = 0;
                    bonds = 0;
                    Bpseudo = 0;
                    pseudoVect = zeros(1,3);
                    contracting = 0;
                    contractilitytime=0;
                    
                end

                
            end
            %=============================================================%  

            
            
            cellAx(i,:,samp) = [a(1),b(1),c(1)];%
            gFiberDir(i,:,samp) =fiber(1,:);
            gPV(i,:,samp) = pseudoVect(1,:);
            gVe(i,:,samp) = ve(1,:);
            gphaseRT(i,samp)=phaseRT(i);

            % fiber_pseudo(i,1,samp)=Lpseudo;
            % fiber_pseudo(i,2,samp)=crossL;
            % fiber_pseudo(i,3,samp)=phaseRT(i);
            % fiber_pseudo(i,4,samp)=dMove;
            % fiber_pseudo(i,5,samp)=newAngle;
            % fiber_pseudo(i,6,samp)=numSearch;
            % fiber_pseudo(i,7,samp)=numElmts;
            % fiber_pseudo(i,8,samp)=bonds;
            % fiber_pseudo(i,9,samp)=angle(min_idx);
            % fiber_pseudo(i,10,samp)=counter;
            % fiber_pseudo(i,11,samp)=Vinst;
            % fiber_pseudo(i,12,samp)=norm(cellPos(pos+1,:)-cellPos(pos,:));
            
            pos = pos+1;

            %prevve(i,:)=ve;
            degredation=0;

            
        end%................... Runtime  ..............................End simulation
        rw = cellPos;

        %figure(1)
        % h = surfl(x1, y1, z1); 
        % set(h, 'FaceAlpha', 0.2)
        % shading interp
        %plot3(rw(:,1),rw(:,2),rw(:,3), linewidth=3)
       
        randomWalk(:,:,par,samp) = rw;
        % set(gcf,"color","w")
        % set(gca,'fontsize',14)
        %axis([-15 15 -15 15 -15 15])
        %axis([-30 30 -30 30 -30 30])
        %axis([-100 100 -100 100 -100 100])
        %hold on
        %view(3);
        

        cellinvasive(par,samp)=norm(cellPos(time(end),:)-cellPos(1,:));

        %avegLpseudo(samp)=mean(nonzeros(gLpseudo(:,samp)))
        %Lpseudonumber(samp)=length(nonzeros(gLpseudo(:,samp)))
        %avecrossL(samp)=mean(nonzeros(g_crossL(:,samp)));
       
        %=================================================================%
        % =============================================================================
% 10. TRAJECTORY-BASED MIGRATION METRICS
% =============================================================================
% Quantify the simulated trajectory after each stochastic realization.
% Metrics include net invasion distance, average migration velocity, and
% persistence length.
%
% Cell Velocity
        %=================================================================%
        w = 1;
        tm = (1:runTime/10:runTime)./60;%..10 time points to get distance
        
        dst = zeros;
        pos = 1;
        for q=1:(size(rw,1))/10:size(rw,1)%........Gets distance at time tm
            Mag = zeros;
            for j=1:q
                if cellPos(j+1,:) ~= cellPos(j,:)
                    Mag(pos,1) = norm(cellPos(j+1,:)-cellPos(j,:));
                    pos = pos+1;
                end
            end
            dst(w,1) = sum(Mag);%.Gets distance traveled at each time point
            w = w+1;
        end
        try%....................Curve fit to get slope for average velocity
            [f,g] = fit(tm',dst,'poly1');
            c = coeffvalues(f);
            avgVel(par,samp) = c(1);%............................Cell Speed
            R1(par,samp) = g.rsquare;%......................Goodness of fit
        catch
            disp('Velocity Fitting Error')
            avgVel(par,samp) = 0;
            R1(par,samp) = 0;
        end
        
        
        %=================================================================%
        % -------------------------------------------------------------------------
% PERSISTENCE LENGTH
% -------------------------------------------------------------------------
% Estimate how long the cell maintains directional memory along its path.
% The persistence length is obtained by fitting mean squared displacement
% along contour length to a persistent-random-walk relation.
%
%Persistence Length 
        %=================================================================%
        %Creates indices and variables
        nPos = zeros;
        nVect = zeros;
        nMag = zeros;
        avgCos = zeros;
        avgRsquared = zeros;
        contourLength = zeros;
        pos = 1;


        %Position, vector, and magnitude for each change in position
        for j=1:size(cellPos,1)-1
            if cellPos(j+1,:) ~= cellPos(j,:)
                nPos(pos,1:3) = cellPos(j,:);
                nVect(pos,1:3) = cellPos(j+1,:)-cellPos(j,:);
                nMag(pos,1) = norm(cellPos(j+1,:)-cellPos(j,:));
                pos = pos+1;
            end
        end
        nPos(pos,1:3) = cellPos(size(cellPos,1),1:3);
        cMag(par,samp) = mean(nMag);
        cellDist(par,samp) = sum(nMag);
        steps(par,samp) = runTime/dt;

        xl =floor(sum(nMag)/4); %floor(sum(nMag)/4);%.....Contour length limit for Lp plots
        n = 0;%.........................................Tracks steps for L0

        %Loop for incresing L0
        for L0=1:xl/300:xl

            n = n + 1;%.............................Increments each L0 step
            endPos = 1;%.............................First direction vector
            m = 0;%....................Tracks cos(theta)/R for each L0 step
            stop = 0;%Breaks while loop if pos equals size of vector matrix
            totalCos = 0;
            totalR = 0;

            %Calculates average cos(theta)/R^2 for given L0
            while endPos < size(nVect,1)

                Length = 0;%........Resets L after moving step length of L0
                m = m + 1;%...............Increments steps for cos(theta)/R
                strt_pos = endPos;%....Gets starting position of L0 contour

                %Gets the contour position at end of length L0
                while Length <= L0
                    %Breaks loop at end of contour 
                    if endPos == size(nVect,1)
                        stop = 1;
                        break;
                    end

                    Length = Length + nMag(endPos);%......Length of contour
                    endPos = endPos + 1;%...Gets end position of L0 contour
                end

                if stop == 0
                    %cos(theta) between vectors at start and end of contour
%                     cos_theta = dot(nVect(endPos-1,:),nVect(strt_pos,:))/...
%                         (nMag(endPos-1)*nMag(strt_pos));
%                     totalCos = totalCos + cos_theta;

                    %Distance R between start and end points of contour
                    Rsquared = (norm(nPos(endPos,:) - nPos(strt_pos,:)))^2;
                    totalR = totalR + Rsquared;
                end
            end

            %avgCosStep = totalCos/m;
            avgRstep = totalR/m;

            %avgCos(n,1) = avgCosStep;%...................Average cos(theta) 
            avgRsquared(n,1) = avgRstep;%.................Average R squared
            contourLength(n,1) = L0;%.......................Contour lengths
        end

        ft = fittype('2*(b1^2)*(exp(-(x/b1))-1+x/b1)');
        try%............................Curve fit to get persistence length
            [f,g] = fit(contourLength,avgRsquared,ft,'StartPoint',1);
            c = coeffvalues(f);
            LpR(par,samp) = c(1);%.......................Persistence length
            R3(par,samp) = g.rsquare;%......................Goodness of fit
        catch
            LpR(par,samp) = 0;
            R3(par,samp) = 0;
        end

        retr(par,samp) = sum(phaseRT == 1)
        outg(par,samp) = sum(phaseRT == 2)
        cont(par,samp) = sum(phaseRT == 3)
    end %........................................................End for sample loop

end %end of par loop





%% 

head = gobjects(samples,1);
vHead = gobjects(samples,1);
curve = gobjects(samples,1);

for samp=1:samples
    curve(samp) = animatedline('Color',abs(rand(3,1)),'LineWidth',2);
end

% =============================================================================
% 11. 3D TRAJECTORY VISUALIZATION
% =============================================================================
% Animate the simulated cells and their trajectories.
%
% Cell colors indicate migration state:
%       Red    = retracting
%       Yellow = pseudopod outgrowth
%       Green  = contraction / productive movement
%
% 3D cell animation
%figure(2)




grid on
%daspect([0 0 0])
xlabel('x[\mum]')
ylabel('y[\mum]')
zlabel('z[\mum]')
%axis([-20 20 -20 20 -20 20])
%axis([-50 50 -50 50 -50 50])
axis([-100 100 -100 100 -100 100])
%axis([-15 15 -15 15 -15 15])

view(3);
hold on

animation=1;
if animation == 1
    myVideo = VideoWriter(filename,"MPEG-4");
    myVideo.FrameRate =80;
    open(myVideo);

    for j=1:120:runTime-1 %3D Animation plot
        %for samp=[4,7,15,19]
        %for samp=[4,7,15,20]
        for samp=1:samples
        %for samp=1
            if gphaseRT(j,samp) == 1
                C = [1 0 0];%...........................................Red
            elseif gphaseRT(j,samp) == 2
                C = [1 1 0];%........................................Yellow
            elseif gphaseRT(j,samp) == 3
                C = [0 1 0];%.........................................Green
            end
            [x,y,z] = ellipsoid(0,0,0,cellAx(j,3,samp),cellAx(j,2,samp),cellAx(j,1,samp),15);
            set(gcf,"color","w")

            addpoints(curve(samp),randomWalk(j,1,samp),randomWalk(j,2,samp),randomWalk(j,3,samp))
            head(samp) = surf(x+randomWalk(j,1,samp),y+randomWalk(j,2,samp),z+randomWalk(j,3,samp),'FaceColor',C);
            %pause(0.1)
            if gphaseRT(j,samp) == 3 || j == 1
                if j == 1
                    if norm(gPV(j,:,samp)) ~= 0
                        %pV2(samp,:) = gpd(j,:,samp);
                        pV2(samp,:) = gPV(j,:,samp);
                    else
                        pV2(samp,:) = gPV(j,:,samp);
                        pV2(samp,:) = gFiberDir(j,:,samp);
                        pV2(samp,:) = gVe(j,:,samp);
                    end
                end
                %rotate(head(samp),[0 0 1],atand(pV2(samp,2)/pV2(samp,1)),randomWalk(j,:,samp))
                rotate(head(samp),cross([0 0 1],pV2(samp,:)),acosd(pV2(samp,3)/norm(pV2(samp,:))),randomWalk(j,:,samp))
                alpha(head(samp),0.25);
            end
        end
        frame = getframe(gcf); %get frame
        writeVideo(myVideo, frame);
        drawnow
        %delete(vHead(:));
        delete(head(:));

    end   
    close(myVideo)
end

% =============================================================================
% 12. SAVE SIMULATION OUTPUTS
% =============================================================================
% Save the final trajectory and migration metrics for downstream analysis.
%
% cellinvasive : net displacement from starting position
% randomWalk   : complete 3D trajectory
% MSD          : mean squared displacement values
% LpR          : persistence length
% avgVel       : average migration velocity
%
%save(filename,'avgVel','cellDist','v1','v2','LpR','p1','p2','D','d1','d2','alp','a1','a2','R1','R2','R3','outg','cont','retr')
save(filename(1),'cellinvasive','randomWalk','MSD',"LpR","avgVel")
%save(filename(1),'randomWalk',"ave_cellinv")
%delete(gcp('nocreate'))