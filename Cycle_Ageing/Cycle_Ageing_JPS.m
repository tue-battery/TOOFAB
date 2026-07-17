%% Test setup
load Data.mat
load DFNparams.mat
%%
p=p_Ageing; ph=ph_Ageing;
p.Cbat=Equil2.Cbat;
    
    %Setting ageing 
    p.ageing = 1;
    p.SaveDataTskip = 1;
    p.verbose=2; %2 shows cycling progress, 0 shows no output
    p.set_simp = [1,1,1,1,1,0]; %model simplifications 
    p.ageingMech.Oxidation = 1; %cathode oxidation
    p.ageingMech.SEI = 1; %SEI formation
    p.ageingMech.LPL = 1; %lithium plating
    p.LimitedSEI=1; %SEI limitation
    
    %Setting Ageing Parameters
    %SEI formation
    p.i02 = 0.016;
    p.B1_s=25;
    p.B2_s=5;
    p.lamda_sei = 15;
    p.alpha_c2=0.6;
    %lithium plating
    p.i0lpl=4.1e-4;
    p.B1_l=15;
    p.B1_l=1.5;
    p.alphaa_lpl=0.3;
    p.alphac_lpl=0.7;
    %Cathode Oxidation
    p.i0C=7e-42;
    p.alpha_a2=0.5;
    p.BO=20;

    %Timestamp variables
    p.count=1; 
    p.DCCount=1;
    p.dt_charge =5;
    p.dt_discharge = 5;
    p.dt_CV = 1;
    p.dt=1;
    p.DriveCycle=Data;   
  
    %setting temperature dependent DFN parameters, based on previously
    %estimated grouped parameters. 
    p.Ds_neg=@(T)(ph.Ds_neg*exp((ph.Ds_neg_Ea/(8.3145))*((1/298)-(1/T))))*p.R_neg^2;
    p.Ds_pos=@(T)(ph.Ds_pos*exp((ph.Ds_pos_Ea/(8.3145))*((1/298)-(1/T))))*p.R_pos^2;
    p.De= @(T)(ph.De*exp((ph.DEa/(8.3145))*((1/298)-(1/T))))*(1-p.t_plus)/(p.F*p.A_surf*p.ce0);
    p.k0_neg= @(T)(ph.k0_neg*exp((ph.knEa/(8.3145))*((1/298)-(1/T))))*p.R_neg*p.F/p.ce0^p.alpha_a;
    p.k0_pos= @(T)(ph.k0_pos*exp((ph.kpEa/(8.3145))*((1/298)-(1/T))))*p.R_pos*p.F/p.ce0^p.alpha_a;
    p.TotCycle=301;

    p.VOffset=0;
    p.T_amb=273.15+20;
    %Capacity loss outputs can be found in the results Struct
    results = DFN_Ageing_JPS(@cycle2,2880000,1,p);



%% Cycle Function
function [i_app_next,mem,end_simulation,dt] = cycle2(k,i_app,V,mem,p)
V_max = 4.2;   %Voltage at which the battery switches to CV mode
V_min = 2.5;   %Voltage after which the battery stops discharging
i_min = 0.005;  %Current level to terminate the CV mode with
Kp1 =15;   %Proportional gain of the PI controller for the CV stage
Ki1 = 5;   %Integral gain for the PI controller
Kd1 = 0;
end_simulation=0; %flag to indicate whether the simulation should end after the iteration. Has to be set at every iteration

% Initialize the various parameters. Note how the mem struct is used to
% store the states of the PI controller 
if k==1
    mem.status = 0;
    mem.CheckUp=1;
    mem.e = 0;
    mem.e2 = 0;
    mem.DriveCycleCount=1; %necessary for ageing as we will apply the checkup drivecycle every 25/50 cycles
    mem.T_amb=p.T_amb;
    mem.count=1;
    mem.CCcount=1;
    mem.CheckUpCount=1;
end

% Condition to terminate the initial discharge stage and start the CC
% charge stage
if (V(k)<=V_min && mem.status==0) || mem.status==3
    
    if mem.status==3
        if (V(k)<=V_min) || mem.count-1>=p.DriveCycle{1}.t(end)
            mem.count=1;

                if length(mem.V_resp) < length(p.DriveCycle{1}.I)
                    mem.Cycles{mem.DriveCycleCount}.V=mem.V_resp;
                    mem.Cycles{mem.DriveCycleCount}.I=mem.i_resp;
                    mem.Cycles{mem.DriveCycleCount}.t=mem.t;
                else
                    mem.Cycles{mem.DriveCycleCount}.V=mem.V_resp(1:length(p.DriveCycle{1}.I));
                    mem.Cycles{mem.DriveCycleCount}.I=mem.i_resp(1:length(p.DriveCycle{1}.I));
                    mem.Cycles{mem.DriveCycleCount}.t=mem.t(1:length(p.DriveCycle{1}.I));
                end
                mem.V_resp=[];
                mem.i_resp=[];
                mem.t=[];
                mem.DriveCycleCount=mem.DriveCycleCount+1;
                mem.T_amb=p.T_amb;
                mem.status = 1;  


        end
    else
        if i_app(k)==-1
            mem.CCcount=1;
            mem.Checkup{mem.CheckUpCount}.I=mem.I_checkup;
            mem.Checkup{mem.CheckUpCount}.t=mem.t_checkup;
            mem.I_checkup=[];
            mem.t_checkup=[];
            mem.CheckUpCount=mem.CheckUpCount+1;
            mem.status = 1; 
        else
            mem.status = 1; 
        end
    end
end

% Condition to terminate the CC charge stage and start the CV stage
if V(k)>=4.2 &&mem.status==1
    mem.status = 2; 
    mem.CVcount=1;
end

%Condition to terminate the simulation
if i_app(k)<i_min && mem.status==2 && mem.CycleNo == mem.TotCycle
    end_simulation=1; 
end
% if  mem.CycleNo == mem.TotCycle
%     end_simulation=1; 
% end


%Condition to reset CC-CV Charge Cycle 
if i_app(k)<i_min && mem.status==2
% if mem.CVcount>=length(p.CV_I) && mem.status==2
    mem.Cbat_aged(mem.CycleNo) = mem.Cbat_aged_current;
    mem.CycleTimeStamp(mem.CycleNo) = k;
    mem.CycleNo = mem.CycleNo + 1;
    if  rem(mem.CycleNo,25)==0
        mem.CheckUp=1+mem.CycleNo/25;
        mem.status = 0; 
        mem.count=1;

    else
        mem.status = 0;
    end
    mem.e = 0;
    mem.e2 = 0;
end

%Switch between the different stages
switch mem.status
    case 0 %CC discharge
        % i_app_next = -p.Crate_discharge*Cbat; %value of i_app for the next iteration
        mem.T_amb=p.T_amb;
        mem.t_checkup(1)=1;
        if rem(mem.CycleNo,25)==0||mem.CycleNo==1
            if mem.CCcount>1
                mem.t_checkup(mem.CCcount)=mem.t_checkup(mem.CCcount-1)+5;
            end
            i_app_next=-1;
            mem.I_checkup(mem.CCcount)=i_app_next;
            mem.CCcount=mem.CCcount+1;
            
        else
            i_app_next = -5; %value of i_app for the next iteration
        end
        dt = p.dt_discharge; %assign the chosen step size
    
    case 1 %CC charge
        mem.T_amb=p.T_amb;
        i_app_next = 5; %value of i_app for the next iteration
        dt = p.dt_charge; %assign the chosen step size
      
    case 2 %CV
        %Solve the PI controller system (in this case with a forward Euler
        %discretization
        mem.T_amb=p.T_amb;
        mem.e(k+1) = V_max-V(k); 
        mem.e2(k+1) = (mem.e(k+1)-mem.e(k))/p.dt;
        mem.e3(k+1) = (mem.e2(k+1) - mem.e2(k))/p.dt;
        i_app_next = i_app(k)+(p.dt*(Kp1*mem.e2(k+1)+Ki1*mem.e(k+1) + Kd1*mem.e3(k+1)));
      
        dt = p.dt_CV; %assign the chosen step size

        
    case 3 %DriveCycle
        i_app_next=p.DriveCycle{1}.I(mem.count); 
        dt = 1;
        mem.V_resp(mem.count)=V(k);
        mem.i_resp(mem.count)=i_app(k);
        mem.t(mem.count)=mem.count;
        mem.T_amb=p.DriveCycle{1}.T(mem.count);
        mem.count=mem.count+1;
end
end
