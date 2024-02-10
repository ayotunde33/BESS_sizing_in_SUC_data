u=1;
test_case=5; %This select the test case that is being calculated
load_selx=load_sel_scen; %selected load profiles
wind_selx=wind_sel_scen; %selected wind profiles
solar_selx=solar_sel_scen; %selected solar profiles
%Note that you may have to manually edit the selected solar profiles as it
%may contain negative values which have to be set to zero as solar
%irradiation cannot have negative values. This problem exist in the
%generation of the profiles for solar irradiation. Just check the values in
%solar_sel_scen and use a for loop to set all negative values to zero
for ia=0:2.5:50 %ia is the BESS size    

%Get the Number of Hours to perform optimization operations
nHours = 24;
Time = (1:nHours);

%Enter the Number of Generating Units in the System
num_of_gen =4;

%Create Matrix to store generator data
gen_data = zeros(num_of_gen, 8);

minup=3; %generator minimum up time
mindown=3; %generator minimum down time
startup_cost=360; %startup cost
fuel_c1=-0.0073; %fuel cost coefficient 1
fuel_c2=103.7806; %fuel cost coefficient 2
fuel_c3=593.9206; %fuel cost coefficient 3
max_p=20.2; %maximum generator power
min_p=6.06; %minimum generator power
gen_data=[fuel_c1	fuel_c2	fuel_c3	max_p min_p	minup	mindown	startup_cost;
fuel_c1	fuel_c2	fuel_c3	max_p min_p	minup	mindown	startup_cost;
fuel_c1	fuel_c2	fuel_c3	max_p min_p	minup	mindown	startup_cost;
fuel_c1	fuel_c2	fuel_c3	max_p min_p	minup	mindown	startup_cost;
];

%Create matrix to store wind turbine and solar power for all the 50 selected scenarios
% Each scenarios has 24 elements for the 24 hours in a day, hence the
% dimension is 24 by 50
WT_power=zeros(24,50);
solar_power=zeros(24,50);
for l=1:50
Wind_speed=wind_selx(:,l); %wind speed for each scenario is selected

%wind turbine power in MW for each hour in a scenario is calculated
WT_power(:,l)=0.5*2*1.225*(23545*1)*0.3451*(Wind_speed.^3)/1e6;
 
%Create for loop for setting the wind turbine power to zero for values 
% outside the operating range (between cut-in and cut-off speed)
for i=1:nHours
    if Wind_speed(i)>=12 && Wind_speed(i)<=25 %rated speed is 12 m/s,
        % cut-off speed is 25 m/s
        WT_power(i,l)=8.6*2; %rated wind turbine power is 8.6 MW,
        % there are two wind turbines
    elseif Wind_speed(i)>25 || Wind_speed(i)<3 %cut-in speed is 3 m/s
        WT_power(i,l)=0;
    end
end

solar_irr=solar_selx(:,l); %solar irradiation for each scenario is selected

%create for loop to set negative values in the selected scenarios to zero,
%the scenario generation approach sometimes generate negative values for
%the solar irradiation and we know that solar irradiation cannot be
%negative so the negative entries have to be set to zero
for i=1:24
    if solar_irr(i)<0
        solar_irr(i)=0;
    end
end

n_pv=23007; %number of solar PV panels
inv_eff=0.96; %efficiency of inverter
PV_nom=inv_eff*390*n_pv/1e6;  %rated power of the solar PV array in MW

solar_std_rad=1000;  %standard solar irradiation in W/m2
NOCT=47; %Nominal operating cell temperature
k_temp=-0.0029; %power temperature coefficient
T_stc=25; %standard temperature
MPPT_eff=0.96;
T_amb=zeros(1,24);
%The for loop below calculate the average hourly ambient temperature
for i=1:8760
    ks=mod(i,24);
    if ks==0
        ks=24;
    end
    T_amb(ks)=T_amb(ks)+RES.Temp(i);
end
T_amb=T_amb./365;
%The for loop below calculate the solar PV array power for each hour
for i=1:24
    T_cell(i)=T_amb(i)+((NOCT-20)/800)*solar_irr(i);
solar_power(i,l)=PV_nom*(solar_irr(i)/solar_std_rad)*(1+k_temp*(T_cell(i)-T_stc))*MPPT_eff;
end
end
    
    BESS_nom=ia; %BESS_nom is BESS capacity
if BESS_nom==0
    BESS_nom=1e-6; %If BESS capacity is zero, set it to a very low value to
    % avoid error when the programme runs
end
%BESS_power1 is the optimization variable for BESS discharge power
%BESS_power2 is the optimization variable for BESS charge power
%BESS_cap is the optimization variable for BESS available energy (state of charge)
%BESS_isOn1 is the optimization variable for discharge status of BESS
%BESS_isOn2 is the optimization variable for charge status of BESS
%flex_load1 is the optimization variable for flexible load drawing more power
%flex_load2 is the optimization variable for flexible load drawing less power

BESS_power1 = optimvar('BESS_power1',nHours,50,'LowerBound',0,'UpperBound',BESS_nom*0.96);
BESS_power2 = optimvar('BESS_power2',nHours,50,'LowerBound',-BESS_nom,'UpperBound',0);
BESS_cap = optimvar('BESS_cap',nHours,50,'LowerBound',0,'UpperBound',BESS_nom);
BESS_isOn1 = optimvar('BESS_isOn1',nHours,50,'Type','integer','LowerBound',0,'UpperBound',1);
BESS_isOn2 = optimvar('BESS_isOn2',nHours,50,'Type','integer','LowerBound',0,'UpperBound',1);
flex_load1=optimvar('flex_load1',nHours,50,'LowerBound',0,'UpperBound',6);
flex_load2=optimvar('flex_load2',nHours,50,'LowerBound',-6,'UpperBound',0);
flex_isOn1=optimvar('flex_isOn1',nHours,50,'Type','integer','LowerBound',0,'UpperBound',1);
flex_isOn2=optimvar('flex_isOn2',nHours,50,'Type','integer','LowerBound',0,'UpperBound',1);

%The for loop below set the flexible load to 50% power for all hours for a
%6 MW water injection system on an offshore oil and gas platform

if test_case==3
wis_cap=6;
else 
wis_cap=3;
end

for i=1:24
    for j=1:50
       wis(i,j)=wis_cap;
    end
end

%Pull apart generator properties

%Get Start Up Cost data
StartupCost  = (gen_data(:, 8))';

%The for loop below calculates cost at minimum generation power
OperatingCost = zeros(1, num_of_gen); %OperatingCost is the cost at minimum generation power for each generator

for num = 1:1:num_of_gen
    %Operating Cost = a.(P_min)^2 + b.(P_min) + c, 
    OperatingCost(1, num) = (gen_data(num, 1) * gen_data(num, 5) * gen_data(num, 5)) + (gen_data(num, 2) * gen_data(num, 5)) + gen_data(num, 3);
end

%co2_c1, co2_c2 and co2_c3 are the carbon emission tax coefficients
co2_c1=-0.0325;
co2_c2=461.91;
co2_c3=2643.4;
co2_tax=0.069; %co2_tax is the carbon emission tax in $/kgco2

%The for loop below calculates carbon emission tax at minimum generation power
OperatingCO2Cost = zeros(1, num_of_gen); %OperatingCO2Cost is the carbon emission tax at minimum generation power for each generator

for num = 1:1:num_of_gen    
    %Operating CO2 Cost = x.(P_min)^2 + y.(P_min) + z
    OperatingCO2Cost(1, num) = (co2_c1 * gen_data(num, 5) * gen_data(num, 5)) + (co2_c2 * gen_data(num, 5)) + co2_c3;
end

%Get the Minimum up and down times data
MinimumUpTime = (gen_data(:, 6))';
MinimumDownTime = (gen_data(:, 7))';

%==========PIECEWISE LINEARIZATION OF COST FUNCTIONS==========

%Find Segment Width
%Specify the number of Segments
%Increasing the number of segments improves the results
K = 4; %This is the number of segments
W = zeros(num_of_gen, 1); %This is the width of each segment

for num = 1:1:num_of_gen    
    %W = (P_max - P_min)/K
    W(num, 1) = (gen_data(num, 4) - gen_data(num, 5)) / K;
end

%Get Maximum generation level for each segment of linearized cost curve
MaxGenerationLevel = zeros(1, (num_of_gen*K));

%The for loop below set the maximum generation level for each segment to
%the width of the segment
for num = 1:1:num_of_gen
    for k = 1:1:K
        MaxGenerationLevel(1, ((num - 1)*K)+k) = W(num, 1);
    end
end

%Get Maximum generation level for the actual units
MinPowerGenLevel = gen_data(:, 5)';

%Power_Segment is used to store the end point for each power segment
Power_Segment = zeros(num_of_gen, (K + 1));

%The for loop below find end point for each Power Segment
for num = 1:1:num_of_gen    
    for k = 1:1:(K + 1)
        if k == 1
            Power_Segment(num, k)  = gen_data(num, 5);
        elseif k == (K + 1)
            Power_Segment(num, k) = gen_data(num, 4);
        else
            Power_Segment(num, k) = gen_data(num, 5) + ((k - 1) * W(num, 1));
        end
    end
end

%Cost is a variable to hold the cost for each end point of each power segment
Cost = zeros(num_of_gen, (K + 1));

%The for loop below is used to find respective Cost for each end point of each Power Segment
for num = 1:1:num_of_gen
    for k = 1:1:(K + 1)
        Cost(num, k) = (gen_data(num, 1) * Power_Segment(num, k) * Power_Segment(num, k)) + (gen_data(num, 2) * Power_Segment(num, k)) + gen_data(num, 3);
    end
end

%CO2_Cost is a variable to hold the carbon emission tax for each end point of each power segment
CO2_Cost = zeros(num_of_gen, (K + 1));

%The for loop below is used to find respective carbon emission tax for each end point of each Power Segment
for num = 1:1:num_of_gen
    for k = 1:1:(K + 1)
        CO2_Cost(num, k) = (co2_c1 * Power_Segment(num, k) * Power_Segment(num, k)) + (co2_c2 * Power_Segment(num, k)) + co2_c3;
    end
end

%Slope holds the respective fuel cost Slope for each Power Segment
Slope = zeros(num_of_gen, K);

%The for loop below is used to calculate the fuel cost slope for each power segment 
for num = 1:1:num_of_gen
    for k = 1:1:K
        Slope(num, k) = (Cost(num, (k+1)) - Cost(num, k)) / W(num, 1);
    end
end

%CO2_Slope holds the respective carbon emission Slope for each Power Segment
CO2_Slope = zeros(num_of_gen, K);

%The for loop below is used to calculate the carbon emission tax slope for each power segment
for num = 1:1:num_of_gen
    for k = 1:1:K
        CO2_Slope(num, k) = (CO2_Cost(num, (k+1)) - CO2_Cost(num, k)) / W(num, 1);
    end
end

%===================END of LINEARIZATION======================

%======FORMULATION FOR USE OF 'intlinprog' solver FOR OPTIMIZATION====

%f is used to hold the fuel cost slope for all power segments and all generators
f = zeros((num_of_gen * K), 1);

%Create sub matrices to easily formulate f
sub_f = cell(1, 1);

%The for loop below is used to assign fuel cost slope values to the sub matrices
for num = 1:1:num_of_gen    
    for k = 1:1:K
            sub_f{num, 1}(k, 1) = Slope(num, k);
    end
end

%Convert cell to matrix to formulate f
f = cell2mat(sub_f);

%CO2_f is used to hold the carbon emission tax slope for all power segments and all generators
CO2_f = zeros((num_of_gen * K), 1);

%Create sub matrices to easily formulate CO2_f
CO2_sub_f = cell(1, 1);

%The for loop below is used to assign carbon emission tax slope values to the sub matrices
for num = 1:1:num_of_gen   
    for k = 1:1:K
            CO2_sub_f{num, 1}(k, 1) = CO2_Slope(num, k);
    end
end

%Convert cell to matrix to formulate CO2_f
CO2_f = cell2mat(CO2_sub_f);

%Create maxGenConst, minGenConst and minPowConst to be used in maximum
%generation constraint, minimum generation constraint and minimum power
%constraint later. maxGenConst is the maximum power in each power segment
%for each generator, minGenConst is the minimum power in each power segment
%for each generator which is zero and minPowConst is the minimum power for
%each generator and is given in the generator data which is 6.06 MW

idxHr2ToEnd=2:nHours; %this is an index used in a constraint later
maxGenConst = repmat(MaxGenerationLevel,nHours,50);
minPowConst = repmat(MinPowerGenLevel,nHours,50);

%Create an optimization problem
powerprob = optimproblem;

%Create Optimization Variables
L=48; %This is the number of segments for the battery degradation cost piecewise linearization
%Amount of power generated in an hour by a plant
power = optimvar('power',nHours,800,'LowerBound',0,'UpperBound',maxGenConst);
%Indicator if power segment is operating during an hour 
isOn = optimvar('isOn',nHours,800,'Type','integer','LowerBound',0,'UpperBound',1);
%Indicator if plant is starting up during an hour
startup = optimvar('startup',nHours,200,'Type','integer','LowerBound',0,'UpperBound',1);
%This is an optimization variable to store the amount of discharge power the BESS produces
BESS_thru = optimvar('BESS_thru',L,50,'LowerBound',0,'UpperBound',BESS_nom*12/L);
%This is an indicator if the value in a segment of BESS_thru is greater than zero
isOn_thru = optimvar('isOn_thru',L,50,'Type','integer','LowerBound',0,'UpperBound',1);
powerCO2=optimvar('powerCO2',1,1,'LowerBound',0,'UpperBound',1e7);
isOnCO2=optimvar('isOnCO2',1,1,'LowerBound',0,'UpperBound',1e7);
%Create an optimization variable new_maxGenConst to be used in a constraint later
new_maxGenConst = maxGenConst.*isOn;
isOn2=repmat(isOn(:,1),1,4);
isOn2=[isOn2 repmat(isOn(:,5),1,4)];
isOn2=[isOn2 repmat(isOn(:,9),1,4)];
isOn2=[isOn2 repmat(isOn(:,13),1,4)];
isOn2=[isOn2 repmat(isOn(:,17),1,4)];
isOn2=[isOn2 repmat(isOn(:,21),1,4)];
isOn2=[isOn2 repmat(isOn(:,25),1,4)];
isOn2=[isOn2 repmat(isOn(:,29),1,4)];
isOn2=[isOn2 repmat(isOn(:,33),1,4)];
isOn2=[isOn2 repmat(isOn(:,37),1,4)];
isOn2=[isOn2 repmat(isOn(:,41),1,4)];
isOn2=[isOn2 repmat(isOn(:,45),1,4)];
isOn2=[isOn2 repmat(isOn(:,49),1,4)];
isOn2=[isOn2 repmat(isOn(:,53),1,4)];
isOn2=[isOn2 repmat(isOn(:,57),1,4)];
isOn2=[isOn2 repmat(isOn(:,61),1,4)];
isOn2=[isOn2 repmat(isOn(:,65),1,4)];
isOn2=[isOn2 repmat(isOn(:,69),1,4)];
isOn2=[isOn2 repmat(isOn(:,73),1,4)];
isOn2=[isOn2 repmat(isOn(:,77),1,4)];
isOn2=[isOn2 repmat(isOn(:,81),1,4)];
isOn2=[isOn2 repmat(isOn(:,85),1,4)];
isOn2=[isOn2 repmat(isOn(:,89),1,4)];
isOn2=[isOn2 repmat(isOn(:,93),1,4)];
isOn2=[isOn2 repmat(isOn(:,97),1,4)];
isOn2=[isOn2 repmat(isOn(:,101),1,4)];
isOn2=[isOn2 repmat(isOn(:,105),1,4)];
isOn2=[isOn2 repmat(isOn(:,109),1,4)];
isOn2=[isOn2 repmat(isOn(:,113),1,4)];
isOn2=[isOn2 repmat(isOn(:,117),1,4)];
isOn2=[isOn2 repmat(isOn(:,121),1,4)];
isOn2=[isOn2 repmat(isOn(:,125),1,4)];
isOn2=[isOn2 repmat(isOn(:,129),1,4)];
isOn2=[isOn2 repmat(isOn(:,133),1,4)];
isOn2=[isOn2 repmat(isOn(:,137),1,4)];
isOn2=[isOn2 repmat(isOn(:,141),1,4)];
isOn2=[isOn2 repmat(isOn(:,145),1,4)];
isOn2=[isOn2 repmat(isOn(:,149),1,4)];
isOn2=[isOn2 repmat(isOn(:,153),1,4)];
isOn2=[isOn2 repmat(isOn(:,157),1,4)];
isOn2=[isOn2 repmat(isOn(:,161),1,4)];
isOn2=[isOn2 repmat(isOn(:,165),1,4)];
isOn2=[isOn2 repmat(isOn(:,169),1,4)];
isOn2=[isOn2 repmat(isOn(:,173),1,4)];
isOn2=[isOn2 repmat(isOn(:,177),1,4)];
isOn2=[isOn2 repmat(isOn(:,181),1,4)];
isOn2=[isOn2 repmat(isOn(:,185),1,4)];
isOn2=[isOn2 repmat(isOn(:,189),1,4)];
isOn2=[isOn2 repmat(isOn(:,193),1,4)];
isOn2=[isOn2 repmat(isOn(:,197),1,4)];

%Create an optimization variable isOn1 to hold the status of each generator in each hour
isOn1=[isOn(:,1) isOn(:,5) isOn(:,9) isOn(:,13) isOn(:,17) isOn(:,21) isOn(:,25) isOn(:,29)...
    isOn(:,33) isOn(:,37) isOn(:,41) isOn(:,45) isOn(:,49) isOn(:,53) isOn(:,57) isOn(:,61)...
    isOn(:,65) isOn(:,69) isOn(:,73) isOn(:,77) isOn(:,81) isOn(:,85) isOn(:,89) isOn(:,93)...
    isOn(:,97) isOn(:,101) isOn(:,105) isOn(:,109) isOn(:,113) isOn(:,117) isOn(:,121) isOn(:,125)...
    isOn(:,129) isOn(:,133) isOn(:,137) isOn(:,141) isOn(:,145) isOn(:,149) isOn(:,153) isOn(:,157)...
    isOn(:,161) isOn(:,165) isOn(:,169) isOn(:,173) isOn(:,177) isOn(:,181) isOn(:,185) isOn(:,189)...
    isOn(:,193) isOn(:,197) isOn(:,201) isOn(:,205) isOn(:,209) isOn(:,213) isOn(:,217) isOn(:,221)...
    isOn(:,225) isOn(:,229) isOn(:,233) isOn(:,237) isOn(:,241) isOn(:,245) isOn(:,249) isOn(:,253)...
    isOn(:,257) isOn(:,261) isOn(:,265) isOn(:,269) isOn(:,273) isOn(:,277) isOn(:,281) isOn(:,285)...
    isOn(:,289) isOn(:,293) isOn(:,297) isOn(:,301) isOn(:,305) isOn(:,309) isOn(:,313) isOn(:,317)...
    isOn(:,321) isOn(:,325) isOn(:,329) isOn(:,333) isOn(:,337) isOn(:,341) isOn(:,345) isOn(:,349)...
    isOn(:,353) isOn(:,357) isOn(:,361) isOn(:,365) isOn(:,369) isOn(:,373) isOn(:,377) isOn(:,381)...
    isOn(:,385) isOn(:,389) isOn(:,393) isOn(:,397) isOn(:,401) isOn(:,405) isOn(:,409) isOn(:,413)...
    isOn(:,417) isOn(:,421) isOn(:,425) isOn(:,429) isOn(:,433) isOn(:,437) isOn(:,441) isOn(:,445)...
    isOn(:,449) isOn(:,453) isOn(:,457) isOn(:,461) isOn(:,465) isOn(:,469) isOn(:,473) isOn(:,477)...
    isOn(:,481) isOn(:,485) isOn(:,489) isOn(:,493) isOn(:,497) isOn(:,501) isOn(:,505) isOn(:,509)...
    isOn(:,513) isOn(:,517) isOn(:,521) isOn(:,525) isOn(:,529) isOn(:,533) isOn(:,537) isOn(:,541)...
    isOn(:,545) isOn(:,549) isOn(:,553) isOn(:,557) isOn(:,561) isOn(:,565) isOn(:,569) isOn(:,573)...
    isOn(:,577) isOn(:,581) isOn(:,585) isOn(:,589) isOn(:,593) isOn(:,597) isOn(:,601) isOn(:,605)...
    isOn(:,609) isOn(:,613) isOn(:,617) isOn(:,621) isOn(:,625) isOn(:,629) isOn(:,633) isOn(:,637)...
    isOn(:,641) isOn(:,645) isOn(:,649) isOn(:,653) isOn(:,657) isOn(:,661) isOn(:,665) isOn(:,669)...
    isOn(:,673) isOn(:,677) isOn(:,681) isOn(:,685) isOn(:,689) isOn(:,693) isOn(:,697) isOn(:,701)...
    isOn(:,705) isOn(:,709) isOn(:,713) isOn(:,717) isOn(:,721) isOn(:,725) isOn(:,729) isOn(:,733)...
    isOn(:,737) isOn(:,741) isOn(:,745) isOn(:,749) isOn(:,753) isOn(:,757) isOn(:,761) isOn(:,765)...
    isOn(:,769) isOn(:,773) isOn(:,777) isOn(:,781) isOn(:,785) isOn(:,789) isOn(:,793) isOn(:,797)];

%powerCost holds the value for the additional fuel cost for additional
%generation above the minimum power level. Note that this is calculated using the
%50 selected scenarios with probVecx holding the scenario probability
powerCost = sum(power(:,1:16)*f,1)*probVecx(1)+sum(power(:,17:32)*f,1)*probVecx(2)+...
  sum(power(:,33:48)*f,1)*probVecx(3)+sum(power(:,49:64)*f,1)*probVecx(4)+...
  sum(power(:,65:80)*f,1)*probVecx(5)+sum(power(:,81:96)*f,1)*probVecx(6)+...
  sum(power(:,97:112)*f,1)*probVecx(7)+sum(power(:,113:128)*f,1)*probVecx(8)+...
  sum(power(:,129:144)*f,1)*probVecx(9)+sum(power(:,145:160)*f,1)*probVecx(10)+...
  sum(power(:,161:176)*f,1)*probVecx(11)+sum(power(:,177:192)*f,1)*probVecx(12)+...
  sum(power(:,193:208)*f,1)*probVecx(13)+sum(power(:,209:224)*f,1)*probVecx(14)+...
  sum(power(:,225:240)*f,1)*probVecx(15)+sum(power(:,241:256)*f,1)*probVecx(16)+...
  sum(power(:,257:272)*f,1)*probVecx(17)+sum(power(:,273:288)*f,1)*probVecx(18)+...
  sum(power(:,289:304)*f,1)*probVecx(19)+sum(power(:,305:320)*f,1)*probVecx(20)+...
  sum(power(:,321:336)*f,1)*probVecx(21)+sum(power(:,337:352)*f,1)*probVecx(22)+...
  sum(power(:,353:368)*f,1)*probVecx(23)+sum(power(:,369:384)*f,1)*probVecx(24)+...
  sum(power(:,385:400)*f,1)*probVecx(25)+sum(power(:,401:416)*f,1)*probVecx(26)+...
  sum(power(:,417:432)*f,1)*probVecx(27)+sum(power(:,433:448)*f,1)*probVecx(28)+...
  sum(power(:,449:464)*f,1)*probVecx(29)+sum(power(:,465:480)*f,1)*probVecx(30)+...
  sum(power(:,481:496)*f,1)*probVecx(31)+sum(power(:,497:512)*f,1)*probVecx(32)+...
  sum(power(:,513:528)*f,1)*probVecx(33)+sum(power(:,529:544)*f,1)*probVecx(34)+...
  sum(power(:,545:560)*f,1)*probVecx(35)+sum(power(:,561:576)*f,1)*probVecx(36)+...
  sum(power(:,577:592)*f,1)*probVecx(37)+sum(power(:,593:608)*f,1)*probVecx(38)+...
  sum(power(:,609:624)*f,1)*probVecx(39)+sum(power(:,625:640)*f,1)*probVecx(40)+...
  sum(power(:,641:656)*f,1)*probVecx(41)+sum(power(:,657:672)*f,1)*probVecx(42)+...
  sum(power(:,673:688)*f,1)*probVecx(43)+sum(power(:,689:704)*f,1)*probVecx(44)+...
  sum(power(:,705:720)*f,1)*probVecx(45)+sum(power(:,721:736)*f,1)*probVecx(46)+...
  sum(power(:,737:752)*f,1)*probVecx(47)+sum(power(:,753:768)*f,1)*probVecx(48)+...
  sum(power(:,769:784)*f,1)*probVecx(49)+sum(power(:,785:800)*f,1)*probVecx(50);

%isOnCost holds the value for the fuel cost at minimum power level. Note that this
% is calculated using the selected 50 scenarios with probVecx holding the scenario probability
isOnCost = sum(isOn1(:,1:4)*OperatingCost',1)*probVecx(1)+sum(isOn1(:,5:8)*OperatingCost',1)*probVecx(2)+...
sum(isOn1(:,9:12)*OperatingCost',1)*probVecx(3)+sum(isOn1(:,13:16)*OperatingCost',1)*probVecx(4)+...
sum(isOn1(:,17:20)*OperatingCost',1)*probVecx(5)+sum(isOn1(:,21:24)*OperatingCost',1)*probVecx(6)+...
sum(isOn1(:,25:28)*OperatingCost',1)*probVecx(7)+sum(isOn1(:,29:32)*OperatingCost',1)*probVecx(8)+...
sum(isOn1(:,33:36)*OperatingCost',1)*probVecx(9)+sum(isOn1(:,37:40)*OperatingCost',1)*probVecx(10)+...
sum(isOn1(:,41:44)*OperatingCost',1)*probVecx(11)+sum(isOn1(:,45:48)*OperatingCost',1)*probVecx(12)+...
sum(isOn1(:,49:52)*OperatingCost',1)*probVecx(13)+sum(isOn1(:,53:56)*OperatingCost',1)*probVecx(14)+...
sum(isOn1(:,57:60)*OperatingCost',1)*probVecx(15)+sum(isOn1(:,61:64)*OperatingCost',1)*probVecx(16)+...
sum(isOn1(:,65:68)*OperatingCost',1)*probVecx(17)+sum(isOn1(:,69:72)*OperatingCost',1)*probVecx(18)+...
sum(isOn1(:,73:76)*OperatingCost',1)*probVecx(19)+sum(isOn1(:,77:80)*OperatingCost',1)*probVecx(20)+...
sum(isOn1(:,81:84)*OperatingCost',1)*probVecx(21)+sum(isOn1(:,85:88)*OperatingCost',1)*probVecx(22)+...
sum(isOn1(:,89:92)*OperatingCost',1)*probVecx(23)+sum(isOn1(:,93:96)*OperatingCost',1)*probVecx(24)+...
sum(isOn1(:,97:100)*OperatingCost',1)*probVecx(25)+sum(isOn1(:,101:104)*OperatingCost',1)*probVecx(26)+...
sum(isOn1(:,105:108)*OperatingCost',1)*probVecx(27)+sum(isOn1(:,109:112)*OperatingCost',1)*probVecx(28)+...
sum(isOn1(:,113:116)*OperatingCost',1)*probVecx(29)+sum(isOn1(:,117:120)*OperatingCost',1)*probVecx(30)+...
sum(isOn1(:,121:124)*OperatingCost',1)*probVecx(31)+sum(isOn1(:,125:128)*OperatingCost',1)*probVecx(32)+...
sum(isOn1(:,129:132)*OperatingCost',1)*probVecx(33)+sum(isOn1(:,133:136)*OperatingCost',1)*probVecx(34)+...
sum(isOn1(:,137:140)*OperatingCost',1)*probVecx(35)+sum(isOn1(:,141:144)*OperatingCost',1)*probVecx(36)+...
sum(isOn1(:,145:148)*OperatingCost',1)*probVecx(37)+sum(isOn1(:,149:152)*OperatingCost',1)*probVecx(38)+...
sum(isOn1(:,153:156)*OperatingCost',1)*probVecx(39)+sum(isOn1(:,157:160)*OperatingCost',1)*probVecx(40)+...
sum(isOn1(:,161:164)*OperatingCost',1)*probVecx(41)+sum(isOn1(:,165:168)*OperatingCost',1)*probVecx(42)+...
sum(isOn1(:,169:172)*OperatingCost',1)*probVecx(43)+sum(isOn1(:,173:176)*OperatingCost',1)*probVecx(44)+...
sum(isOn1(:,177:180)*OperatingCost',1)*probVecx(45)+sum(isOn1(:,181:184)*OperatingCost',1)*probVecx(46)+...
sum(isOn1(:,185:188)*OperatingCost',1)*probVecx(47)+sum(isOn1(:,189:192)*OperatingCost',1)*probVecx(48)+...
sum(isOn1(:,193:196)*OperatingCost',1)*probVecx(49)+sum(isOn1(:,197:200)*OperatingCost',1)*probVecx(50);

%powerCO2Cost holds the value for the additional carbon emission tax for 
%power level above minimum power level. Note that this is calculated using
%the selected 50 scenarios with probVecx holding the scenario probability
powerprob.Constraints.CO2_p=powerCO2 >= sum(power(:,1:16)*CO2_f,1)*probVecx(1)+sum(power(:,17:32)*CO2_f,1)*probVecx(2)+...
  sum(power(:,33:48)*CO2_f,1)*probVecx(3)+sum(power(:,49:64)*CO2_f,1)*probVecx(4)+...
  sum(power(:,65:80)*CO2_f,1)*probVecx(5)+sum(power(:,81:96)*CO2_f,1)*probVecx(6)+...
  sum(power(:,97:112)*CO2_f,1)*probVecx(7)+sum(power(:,113:128)*CO2_f,1)*probVecx(8)+...
  sum(power(:,129:144)*CO2_f,1)*probVecx(9)+sum(power(:,145:160)*CO2_f,1)*probVecx(10)+...
  sum(power(:,161:176)*CO2_f,1)*probVecx(11)+sum(power(:,177:192)*CO2_f,1)*probVecx(12)+...
  sum(power(:,193:208)*CO2_f,1)*probVecx(13)+sum(power(:,209:224)*CO2_f,1)*probVecx(14)+...
  sum(power(:,225:240)*CO2_f,1)*probVecx(15)+sum(power(:,241:256)*CO2_f,1)*probVecx(16)+...
  sum(power(:,257:272)*CO2_f,1)*probVecx(17)+sum(power(:,273:288)*CO2_f,1)*probVecx(18)+...
  sum(power(:,289:304)*CO2_f,1)*probVecx(19)+sum(power(:,305:320)*CO2_f,1)*probVecx(20)+...
  sum(power(:,321:336)*CO2_f,1)*probVecx(21)+sum(power(:,337:352)*CO2_f,1)*probVecx(22)+...
  sum(power(:,353:368)*CO2_f,1)*probVecx(23)+sum(power(:,369:384)*CO2_f,1)*probVecx(24)+...
  sum(power(:,385:400)*CO2_f,1)*probVecx(25)+sum(power(:,401:416)*CO2_f,1)*probVecx(26)+...
  sum(power(:,417:432)*CO2_f,1)*probVecx(27)+sum(power(:,433:448)*CO2_f,1)*probVecx(28)+...
  sum(power(:,449:464)*CO2_f,1)*probVecx(29)+sum(power(:,465:480)*CO2_f,1)*probVecx(30)+...
  sum(power(:,481:496)*CO2_f,1)*probVecx(31)+sum(power(:,497:512)*CO2_f,1)*probVecx(32)+...
  sum(power(:,513:528)*CO2_f,1)*probVecx(33)+sum(power(:,529:544)*CO2_f,1)*probVecx(34)+...
  sum(power(:,545:560)*CO2_f,1)*probVecx(35)+sum(power(:,561:576)*CO2_f,1)*probVecx(36)+...
  sum(power(:,577:592)*CO2_f,1)*probVecx(37)+sum(power(:,593:608)*CO2_f,1)*probVecx(38)+...
  sum(power(:,609:624)*CO2_f,1)*probVecx(39)+sum(power(:,625:640)*CO2_f,1)*probVecx(40)+...
  sum(power(:,641:656)*CO2_f,1)*probVecx(41)+sum(power(:,657:672)*CO2_f,1)*probVecx(42)+...
  sum(power(:,673:688)*CO2_f,1)*probVecx(43)+sum(power(:,689:704)*CO2_f,1)*probVecx(44)+...
  sum(power(:,705:720)*CO2_f,1)*probVecx(45)+sum(power(:,721:736)*CO2_f,1)*probVecx(46)+...
  sum(power(:,737:752)*CO2_f,1)*probVecx(47)+sum(power(:,753:768)*CO2_f,1)*probVecx(48)+...
  sum(power(:,769:784)*CO2_f,1)*probVecx(49)+sum(power(:,785:800)*CO2_f,1)*probVecx(50);

%isOnCO2Cost holds the carbon emission tax at minimum power level. Note that this
% is calculated using the selected 50 scenarios with probVecx holding the scenario probability
powerprob.Constraints.CO2_o=isOnCO2 >= sum(isOn1(:,1:4)*OperatingCO2Cost',1)*probVecx(1)+sum(isOn1(:,5:8)*OperatingCO2Cost',1)*probVecx(2)+...
sum(isOn1(:,9:12)*OperatingCO2Cost',1)*probVecx(3)+sum(isOn1(:,13:16)*OperatingCO2Cost',1)*probVecx(4)+...
sum(isOn1(:,17:20)*OperatingCO2Cost',1)*probVecx(5)+sum(isOn1(:,21:24)*OperatingCO2Cost',1)*probVecx(6)+...
sum(isOn1(:,25:28)*OperatingCO2Cost',1)*probVecx(7)+sum(isOn1(:,29:32)*OperatingCO2Cost',1)*probVecx(8)+...
sum(isOn1(:,33:36)*OperatingCO2Cost',1)*probVecx(9)+sum(isOn1(:,37:40)*OperatingCO2Cost',1)*probVecx(10)+...
sum(isOn1(:,41:44)*OperatingCO2Cost',1)*probVecx(11)+sum(isOn1(:,45:48)*OperatingCO2Cost',1)*probVecx(12)+...
sum(isOn1(:,49:52)*OperatingCO2Cost',1)*probVecx(13)+sum(isOn1(:,53:56)*OperatingCO2Cost',1)*probVecx(14)+...
sum(isOn1(:,57:60)*OperatingCO2Cost',1)*probVecx(15)+sum(isOn1(:,61:64)*OperatingCO2Cost',1)*probVecx(16)+...
sum(isOn1(:,65:68)*OperatingCO2Cost',1)*probVecx(17)+sum(isOn1(:,69:72)*OperatingCO2Cost',1)*probVecx(18)+...
sum(isOn1(:,73:76)*OperatingCO2Cost',1)*probVecx(19)+sum(isOn1(:,77:80)*OperatingCO2Cost',1)*probVecx(20)+...
sum(isOn1(:,81:84)*OperatingCO2Cost',1)*probVecx(21)+sum(isOn1(:,85:88)*OperatingCO2Cost',1)*probVecx(22)+...
sum(isOn1(:,89:92)*OperatingCO2Cost',1)*probVecx(23)+sum(isOn1(:,93:96)*OperatingCO2Cost',1)*probVecx(24)+...
sum(isOn1(:,97:100)*OperatingCO2Cost',1)*probVecx(25)+sum(isOn1(:,101:104)*OperatingCO2Cost',1)*probVecx(26)+...
sum(isOn1(:,105:108)*OperatingCO2Cost',1)*probVecx(27)+sum(isOn1(:,109:112)*OperatingCO2Cost',1)*probVecx(28)+...
sum(isOn1(:,113:116)*OperatingCO2Cost',1)*probVecx(29)+sum(isOn1(:,117:120)*OperatingCO2Cost',1)*probVecx(30)+...
sum(isOn1(:,121:124)*OperatingCO2Cost',1)*probVecx(31)+sum(isOn1(:,125:128)*OperatingCO2Cost',1)*probVecx(32)+...
sum(isOn1(:,129:132)*OperatingCO2Cost',1)*probVecx(33)+sum(isOn1(:,133:136)*OperatingCO2Cost',1)*probVecx(34)+...
sum(isOn1(:,137:140)*OperatingCO2Cost',1)*probVecx(35)+sum(isOn1(:,141:144)*OperatingCO2Cost',1)*probVecx(36)+...
sum(isOn1(:,145:148)*OperatingCO2Cost',1)*probVecx(37)+sum(isOn1(:,149:152)*OperatingCO2Cost',1)*probVecx(38)+...
sum(isOn1(:,153:156)*OperatingCO2Cost',1)*probVecx(39)+sum(isOn1(:,157:160)*OperatingCO2Cost',1)*probVecx(40)+...
sum(isOn1(:,161:164)*OperatingCO2Cost',1)*probVecx(41)+sum(isOn1(:,165:168)*OperatingCO2Cost',1)*probVecx(42)+...
sum(isOn1(:,169:172)*OperatingCO2Cost',1)*probVecx(43)+sum(isOn1(:,173:176)*OperatingCO2Cost',1)*probVecx(44)+...
sum(isOn1(:,177:180)*OperatingCO2Cost',1)*probVecx(45)+sum(isOn1(:,181:184)*OperatingCO2Cost',1)*probVecx(46)+...
sum(isOn1(:,185:188)*OperatingCO2Cost',1)*probVecx(47)+sum(isOn1(:,189:192)*OperatingCO2Cost',1)*probVecx(48)+...
sum(isOn1(:,193:196)*OperatingCO2Cost',1)*probVecx(49)+sum(isOn1(:,197:200)*OperatingCO2Cost',1)*probVecx(50);

powerCO2Cost=powerCO2*co2_tax;
isOnCO2Cost=isOnCO2*co2_tax;
%startupCost holds the value for the startupCost. Note that this is calculated 
%using the selected 50 scenarios with probVecx holding the scenario probability
%Note that startup is the same for all scenarios
startupCost = sum(startup(:,1:4)*StartupCost',1)*probVecx(1)+sum(startup(:,5:8)*StartupCost',1)*probVecx(2)+...
sum(startup(:,9:12)*StartupCost',1)*probVecx(3)+sum(startup(:,13:16)*StartupCost',1)*probVecx(4)+...
sum(startup(:,17:20)*StartupCost',1)*probVecx(5)+sum(startup(:,21:24)*StartupCost',1)*probVecx(6)+...
sum(startup(:,25:28)*StartupCost',1)*probVecx(7)+sum(startup(:,29:32)*StartupCost',1)*probVecx(8)+...
sum(startup(:,33:36)*StartupCost',1)*probVecx(9)+sum(startup(:,37:40)*StartupCost',1)*probVecx(10)+...
sum(startup(:,41:44)*StartupCost',1)*probVecx(11)+sum(startup(:,45:48)*StartupCost',1)*probVecx(12)+...
sum(startup(:,49:52)*StartupCost',1)*probVecx(13)+sum(startup(:,53:56)*StartupCost',1)*probVecx(14)+...
sum(startup(:,57:60)*StartupCost',1)*probVecx(15)+sum(startup(:,61:64)*StartupCost',1)*probVecx(16)+...
sum(startup(:,65:68)*StartupCost',1)*probVecx(17)+sum(startup(:,69:72)*StartupCost',1)*probVecx(18)+...
sum(startup(:,73:76)*StartupCost',1)*probVecx(19)+sum(startup(:,77:80)*StartupCost',1)*probVecx(20)+...
sum(startup(:,81:84)*StartupCost',1)*probVecx(21)+sum(startup(:,85:88)*StartupCost',1)*probVecx(22)+...
sum(startup(:,89:92)*StartupCost',1)*probVecx(23)+sum(startup(:,93:96)*StartupCost',1)*probVecx(24)+...
sum(startup(:,97:100)*StartupCost',1)*probVecx(25)+sum(startup(:,101:104)*StartupCost',1)*probVecx(26)+...
sum(startup(:,105:108)*StartupCost',1)*probVecx(27)+sum(startup(:,109:112)*StartupCost',1)*probVecx(28)+...
sum(startup(:,113:116)*StartupCost',1)*probVecx(29)+sum(startup(:,117:120)*StartupCost',1)*probVecx(30)+...
sum(startup(:,121:124)*StartupCost',1)*probVecx(31)+sum(startup(:,125:128)*StartupCost',1)*probVecx(32)+...
sum(startup(:,129:132)*StartupCost',1)*probVecx(33)+sum(startup(:,133:136)*StartupCost',1)*probVecx(34)+...
sum(startup(:,137:140)*StartupCost',1)*probVecx(35)+sum(startup(:,141:144)*StartupCost',1)*probVecx(36)+...
sum(startup(:,145:148)*StartupCost',1)*probVecx(37)+sum(startup(:,149:152)*StartupCost',1)*probVecx(38)+...
sum(startup(:,153:156)*StartupCost',1)*probVecx(39)+sum(startup(:,157:160)*StartupCost',1)*probVecx(40)+...
sum(startup(:,161:164)*StartupCost',1)*probVecx(41)+sum(startup(:,165:168)*StartupCost',1)*probVecx(42)+...
sum(startup(:,169:172)*StartupCost',1)*probVecx(43)+sum(startup(:,173:176)*StartupCost',1)*probVecx(44)+...
sum(startup(:,177:180)*StartupCost',1)*probVecx(45)+sum(startup(:,181:184)*StartupCost',1)*probVecx(46)+...
sum(startup(:,185:188)*StartupCost',1)*probVecx(47)+sum(startup(:,189:192)*StartupCost',1)*probVecx(48)+...
sum(startup(:,193:196)*StartupCost',1)*probVecx(49)+sum(startup(:,197:200)*StartupCost',1)*probVecx(50);

%=====================BESS DEGRADATION COST MODELLING=====================
%alpha_cyc and beta_cyc are parameters that model the cycling aging, alpha_cal
%and beta_cal are parameters that model calendar aging
alpha_cyc=4.42e-5;
beta_cyc=0.02676;
alpha_cal=1.985e-7;
beta_cal=0.051;
BESS_usable=0.88;
BESS_cost=500000*BESS_nom; %This is the capital cost of the BESS

%Power_Segment2 defines the segment for BESS usage. The BESS discharge
%power is used to calculate BESS degradation cost
Power_Segment2 = zeros(1, (K + 1));

for k = 1:1:(L + 1) %This for loop defines the segment for BESS usage
%BESS_uable is used in the expressions below because not 100% of the BESS is usable
        if k == 1
            Power_Segment2(k)  = 0;
        elseif k == (L + 1)
            Power_Segment2(k) = BESS_nom*BESS_usable*12; 
        else
            Power_Segment2(k) = (k - 1) * BESS_nom*BESS_usable*12/L;
        end
end

%The for loop below is used to calculate the BESS degradation cost at the
%end of each BESS usage segment
    for k = 1:1:(L + 1) 
         Cost1(k) = BESS_cost/(20/(alpha_cyc*exp(beta_cyc*298)*...
             (Power_Segment2(k)/(BESS_nom))^0.5+alpha_cal*exp(beta_cal*298)*(1/30)^0.5))^2;
    end

%The for loop below is used to calcul.ate the BESS degradation cost slope
%for each BESS usage segment
    for k = 1:1:L
%0.95 is used in the expression below because only 95% of the BESS is usable
        Slope1(k) = (Cost1(k+1) - Cost1(k)) / (BESS_nom*BESS_usable*12/L);
    end
    
%f21 is a variable that holds the BESS degradation cost slope for all BESS usage segments
    f21 = zeros((L), 1);

%Create sub matrices to easily formulate f21
sub_f21 = cell(1, 1);

%The for loop below is used to assign BESS degradation cost slope values to the sub matrices
    for k = 1:1:L
            sub_f21{1, 1}(k, 1) = Slope1(k);    
    end
 
    %Convert cell to matrix to formulate f21
f21 = cell2mat(sub_f21);

%bessCost is used to hold the BESS degradation cost. Note that this is calculated
% using the 50 selected scenarios where probVecx holds the scenario probability.
% Note that at the minimum the BESS has a degradation cost equal to Cost1(1)
bessCost = (BESS_thru(:,1)'*f21+Cost1(1))*probVecx(1)+(BESS_thru(:,2)'*f21+Cost1(1))*probVecx(2)+...
(BESS_thru(:,3)'*f21+Cost1(1))*probVecx(3)+(BESS_thru(:,4)'*f21+Cost1(1))*probVecx(4)+...
(BESS_thru(:,5)'*f21+Cost1(1))*probVecx(5)+(BESS_thru(:,6)'*f21+Cost1(1))*probVecx(6)+...
(BESS_thru(:,7)'*f21+Cost1(1))*probVecx(7)+(BESS_thru(:,8)'*f21+Cost1(1))*probVecx(8)+...
(BESS_thru(:,9)'*f21+Cost1(1))*probVecx(9)+(BESS_thru(:,10)'*f21+Cost1(1))*probVecx(10)+...
(BESS_thru(:,11)'*f21+Cost1(1))*probVecx(11)+(BESS_thru(:,12)'*f21+Cost1(1))*probVecx(12)+...
(BESS_thru(:,13)'*f21+Cost1(1))*probVecx(13)+(BESS_thru(:,14)'*f21+Cost1(1))*probVecx(14)+...
(BESS_thru(:,15)'*f21+Cost1(1))*probVecx(15)+(BESS_thru(:,16)'*f21+Cost1(1))*probVecx(16)+...
(BESS_thru(:,17)'*f21+Cost1(1))*probVecx(17)+(BESS_thru(:,18)'*f21+Cost1(1))*probVecx(18)+...
(BESS_thru(:,19)'*f21+Cost1(1))*probVecx(19)+(BESS_thru(:,20)'*f21+Cost1(1))*probVecx(20)+...
(BESS_thru(:,21)'*f21+Cost1(1))*probVecx(21)+(BESS_thru(:,22)'*f21+Cost1(1))*probVecx(22)+...
(BESS_thru(:,23)'*f21+Cost1(1))*probVecx(23)+(BESS_thru(:,24)'*f21+Cost1(1))*probVecx(24)+...
(BESS_thru(:,25)'*f21+Cost1(1))*probVecx(25)+(BESS_thru(:,26)'*f21+Cost1(1))*probVecx(26)+...
(BESS_thru(:,27)'*f21+Cost1(1))*probVecx(27)+(BESS_thru(:,28)'*f21+Cost1(1))*probVecx(28)+...
(BESS_thru(:,29)'*f21+Cost1(1))*probVecx(29)+(BESS_thru(:,30)'*f21+Cost1(1))*probVecx(30)+...
(BESS_thru(:,31)'*f21+Cost1(1))*probVecx(31)+(BESS_thru(:,32)'*f21+Cost1(1))*probVecx(32)+...
(BESS_thru(:,33)'*f21+Cost1(1))*probVecx(33)+(BESS_thru(:,34)'*f21+Cost1(1))*probVecx(34)+...
(BESS_thru(:,35)'*f21+Cost1(1))*probVecx(35)+(BESS_thru(:,36)'*f21+Cost1(1))*probVecx(36)+...
(BESS_thru(:,37)'*f21+Cost1(1))*probVecx(37)+(BESS_thru(:,38)'*f21+Cost1(1))*probVecx(38)+...
(BESS_thru(:,39)'*f21+Cost1(1))*probVecx(39)+(BESS_thru(:,40)'*f21+Cost1(1))*probVecx(40)+...
(BESS_thru(:,41)'*f21+Cost1(1))*probVecx(41)+(BESS_thru(:,42)'*f21+Cost1(1))*probVecx(42)+...
(BESS_thru(:,43)'*f21+Cost1(1))*probVecx(43)+(BESS_thru(:,44)'*f21+Cost1(1))*probVecx(44)+...
(BESS_thru(:,45)'*f21+Cost1(1))*probVecx(45)+(BESS_thru(:,46)'*f21+Cost1(1))*probVecx(46)+...
(BESS_thru(:,47)'*f21+Cost1(1))*probVecx(47)+(BESS_thru(:,48)'*f21+Cost1(1))*probVecx(48)+...
(BESS_thru(:,49)'*f21+Cost1(1))*probVecx(49)+(BESS_thru(:,50)'*f21+Cost1(1))*probVecx(50);

%====================END BESS DEGRADATION COST MODELLING===================

src=50; %src is the spinning reserve cost per available MW for the generators

%{
new_maxGenConst2 = maxGenConst.*isOn2; %This variable is used for computing the spinning reserve cost

sr_cost=(sum(sum(new_maxGenConst2(:,1:16)-power(:,1:16),2))*src)*probVecx(1)+(sum(sum(new_maxGenConst2(:,17:32)-power(:,17:32),2))*src)*probVecx(2)+...
(sum(sum(new_maxGenConst2(:,33:48)-power(:,33:48),2))*src)*probVecx(3)+(sum(sum(new_maxGenConst2(:,49:64)-power(:,49:64),2))*src)*probVecx(4)+...
(sum(sum(new_maxGenConst2(:,65:80)-power(:,65:80),2))*src)*probVecx(5)+(sum(sum(new_maxGenConst2(:,81:96)-power(:,81:96),2))*src)*probVecx(6)+...
(sum(sum(new_maxGenConst2(:,97:112)-power(:,97:112),2))*src)*probVecx(7)+(sum(sum(new_maxGenConst2(:,113:128)-power(:,113:128),2))*src)*probVecx(8)+...
(sum(sum(new_maxGenConst2(:,129:144)-power(:,129:144),2))*src)*probVecx(9)+(sum(sum(new_maxGenConst2(:,145:160)-power(:,145:160),2))*src)*probVecx(10)+...
(sum(sum(new_maxGenConst2(:,161:176)-power(:,161:176),2))*src)*probVecx(11)+(sum(sum(new_maxGenConst2(:,177:192)-power(:,177:192),2))*src)*probVecx(12)+...
(sum(sum(new_maxGenConst2(:,193:208)-power(:,193:208),2))*src)*probVecx(13)+(sum(sum(new_maxGenConst2(:,209:224)-power(:,209:224),2))*src)*probVecx(14)+...
(sum(sum(new_maxGenConst2(:,225:240)-power(:,225:240),2))*src)*probVecx(15)+(sum(sum(new_maxGenConst2(:,241:256)-power(:,241:256),2))*src)*probVecx(16)+...
(sum(sum(new_maxGenConst2(:,257:272)-power(:,257:272),2))*src)*probVecx(17)+(sum(sum(new_maxGenConst2(:,273:288)-power(:,273:288),2))*src)*probVecx(18)+...
(sum(sum(new_maxGenConst2(:,289:304)-power(:,289:304),2))*src)*probVecx(19)+(sum(sum(new_maxGenConst2(:,305:320)-power(:,305:320),2))*src)*probVecx(20)+...
(sum(sum(new_maxGenConst2(:,321:336)-power(:,321:336),2))*src)*probVecx(21)+(sum(sum(new_maxGenConst2(:,337:352)-power(:,337:352),2))*src)*probVecx(22)+...
(sum(sum(new_maxGenConst2(:,353:368)-power(:,353:368),2))*src)*probVecx(23)+(sum(sum(new_maxGenConst2(:,369:384)-power(:,369:384),2))*src)*probVecx(24)+...
(sum(sum(new_maxGenConst2(:,385:400)-power(:,385:400),2))*src)*probVecx(25)+(sum(sum(new_maxGenConst2(:,401:416)-power(:,401:416),2))*src)*probVecx(26)+...
(sum(sum(new_maxGenConst2(:,417:432)-power(:,417:432),2))*src)*probVecx(27)+(sum(sum(new_maxGenConst2(:,433:448)-power(:,433:448),2))*src)*probVecx(28)+...
(sum(sum(new_maxGenConst2(:,449:464)-power(:,449:464),2))*src)*probVecx(29)+(sum(sum(new_maxGenConst2(:,465:480)-power(:,465:480),2))*src)*probVecx(30)+...
(sum(sum(new_maxGenConst2(:,481:496)-power(:,481:496),2))*src)*probVecx(31)+(sum(sum(new_maxGenConst2(:,497:512)-power(:,497:512),2))*src)*probVecx(32)+...
(sum(sum(new_maxGenConst2(:,513:528)-power(:,513:528),2))*src)*probVecx(33)+(sum(sum(new_maxGenConst2(:,529:544)-power(:,529:544),2))*src)*probVecx(34)+...
(sum(sum(new_maxGenConst2(:,545:560)-power(:,545:560),2))*src)*probVecx(35)+(sum(sum(new_maxGenConst2(:,561:576)-power(:,561:576),2))*src)*probVecx(36)+...
(sum(sum(new_maxGenConst2(:,577:592)-power(:,577:592),2))*src)*probVecx(37)+(sum(sum(new_maxGenConst2(:,593:608)-power(:,593:608),2))*src)*probVecx(38)+...
(sum(sum(new_maxGenConst2(:,609:624)-power(:,609:624),2))*src)*probVecx(39)+(sum(sum(new_maxGenConst2(:,625:640)-power(:,625:640),2))*src)*probVecx(40)+...
(sum(sum(new_maxGenConst2(:,641:656)-power(:,641:656),2))*src)*probVecx(41)+(sum(sum(new_maxGenConst2(:,657:672)-power(:,657:672),2))*src)*probVecx(42)+...
(sum(sum(new_maxGenConst2(:,673:688)-power(:,673:688),2))*src)*probVecx(43)+(sum(sum(new_maxGenConst2(:,689:704)-power(:,689:704),2))*src)*probVecx(44)+...
(sum(sum(new_maxGenConst2(:,705:720)-power(:,705:720),2))*src)*probVecx(45)+(sum(sum(new_maxGenConst2(:,721:736)-power(:,721:736),2))*src)*probVecx(46)+...
(sum(sum(new_maxGenConst2(:,737:752)-power(:,737:752),2))*src)*probVecx(47)+(sum(sum(new_maxGenConst2(:,753:768)-power(:,753:768),2))*src)*probVecx(48)+...
(sum(sum(new_maxGenConst2(:,769:784)-power(:,769:784),2))*src)*probVecx(49)+(sum(sum(new_maxGenConst2(:,785:800)-power(:,785:800),2))*src)*probVecx(50);

%sr_cost holds the value for the cost of spinning reserve for all
%generators. Note that this is calculated using the 50 selected scenarios
%where probVecx is the scenario probability
%}

%The objective function is formulated here
if test_case == 4
powerprob.Objective = (powerCost + isOnCost + startupCost...
    +BESS_cost/(12.5*365)+powerCO2Cost+isOnCO2Cost+sum(startup(1:96))*(110.56));%79.85 BESS_cost/(20*365)+bessCost+BESS_cost*sum(BESS_power1,1)/(BESS_nom*9400);
else
  powerprob.Objective = (powerCost + isOnCost + startupCost...
    +bessCost+powerCO2Cost+isOnCO2Cost+sum(startup(1:96))*(110.56));%BESS_cost/(20*365)+bessCost+BESS_cost*sum(BESS_power1,1)/(BESS_nom*9400);
end  
%+sr_cost
%==========================BEGIN CONSTRAINTS===============================
powerprob.Constraints.isDemandMet = optimconstr(24,50);
powerprob.Constraints.WTpower = optimconstr(24,50);
powerprob.Constraints.solar_power = optimconstr(24,50);
powerprob.Constraints.flex_load1 = optimconstr(24,50);
powerprob.Constraints.flex_load2 = optimconstr(24,50);
powerprob.Constraints.seg1a=optimconstr(24,50);
powerprob.Constraints.seg1b=optimconstr(24,50);
powerprob.Constraints.seg1c=optimconstr(24,50);
powerprob.Constraints.seg1d=optimconstr(24,50);
powerprob.Constraints.seg2a=optimconstr(24,50);
powerprob.Constraints.seg2b=optimconstr(24,50);
powerprob.Constraints.seg2c=optimconstr(24,50);
powerprob.Constraints.seg2d=optimconstr(24,50);
powerprob.Constraints.seg3a=optimconstr(24,50);
powerprob.Constraints.seg3b=optimconstr(24,50);
powerprob.Constraints.seg3c=optimconstr(24,50);
powerprob.Constraints.seg3d=optimconstr(24,50);

powerprob.Constraints.powerOnlyWhenOn = power <= new_maxGenConst; 
BESS_eff=0.93;
for i=1:50
    Load(:,i)=load_selx(:,i)+(4.5*ia/(1000)); %4.5 kW/MWh is the BESS HVAC power
    
    %load demand constraints
powerprob.Constraints.isDemandMet(:,i) = sum(power(:,1+16*(i-1):16*i),2) +BESS_power1(:,i)+BESS_power2(:,i)+WT_power(:,i)...
   +solar_power(:,i) >= Load(:,i) - sum((minPowConst(:,1+4*(i-1):4*i).*isOn1(:,1+4*(i-1):4*i)),2)+flex_load1(:,i)+flex_load2(:,i); %

%flexible load constraints
powerprob.Constraints.flex_load1(:,i) =flex_load1(:,i)<=6-wis(:,i);
powerprob.Constraints.flex_load2(:,i) =flex_load2(:,i)>=wis(:,i)-6;

%Constraints that enforces the power segments being used from the
%costliest to the cheapest
powerprob.Constraints.seg1a(:,i)=power(:,1+16*(i-1))>=isOn(:,2+16*(i-1)).*W(1,1);
powerprob.Constraints.seg1b(:,i)=power(:,5+16*(i-1))>=isOn(:,6+16*(i-1)).*W(2,1);
powerprob.Constraints.seg1c(:,i)=power(:,9+16*(i-1))>=isOn(:,10+16*(i-1)).*W(3,1);
powerprob.Constraints.seg1d(:,i)=power(:,13+16*(i-1))>=isOn(:,14+16*(i-1)).*W(4,1);
powerprob.Constraints.seg2a(:,i)=power(:,2+16*(i-1))>=isOn(:,3+16*(i-1)).*W(1,1);
powerprob.Constraints.seg2b(:,i)=power(:,6+16*(i-1))>=isOn(:,7+16*(i-1)).*W(2,1);
powerprob.Constraints.seg2c(:,i)=power(:,10+16*(i-1))>=isOn(:,11+16*(i-1)).*W(3,1);
powerprob.Constraints.seg2d(:,i)=power(:,14+16*(i-1))>=isOn(:,15+16*(i-1)).*W(4,1);
powerprob.Constraints.seg3a(:,i)=power(:,3+16*(i-1))>=isOn(:,4+16*(i-1)).*W(1,1);
powerprob.Constraints.seg3b(:,i)=power(:,7+16*(i-1))>=isOn(:,8+16*(i-1)).*W(2,1);
powerprob.Constraints.seg3c(:,i)=power(:,11+16*(i-1))>=isOn(:,12+16*(i-1)).*W(3,1);
powerprob.Constraints.seg3d(:,i)=power(:,15+16*(i-1))>=isOn(:,16+16*(i-1)).*W(4,1);
end

%More flexible load constraints
powerprob.Constraints.flex_load3 = flex_load1 <= 6.*flex_isOn1;
powerprob.Constraints.flex_load4 = flex_load2 >= -6.*(flex_isOn2);
powerprob.Constraints.flex_load5 = flex_isOn1+flex_isOn2 == 1; 
powerprob.Constraints.flex_load6 = sum(flex_load1,1)+sum(flex_load2,1) == 0; 

oneOn=1:4;
%The constraint below is used to enforce the same unit commitment for all scenarios
powerprob.Constraints.isOnOne = optimconstr(24,4,49);
for l=1:49
powerprob.Constraints.isOnOne(:,:,l)=isOn1(:,oneOn+(l-1)*4)==isOn1(:,oneOn+(l)*4);
end

%===BESS constraints===
%0.95 is used in some of the constraints in this section because the usable energy is 95%
%The constraint below set the initial state of charge of the BESS
powerprob.Constraints.BESS1 = BESS_cap(1,:) == (0.605-0.05)*BESS_nom;

%The constraints below limit the discharge and charge power based on the
%state of charge
powerprob.Constraints.BESS2 = BESS_power1 <= BESS_cap;
powerprob.Constraints.BESS3 = BESS_power2 >= (BESS_cap-BESS_nom*BESS_usable)/BESS_eff;

%The constraint below calculates the state of charge based on the discharge and charge power
powerprob.Constraints.BESS4 = optimconstr(nHours-1,50);
for l=1:50
for i=2:nHours
powerprob.Constraints.BESS4(i,l) = BESS_cap(i,l) == BESS_cap(i-1,l)-BESS_power1(i-1,l)-BESS_power2(i-1,l)*BESS_eff;
end
end

%BESS discharge/charge status constraints
powerprob.Constraints.BESS5 = BESS_power1 <= (BESS_nom*BESS_usable).*BESS_isOn1;
powerprob.Constraints.BESS6 = BESS_power2 >= (-BESS_nom*BESS_usable/BESS_eff).*(BESS_isOn2);
powerprob.Constraints.BESS7 = BESS_isOn1+BESS_isOn2 == 1; 

%BESS usage constraints
powerprob.Constraints.BESS8 = optimconstr(1,50);
for l=1:50
powerprob.Constraints.BESS8(1,l) = sum(BESS_thru(:,l))==sum(BESS_power1(:,l));
end
powerprob.Constraints.BESS9 = BESS_thru<=isOn_thru.*(BESS_nom*BESS_usable*12/(L));

%Constraints setting initial state of charge equal to final state of charge
powerprob.Constraints.BESS10 = BESS_cap(1,:)==BESS_cap(24,:)-BESS_power1(24,:)-BESS_power2(24,:)*BESS_eff;

%Constraint ensuring that the BESS usage segments are used starting from 
%the costliest to the cheapest
powerprob.Constraints.BESS11 = optimconstr(L-1,50);
for l=1:50
for i=1:L-1
powerprob.Constraints.BESS11(i,l)=BESS_thru(i,l)>=isOn_thru(i+1,l)*(BESS_nom*BESS_usable*12/(L));
end
end

%{
powerprob.Constraints.isOna19=sum(isOn(1:24))+sum(isOn(97:120))+...
    sum(isOn(193:216))+sum(isOn(289:312))<=65;
%powerprob.Constraints.isOna20=sum(startup(1:96))<=5;
%
powerprob.Constraints.isOna1a=isOn(:,1)==sol5c_25(:,1);
powerprob.Constraints.isOna1b=isOn(:,5)==sol5c_25(:,2);
powerprob.Constraints.isOna1c=isOn(:,9)==sol5c_25(:,3);
powerprob.Constraints.isOna1d=isOn(:,13)==sol5c_25(:,4);
%}
%{
powerprob.Constraints.isOna1=isOn(9,1)==0;
powerprob.Constraints.isOna2=isOn(9,5)==0;
%
powerprob.Constraints.isOna3=isOn(10,1)==0;
powerprob.Constraints.isOna4=isOn(10,5)==0;
%
%powerprob.Constraints.isOna17=isOn(1,1)==0;
%powerprob.Constraints.isOna18=isOn(1,5)==0;
%{
powerprob.Constraints.isOna5=isOn(15,1)==0;
powerprob.Constraints.isOna6=isOn(15,5)==0;
%
powerprob.Constraints.isOna7=isOn(16,1)==0;
powerprob.Constraints.isOna8=isOn(16,5)==0;
%}
%powerprob.Constraints.isOna9=isOn(21,1)==0;
%powerprob.Constraints.isOna10=isOn(21,5)==0;

powerprob.Constraints.isOna13=isOn(22,1)==0;
powerprob.Constraints.isOna14=isOn(22,5)==0;

powerprob.Constraints.isOna15=isOn(23,1)==0;
powerprob.Constraints.isOna16=isOn(23,5)==0;

%powerprob.Constraints.isOna11=isOn(3,1)==0;
%powerprob.Constraints.isOna12=isOn(3,5)==0;
%}

%enforce startup=1 when moving from off to on
% no need to enforce startup=0 at other times since minimizing cost forces it
powerprob.Constraints.startupConst = -isOn1(idxHr2ToEnd-1,:) + isOn1(idxHr2ToEnd,:) - startup(idxHr2ToEnd,:) <= 0;

%Min uptime Constraint
powerprob.Constraints.minUptimeConst = optimconstr(nHours,actual_plant,50);
for l=1:50
for jj = 1:num_of_gen
    for kk = 1:nHours % based on possible startups; no penalty at end for running over
        if kk > nHours-MinimumUpTime(jj)
            sumidx = kk:nHours;
        else
            sumidx = kk:kk+MinimumUpTime(jj)-1;
        end
        powerprob.Constraints.minUptimeConst(kk,jj,l) = ...
            startup(kk,jj+(l-1)*4) - sum(isOn1(sumidx,jj+(l-1)*4)/length(sumidx)) <= 0;
    end
end

%Min downtime Constraint
powerprob.Constraints.minDowntimeConst = optimconstr(nHours,actual_plant,50);
for jj = 1:num_of_gen
    for kk = 2:nHours % based on possible startups; no penalty at beginning
        if kk <= MinimumDownTime(jj)
            sumidx = 1:kk-1;
        else
            sumidx = kk-MinimumDownTime(jj):kk-1;
        end
        powerprob.Constraints.minDowntimeConst(kk,jj,l) = ...
            startup(kk,jj+(l-1)*4) + sum(isOn1(sumidx,jj+(l-1)*4)/length(sumidx)) <= 1;
    end
end
end
%==========================END OF CONSTRAINTS==============================
%==========================END of FORMULATION==============================

% options for the optimization algorithm, here we set the max time it can run for
options = optimoptions('intlinprog','MaxTime',60000);
% call the optimization solver to find the best solution
[sol,TotalCost,exitflag,output] = solve(powerprob,'Options',options);

if test_case==3
TUCC_s3(u)=TotalCost; %Total unit commitment cost
co2_s3(u)=sol.powerCO2+sol.isOnCO2+sum(sol.startup(1:96))*1602.3; %CO2 emissions
elseif test_case==4
TUCC_s4(u)=TotalCost; %Total unit commitment cost
co2_s4(u)=sol.powerCO2+sol.isOnCO2+sum(sol.startup(1:96))*1602.3; %CO2 emissions
else
TUCC_s5(u)=TotalCost; %Total unit commitment cost
co2_s5(u)=sol.powerCO2+sol.isOnCO2+sum(sol.startup(1:96))*1602.3; %CO2 emissions
end
u=u+1;

end