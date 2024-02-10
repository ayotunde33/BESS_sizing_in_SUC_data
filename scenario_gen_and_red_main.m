%load wind_w, solar_w and load_w into workspace
wind_speed_grouped=grouping_function(wind_w,24);
solar_irr_grouped=grouping_function(solar_w,24);
load_grouped=grouping_function(load_w,24);

generated_scen_load=scenario_generation(load_grouped,1000);
generated_scen_wind=scenario_generation(wind_speed_grouped,1000);
generated_scen_solar=scenario_generation(solar_irr_grouped,1000);
ranked_scenarios_load=kantorovich_ranking(generated_scen_load,1000);
ranked_scenarios_wind=kantorovich_ranking(generated_scen_wind,1000);
ranked_scenarios_solar=kantorovich_ranking(generated_scen_solar,1000);
%

ranked_scenarios_load_new=ranked_scenarios_load;
ranked_scenarios_wind_new=ranked_scenarios_wind;
ranked_scenarios_solar_new=ranked_scenarios_solar;
lx=1;
gx=1;
while lx<=37
    for ix=1:27
lh=rand(1)/3;
wh=rand(1)/3;
sh=rand(1)/3;
lm=0.333+rand(1)/3;
wm=0.333+rand(1)/3;
sm=0.333+rand(1)/3;
ll=0.666+rand(1)/3;
wl=0.666+rand(1)/3;
sl=0.666+rand(1)/3;
load_high=(abs(lh-ranked_scenarios_load_new));
wind_high=(abs(wh-ranked_scenarios_wind_new));
solar_high=(abs(sh-ranked_scenarios_solar_new));
load_medium=(abs(lm-ranked_scenarios_load_new));
wind_medium=(abs(wm-ranked_scenarios_wind_new));
solar_medium=(abs(sm-ranked_scenarios_solar_new));
load_low=(abs(ll-ranked_scenarios_load_new));
wind_low=(abs(wl-ranked_scenarios_wind_new));
solar_low=(abs(sl-ranked_scenarios_solar_new));

if ix==1 %load high wind high solar high
lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==2 %load high wind high solar medium
lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==3 %load high wind high solar low
lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;    
elseif ix==4 %load high wind medium solar high
lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==5 %load high wind medium solar medium
    lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==6 %load high wind medium solar low
    lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==7 %load high wind low solar high
lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==8 %load high wind low solar medium 
lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==9 %load high wind low solar low
    lh_a=find(load_high==min(load_high));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==10 %load medium wind high solar high
    lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==11 %load medium wind high solar medium
     lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==12 %load medium wind high solar low
     lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==13 %load medium wind medium solar high
     lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==14 %load medium wind medium solar medium
 lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==15 %load medium wind medium solar low
 lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==16 %load medium wind low solar high
 lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==17 %load medium wind low solar medium
     lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==18 %load medium wind low solar low
 lh_a=find(load_medium==min(load_medium));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==19 %load low wind high solar high
     lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==20 %load low wind high solar medium
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==21 %load low wind high solar low
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_high==min(wind_high));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==22 %load low wind medium solar high
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==23 %load low wind medium solar medium
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==24 %load low wind medium solar low
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_medium==min(wind_medium));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==25 %load low wind low solar high
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_high==min(solar_high));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
elseif ix==26 %load low wind low solar medium
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_medium==min(solar_medium));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
else %load low wind low solar low
    lh_a=find(load_low==min(load_low));
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(lh_a));
ranked_scenarios_load_new(lh_a)=[];
wh_a=find(wind_low==min(wind_low));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(wh_a));
ranked_scenarios_wind_new(wh_a)=[];
sh_a=find(solar_low==min(solar_low));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(sh_a));
ranked_scenarios_solar_new(sh_a)=[];
gx=gx+1;
end

%{
if size(find(lh_i==lh_i(l)),2)>0

lh_rank(l)=ranked_scenarios_load(lh_i(l));
end
%}
    end
    lx=lx+1;

end
load_i(gx)=find(ranked_scenarios_load==ranked_scenarios_load_new(1));
wind_i(gx)=find(ranked_scenarios_wind==ranked_scenarios_wind_new(1));
solar_i(gx)=find(ranked_scenarios_solar==ranked_scenarios_solar_new(1));
gx=gx+1;

[load_sel_scen, wind_sel_scen, solar_sel_scen, probVeci]=scenario_reduction(ranked_scenarios_load,...
    ranked_scenarios_wind,ranked_scenarios_solar,load_i,wind_i,solar_i,...
    generated_scen_load,generated_scen_wind,generated_scen_solar,50);
%{
load high, wind high, solar high
load high, wind high, solar medium
load high, wind high, solar low
load high, wind medium, solar high
load high, wind medium, solar medium
load high, wind medium, solar low
load high, wind low, solar high
load high, wind low, solar medium
load high, wind low, solar low

load medium, wind high, solar high
load medium, wind high, solar medium
load medium, wind high, solar low
load medium, wind medium, solar high
load medium, wind medium, solar medium
load medium, wind medium, solar low
load medium, wind low, solar high
load medium, wind low, solar medium
load medium, wind low, solar low

load low, wind high, solar high
load low, wind high, solar medium
load low, wind high, solar low
load low, wind medium, solar high
load low, wind medium, solar medium
load low, wind medium, solar low
load low, wind low, solar high
load low, wind low, solar medium
load low, wind low, solar low

%}