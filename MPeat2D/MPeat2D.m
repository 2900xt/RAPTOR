%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% MPeat2D Program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
%MPeat2D is a fully coupled mechanical, ecological, and hydrological model
%of long-term peatland growth in two dimensions, which takes spatial 
%variability and structure into consideration. MPeat2D is developed based on 
%poroelasticity that couples fluid flow and solid deformation. 
%The formulation of MPeat2D is divided into mechanical,
%ecological, and hydrological submodels.

%First version 2023

%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
t_extent = 5000; %Simulation time (year)
oxic_decay = 5e-2; %Unsaturated zone decay rate(yr^{-1})
anoxic_decay = 8e-5; %Saturated zone decay rate(yr^{-1})
density_initial = 50; %Initial peat bulk density (kg m^{-3})
density_param = 3; %Bulk density parameter (-)
porosity_initial = 0.8; %Initial active porosity (-) 
porosity_param = 2; %Active porosity parameter (-) 
hydraulic_initial = 1e-2/3.171e-8;  %Hydraulic conductivity initial value (m yr{-1}) 
hydraulic_param = 15; %Hydraulic conductivity parameter (-)
modulus_param_1 = 4e5; %Young’s modulus parameter 1 (Pa) 
modulus_param_2 = 0.1; %Young’s modulus parameter 2 (-)
nu =0.2; %Poisson's ratio (-) 
ss_sat_initial = 1.4e-2; %Initial specific storage (m^{-1})
carbon_content = 0.47; %Carbon content (-)
sw_un = 0.4; %Degree of saturation of water in the unsaturated zone(-)
sw_sat = 1; %Degree of saturation of water in the saturated zone (-)
weight_water = 9800; %Specific weight of water (N m^{-3})
water_retention_1 = 0.5; %Water retention empirical constant 1 (-)
water_retention_2 = 0.4; %Water retention empirical constant 2 (m^{-1})
water_density = 997; %Water density (kg m^{-3})
shrub_wet_constant = 0.4; %Shrub wet condition constant (-)
sedge_wet_constant = 0.4; %Sedge or herb wet condition constant (-)
spagnum_wet_constant = 20; %Sphagnum wet condition constant (-)
a1 = 1.25; %Shrub-Young's modulus parameter
a2 = 1; %Sedge- or herb-Young's modulus parameter
a3 = 0.75; %Sphagnum-Young's modulus parameter
percentage_stress = 0.05; %Proportion of residual stress (Pa)
gravity = 9.8; %Gravitational acceleration (m s^{-2})
outside_density = 200; %Bulk density outside the peatland domain (kg m^{-3})
outside_porosity = 0.2; %Active porosity outside the peatland domain (-)
outside_hydraulic = 1e-8/3.171e-8; %Hydraulic conductivity outside the peatland domain (m yr{-1})
outside_modulus = 1e6; %Young's modulus outside peatland domain (Pa)
outside_ss = 1.4e-2; %Specific storage outside peatland domain (m^{-1})
%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Domain and Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
L = 500; %Domain length (m)
W = 500; %Domain width (m)
outside_length = 50; %Outside doamain length (m)
total_length = L+outside_length; %Total length (m)
grid_spatial = 50; %Number of grid
nx = grid_spatial; %Grid point in the horizontal direction
ny = grid_spatial; %Grid point in the vertical direction
za = 1; %Annual time counter
zb = 0; %Layer index counter
time_counter = 0;
hydrological_iteration = 1000; 
mechanic_time_average = 100;
annual_timestep_mechanic = 10; 
dx = L/nx;
dy = W/ny;
U = dx:dx:L;
v = dy:dy:W;
[x,y] = meshgrid(U,v);
N = nx;
annual_tsteps = 10;
timestep = 1/(annual_tsteps);
mechanic_timestep = 1/annual_timestep_mechanic;
%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Initial Array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
H_total = zeros(nx,ny); %Initial array for water table height
layer_mass = ones(grid_spatial,1)*(9.3^2)*0.001; %Initial array for layer mass
layer_initial_mass =layer_mass; %Initial array for layer initial mass
layer_remaining_mass = ones(grid_spatial,1); %Initial array for layer remaining mass
wet_proportion = ones(grid_spatial,1); %Initial array for wet condition
layer_thickness = layer_mass/density_initial; %Initial array for layer thickness
layer_elevation = layer_thickness; %Initial array for layer elevation
%-------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% Main Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------------
for ma = 1 : t_extent/mechanic_time_average
    mechanic_time(ma) = ma*mechanic_time_average;
end

while (1)

while (1)
time_counter = time_counter + 1;
if(time_counter > annual_tsteps)
  break
end
bog_height = 0;
%The calculation of layer mass, layer initial mass, layer remaining mass,
%layer thickness, layer elevation, and wet proportion of the peat.
%The loop for the vertical direction is nested inside the loop for lateral or horizontal direction 
for sa = 1:grid_spatial 
 if za ==1
 layer_mass(sa,1) = (9.3^2)*0.001;
 layer_initial_mass(sa,1) = (9.3^2)*0.001;
 layer_remaining_mass(sa,1) = 1;
 wet_proportion(sa,1) = 1;
 layer_thickness(sa,1) = layer_mass(sa,1)/density_initial;
 layer_elevation(sa,1) = layer_thickness(sa,1);
 else
 layer_mass(sa,za) = layer_mass_spatial(sa);
 layer_initial_mass(sa,za) = layer_initial_mass_spatial(sa);
 layer_remaining_mass(sa,za) = layer_remaining_mass_spatial(sa);
 wet_proportion(sa,za) = wet_proportion_spatial(sa);
 layer_thickness(sa,za) = layer_thickness_spatial(sa);
 layer_elevation(sa,za) = layer_elevation_spatial(sa);    
 end
 
%The changes of layer mass  with constant decay rate for each zone 
for zb = 1: za 
layer_mass(sa,zb) =(layer_mass(sa,zb) .* wet_proportion(sa,zb).*(exp(- anoxic_decay .* timestep)))+(layer_mass(sa,zb) .*(1 - wet_proportion(sa,zb)).*(exp(- oxic_decay .* timestep )));
layer_remaining_mass(sa,zb) = layer_mass(sa,zb) ./ max(layer_mass,[],'all');
if za > mechanic_time_average
    layer_thickness(sa, zb) = layer_mass(sa,zb) ./ mean(density); %Layer thickness calculation
else
    layer_thickness(sa,zb) = layer_mass(sa,zb) ./ density_initial;
end
if(zb == 1)
    layer_elevation(sa,zb) = layer_thickness(sa,zb); 
else
    layer_elevation(sa,zb) = layer_thickness(sa,zb) + layer_elevation(sa,zb - 1); %Layer elevation calculatio
end

end
%The peatland height based on the layer elevation value 
    bog_height = layer_elevation(sa,za);
    peatland_height(sa,1) = bog_height; 
    peatland_height_time(za,sa) = bog_height;
end

end
time_counter = 0;

%The changes in peat physical properties that affect water table position
 if za > mechanic_time_average
     hydraulic_vector_wt = hydraulic_update_wt'; %The variations of hydraulic conductivity due to mechanical compaction
     porosity_vector_wt = porosity_update_wt'; %The variations of active porosity due to mechanical compaction
     trans_factor = mean(density)/density_initial;
 else
     hydraulic_vector_wt = ones(N*N,1).*(hydraulic_initial);
     porosity_vector_wt = ones(N*N,1).*(porosity_initial);
     trans_factor = 1;
 end

%%Hydrological Submodel Main Calculation%% 
 M_vector = ss_sat_initial.*(dx.^2)./(hydraulic_vector_wt);
 d = [-N,-1,0,1,N];
 v1 = [ones(1,(N^2)-N), zeros(1,N)];
 v21 = ones(1,N*N-1);
 for i =1:N-1
     v21(i*N)=0;
 end
 v2 = [v21,0];
 
 v3 = ones(1,N*N);
 v3(1) = -(2+M_vector(1));
 v3(N*N) = -(6+M_vector(N*N));
 v3(N) = -(4+M_vector(N));
 v3(N*N-N+1) = -(4+M_vector(N*N-N+1));
 for i = 2:N-1
    v3(i) = -(3+M_vector(i));
    v3((i-1)*N+1)=-(3+M_vector((i-1)*N+1));
    v3(N*N-N+i)=-(5+M_vector(N*N-N+i));
    v3(i*N)=-(5+M_vector(i*N));
 end
 
 for i=1:N-2
    for j=2:N-1
        v3(i*N+j)=-(4+M_vector(i*N+j));
    end
 end
 
 v4=[0,v21];
 v5=[zeros(1,N),ones(1,(N^2)-N)];
 
 B=[v1',v2',v3',v4',v5'];
 A=spdiags(B,d,(N^2),(N^2));
 
 %Net rainfall as a source term for the water table 
 rain_source = ones(N*N,1)*(net_rainfall(za));
 rain_source(N*N)=0;
 rain_source(1)= (net_rainfall(za));
 rain_source(N)= 0;
 rain_source(N*N-N+1)= 0;
 for i=2:(N-1)
     rain_source(i)= (net_rainfall(za));
     rain_source(i*N)= 0;
     rain_source(N*N-N+i)= 0;
     rain_source((i-1)*N+1)= (net_rainfall(za));
 end
 
 h=zeros(N*N,1);

 for i=1:hydrological_iteration
     b=(-M_vector.*h-(rain_source/trans_factor))/hydrological_iteration;
     h_new=A\b;
     h=h_new/hydrological_iteration;
     H=reshape(h, N, []).';
 end
 %The position of water table
 H_new = H;
 H_total = H_total+H_new;
 wt_height = H_total(:,1);
 
%Check position of water table relative to the peatland surface
for sb = 1:grid_spatial
if wt_height(sb) > peatland_height(sb) 
    wt_height(sb) = peatland_height(sb);
end

%Calculate the water-table depth.
wt_depth(sb) = peatland_height(sb) - wt_height(sb);
wt_depth_time (za,sb) = wt_depth(sb);

%Plant functional type estimation consiting of shrub, sedge, and sphagnum
    shrub_proportion(sb) = (2.23*wt_depth(sb))-0.28; %shrub linear regression
    sedge_proportion(sb) = (-1.42*wt_depth(sb))+0.63; %sedge linear regression
    spagnum_proportion(sb) = (-0.81*wt_depth(sb))+0.64; %sphagnum linear regression

%Check the PFT proportion    
if shrub_proportion(sb) < 0
    shrub_proportion_transition(sb) = 0;
else
    shrub_proportion_transition(sb) = shrub_proportion(sb);
end
if sedge_proportion(sb) < 0
    sedge_proportion_transition(sb) = 0;    
else
    sedge_proportion_transition(sb) = sedge_proportion(sb);   
end
if spagnum_proportion(sb) < 0
    spagnum_proportion_transition(sb) = 0; 
else
    spagnum_proportion_transition(sb) = spagnum_proportion(sb); 
end
PFT_total_proportion(sb) = shrub_proportion_transition(sb)+sedge_proportion_transition(sb)+spagnum_proportion_transition(sb);
shrub_proportion_final(sb) = shrub_proportion_transition(sb)/PFT_total_proportion(sb);
sedge_proportion_final(sb) = sedge_proportion_transition(sb)/PFT_total_proportion(sb);
spagnum_proportion_final(sb) = spagnum_proportion_transition(sb)/PFT_total_proportion(sb);

shrub_proportion_final_time(za,sb) = shrub_proportion_final(sb);
sedge_proportion_final_time(za,sb) = sedge_proportion_final(sb);
spagnum_proportion_final_time(za,sb) = spagnum_proportion_final(sb);

% The wet proportion calculation
for zb = 1: za
if zb == 1
    if(wt_height(sb) >= layer_elevation(sb,zb))
    wet_proportion(sb,1) = 1;
    else
    wet_proportion(sb,1) = wt_height(sb) ./ layer_thickness(sb,1);
    end
elseif wt_height(sb) >= layer_elevation(sb,zb) 
    wet_proportion(sb,zb) = 1;
elseif wt_height(sb) <= layer_elevation(sb,zb - 1) 
    wet_proportion(sb,zb) = 0;
else
    wet_proportion(sb,zb) =(wt_height(sb) - layer_elevation(sb,zb - 1))./ layer_thickness(sb,zb);
end

end
end
 
%Check the total time 
za = za + 1;
if za > t_extent;
  break;
end
%New layer initial properties and plant weight
for sc=1:grid_spatial
    if wt_depth(sc) > 0.668
        layer_mass(za) = 0.0000001;
    else
        layer_mass(za) =((9.3 +(133 * wt_depth(sc)) -(0.022*((wt_depth(sc)*100)^2)))^2)*0.001*((0.1575 * temp(za)) + 0.009132);
    end

layer_mass_spatial(sc) = layer_mass(za);
layer_initial_mass_spatial(sc) = layer_mass(za);
layer_remaining_mass_spatial(sc) = 1;
wet_proportion_spatial(sc) = 0;
layer_addition_weight(sc) = layer_initial_mass_spatial(sc)*gravity;
plant_weight(sc) = (((((10^(((log10(layer_initial_mass_spatial(sc)))+0.409)/0.985))*shrub_proportion_final(sc)*(1+shrub_wet_constant))+((10^(((log10(layer_initial_mass_spatial(sc)))+0.001)/1))*sedge_proportion_final(sc)*(1+sedge_wet_constant))))+((0.144*spagnum_proportion_final(sc))*(1+spagnum_wet_constant)))*gravity;
plant_weight_time(za,sc) = plant_weight(sc);
layer_thickness_spatial(sc) = layer_mass_spatial(sc) ./ density_initial;
layer_elevation_spatial(sc) = peatland_height(sc,1) + layer_thickness_spatial(sc);
end


%%Mechanical Submodel Main Calculation%% 
mb = ismember(za,mechanic_time);
if mb == 1
%Mesh for finite element calculation
substrate_node = H_total(:,grid_spatial);
x_substrate = U';
y_substrate = zeros(ny,1);
x_saturated = reshape(x',[],1);
y_saturated = reshape(H_total,[],1);
x_unsaturated = U';
y_unsaturated = peatland_height;
x_outside = L+dx : dx : total_length;
y_outside = linspace(0,min(y_unsaturated), ny);
x_centre = zeros(ny,1);
y_centre = linspace(max(substrate_node),max(y_unsaturated), ny);
[x_outside_combine,y_outside_combine] = meshgrid(x_outside,y_outside);
x_outside_final = reshape(x_outside_combine,[],1);
y_outside_final = reshape(y_outside_combine,[],1);
x_mesh_final = [x_saturated;x_unsaturated; x_outside_final; x_centre];
y_mesh_final = [y_saturated;y_unsaturated; y_outside_final; y_centre'];
mesh_grid = delaunay(x_mesh_final,y_mesh_final);
node_number = size(x_mesh_final,1);
element_number = size(mesh_grid,1);
% Finite element input 
input.constant      = [node_number element_number 3 3 2 2];
nodes               = node_number; %Number of nodes
numel               = element_number; %Number of element
nodaldof            = 3; %Number of degree of freedom (vertical displacement, horizontal displacement, and pore pressure)
nodesperelement     = 3; %Number of nodes per element
numat               = element_number; %Number of material for nonlinear poroelasticity
ndimensions         = 2; %Number of dimension
numeq               = nodes*nodaldof;
totalelementdof     = nodaldof*nodesperelement;

Kg  = zeros(numeq,numeq);
Cg  = zeros(numeq,numeq);
Mg  = zeros(numeq,numeq);
f   = zeros(numeq,1);
u   = zeros(numeq,1);

%Bottom, outter, and outside node
outter_node_initial = ismember(y_mesh_final,y_unsaturated);
outside_node_initial = ismember(y_mesh_final,y_outside_final);
buttom_node_initial = ismember(y_mesh_final,substrate_node);
for aa = 1:node_number
    if buttom_node_initial(aa) == 1 
       bottom_layer_sign(aa) = 1 ;
    else
       bottom_layer_sign(aa) = 0 ;
    end
    if outter_node_initial(aa) ==1 
        outter_layer_sign(aa) = 1;
    else
        outter_layer_sign(aa) = 0;
    end
    if  outter_layer_sign(aa) == 1 
        outter_node(aa) = aa;
    end
    if outside_node_initial(aa) ==1
        outside_layer_sign(aa) = 1;
    else
        outside_layer_sign(aa) = 0;
    end
    if x_mesh_final(aa) == 0
        centre_layer_sign(aa) = 1;
    else
        centre_layer_sign(aa) = 0;
    end
    if  outter_layer_sign(aa) == 1 | outside_layer_sign(aa) == 1 |centre_layer_sign(aa) == 1;
        node_bottom_outter(aa) =1;
    else
        node_bottom_outter(aa) =0;
    end
end
total_node_un = sum(outter_layer_sign);

%Node element
node_element_1 = ismember(mesh_grid(:,1),outter_node);
node_element_2 = ismember(mesh_grid(:,2),outter_node);
node_element_3 = ismember(mesh_grid(:,3),outter_node);

for gg = 1:element_number
    if  node_element_1(gg) == 1 
        node_combine(gg) = 1;
    elseif node_element_2(gg) == 1
        node_combine(gg) = 1;
    elseif node_element_3(gg) == 1
        node_combine(gg) = 1;
    else
        node_combine(gg) = 2;
    end
end
%Node coordinate 
for bb = 1:node_number    
    first_column_ne(bb) = bb;  
    second_column_ne(bb) = x_mesh_final(bb,1);
    third_column_ne(bb) =  y_mesh_final(bb,1);
    node_data = [first_column_ne' second_column_ne' third_column_ne' ];
end
%Element coordinate 
for cc = 1:element_number
    first_column_et(cc) = cc;
    second_column_et(cc) = cc; 
    third_column_et(cc) = mesh_grid(cc,1);
    fourth_column_et(cc) = mesh_grid(cc,2);
    fifth_column_et(cc) = mesh_grid(cc,3);
    element_data = [first_column_et' second_column_et' third_column_et' fourth_column_et' fifth_column_et'];
end 

%Constrain or boundary condition [node_number x_displacement y_displacement pore_water_pressure]
% 1 fixed displacement or pore water pressure
% 0 non-fixed displacement or pore water pressure
for dd = 1:node_number  
    first_column_cn(dd) = dd;
    if x_mesh_final(dd) == min(x_mesh_final)
        second_column_cn(dd) = 1;
    elseif x_mesh_final(dd) == max(x_mesh_final)
        second_column_cn(dd) = 1;
    else
        second_column_cn(dd) = 0;   
    end
    if bottom_layer_sign(dd) == 1
        third_column_cn(dd) = 1;
    elseif x_mesh_final(dd) == min(x_mesh_final) & y_mesh_final(dd) == max(substrate_node)
        third_column_cn(dd) = 1;
    else 
        third_column_cn(dd) = 0;
    end
    if x_mesh_final(dd) == max(x_mesh_final)
        fourth_column_cn(dd) = 1;
    else
        fourth_column_cn(dd) = 0;
    end
    constrain_data = [first_column_cn' second_column_cn' third_column_cn' fourth_column_cn'];
end
%Load data 
%Magnitude of prescribed 
%load_solid : load or displacement
%load_fluid : pore_water_pressure or flux
sd = 1;
for ee = 1:node_number
            first_column_ld(ee) = ee;
            second_column_ld(ee) = 0;
            if x_mesh_final(ee) > L
                third_column_ld(ee) = 0;
            elseif x_mesh_final(ee) == 0
                third_column_ld(ee) = 0;
            elseif bottom_layer_sign(ee) == 1
                third_column_ld(ee) = 0;
            elseif outter_layer_sign(ee) == 1                
                third_column_ld(ee) = -((plant_weight(sd) +layer_addition_weight(sd))+((((1-porosity_initial)*density_initial) + (porosity_initial*water_density*sw_un))*gravity));
                sd = sd+1;
            elseif za == mechanic_time_average              
               third_column_ld(ee) = -((((1-porosity_initial)*density_initial) + (porosity_initial*water_density))*gravity);
            else
                third_column_ld(ee) = -((((1-porosity_update(ee))*density_update(ee)) + (porosity_update(ee)*water_density))*gravity);
            end
            if x_mesh_final(ee) == max(x_mesh_final)
                fourth_column_ld(ee) = 0;
            elseif outter_layer_sign(ee) == 1
                fourth_column_ld(ee) = net_rainfall(za);
            else
                fourth_column_ld(ee) = 0;
            end            
            load_data = [first_column_ld' second_column_ld' third_column_ld' fourth_column_ld'];          
end

%Peat physical properties calculation 
 for nn = 1:node_number
     if x_mesh_final(nn) > L
        modulus_update(nn) = outside_modulus;
        porosity_update(nn) = outside_porosity;
        density_update(nn) = outside_density;
        hydraulic_update(nn) = outside_hydraulic; 
        ss_update(nn) = outside_ss;
     elseif za == mechanic_time_average
        modulus_update(nn) = modulus_param_1*(1+(1^modulus_param_2))*(a1*mean(shrub_proportion_final_time, 'all')+a2*mean(sedge_proportion_final_time,'all')+a3*mean(spagnum_proportion_final_time,'all'));
        porosity_update(nn) = porosity_initial;
        density_update(nn) = density_initial;
        hydraulic_update(nn) = hydraulic_initial; 
        ss_update(nn) = ss_sat_initial;       
     elseif outter_layer_sign(nn) ==1
        modulus_update(nn) = modulus_param_1*(1+(1^modulus_param_2))*(a1*mean(shrub_proportion_final)+a2*mean(sedge_proportion_final)+a3*mean(spagnum_proportion_final));
        porosity_update(nn) = porosity_initial;
        density_update(nn) = density_initial;
        hydraulic_update(nn) = hydraulic_initial; 
        ss_update(nn) = porosity_update(nn)/(((weight_water*(1-water_retention_1))/(water_retention_1*water_retention_2))*(sw_un^(-1/water_retention_1))*((1-sw_un^(1/water_retention_1))^water_retention_1));
     else
        modulus_update(nn) = modulus_param_1*(1+(mean(layer_remaining_mass, 'all'))^modulus_param_2)*(a1*mean(shrub_proportion_final_time, 'all')+a2*mean(sedge_proportion_final_time,'all')+a3*mean(spagnum_proportion_final_time,'all'));
        porosity_update(nn) = (porosity_update(nn)+(porosity_param*(grad_total(nn))))/(1+(grad_total(nn)));
        density_update(nn) = density_update(nn)/(1+density_param*(grad_total(nn)));
        hydraulic_update(nn) = hydraulic_initial*((porosity_update(nn)/porosity_initial)^hydraulic_param);
        ss_update(nn) = ss_sat_initial;
     end
 end
 %Peat physical properties update for hydrological submodel
 hydraulic_update_wt = hydraulic_update;
 hydraulic_update_wt(node_bottom_outter==1)=[];
 porosity_update_wt = porosity_update;
 porosity_update_wt(node_bottom_outter==1)=[];
 
%Material array
 for kk=1:3
    for hh = 1:numat
        modulus_element(hh,kk) = modulus_update(mesh_grid(hh,kk));
        hydraulic_element(hh,kk) = hydraulic_update(mesh_grid(hh,kk));
        porosity_element(hh,kk) = porosity_update(mesh_grid(hh,kk));
        ss_element(hh,kk) = ss_update(mesh_grid(hh,kk));
    end
 end
  
 for mm=1:numat
     first_column_ml(mm) = mean(modulus_element(mm,:));
     second_column_ml(mm) = nu;
     third_column_ml(mm) = mean(hydraulic_element(mm,:));
     if node_combine(mm) == 1
         fourth_column_ml(mm) = sw_un;
     else
         fourth_column_ml(mm) = sw_sat;
     end
     fifth_column_ml(mm) = mean(ss_element(mm,:));
     sixth_column_ml(mm) = mean(porosity_element(mm,:));
 end
material =[first_column_ml' second_column_ml'  third_column_ml' fourth_column_ml'  fifth_column_ml' sixth_column_ml'];

%Global stifness matrix
for iel=1:numel
    [local_data] = localize (iel, node_data, element_data, material, input);
            ElementNo         = iel;
            MaterialArray     = local_data.mater;
            CoordinateArray   = local_data.coords;
            dofAddressArray   = local_data.dofs;
%Create local stiffness matrix             
    MaterialType = element_data(iel,2);   
    [Kl,Cl,Ml] = element(local_data);
%Assemble global stiffness matrix from local stiffness matrix   
    [Kg] = assemble (Kg, Kl, local_data, totalelementdof);
    [Cg] = assemble (Cg, Cl, local_data, totalelementdof); 
    [Mg] = assemble (Mg, Ml, local_data, totalelementdof);
    
end
%Assemble load and constrain  vectors for solid and fluid
[f]     = vectorstretch (load_data, nodes, 1+nodaldof);  
[con]   = vectorstretch (constrain_data,  nodes, 1+nodaldof) ; 
u  = f.* (con);                                        
f  = f.* (1-con);                                      
f  = f - Kg*u;                                         
Ktemp = Kg; 
Ctemp = Cg; 
Mtemp = Mg;       
ftemp = f;  
% Residual stress
if za == mechanic_time_average
      utemp =   zeros(numeq,1); 
else
      utemp_initial =  zeros(numeq,1);
          for tt = 1:size(utemp,1)
                utemp_initial(tt) = utemp(tt);
          end
      utemp = (utemp_initial/za)*percentage_stress;
end
%Reorder the mtrix       
for ieq = numeq: -1: 1
    if (con(ieq,1) == 1) Ktemp(:,ieq)=[]; Ktemp(ieq,:)=[]; 
                         Ctemp(:,ieq)=[]; Ctemp(ieq,:)=[]; 
                         Mtemp(:,ieq)=[]; Mtemp(ieq,:)=[]; 
                         ftemp(ieq,:)=[]; utemp(ieq,:)=[]; 
    end
end
%Solve system equations    
for mt = 1:annual_timestep_mechanic                
    utemp = (Ktemp + (1/mechanic_timestep)*Ctemp)\(ftemp + (1/mechanic_timestep)*(Ctemp*utemp));
    ftemp = ftemp*0;                
    warning('off')
end
icount=0;                                                   
for ieq = 1:numeq
    if (con(ieq,1) == 0) icount = icount + 1; u(ieq,1) = utemp(icount,1);
    end
end
u;    
%Displacement and pore water pressure
for ff=1:node_number 
    dis_x(ff)=u((3*ff)-2);
    dis_y(ff)=u((3*ff)-1);
    press(ff)=u(ff*3);
end   
dis_x_new = dis_x;
dis_y_new = dis_y;   
if za == mechanic_time_average
   dis_x_total = dis_x_new;
   dis_y_total = dis_y_new;
else
   dis_x_total = dis_x_total + dis_x_new;
   dis_y_total = dis_y_total + dis_y_new; 
end  
%Strain calculation
[strainxx, strainxy] = trigradient((x_mesh_final),(y_mesh_final),(dis_x_new'));
[strainyx, strainyy] = trigradient((x_mesh_final),(y_mesh_final),(dis_y_new'));
shear_strain = 0.5*((strainxy)+(strainyx));
grad_total = ((strainxx) + (strainyy));
      
for qq =1:node_number
    if x_mesh_final(qq) < dx | x_mesh_final(qq) > L
        remove_node(qq) =1;
    else
        remove_node(qq) =0;
    end
end
%Remove outside node
x_coordinate_final = x_mesh_final';
x_coordinate_final(remove_node==1)=[];
y_coordinate_final = y_mesh_final';
y_coordinate_final(remove_node==1)=[];
porosity = porosity_update;
porosity(remove_node==1)=[];
density = density_update;
density(remove_node==1)=[];
hydraulic = hydraulic_update;
hydraulic(remove_node==1)=[];
modulus = modulus_update;
modulus(remove_node==1)=[]; 

for ss = 1:grid_spatial
    porosity(N*N-N+ss) = porosity(N*N-2*N+ss);
    density(N*N-N+ss) = density(N*N-2*N+ss);
    hydraulic(N*N-N+ss) = hydraulic(N*N-2*N+ss);
    modulus(N*N-N+ss) = modulus(N*N-2*N+ss);
    substrate_porosity(ss) = porosity(N*N-N+ss);
    substrate_density(ss) = density(N*N-N+ss);
    substrate_hydraulic(ss) = hydraulic(N*N-N+ss);
    substrate_modulus(ss) = modulus(N*N-N+ss); 
end
%Peat physical properties plot
x_plot = [x_coordinate_final x_substrate'];
y_plot = [y_coordinate_final y_substrate'];
porosity_plot = [porosity substrate_porosity];
density_plot = [density substrate_density];
hydraulic_plot = [hydraulic substrate_hydraulic];
modulus_plot = [modulus substrate_modulus];
end 

end
%Carbon calculation at the centre
peatland_height_centre = peatland_height_time(:,1);
cumulative_carbon_centre = peatland_height_centre.*mean(density_plot).*carbon_content;



