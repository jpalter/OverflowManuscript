%plot_outflow_bathymetry

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run this code with the following choices in order to: 
% 
% a) Create the trajectory maps in Figure 3 for the streamtube model.  Do 
% not loop over all temperatures and salinities.  Rather, choose default 
% source water temperatures and salinities (i.e. those from Table S1) and 
% run the solely the inner content of the loops without looping.  Mapping 
% scripts provided.
%
% b) Create the colors representing the streamtube model transports for the
% Denmark Strait in Figure 4a.  Run the first set of nested loops.  Save
% output  Si CTi Ti, which are, respectively, the Absolute Salinity,
% Conservative Temperature and Transports (length = ntimesteps, requires
% averaging).  These will be plotted in the Python code, which deals with
% the CMIP6 output.
%
% c) Create the colors representing the streamtube model transports for the
% Faroe Bank Channel.  As with the previous step, save outputs to use in
% CMIP6 script
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bathy = netcdf_load('gebco_2025_n70.0_s61.0_w-40.0_e-3.0.nc');
D_map = double(bathy.elevation); %take depth as positive number (multiply elevation by -1)

%Smooth bathymetry to stop runaway entrainment
% Simple 3x3 moving average (replace with Gaussian if desired)
D_map_smooth = conv2(D_map, ones(10)/10^2, 'same');

%Beautiful bathymetry map
figure
axesm('MapProjection', 'miller', 'MapLatLimit', [61 69], ...
'MapLonLimit', [-40 -3])
pcolorm(bathy.lat, bathy.lon, D_map_smooth')
shading flat

setm(gca, 'MLabelLocation', 5, 'PLabelLocation', 2)
setm(gca, 'MLineLocation', 5, 'PLineLocation', 2)
cmocean('topo','pivot',0);colorbar
set(gcf, 'Color', 'w')
setm(gca,'fontsize',14)
colorbar

gridm on
mlabel on
plabel on
tightmap

%Now make bathymetry transect along Greenland-Scotland Ridge between Faroe
%Islands -> Iceland -> Greenland

% Define transect waypoints
E = -3.1;
midlo = -18;
W=-39;
S = 61.1;
midla = 65.5;
N = 67;

% Create fine-resolution path along transect
n_pts = 200;
seg1_lat = linspace(S, midla, n_pts);
seg1_lon = linspace(E, midlo, n_pts);
seg2_lat = linspace(midla, N, n_pts);
seg2_lon = linspace(midlo, W, n_pts);

transect_lat_fine = [seg1_lat seg2_lat(2:end)];
transect_lon_fine = [seg1_lon seg2_lon(2:end)];

hold on;
%plotm(seg1_lat,seg1_lon,'k')
%plotm(seg2_lat,seg2_lon,'k')


%Plot bathymetry along transect
[loni, lati] = meshgrid(double(bathy.lon), double(bathy.lat));  % lon first, then lat


depth_transect = interp2(loni, lati, D_map', ...
                 transect_lon_fine, transect_lat_fine, ...
                 'nearest');
x = transect_lon_fine(:).';      % longitude
z = depth_transect(:).';         % negative depths

% Choose a bottom level a bit deeper than the deepest point
z_bottom = min(z) * 1.05;        % e.g. -4200 if min(z) ~ -4000

figure('Color','w')

% Polygon: along seafloor, then down to bottom, back along bottom, up to start
x_poly = [x,                fliplr(x)];
z_poly = [z, fliplr(repmat(z_bottom, size(z)))];

fill(x_poly, z_poly, [0.7 0.7 0.7], 'EdgeColor','none')  % gray seafloor
hold on

xlabel('Longitude (°)')
ylabel('Depth (m)')

xlim([min(x) max(x)])
ylim([z_bottom 0])           % 0 (surface) at top, z_bottom at bottom

ax = gca;
ax.Color = 'w';              % water is white
box on


%Now run default Price and Barringer
n_cts = 20;
n_sts = 22;
CT_range = linspace(-.1,5.5,n_cts);
S_range = linspace(34.5,35.25,n_sts);
n_time_steps = 701;

 Us= NaN(n_cts,n_sts,n_time_steps);
 Hs= NaN(n_cts,n_sts,n_time_steps);
 Ws= NaN(n_cts,n_sts,n_time_steps);
 p_mid = NaN(n_cts,n_sts,n_time_steps);
 distance= NaN(n_cts,n_sts,n_time_steps);
 TTTs= NaN(n_cts,n_sts,n_time_steps);
 SSSs= NaN(n_cts,n_sts,n_time_steps);


  %load bathymetry for the region from GEBCO2025; %Note: Hiroki had to do
        %massive smoothing to deal with numerical noise. We will see if that's
        %necessary
      
        
for CTs =1:n_cts  
    for Ss =1:n_sts  
        CTs
        %Define and initialize variableso

        dt = 60*60*24/50; %400 time steps per day
        %Ambient water properties
        theta_ambient = 3.7; %User defined ambient potential temperature [°C]
        S_ambient = 34.9; %User defined ambient absolute salinity [g/kg]

        %Source water properties
        theta_source =CT_range(CTs);%0;% 0.5; %User defined source water potential temperature [°C]; 0.3°C from Eldevick 2009- similar in Sun et al 2025; set to 0°C in P&B94
        S_source = S_range(Ss);%34.9 %User defined source water practical salinity [g/kg]; 34.9 from PB94 (similar in Eldevick et al 2009 and Sun et al 2025)

        %constants
        rho_not = 1030; %user-defined reference density, a constant [kg/m3]
        Cd = 4e-3;%3e-3; % Drag coefficient (~50% uncertainty) [unitless]
        g = 9.8; %gravitational constant [m/s2]

        %site parameters and initial conditions
        H = NaN(n_time_steps,1);
        W = NaN(n_time_steps,1);
        H(1) = 150; %initial thickness of overflow, 200 from Hiroki; 150 m from Table 1 in P&B94 [m]
        W(1) = 50e3; %initial width of overflow, 80km from Hiroki; 50 km from Table 1 in P&B94 [m]

        lon_initial = -27.90; %define grid cell in bathymetry? 65.93°N, 27.90°W - Hiroki
        lat_initial = 65.93;
        initial_direction = 225; %an angle representing the dirction of flow, taken from Hiroki and P&B94
        initial_transport = 3e6; % [m3/s]
        initial_speed = 0.01;%0.35;%initial_transport/(H(1)*W(1)); %[m/s] Taken

        %pre-allocate arrays, so they don't grow in loop
        u = NaN(n_time_steps,1);
        v = NaN(n_time_steps,1);
        x = NaN(n_time_steps,1);
        y = NaN(n_time_steps,1);
        u(1) = initial_speed*cosd(initial_direction);
        v(1) = initial_speed*sind(initial_direction);
        x(1) =  lon_initial;
        y(1) =  lat_initial;

        f = gsw_f(y(1));
       

        D = NaN(n_time_steps,1);
        D(1) = interp2(bathy.lon,bathy.lat,D_map_smooth',x(1),y(1));

        %calculate in situ density from user-defined theta and S for ambient and
        %source waters:

        % Mid-depth (in meters)
        depth_mid = D(1) + H(1)/2;
        % Convert to pressure (dbar)
        p_mid = gsw_p_from_z(depth_mid, y(1));
        depth_top = D(1)+H(1);
        p_top=gsw_p_from_z(depth_top, y(1));

        CT_ambient = gsw_CT_from_pt(S_ambient,theta_ambient);
        S_ambient = gsw_SA_from_SP(S_ambient,p_top,lon_initial,lat_initial);

        rho_ambient = gsw_rho(S_ambient, CT_ambient, p_top);
        
        S_source = gsw_SA_from_SP(S_source,p_mid,lon_initial,lat_initial);
        CT_source = gsw_CT_from_pt(S_source,theta_source);

        rho_source = gsw_rho(S_source,CT_source,p_top);
        rho_dense(1) = rho_source;
        rho_dense(2:n_time_steps) = NaN;

        %initial conditions of the overflow are set to the source water properties
        CT_dense = NaN(n_time_steps,1);
        S_dense = NaN(n_time_steps,1);

        CT_dense(1) = CT_source;
        S_dense(1) = S_source;
        E = NaN(n_time_steps,1);

        %The delrho term is the difference between the in situ density between the
        %outflow water and the ambient oceanic water.  This is defined in 2 ways,
        %but the one relevant for entrainment, delrho+, is taken as the difference between
        %the water at the top of the outflow (D-H, i.e. ocean depth minus the
        %thickness of the outflow).  The other delrho (no plus sign) is the density
        %difference between the outflow and oceanic water at mid-depth of the
        %outflow. I'm setting a single constant value for
        %the ambient seawater, rather than defining an entire profile.
        %This a simplification, but not a big one, given week
        %stratification at the depth of the outflow.

        %in the first time step, the delta is between source water and ambient
        %water.  In the loop, the source water will be replaced by the dense water
        %characteristics, which change vian entrainment/mixing

        delrho_plus = NaN(n_time_steps,1);
        delrho_plus(1) = rho_source - rho_ambient;

        deltheta_plus = NaN(n_time_steps,1);
        deltheta_plus(1) = CT_source - CT_ambient;
        delS_plus = NaN(n_time_steps,1);
        delS_plus(1) = S_source - S_ambient;

        %Here define the delrho used in the momentum equation (these are different
        %only due to the different pressures used in the in situ density
        %calculation)
        delrho_momentum = NaN(n_time_steps,1);

        rho_dense_mid = gsw_rho(S_source, CT_source, p_mid);
        rho_ambient_mid = gsw_rho(S_ambient,CT_ambient,p_mid);
        delrho_momentum(1) = rho_dense_mid - rho_ambient_mid;

        delta_CT = 0;
        delta_S = 0;
        distance(1) = 0;
        distance(2:n_time_steps) = NaN(n_time_steps-1,1);
        for t = 1:n_time_steps
            if delrho_plus(t) <0 
                u(t) = 0;
                v(t) = 0;
                H(t) = 0;
                x(t) = 0;
                y(t) = 0;
                rho_dense(t) = 0;
                deltheta_plus(t) = 0;
                continue
            end


            %Step 1: time step the differential equation for new position, water depth, and
            %entrainment based on initial conditions, then t=t1 (the loop will solve
            %for theta, S, density, velocity, H and W for t=t + 1)
            delta_x_m = u(t)*dt;
            delta_y_m = v(t)*dt;
            delta_x_lon = delta_x_m / (111e3 * cosd(y(t)));
            delta_y_lat = delta_y_m / 111e3;

            distance(CTs,Ss,t+1) =distance(CTs,Ss,t) + (abs(delta_x_m) + abs(delta_y_m))/1000;
            

            x(t+1) = delta_x_lon + x(t);
            y(t+1) = delta_y_lat + y(t);

            D(t+1) = interp2(bathy.lon,bathy.lat,D_map_smooth',x(t+1),y(t+1));

            %Solve for entrainment, E
            Ri(t) = (g*delrho_plus(t)*H(t))/(rho_not*(u(t)^2 + v(t)^2)); %layer Richardson number

            if Ri(t) <= 0.8
                E(t) = sqrt(u(t)^2 + v(t)^2)*(0.08-(0.1*Ri(t)))/(1+5*Ri(t));
            else
                E(t)=0;
            end

            %all previous solved using the time step t; before the end of the loop

            %Step 2: Solve for velocity in time step t+1

            %(2a) calculate velocity tendency due to bottom friction
            tau_b_x =  rho_not*Cd*u(t)*sqrt(u(t)^2 + v(t)^2); %zonal bottom friction
            tau_b_y = rho_not*Cd*v(t)*sqrt(u(t)^2 + v(t)^2);  %meridional bottom friction
            Bottom_stress_x = tau_b_x/(rho_not*H(t));
            Bottom_stress_y = tau_b_y/(rho_not*H(t));

            %(2b) calculate velocity tendency due to pressure term

            % (2b.1) find local slope of bathymetry
            Dxxx = interp2(bathy.lon,bathy.lat,D_map_smooth',[x(t) - 0.01 x(t) x(t)+0.01],[y(t) y(t) y(t)]);
            Dyyy = interp2(bathy.lon,bathy.lat,D_map_smooth',[x(t) x(t) x(t)],[y(t) - 0.01 y(t) y(t)+0.01]);
            dDdx = (Dxxx(3) - Dxxx(1))/(sw_dist([y(t), y(t)],[x(t) - 0.01, x(t)+0.01],'km')*1e3);
            dDdy = (Dyyy(3) - Dyyy(1))/(sw_dist([y(t)-0.01, y(t)+0.01],[x(t), x(t)],'km')*1e3);

            % (2b.2) use local bathymetry slope and the density difference from
            % t=1 to solve for pressure tendency term
            Pressure_tendency_x = -(g*delrho_momentum(t) * dDdx)/rho_not;
            Pressure_tendency_y = -(g*delrho_momentum(t) * dDdy)/rho_not;

            %(2c) Calculate velocity tendency due to entrainment term
            Entrainment_x = (E(t) * u(t)) / H(t);
            Entrainment_y = (E(t) * v(t)) / H(t);

            %(2d) Put it all together to solve for u, v at t=t+1
            RHS_x =  Pressure_tendency_x - Bottom_stress_x - Entrainment_x; %tendency terms for zonal velocity
            RHS_y = Pressure_tendency_y - Bottom_stress_y - Entrainment_y;

            %Forward difference
            %    u(t+1) = u(t) + (RHS_x + f*v(t))*dt;  %linear, zonal momentum equation with Coriolis term giving the velocity
            %    v(t+1) = v(t) + (RHS_y - f*u(t))*dt;  %linear, meridional momentum equation with Coriolis term giving the velocity

            %Backward (implicit) difference, help from chatgpt
            A = [1, -dt*f; dt*f, 1];
            b = [u(t) + dt*RHS_x; v(t) + dt*RHS_y];
            uv_new = A \ b;
            u(t+1) = uv_new(1);
            v(t+1) = uv_new(2);

            %Step 3:  Solve for theta, S, rho, H and W for overflow at t=t+1

            % (3a) Theta and S from conservation equations
            %Forward Euler - can create dense water lighter than ambient because of
            %"jumps"
            %CT_dense(t+1) = ((-E(t) * (CT_dense(t) - CT_ambient))/H(t))*dt + CT_dense(t);
            %S_dense(t+1) = ((-E(t) * (S_dense(t) - S_ambient))/H(t))*dt + S_dense(t);
            %Replace with implicit (or backwards) Eueler
            CT_dense(t+1) = (CT_dense(t) + (E(t)*dt/H(t))*CT_ambient) / (1 + E(t)*dt/H(t));
            S_dense(t+1)  = (S_dense(t)  + (E(t)*dt/H(t))*S_ambient)  / (1 + E(t)*dt/H(t));

            %(3b) W from parameterization and H from continuity

            %calculate W at t+1 from parameterization, depends only on thickness and
            %bottom friction (equations 7 and 8 from P&B1994)

            %coded it as in the commented script , which seemed to work,
            %but I think it misread equation 7 which has the vector U and
            %the scalar U (for speed
            % dWdx = 2*tau_b_x/(u(t)*rho_not*H(t)*f);
            % dWdy = 2*tau_b_y/(v(t)*rho_not*H(t)*f);
           speed = sqrt(u(t)^2 + v(t)^2);
           dWdx = 2*tau_b_x*speed/(f*u(t)^2*rho_not*H(t));
           dWdy = 2*tau_b_y*speed/(f*v(t)^2*rho_not*H(t));

            dWdt = u(t)*dWdx + v(t)*dWdy;
            W(t+1) = W(t)+(dWdt*dt);


            %calculate H from continuity, based on equation 6 P&B but in
            %flux form to maintain numerical stability
           
            % A morning of idiotic circular work to calculate the Hdudx
            % term :(
           % dudx = (u(t+1) - u(t))/delta_x_m
           % dvdy = (v(t+1)-v(t))/delta_y_m
           % udHdx =(E(t) - H(t)*u(t)*dWdx/W(t) - H(t)*dudx) %
           % vdHdy = (E(t) - H(t)*v(t)*dWdy/W(t)- H(t)*dvdy)
          DuHdxplusdvHdy = E(t) - H(t)*v(t)*dWdy/W(t) - H(t)*u(t)*dWdx/W(t);
          DHDt = DuHdxplusdvHdy;


            H(t+1) = H(t) + (DHDt*dt);
   


            %update coriolis
            f = gsw_f(y(t+1));

            %finally, have what we need to calculate the evolving density of the
            %overflow, which becomes essential for the next time step
            depth_mid = D(t+1) + H(t+1)/2;
            depth_top = D(t+1) + H(t+1);
            % Convert to pressure (dbar)
            p_mid = gsw_p_from_z(depth_mid, y(t+1));
            p_top = gsw_p_from_z(depth_top, y(t+1));

            rho_dense(t+1) = gsw_rho(S_dense(t+1) ,CT_dense(t+1),p_top); %Need to solve for depth of overflow ...
            rho_ambient_top = gsw_rho(S_ambient, CT_ambient, p_top);
            rho_ambient_mid(t+1) = gsw_rho(S_ambient, CT_ambient, p_mid);
            delrho_plus(t+1) = rho_dense(t) - rho_ambient_top;
            rho_dense_mid(t+1) = gsw_rho(S_dense(t+1), CT_dense(t+1),p_mid);
            delrho_momentum(t+1) = rho_dense_mid(t+1) - rho_ambient_mid(t+1);
            gradD(t) = sqrt(dDdx^2 + dDdy^2);

            %disp([t, H(t), delrho_plus(t), u(t)^2 + v(t)^2, Ri(t), dDdx])
            % --- debug print for a few steps ---------------------------------------
%             if mod(t,10)==0 || t<=5
%                 fprintf('t=%d Ptx=%.4e Btx=%.4e Entx=%.4e RHSx=%.4e  U=%.4f\n', ...
%                     t, Pressure_tendency_x, -Bottom_stress_x, Entrainment_x, RHS_x, sqrt(u(t)^2+v(t)^2));
%             end
        end  %stepping through time

        timeax = linspace(1,round(dt*n_time_steps)/(60*60*24),n_time_steps+1);
        
        %figure; subplot(3,1,1);plot(distance,u); 
        %title(strcat('S-initial',num2str(S_range(Ss)),'CT-initial =', num2str(CT_range(CTs))));
        %hold on;plot(distance,v);legend('u','v');grid on;title('velocity');
        %subplot(3,1,2); plot(distance(1:end-1),gradD);title('Grad(D)');grid on;
        %subplot(3,1,3);plot(distance(1:end-1),E); title('E');grid on;
        %xlabel('time (days)')

        %figure;plot(distance,delrho_momentum);title('delrho_momentum');grid on

        %figure;subplot(2,1,1); plot(distance,W);title('W (km)');subplot(2,1,2);plot(distance,H);title('H (m)');

        U = sqrt(u.^2 + v.^2);
        if Ss == 1
        figure;plot(U.*H.*W/1e6)
        end
        Us(CTs,Ss,1:n_time_steps+1) = U;
        Hs(CTs,Ss,1:n_time_steps+1) = H;
        Ws(CTs,Ss,1:n_time_steps+1) = W;
        p_mids(CTs,Ss,1:n_time_steps+1) = p_mid;
        TTTs(CTs,Ss,1:n_time_steps+1) = S_dense; 
        SSSs(CTs,Ss,1:n_time_steps+1) = CT_dense;
        Ss
    end
    CTs 
end

figure(1)
hold on;plotm(y,x,'k','linewidth',2);
plotm(y(1:50:end),y(1:50:end),'k.','markersize',10)

[CTi Si] = meshgrid(CT_range,SA_range);
SA_range = gsw_SA_from_SP(S_range,650,lon_initial,lat_initial);
T = Us.*Hs.*Ws;
%save DS_transport_dataSAiCTi.mat CTi Si T

clear all

n_cts = 20;
n_sts = 22;
CT_range = linspace(-.1,5.5,n_cts);  %when I run loops I need this 
S_range = linspace(34.6,35.45,n_sts);%when I run loops I need this 
  
n_time_steps = 700;
dt = 60*60*24/50; %time step in seconds (denominator is time steps per day)

 Us= NaN(n_cts,n_sts,n_time_steps);
 Hs= NaN(n_cts,n_sts,n_time_steps);
 Ws= NaN(n_cts,n_sts,n_time_steps);
 p_mid = NaN(n_cts,n_sts,n_time_steps);
 distance= NaN(n_cts,n_sts,n_time_steps);
 TTTs= NaN(n_cts,n_sts,n_time_steps);
 SSSs= NaN(n_cts,n_sts,n_time_steps);

bathy = netcdf_load('gebco_2025_n64.0_s59.0_w-30.0_e-3.0.nc');
D_map = double(bathy.elevation); %take depth as positive number (multiply elevation by -1)

%Smooth bathymetry to stop runaway entrainment
% Simple  moving average (replace with Gaussian if desired)
D_map_smooth = conv2(D_map, ones(30)/30^2, 'same');
        
for CTs =1:n_cts
    for Ss = 1:n_sts

        CTs
        Ss
        %Define and initialize variables
    
        %Ambient water properties
        theta_ambient = 7; %User defined ambient potential temperature [°C] - 7°C taken from P&B Table 1
        S_ambient = 35.15; %User defined ambient practical salinity [g/kg] -35.15 taken from P&B Table 1

        %Source water properties
        theta_source = CT_range(CTs);%0;% 0.5; %User defined source water potential temperature [°C]; 0 from PB94 but 0.3 from Eldevik
        S_source = S_range(Ss);%34.9;%34.87; %User defined source water practical salinity [g/kg]; Sp= 34.92 from PB94, very close in Eldevic

        %constants
        rho_not = 1030; %user-defined reference density, a constant [kg/m3]
        Cd = 4e-3;%3e-3; % Drag coefficient (~50% uncertainty) [unitless]
        g = 9.8; %gravitational constant [m/s2]

        %site parameters and initial conditions
        H = NaN(n_time_steps,1);
        W = NaN(n_time_steps,1);
        H(1) = 150; %initial thickness of overflow, 200 from Hiroki; 150 m from Table 1 in P&B94 [m]
        W(1) = 15e3;%50e3; %initial width of overflow, 80km from Hiroki; 50 km from Table 1 in P&B94 [m]

        lon_initial = -12.4; %define grid cell in bathymetry? 65.93°N, 27.90°W - Hiroki
        lat_initial = 61.73;
        initial_direction = -200; %an angle representing the initial dirction of flow, in P&B - northwestward
        initial_transport = 1.7e6; % [m3/s]
        initial_speed = 0.08;%P&B describe 0.8m/s for initial speed (though in the Den Strait code, I had to really start slow)

        %pre-allocate arrays, so they don't grow in loop
        u = NaN(n_time_steps,1);
        v = NaN(n_time_steps,1);
        x = NaN(n_time_steps,1);
        y = NaN(n_time_steps,1);
        u(1) = initial_speed*cosd(initial_direction);
        v(1) = initial_speed*sind(initial_direction);
        x(1) =  lon_initial;
        y(1) =  lat_initial;

        f = gsw_f(y(1));
       

        D = NaN(n_time_steps,1);
        D(1) = interp2(bathy.lon,bathy.lat,D_map_smooth',x(1),y(1));

        %calculate in situ density from user-defined theta and S for ambient and
        %source waters:

        % Mid-depth (in meters)
        depth_mid = D(1) + H(1)/2;
        % Convert to pressure (dbar)
        p_mid = gsw_p_from_z(depth_mid, y(1));
        depth_top = D(1)+H(1);
        p_top=gsw_p_from_z(depth_top, y(1));

        CT_ambient = gsw_CT_from_pt(S_ambient,theta_ambient);
        S_ambient = gsw_SA_from_SP(S_ambient,p_top,lon_initial,lat_initial);

        rho_ambient = gsw_rho(S_ambient, CT_ambient, p_top);
        S_source = gsw_SA_from_SP(S_source,p_mid,lon_initial,lat_initial);
        CT_source = gsw_CT_from_pt(S_source,theta_source);

        rho_source = gsw_rho(S_source, CT_source,p_top);
        rho_dense(1) = rho_source;
        rho_dense(2:n_time_steps) = NaN;

        %initial conditions of the overflow are set to the source water properties
        CT_dense = NaN(n_time_steps,1);
        S_dense = NaN(n_time_steps,1);

        CT_dense(1) = CT_source;
        S_dense(1) = S_source;
        E = NaN(n_time_steps,1);

        %The delrho term is the difference between the in situ density between the
        %outflow water and the ambient oceanic water.  This is defined in 2 ways,
        %but the one relevant for entrainment, delrho+, is taken as the difference between
        %the water at the top of the outflow (D-H, i.e. ocean depth minus the
        %thickness of the outflow).  The other delrho (no plus sign) is the density
        %difference between the outflow and oceanic water at mid-depth of the
        %outflow. In the first time through running the code, I'm setting a single constant value for
        %the ambient seawater, rather than defining an entire profile.  This is slightly  too
        %simplistic, given a bit of stratification at the depth range of
        %entrainment (see Hiroki's Figure 4), but a good first step...

        %in the first time step, the delta is between source water and ambient
        %water.  In the loop, the source water will be replaced by the dense water
        %characteristics, which change vian entrainment/mixing

        delrho_plus = NaN(n_time_steps,1);
        delrho_plus(1) = rho_source - rho_ambient;

        deltheta_plus = NaN(n_time_steps,1);
        deltheta_plus(1) = CT_source - CT_ambient;
        delS_plus = NaN(n_time_steps,1);
        delS_plus(1) = S_source - S_ambient;

        %Here define the delrho used in the momentum equation (these are different
        %only due to the different pressures used in the in situ density
        %calculation)
        delrho_momentum = NaN(n_time_steps,1);

        rho_dense_mid = gsw_rho(S_source, CT_source, p_mid);
        rho_ambient_mid = gsw_rho(S_ambient,CT_ambient,p_mid);
        delrho_momentum(1) = rho_dense_mid - rho_ambient_mid;

        delta_CT = 0;
        delta_S = 0;
        distance(1) = 0;
        distance(2:n_time_steps) = NaN(n_time_steps-1,1);
        for t = 1:n_time_steps
            if delrho_plus(t) <0 
                u(t) = 0;
                v(t) = 0;
                H(t) = 0;
                x(t) = 0;
                y(t) = 0;
                rho_dense(t) = 0;
                deltheta_plus(t) = 0;
                continue
            end


            %Step 1: time step the differential equation for new position, water depth, and
            %entrainment based on initial conditions, then t=t1 (the loop will solve
            %for theta, S, density, velocity, H and W for t=t + 1)
            delta_x_m = u(t)*dt;
            delta_y_m = v(t)*dt;
            delta_x_lon = delta_x_m / (111e3 * cosd(y(t)));
            delta_y_lat = delta_y_m / 111e3;

            distance(CTs,Ss,t+1) =distance(CTs,Ss,t) + (abs(delta_x_m) + abs(delta_y_m))/1000;
            

            x(t+1) = delta_x_lon + x(t);
            y(t+1) = delta_y_lat + y(t);

            D(t+1) = interp2(bathy.lon,bathy.lat,D_map_smooth',x(t+1),y(t+1));

            %Solve for entrainment, E
            Ri(t) = (g*delrho_plus(t)*H(t))/(rho_not*(u(t)^2 + v(t)^2)); %layer Richardson number

            if Ri(t) <= 0.8
                E(t) = sqrt(u(t)^2 + v(t)^2)*(0.08-(0.1*Ri(t)))/(1+5*Ri(t));
            else
                E(t)=0;
            end

            %all previous solved using the time step t; before the end of the loop

            %Step 2: Solve for velocity in time step t+1

            %(2a) calculate velocity tendency due to bottom friction
            tau_b_x =  rho_not*Cd*u(t)*sqrt(u(t)^2 + v(t)^2); %zonal bottom friction
            tau_b_y = rho_not*Cd*v(t)*sqrt(u(t)^2 + v(t)^2);  %meridional bottom friction
            Bottom_stress_x = tau_b_x/(rho_not*H(t));
            Bottom_stress_y = tau_b_y/(rho_not*H(t));

            %(2b) calculate velocity tendency due to pressure term

            % (2b.1) find local slope of bathymetry
            Dxxx = interp2(bathy.lon,bathy.lat,D_map_smooth',[x(t) - 0.01 x(t) x(t)+0.01],[y(t) y(t) y(t)]);
            Dyyy = interp2(bathy.lon,bathy.lat,D_map_smooth',[x(t) x(t) x(t)],[y(t) - 0.01 y(t) y(t)+0.01]);
            dDdx = (Dxxx(3) - Dxxx(1))/(sw_dist([y(t), y(t)],[x(t) - 0.01, x(t)+0.01],'km')*1e3);
            dDdy = (Dyyy(3) - Dyyy(1))/(sw_dist([y(t)-0.01, y(t)+0.01],[x(t), x(t)],'km')*1e3);

            %    dDdx(dDdx > 0.012) = 0.012;
            %    dDdy(dDdy > 0.012) = 0.012;
            %    dDdx(dDdx<-0.012) = -0.012;
            %    dDdy(dDdy<-0.012) = -0.012;

            % (2b.2) use local bathymetry slope and the density difference from
            % t=1 to solve for pressure tendency term
            Pressure_tendency_x = -(g*delrho_momentum(t) * dDdx)/rho_not;
            Pressure_tendency_y = -(g*delrho_momentum(t) * dDdy)/rho_not;

            %(2c) Calculate velocity tendency due to entrainment term
            Entrainment_x = (E(t) * u(t)) / H(t);
            Entrainment_y = (E(t) * v(t)) / H(t);

            %(2d) Put it all together to solve for u, v at t=t+1
            RHS_x =  Pressure_tendency_x - Bottom_stress_x - Entrainment_x; %tendency terms for zonal velocity
            RHS_y = Pressure_tendency_y - Bottom_stress_y - Entrainment_y;

            %Forward difference
            %    u(t+1) = u(t) + (RHS_x + f*v(t))*dt;  %linear, zonal momentum equation with Coriolis term giving the velocity
            %    v(t+1) = v(t) + (RHS_y - f*u(t))*dt;  %linear, meridional momentum equation with Coriolis term giving the velocity

            %Backward (implicit) difference, help from chatgpt
            A = [1, -dt*f; dt*f, 1];
            b = [u(t) + dt*RHS_x; v(t) + dt*RHS_y];
            uv_new = A \ b;
            u(t+1) = uv_new(1);
            v(t+1) = uv_new(2);

            %Step 3:  Solve for theta, S, rho, H and W for overflow at t=t+1

            % (3a) Theta and S from conservation equations
            %Forward Euler - can create dense water lighter than ambient because of
            %"jumps"
            %CT_dense(t+1) = ((-E(t) * (CT_dense(t) - CT_ambient))/H(t))*dt + CT_dense(t);
            %S_dense(t+1) = ((-E(t) * (S_dense(t) - S_ambient))/H(t))*dt + S_dense(t);
            %Replace with implicit (or backwards) Eueler
            CT_dense(t+1) = (CT_dense(t) + (E(t)*dt/H(t))*CT_ambient) / (1 + E(t)*dt/H(t));
            S_dense(t+1)  = (S_dense(t)  + (E(t)*dt/H(t))*S_ambient)  / (1 + E(t)*dt/H(t));

            %(3b) W from parameterization and H from continuity

            %calculate W at t+1 from parameterization, depends only on thickness and
            %bottom friction (equations 7 and 8 from P&B1994)

            %coded it as in the commented script , which seemed to work,
            %but I think it misread equation 7 which has the vector U and
            %the scalar U (for speed

            speed = sqrt(u(t)^2 + v(t)^2);
            DWDt = (2*Cd/H(t)*f) * (u(t) + v(t));
            W(t+1) = W(t) + (DWDt*dt);

            DHDt  = E(t)- (H(t)*DWDt/W(t));

            H(t+1) = H(t) + (DHDt*dt);
           
   

            %update coriolis
            f = gsw_f(y(t+1));

            %finally, have what we need to calculate the evolving density
            %of thef
            %overflow, which becomes essential for the next time step
            depth_mid = D(t+1) + H(t+1)/2;
            depth_top = D(t+1) + H(t+1);
            % Convert to pressure (dbar)
            p_mid = gsw_p_from_z(depth_mid, y(t+1));
            p_top = gsw_p_from_z(depth_top, y(t+1));

            rho_dense(t+1) = gsw_rho(S_dense(t+1) ,CT_dense(t+1),p_top); %Need to solve for depth of overflow ...
            rho_ambient_top = gsw_rho(S_ambient, CT_ambient, p_top);
            rho_ambient_mid(t+1) = gsw_rho(S_ambient, CT_ambient, p_mid);
            delrho_plus(t+1) = rho_dense(t) - rho_ambient_top;
            rho_dense_mid(t+1) = gsw_rho(S_dense(t+1), CT_dense(t+1),p_mid);
            delrho_momentum(t+1) = rho_dense_mid(t+1) - rho_ambient_mid(t+1);
            gradD(t) = sqrt(dDdx^2 + dDdy^2);

            %disp([t, H(t), delrho_plus(t), u(t)^2 + v(t)^2, Ri(t), dDdx])
            % --- debug print for a few steps ---------------------------------------
%             if mod(t,10)==0 || t<=5
%                 fprintf('t=%d Ptx=%.4e Btx=%.4e Entx=%.4e RHSx=%.4e  U=%.4f\n', ...
%                     t, Pressure_tendency_x, -Bottom_stress_x, Entrainment_x, RHS_x, sqrt(u(t)^2+v(t)^2));
%             end
        end  %stepping through time

        timeax = linspace(1,round(dt*n_time_steps)/(60*60*24),n_time_steps+1);
        
        %figure; subplot(3,1,1);plot(distance,u); 
        %title(strcat('S-initial',num2str(S_range(Ss)),'CT-initial =', num2str(CT_range(CTs))));
        %hold on;plot(distance,v);legend('u','v');grid on;title('velocity');
        %subplot(3,1,2); plot(distance(1:end-1),gradD);title('Grad(D)');grid on;
        %subplot(3,1,3);plot(distance(1:end-1),E); title('E');grid on;
        %xlabel('time (days)')

        %figure;plot(distance,delrho_momentum);title('delrho_momentum');grid on

        %figure;subplot(2,1,1); plot(distance,W);title('W (km)');subplot(2,1,2);plot(distance,H);title('H (m)');

        U = sqrt(u.^2 + v.^2);
        if Ss == 1
        figure;plot(U.*H.*W/1e6)
        end
        Us(CTs,Ss,1:n_time_steps+1) = U;
        Hs(CTs,Ss,1:n_time_steps+1) = H;
        Ws(CTs,Ss,1:n_time_steps+1) = W;
        p_mids(CTs,Ss,1:n_time_steps+1) = p_mid;
        TTTs(CTs,Ss,1:n_time_steps+1) = S_dense; 
        SSSs(CTs,Ss,1:n_time_steps+1) = CT_dense;
        Ss
    end
    CTs 
end

hold on;plotm(y,x,'r','linewidth',2);
%plotm(y(1:50:end),x(1:50:end),'r+','markersize',20)
% Center at 67.5°N, 26°W with 1° x 1° box
center_lat = 67.5;
center_lon = -26;

% Box boundaries (0.5° in each direction)
lat_min = center_lat - 0.5;
lat_max = center_lat + 0.5;
lon_min = center_lon - 0.5;
lon_max = center_lon + 0.5;
plotm([lat_min lat_max], [lon_min lon_min], 'k-', 'LineWidth', 2)  % 
plotm([lat_max lat_max], [lon_min lon_max], 'k-', 'LineWidth', 2)  % 
plotm([lat_max lat_min], [lon_max lon_max], 'k-', 'LineWidth', 2)  % 
plotm([lat_min lat_min], [lon_max lon_min], 'k-', 'LineWidth', 2)  %  

center_lat = 64;
center_lon = -8;

% Box boundaries (0.5° in each direction)
lat_min = center_lat - 0.5;
lat_max = center_lat + 0.5;
lon_min = center_lon - 0.5;
lon_max = center_lon + 0.5;
plotm([lat_min lat_max], [lon_min lon_min], 'r-', 'LineWidth', 2)  % 
plotm([lat_max lat_max], [lon_min lon_max], 'r-', 'LineWidth', 2)  % 
plotm([lat_max lat_min], [lon_max lon_max], 'r-', 'LineWidth', 2)  % 
plotm([lat_min lat_min], [lon_max lon_min], 'r-', 'LineWidth', 2)  %  