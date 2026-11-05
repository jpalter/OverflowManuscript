% Price and Baringer (1994) Streamtube model
% Now that I got the model working, loop over the whole thing, for
% different combinations of CT_source and S_source

figure
n_cts = 22;
n_sts = 22;
CT_range = linspace(-1,5.5,n_cts);  %when I run loops I need this
S_range = linspace(34.5,35.25,n_sts);%when I run loops I need this
n_time_steps = 2000; % 40 days with 50 time steps per day (as dt defined below)
dt = 60*60*24/50; %time step in seconds (denominator is time steps per day, about a half hour)

%initalize arrays that will be populated in loop
Us= NaN(n_cts,n_sts,n_time_steps);
Hs= NaN(n_cts,n_sts,n_time_steps);
Ws= NaN(n_cts,n_sts,n_time_steps);
p_mid = NaN(n_cts,n_sts,n_time_steps);
distance= NaN(n_cts,n_sts,n_time_steps);
TTTs= NaN(n_cts,n_sts,n_time_steps);
SSSs= NaN(n_cts,n_sts,n_time_steps);
lonss= NaN(n_cts,n_sts,n_time_steps);
latss= NaN(n_cts,n_sts,n_time_steps);


%load bathymetry for the region from GEBCO2025; 
bathy = netcdf_load('gebco_2025_n64.0_s59.0_w-30.0_e-3.0.nc');
D_map = double(bathy.elevation); %take depth as positive number (multiply elevation by -1)

%Smooth bathymetry to stop runaway entrainment at very rough spots in
%native resolution)
D_map_smooth = conv2(D_map, ones(30)/30^2, 'same'); %30x30 boxcar.


for CTs =1:length(CT_range)
    for Ss =1:length(S_range)

        %Define and initialize variables
        %Ambient water properties
        theta_ambient = 7; %User defined ambient potential temperature [°C] - taken from P&B Table 1
        S_ambient = 35.15; %User defined ambient absolute salinity [g/kg] -taken from P&B Table 1

        %Source water properties
        %Looping over hypothetical range of source water temp and salinity
        theta_source =CT_range(CTs); % [degC]
        S_source = S_range(Ss); %[g/kg]

        %constants
        rho_not = 1030; %user-defined reference density, a constant [kg/m3]
        Cd = 4e-3;% % Drag coefficient [unitless]
        g = 9.8; %gravitational constant [m/s2]

        %site parameters and initial conditions
        H = NaN(n_time_steps,1); %height of streamtube 
        W = NaN(n_time_steps,1); %width of streamtube 
        H(1) = 150; %initial thickness of overflow, 150 m from Table 1 in P&B94 [m]
        W(1) = 15e3;%initial width of overflow [m]

        %define initial grid cell in bathymetry - must be in FBC
        lon_initial = -12.4;
        lat_initial = 61.73;
        initial_direction = -200; %an angle representing the initial dirction of flow -  northwestward
        initial_speed = 0.08;%[m/s]

        %pre-allocate arrays, so they don't grow in loop
        u = NaN(n_time_steps,1);
        v = NaN(n_time_steps,1);
        x = NaN(n_time_steps,1);
        y = NaN(n_time_steps,1);
        u(1) = initial_speed*cosd(initial_direction);
        v(1) = initial_speed*sind(initial_direction);
        x(1) =  lon_initial;
        y(1) =  lat_initial;

        f = gsw_f(y(1)); %Coriolis parameter [1/s]

        D = NaN(n_time_steps,1);
        D(1) = interp2(bathy.lon,bathy.lat,D_map_smooth',x(1),y(1)); %Water depth from bathymetry [m]

        %calculate in situ density from user-defined theta and S for
        %ambient and source waters:
        % Mid-depth (in meters)
        depth_mid = D(1) + H(1)/2;
        % Convert to pressure (dbar)
        p_mid = gsw_p_from_z(depth_mid, y(1));
        depth_top = D(1)+H(1);
        p_top=gsw_p_from_z(depth_top, y(1));

        CT_ambient = gsw_CT_from_pt(S_ambient,theta_ambient);
        rho_ambient = gsw_rho(S_ambient, CT_ambient, p_top);
        CT_source = gsw_CT_from_pt(S_source,theta_source);

        rho_source = gsw_rho(S_range(Ss), CT_range(CTs),p_top);
        rho_dense(1) = rho_source;
        rho_dense(2:n_time_steps) = NaN;

        %initial conditions of the overflow are set to the source water properties
        CT_dense = NaN(n_time_steps,1);
        S_dense = NaN(n_time_steps,1);

        CT_dense(1) = CT_source;
        S_dense(1) = S_source;
        E = NaN(n_time_steps,1);

        %The delrho term is the difference between the in situ density
        %between the outflow water and the ambient oceanic water.  This is
        %defined in 2 ways, but the one relevant for entrainment, delrho+,
        %is taken as the difference between the water at the top of the
        %outflow (D-H, i.e. ocean depth minus the thickness of the
        %outflow).  The other delrho (no plus sign) is the density
        %difference between the outflow and oceanic water at mid-depth of
        %the outflow. I'm setting a single constant value for the ambient
        %seawater, rather than defining an entire profile to match
        %conditions along the depth of the streamtube. Given the weak
        %stratification, I believe this is a reasonable simplification to
        %explore the system.

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
            if delrho_plus(t) <0  %assume no outflow if there is no density gradient at the sill
                u(t) = 0;
                v(t) = 0;
                H(t) = 0;
                x(t) = 0;
                y(t) = 0;
                rho_dense(t) = 0;
                deltheta_plus(t) = 0;
                continue
            end


            %Step 1: time step the differential equation for new position,
            %water depth, and entrainment based on initial conditions, then
            % the loop will solve for theta, S, density, velocity, H
            %and W for t=t + 1)
            delta_x_m = u(t)*dt;
            delta_y_m = v(t)*dt;
            delta_x_lon = delta_x_m / (111e3 * cosd(y(t)));
            delta_y_lat = delta_y_m / 111e3;

            distance(CTs,Ss,t+1) =distance(t)+ (abs(delta_x_m) + abs(delta_y_m))/1000;

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

            % (2b.2) use local bathymetry slope and the density difference
            % from t=1 to solve for pressure tendency term
            Pressure_tendency_x = -(g*delrho_momentum(t) * dDdx)/rho_not;
            Pressure_tendency_y = -(g*delrho_momentum(t) * dDdy)/rho_not;

            %(2c) Calculate velocity tendency due to entrainment term
            Entrainment_x = (E(t) * u(t)) / H(t);
            Entrainment_y = (E(t) * v(t)) / H(t);

            %(2d) Put it all together to solve for total tendency for zonal
            %velocity
            RHS_x =  Pressure_tendency_x - Bottom_stress_x - Entrainment_x; 
            RHS_y = Pressure_tendency_y - Bottom_stress_y - Entrainment_y;

            %Forward difference - very unstable (therefore commented)
            %    u(t+1) = u(t) + (RHS_x + f*v(t))*dt;  %linear, zonal momentum equation with Coriolis term giving the velocity
            %    v(t+1) = v(t) + (RHS_y - f*u(t))*dt;  %linear, meridional momentum equation with Coriolis term giving the velocity

            %Backward (implicit) difference
            A = [1, -dt*f; dt*f, 1];
            b = [u(t) + dt*RHS_x; v(t) + dt*RHS_y];
            uv_new = A \ b;
            u(t+1) = uv_new(1);
            v(t+1) = uv_new(2);

            %Step 3:  Solve for theta, S, rho, H and W for overflow at t=t+1

            % (3a) Theta and S from conservation equations
            %Forward Euler - can create dense water lighter than ambient
            %because of "jumps"
            %CT_dense(t+1) = ((-E(t) * (CT_dense(t) - CT_ambient))/H(t))*dt + CT_dense(t);
            %S_dense(t+1) = ((-E(t) * (S_dense(t) - S_ambient))/H(t))*dt + S_dense(t);
            
            %Replace with implicit (or backwards) Euler
            CT_dense(t+1) = (CT_dense(t) + (E(t)*dt/H(t))*CT_ambient) / (1 + E(t)*dt/H(t));
            S_dense(t+1)  = (S_dense(t)  + (E(t)*dt/H(t))*S_ambient)  / (1 + E(t)*dt/H(t));

            %(3b) W from parameterization and H from continuity

            %calculate W at t+1 from parameterization, depends only on thickness and
            %bottom friction (equations 7 and 8 from P&B1994)

            dWdx = 2*tau_b_x*speed/(f*u(t)^2*rho_not*H(t));
            dWdy = 2*tau_b_y*speed/(f*v(t)^2*rho_not*H(t));

            dWdt = u(t)*dWdx + v(t)*dWdy;
            W(t+1) = W(t)+(dWdt*dt);


            %calculate H from continuity, based on equation 6 P&B but in
            %flux form to maintain numerical stability

            DuHdxplusdvHdy = E(t) - H(t)*v(t)*dWdy/W(t) - H(t)*u(t)*dWdx/W(t);
            DHDt = DuHdxplusdvHdy;

            H(t+1) = H(t) + (DHDt*dt);

            %update coriolis
            f = gsw_f(y(t+1));

            %finally, have what we need to calculate the evolving density
            %of the overflow, which becomes essential for the nexttime step
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
            figure;subplot(3,1,1);plot(U);subplot(3,1,2);plot(H),subplot(3,1,3);plot(W); title({num2str(CTs) 'k',num2str(Ss)});

            %figure;plot(U.*H.*W/1e6)
        end
        Us(CTs,Ss,1:n_time_steps+1) = U;
        Hs(CTs,Ss,1:n_time_steps+1) = H;
        Ws(CTs,Ss,1:n_time_steps+1) = W;
        p_mids(CTs,Ss,1:n_time_steps+1) = p_mid;
        TTTs(CTs,Ss,1:n_time_steps+1) = S_dense;
        SSSs(CTs,Ss,1:n_time_steps+1) = CT_dense;
        lonss(CTs,Ss,1:n_time_steps+1) = x;
        latss(CTs,Ss,1:n_time_steps+1) = y;
        %figure;plot(U.*H.*W/1e6)
        %figure;subplot(3,1,1);plot(U);subplot(3,1,2);plot(H),subplot(3,1,3);plot(W); title({num2str(CTs) 'k',num2str(Ss)});

        Ss
    end
    CTs
end

[Si CTi] = meshgrid(S_range,CT_range);
T = Us.*Hs.*Ws;
figure;pcolor(Si,CTi,mean(T(:,:,end-300:end-50),3)/1e6);colorbar
