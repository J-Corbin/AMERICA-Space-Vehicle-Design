function [Data_per_sol] = CommsAccess(et)
%% FUNCTION DESCRIPTION:
% This function is capable of taking a ephemeris time period, then determining 
% line of sight intervals between a mars lander (with given coordinates) and any 
% mars orbiter. 
% 
% INPUTS:
% et = an nx1 array of ephemeris times (must be within range of spice data period)
% OUTPUTS:
% Data_per_sol = total amount of data, in bits, per sol that can be transferred
% from the lander to all (or any subset) of the orbiters.
% 
% Developed by William Parker, November 2019

%% Load SPICE data
DIR1 = ['C:\Users\Will\Documents\MATLAB\Moon Project\spice kernels\'];

%Fixed Kernels
cspice_furnsh([DIR1 'naif0012.tls']);                      % Generic LSK       (Leapsecond Kernel)
%MRO
cspice_furnsh([DIR1 'mro_v15.tf']);   % Mission FK 27 Nov 2012       (Frames Kernel)
cspice_furnsh([DIR1 'mro_sclkscet_00050_65536.tsc']);      % Spacecraft SCLK 27 Nov 2012   (Spacecraft Clock Kernel)
%Odyssey
cspice_furnsh([DIR1 'm01_v28.tf']);                        % Mission FK 17 March 2006      (Frames Kernel)
cspice_furnsh([DIR1 'orb1_sclkscet_00193.tsc']);      % Spacecraft SCLK 27 Nov 2012   (Spacecraft Clock Kernel)
%MAVEN
cspice_furnsh([DIR1 'maven_v09.tf']);                        % Mission FK 27 Nov 2012       (Frames Kernel)
cspice_furnsh([DIR1 'mvn_sclkscet_00049.tsc']);      % Spacecraft SCLK 27 Nov 2012   (Spacecraft Clock Kernel)

% cspice_furnsh([DIR1 'cas_iss_v10.ti']);                    % Camera IK         (ISS Instrument Kernel)
cspice_furnsh([DIR1,'pck00010.tpc']);  %PCK planetary characteristics kernel (RADIUS)


olddir = cd(DIR1);
files_bsp = dir('*.bsp');
cd(olddir)
for ii = 1:length(files_bsp)
    cspice_furnsh([DIR1 files_bsp(ii).name]);
end



% dates with viable SPICE data based on loaded kernels:
% MRO: 
% START_TIME                   = 2017-07-01T00:00:01
% STOP_TIME                    = 2017-10-01T01:00:00
% 
% Odyssey:
% START_TIME                   = 2017-07-01T00:00:00
% STOP_TIME                    = 2017-10-01T01:00:00
% 
% MAVEN: 
% START_TIME                   = 2017-07-01T00:00:00
% STOP_TIME                    = 2017-10-01T01:00:00


%Assemble current epoch string and covert to et
et = et:60:et+3600*24;

%determine position of lander in IAU-Mars frame
r = 3389.5;
lander_deg = [-137.8, -5.4];
temp = r.*cosd(5.4);
lander = [-temp*sind(47.8),-temp*cosd(47.8), -r*sind(5.4)].';
k = norm(lander);

%get position of satellites relative to Mars center at each ephemeris time
 [MRO_pos,~] = cspice_spkpos('MRO', et,'IAU_MARS', 'NONE', 'MARS');
 [Odyssey_pos,~] = cspice_spkpos('Mars Odyssey', et,'IAU_MARS', 'NONE', 'MARS');
 [Maven_pos,~] = cspice_spkpos('MAVEN', et,'IAU_MARS', 'NONE', 'MARS');
   
 
%%Calculate position of TGO from orbit parameters (No SPICE data available)
%  mu = 42828;
%  a = 3785.7;
%  e = 0.0071;
%  h = sqrt(a*mu*(1-e^2));
%  incl = 73.5;
%  RAAN = 180.6;
%  ArgPer = 266.8;
%  Anom = 63.6;
 
%  TGO_pos = ones(length(et),3);
%  for i = 1:length(et)
%  TGO_pos = MyOrbitWithAnomaly(et,mu, a , e , h, incl, RAAN, ArgPer, Anom); 
%  end
 
% Get position of orbiters in J2000 for plotting purposes
%    [MRO_pos,~] = cspice_spkpos('MRO', et,'J2000', 'NONE', 'MARS');
%   [Odyssey_pos,~] = cspice_spkpos('Mars Odyssey', et,'J2000', 'NONE', 'MARS');
%    [Maven_pos,~] = cspice_spkpos('MAVEN', et,'J2000', 'NONE', 'MARS');
%    
%    
% Get line of sight from Mars to Earth in J2000
%    [Earth_pos,~] = cspice_spkpos('EARTH', et,'J2000', 'NONE', 'MARS');
%       
%       U_MRO = MRO_pos/norm(MRO_pos);
%       U_earth = Earth_pos/norm(Earth_pos);
%       if U_earth*U_MRO' > 0
%           connect_val = 500e3;
%       else 
%           connect_val = 0;
%       end
%       
%       plot(et, connect_val);



%    plot3(MRO_pos(1,:), MRO_pos(2,:), MRO_pos(3,:),'r.');
%    hold on
%    scatter3(Odyssey_pos(1,:), Odyssey_pos(2,:),Odyssey_pos(3,:),'g.');
%    scatter3(Maven_pos(1,:), Maven_pos(2,:), Maven_pos(3,:), 'b.');
%       scatter3(TGO_pos(1,:), TGO_pos(2,:),TGO_pos(3,:),'m.');
% plot3(lander(1),lander(2),lander(3), '*c','MarkerSize',10);
%    legend('MRO','Mars Odyssey','MAVEN','TGO','Landing Site');
%    [sx,sy,sz] = sphere();
%    r = 3389.5;
%    hs1 = surf( r*sx, r*sy, r*sz);
%    set(hs1, 'FaceColor', [0.8500, 0.3250, 0.0980]	);
%    xlabel('[km]');
%    ylabel('[km]');
%    zlabel('[km]');
%    title('Communications Orbits');
%    axis equal
   
%    hold off
   
%initialize variables   
   MRO_connect = zeros(3,length(et));
       Odyssey_connect = zeros(3,length(et));
       Maven_connect = zeros(3,length(et));
%         TGO_connect = zeros(3,length(et));

%Determine unit vector from lander to each orbiter at each ephemeris time
   for i = 1:length(et)
       MRO_connect(:,i) = MRO_pos(:,i)- lander;
       dr_MRO(i) = sqrt(MRO_connect(1,i)^2+MRO_connect(2,i)^2+MRO_connect(3,i)^2);
       u_MRO(:,i) = MRO_connect(:,i)./dr_MRO(i);
%        D_MRO = det([lander(1),MRO_pos(1,i);lander(2), MRO_pos(2,i); lander(3), MRO_pos(3)]);
       Odyssey_connect(:,i) = Odyssey_pos(:,i) -lander;
       dr_Od(i) = sqrt(Odyssey_connect(1,i)^2+Odyssey_connect(2,i)^2+Odyssey_connect(3,i)^2);
      u_Odyssey(:,i) = Odyssey_connect(:,i)./dr_Od(i);
       Maven_connect(:,i) = Maven_pos(:,i) - lander;
       dr_Maven(i) = sqrt(Maven_connect(1,i)^2+Maven_connect(2,i)^2+Maven_connect(3,i)^2);
       u_Maven(:,i) = Maven_connect(:,i)./dr_Maven(i);
%         TGO_connect(:,i) = TGO_pos(:,i) - lander;
%        dr_TGO(i) = sqrt(TGO_connect(1,i)^2+TGO_connect(2,i)^2+TGO_connect(3,i)^2);
%        u_TGO(:,i) = TGO_connect(:,i)./dr_TGO(i);
   end
  
   check = (u_MRO.*(lander-[0,0,0]')).^2-(lander-[0,0,0]').^2 - r.^2;

% for i = 1
% u = 0; v = 0; w = 0; 
% X1 = lander(1); Y1 = lander(2); Z1 = lander(3);
% X2 = MRO_pos(1,i); Y2 = MRO_pos(2,i); Z2 = MRO_pos(3,i);
% a=u^2+v^2+w^2;
% b=-2*(u*(X2-X1)+v*(Y2-Y1)+w*(Z2-Z1));
% c=(X2-X1)^2+(Y2-Y1)^2+(Z2-Z1)^2-r^2;
% d(i)=b^2-4*a*c;
% end

for i = 1:length(et)
et_fig(i) = (et(i)-et(1));
end

%Determine whether mars interferes with line of sight from lander to
%orbiter at each ephemeris time. Use data rates for the orbiters to scale
%the results

MRO_inview = zeros(length(et),1);
Odyssey_inview = zeros(length(et),1);
Maven_inview = zeros(length(et),1);
TGO_inview = zeros(length(et),1);
for i = 1:length(et)
    if norm(lander+u_MRO(:,i))< r
    MRO_inview(i) = 0;
    else 
    MRO_inview(i) = 128e3;
    end
    
    if norm(lander+u_Odyssey(:,i))<r
        Odyssey_inview(i) = 0;
    else 
        Odyssey_inview(i) = 256e3;
    end
    
    if norm(lander+u_Maven(:,i))<r
        Maven_inview(i) = 0;
    else 
        Maven_inview(i) = 32e3;
    end
    if 19.875<= et_fig(i)/3600 && et_fig(i)/3600<=20.125 || 14.875<=et_fig(i)/3600 && et_fig(i)/3600 <= 15.125
        
        TGO_inview(i) = 512e3;
    else
        TGO_inview(i) = 0;
    end

end


% for i = 1:length(et)
% MRO_check(i) = (u_MRO(:,i)'*(lander)).^2 - ((lander'*lander)-r^2);
% Odyssey_check(i) = (u_Odyssey(:,i)'*(lander)).^2 - ((lander'*lander)-r^2);
% Maven_check(i) = (u_Maven(:,i)'*(lander)).^2 - ((lander'*lander)-r^2);
% if MRO_check(i) < 0
%     MRO_check(i) = 2048;
% else 
%     MRO_check(i) = 0;
% end
% if Odyssey_check(i) < 0
%     Odyssey_check(i) = 256;
% else 
%     Odyssey_check(i) = 0;
% end
% if Maven_check(i) < 0
%     Maven_check(i) = 2048;
% else 
%     Maven_check(i) = 0;
% end
% 
%     
% end


% 
% figure
% plot(et_fig./3600, MRO_inview/1e3,'r');
% hold on
% plot(et_fig./3600, Odyssey_inview/1e3,'g');
% plot(et_fig./3600,Maven_inview/1e3,'b');
% plot(et_fig./3600,TGO_inview/1e3,'m');
% legend('MRO','Odyssey','MAVEN','TGO');
% % title('Communications Access Windows per Sol')
% ylabel('Data transfer capability [kbps]');
% xlabel('Time [h]')
% xlim([0,24]);
% hold off
% 
for i = 1:length(et)
    TotCommsCap(i) = sum([MRO_inview(i), TGO_inview(i)./2]); %sum([MRO_inview(i), Odyssey_inview(i), Maven_inview(i),TGO_inview(i)]);
end
% figure
% plot(et_fig./3600, TotCommsCap, '-');
% % title('Total Data Transfer Capacity Per Sol');
% ylabel('Data transfer capability [kbps]');
% xlabel('Time [h]')

% Maven_Data_per_sol = trapz(et, Maven_inview);
% Odyssey_Data_per_sol = trapz(et, Odyssey_inview);
% MRO_Data_per_sol = trapz(et, MRO_inview);
%amount of data able to transfer per day in bits
Data_per_sol = trapz(et, TotCommsCap)
mbval = (Data_per_sol/8)/1e6;
mbval_week = mbval*7;
   cspice_kclear
end


