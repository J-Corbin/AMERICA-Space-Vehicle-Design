clc, clear
close all

%Lander Demensions
DH = 1.1; %deck height meters
S1 = 1.66; %side one length side robotic arm on 
S2 = 1.34; % length of sides attached to the 
S3 = 2.9; %length of side with solar pannel on it

%Required Values
X_min = .5; %distance of robotic arm to edge of lander
Y_min = 10; %Minimum depth RA needs to dig
L_min = .5; %minimum length of segment of robotic arm
L_max = 2.35; %Maximum length of robotic arm
Beta = 130; %Angle of the edges of the deck measured on the inside
Alpha = 133.5;
%Optimization values
D_ang = .01; %incriment between each angle

Ang_1_1 = atan((S1/2)/(X_min));
Hyp1 = sqrt((S1/2)^2+(X_min)^2);
Hyp2 = sqrt((Hyp1*sind(Ang_1_1)+S2*cosd(180-Beta))^2+(Hyp1*cosd(Ang_1_1)+S2*sind(180-Beta))^2);
Ang_2_2 = Beta-Ang_1_1-90;
Ang_2_1 = (S2*Ang_2_2)/(Hyp2);
Ang1 = Ang_1_1;
Ang2 =  Ang_2_1;



L = [1.2,1.2];
i = 0;
for Gamma = 0:D_ang:110
    i = i+1;
   if Gamma <= Ang1;
       min(i) = X_min/cosd(Gamma);
   elseif Gamma <= Ang2;
       min(i) = sind(Beta-(90-Ang1))*Hyp1/sind(180-Ang_2_2-(Gamma-Ang1));
   else 
       Ang_3_2 = Alpha - (180-Ang_2_2-(Gamma-Ang1));
       Ang = Gamma - Ang2;
       min(i) =  (sind(Ang_3_2)*Hyp2)/sind(180-Ang_3_2-Gamma);
   end
        
end
min
hold on
plot(0:D_ang:110,min)
ylim([0,20])
