
function Area = find_range(L1,L2,D_Theta,Depth)
i = 0;%initialize index for counting the max and min values of the arm
for Alpha = 0:D_Theta:90
    i=i+1;%increase first index
    ii= 0;%create a second index to account for depth as the second lengths rotates
    Xmax = 0;%resets the maximum x value form the last loop
    Xmin = 100;%resets the minimum x value form the last loop
   for Beta = 0:D_Theta:180
       ii = ii+1;%increases second index
       X(ii) = L1*cosd(Alpha)+ L2*cosd(180-Alpha-Beta);%find the horizontal length of the arm
       Y(ii) = L1*sind(Alpha) - L2*sind(180-Alpha-Beta)+1.1; %find the veriticle distance of the arm
       if Y(ii) <= Depth; %Checks to see if depth requirement is met
          Xmax = max(X(ii),Xmax);%Finds the max X value at this Alpha
          Xmin = min(X(ii),Xmin);%FInds the minimum X value for this alpha
          Ymax(i) = Y(ii);%Recirds the depth
       end
   end
   Xm(i) = Xmax;% finds the max of all the maximums for all thetas
   Xmi(i)= Xmin;% finds the min of all the maximums for all thetas
end
Area = -Calc_Area(min(Xmi),max(Xm),D_Theta); %uses calc area function to find the sample area

