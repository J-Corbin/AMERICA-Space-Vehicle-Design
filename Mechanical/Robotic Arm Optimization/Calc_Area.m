function A = Calc_Area(X1,X2,D_Theta)
    A_1 = pi*X1^2*(150/360); %find the area of the part of the circle of the inner radii
    A_2 = pi*X2^2*(150/360);% find the area of the part of the circle for the outer radii
    A = A_2 - A_1; %subtract the the inner raddi from the outer to determine sample area
   
end