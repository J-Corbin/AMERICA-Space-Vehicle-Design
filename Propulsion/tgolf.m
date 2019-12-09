function acceleration=tgolf(v, r, g)
    % Analytical solution for t_go from Optimal Guidance Law for Planetary
    % Landing
    
    %% Solve equations 63, 64, and 65
    a=( -2*(dot(v,v)))/(g^2/2);
    b=( -12*(dot(v,r)))/(g^2/2);
    c=( -18*(dot(r,r)))/(g^2/2);
    
    %% Solve for equations 69 and 70
    alpha= 1/3*(3*((a^2)-4*c)-4*(a^2));
    beta=1/27*(16*(a^3)-18*a*((a^2)-4*c)-27*(b^2));
    
    %% Solve for equation 72
    delta=((alpha^3)/27)+((beta^2)/4);
    
    %% Solve for equation 71
    t1=((-beta/2)+(sqrt(delta)));
    t2=(-beta/2)-(sqrt(delta));
    z=sign(t1)*abs(t1)^(1/3)+sign(t2)*abs(t2)^(1/3);
    
    %% Solve for nu from equation 68
    nu=z-(2*a/3);
    
    %% Solve equations 73 and 74
    epsilon=(a+nu)/2;
    phi=(-b/(2*nu));
    
    %% Solve for the 4 solutions of t_go
    tgo1=( sqrt(nu)+sqrt(nu-4*(epsilon-(sqrt(nu)*phi))))/2;
    tgo2=( sqrt(nu)-sqrt(nu-4*(epsilon-(sqrt(nu)*phi))))/2;
    tgo3=(-sqrt(nu)+sqrt(nu-4*(epsilon+(sqrt(nu)*phi))))/2;
    tgo4=(-sqrt(nu)-sqrt(nu-4*(epsilon+(sqrt(nu)*phi))))/2;
    
    %% Combine all solutions into single array and find real solution
    tgos=[tgo1,tgo2,tgo3,tgo4];
    
    % Iterate over entry to find real value
    for i=1:4
        if ~isreal(tgos(i))
            tgos(i)=0;
        end 
    end
    
    % Extract optimal t_go value
    besttgo=max(tgos);
    
    %% Calculate acceleration
    acceleration= -4*v/besttgo -6*r/(besttgo^2)-[0, 0, g];
end



    
    




