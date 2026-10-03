%% Routh-Hurwitz stability criterion

%coefficients of Characteristic equation for linearized sys around
%q1*=q2*=q3*=q4*=0 (martix A) and coefficients of Characteristic equation for linearized sys around
%q1*=q2*=q3*=q4*=pi (matrix(F)

%  in this program you must give your system coefficients and the
%  Routh-Hurwitz table would be shown
%% Initialization (A)
clear ; close all; clc

coeffVectorA =[1  48/5  9099/100  57699/125  170523387/80000   308585457/50000  252673892451/16000000  447464318121/20000000  699588608555601/25600000000];
ceoffLength = length(coeffVectorA);
rhTableColumn = round(ceoffLength/2);

rhTable = zeros(ceoffLength,rhTableColumn);

rhTable(1,:) = coeffVectorA(1,1:2:ceoffLength);

if (rem(ceoffLength,2) ~= 0)
    
    rhTable(2,1:rhTableColumn - 1) = coeffVectorA(1,2:2:ceoffLength);
else
    
    rhTable(2,:) = coeffVectorA(1,2:2:ceoffLength);
end
%% Calculate Routh-Hurwitz table's rows

epss = 0.01;

for i = 3:ceoffLength
   
    
    if rhTable(i-1,:) == 0
        order = (ceoffLength - i);
        cnt1 = 0;
        cnt2 = 1;
        for j = 1:rhTableColumn - 1
            rhTable(i-1,j) = (order - cnt1) * rhTable(i-2,cnt2);
            cnt2 = cnt2 + 1;
            cnt1 = cnt1 + 2;
        end
    end
    
    for j = 1:rhTableColumn - 1
       
        firstElemUpperRow = rhTable(i-1,1);
        
     
        rhTable(i,j) = ((rhTable(i-1,1) * rhTable(i-2,j+1)) - ....
            (rhTable(i-2,1) * rhTable(i-1,j+1))) / firstElemUpperRow;
    end
    
    
   
    if rhTable(i,1) == 0
        rhTable(i,1) = epss;
    end
end
%%  Compute number of right hand side poles(unstable poles)

unstablePoles = 0;

for i = 1:ceoffLength - 1
    if sign(rhTable(i,1)) * sign(rhTable(i+1,1)) == -1
        unstablePoles = unstablePoles + 1;
    end
end

fprintf('\n Routh-Hurwitz Table:\n')
rhTable

if unstablePoles == 0
    fprintf('~~~~~> sys with A is a stable system! <~~~~~\n')
else
    fprintf('~~~~~> sys with A is an unstable system! <~~~~~\n')
end
fprintf('\n Number of right hand side poles =%2.0f\n',unstablePoles)

fprintf('\n Given polynomial coefficients roots of A :\n')
    sysRoots = roots(coeffVectorA)
    
    
    %% Initialization (F)

coeffVectorF =[1 48/5 3213/100 189/5 -6197877/8000 -726327/2000 -3506305563/16000000 8464383279/20000000 11694795124401/25600000000];
ceoffLength = length(coeffVectorF);
rhTableColumn = round(ceoffLength/2);

rhTable = zeros(ceoffLength,rhTableColumn);

rhTable(1,:) = coeffVectorF(1,1:2:ceoffLength);

if (rem(ceoffLength,2) ~= 0)
    
    rhTable(2,1:rhTableColumn - 1) = coeffVectorF(1,2:2:ceoffLength);
else
    
    rhTable(2,:) = coeffVectorF(1,2:2:ceoffLength);
end
%% Calculate Routh-Hurwitz table's rows

epss = 0.01;

for i = 3:ceoffLength
   
    
    if rhTable(i-1,:) == 0
        order = (ceoffLength - i);
        cnt1 = 0;
        cnt2 = 1;
        for j = 1:rhTableColumn - 1
            rhTable(i-1,j) = (order - cnt1) * rhTable(i-2,cnt2);
            cnt2 = cnt2 + 1;
            cnt1 = cnt1 + 2;
        end
    end
    
    for j = 1:rhTableColumn - 1
       
        firstElemUpperRow = rhTable(i-1,1);
        
     
        rhTable(i,j) = ((rhTable(i-1,1) * rhTable(i-2,j+1)) - ....
            (rhTable(i-2,1) * rhTable(i-1,j+1))) / firstElemUpperRow;
    end
    
    
   
    if rhTable(i,1) == 0
        rhTable(i,1) = epss;
    end
end
%%  Compute number of right hand side poles(unstable poles)

unstablePoles = 0;

for i = 1:ceoffLength - 1
    if sign(rhTable(i,1)) * sign(rhTable(i+1,1)) == -1
        unstablePoles = unstablePoles + 1;
    end
end

fprintf('\n Routh-Hurwitz Table:\n')
rhTable

if unstablePoles == 0
    fprintf('~~~~~> sys with F is a stable system! <~~~~~\n')
else
    fprintf('~~~~~> sys with F is an unstable system! <~~~~~\n')
end
fprintf('\n Number of right hand side poles =%2.0f\n',unstablePoles)

fprintf('\n Given polynomial coefficients roots of F :\n')
    sysRoots = roots(coeffVectorF)
