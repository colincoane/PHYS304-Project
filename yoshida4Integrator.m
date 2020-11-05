function [rVec, vVec] = yoshida4Integrator(nBodies, r0, v0, m, t0, tEnd, dt)
%YOSHIDA4INTEGRATOR uses Yoshida's 4th order Integrator to solve for the 
%motion of nBodies acting under gravitational attractions. Input r0 in au, 
%v0 in au/day, times in days.

%Gravitational Constant
Gm = 6.67430*10^-11;     % [m^3 kg^-1 s^-2] Gravitational constant
mToAu = 1.495979*10^11;  % [m/AU]
Gau = Gm / (mToAu^3);    % [au^3 kg^-1 s^-2]
G = Gau * ((3600*24)^2); % [au^3 kg^-1 day^-2]

%Create initial accelerations
a0 = zeros(size(r0));

%Populate initial accelerations
for idx1 = 1:nBodies
    for idx2 = 1:nBodies
        if idx2 ~= idx1
            dx = r0(1,idx1) - r0(1,idx2);
            dy = r0(2,idx1) - r0(2,idx2);
            dz = r0(3,idx1) - r0(3,idx2);
            dist = sqrt(dx^2 + dy^2 + dz^2);
            a0(:,idx1) = a0(:,idx1) - (G*m(idx2)/(dist^3)).*(r0(:,idx1) - r0(:,idx2));
        end
    end
end

%Create position vectors for bodies (USE CELLS)
rVec = cell(1,nBodies); % Each cell is for a single body
vVec = cell(1,nBodies);
for index = 1:nBodies
    rVec{1,index} = [r0(:,index),zeros(3,round((tEnd-t0)/dt))];
    vVec{1,index} = [v0(:,index),zeros(3,round((tEnd-t0)/dt))];
end

a = zeros(size(a0));

%Yoshida Coefficients
w0 = -((2^(1/3))/(2 - (2^(1/3))));
w1 = 1/(2 - (2^(1/3)));
c1 = w1/2;
c4 = w1/2;
c2 = (w0 + w1)/2;
c3 = (w0 + w1)/2;
d1 = w1;
d2 = w0;
d3 = w1;
d4 = 0;

%Create intermediate variables
% x1 = zeros(3,nBodies);
% x2 = x1;
% x3 = x1;
% x4 = x1;
% v1 = x1;
% v2 = x1;
% v3 = x1;
% v4 = x1;

%Propagate
for index = 1:(length(rVec{1,1})-1)
    
    %Get initial values xi = r, vi = v
    xi = zeros(3,nBodies);
    vi = xi;
    for id = 1:nBodies
        xi(:,id) = rVec{1,id}(:,index);
        vi(:,id) = vVec{1,id}(:,index);
    end
    
    %Calculate x1
    x1 = xi + c1.*dt.*vi;
    %Calculate a(x1)
    a = zeros(size(a));
    for idx1 = 1:nBodies
        for idx2 = 1:nBodies
            if idx2 ~= idx1
                dx = x1(1,idx1) - x1(1,idx2);
                dy = x1(2,idx1) - x1(2,idx2);
                dz = x1(3,idx1) - x1(3,idx2);
                dist = sqrt(dx^2 + dy^2 + dz^2);
                a(:,idx1) = a(:,idx1) - (G*m(idx2)/(dist^3)).*(x1(:,idx1) - x1(:,idx2));
            end
        end
    end
    %Calculate v1
    v1 = vi + d1.*dt.*a;
    %Calculate x2
    x2 = x1 + c2.*dt.*v1;
    %Calculate a(x2)
    a = zeros(size(a));
    for idx1 = 1:nBodies
        for idx2 = 1:nBodies
            if idx2 ~= idx1
                dx = x2(1,idx1) - x2(1,idx2);
                dy = x2(2,idx1) - x2(2,idx2);
                dz = x2(3,idx1) - x2(3,idx2);
                dist = sqrt(dx^2 + dy^2 + dz^2);
                a(:,idx1) = a(:,idx1) - (G*m(idx2)/(dist^3)).*(x2(:,idx1) - x2(:,idx2));
            end
        end
    end
    %Calculate v2
    v2 = v1 + d2.*dt.*a;
    %Calculate x3
    x3 = x2 + c3.*dt.*v2;
    %Calculate a(x3)
    a = zeros(size(a));
    for idx1 = 1:nBodies
        for idx2 = 1:nBodies
            if idx2 ~= idx1
                dx = x3(1,idx1) - x3(1,idx2);
                dy = x3(2,idx1) - x3(2,idx2);
                dz = x3(3,idx1) - x3(3,idx2);
                dist = sqrt(dx^2 + dy^2 + dz^2);
                a(:,idx1) = a(:,idx1) - (G*m(idx2)/(dist^3)).*(x3(:,idx1) - x3(:,idx2));
            end
        end
    end
    %Calculate v3
    v3 = v2 + d3.*dt.*a;
    %Calculate x4
    x4 = x3 + c4.*dt.*v3;
    %Calculate v4
    v4 = v3 + d4.*dt.*a;
    
    for id = 1:nBodies
        rVec{1,id}(:,index + 1) = x4(:,id);
        vVec{1,id}(:,index + 1) = v4(:,id);
    end
    
    
end

end