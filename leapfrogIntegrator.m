function [rVec, vVec] = leapfrogIntegrator(nBodies, r0, v0, m, t0, tEnd, dt)
%LEAPFROGINTEGRATOR uses Leapfrog Integration to solve for the motion of nBodies
%acting under gravitational attractions. Input r0 in au, v0 in au/day,
%times in days.

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

ai = a0;
ai_1 = ai;

%Propagate
for index = 1:(length(rVec{1,1})-1)
    
    %Update current acceleration values
    ai = ai_1;
    %Clear a_i + 1
    ai_1 = zeros(size(ai));
    
    %Propagate position
    r = zeros(3,nBodies);
    for id = 1:nBodies
        rVec{1,id}(:,index + 1) = rVec{1,id}(:,index) ...
            + dt.*vVec{1,id}(:,index) + ai(:,id)*(dt^2)/2;
        r(:,id) = rVec{1,id}(:,index + 1);
    end
    
    %Update acceleration a_i+1
    for idx1 = 1:nBodies
        for idx2 = 1:nBodies
            if idx2 ~= idx1
                dx = r(1,idx1) - r(1,idx2);
                dy = r(2,idx1) - r(2,idx2);
                dz = r(3,idx1) - r(3,idx2);
                dist = sqrt(dx^2 + dy^2 + dz^2);
                ai_1(:,idx1) = ai_1(:,idx1) - (G*m(idx2)/(dist^3)).*(r(:,idx1) - r(:,idx2));
            end
        end
    end
    
    %Update velocity
    for id = 1:nBodies
        vVec{1,id}(:,index + 1) = vVec{1,id}(:,index) ...
            + (dt/2).*(ai(:,id) + ai_1(:,id));
    end
    
    
end

end