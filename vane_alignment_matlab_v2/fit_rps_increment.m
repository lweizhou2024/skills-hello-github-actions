function [R,t,info] = fit_rps_increment(measured,target,normals,weights)
%FIT_RPS_INCREMENT Minimize normal-direction RPS residuals.
%
% residual_i = n_i dot (R*q_i + t - p_i)

    measured = double(measured);
    target = double(target);
    normals = normalizeRows(double(normals));

    if nargin < 4 || isempty(weights)
        weights = ones(size(measured,1),1);
    end
    weights = double(weights(:));

    % Linearized starting point:
    % n dot (omega x q + t) = n dot (p-q)
    J = [cross(measured,normals,2),normals];
    rightSide = sum((target-measured).*normals,2);

    weightedJ = sqrt(weights).*J;
    weightedRightSide = sqrt(weights).*rightSide;

    singularValues = svd(weightedJ,0);
    rankTolerance = max(size(weightedJ))*eps(max(singularValues));
    systemRank = sum(singularValues > rankTolerance);

    xLinear = pinv(weightedJ)*weightedRightSide;

    characteristicLength = sqrt(mean(sum( ...
        (measured-mean(measured,1)).^2,2)));
    if characteristicLength < eps
        characteristicLength = 1;
    end

    x0 = [xLinear(1:3);xLinear(4:6)/characteristicLength];

    objective = @(x) objectiveFunction( ...
        x,measured,target,normals,weights,characteristicLength);

    options = optimset( ...
        'Display','off', ...
        'TolX',1e-11, ...
        'TolFun',1e-16, ...
        'MaxIter',4000, ...
        'MaxFunEvals',12000);

    [bestX,objectiveValue,exitFlag,optimizationOutput] = ...
        fminsearch(objective,x0,options);

    R = rotationVectorToMatrix(bestX(1:3));
    t = characteristicLength*bestX(4:6);

    fitted = apply_rigid(measured,R,t);
    residual = sum((fitted-target).*normals,2);

    if numel(singularValues) >= 6 && singularValues(end) > 0
        conditionNumber = singularValues(1)/singularValues(end);
    else
        conditionNumber = Inf;
    end

    info.objective = objectiveValue;
    info.normalResiduals = residual;
    info.normalRMS = sqrt(sum(weights.*residual.^2)/sum(weights));
    info.rank = systemRank;
    info.singularValues = singularValues;
    info.conditionNumber = conditionNumber;
    info.exitFlag = exitFlag;
    info.optimizationOutput = optimizationOutput;
end


function value = objectiveFunction(x,measured,target,normals,weights,L)
    R = rotationVectorToMatrix(x(1:3));
    t = L*x(4:6);
    fitted = apply_rigid(measured,R,t);
    residual = sum((fitted-target).*normals,2);
    value = sum(weights.*residual.^2);
end


function R = rotationVectorToMatrix(rotationVector)
    angle = norm(rotationVector);

    if angle < 1e-14
        K = skew(rotationVector);
        R = eye(3)+K+0.5*K*K;
        return;
    end

    axis = rotationVector/angle;
    K = skew(axis);
    R = eye(3)+sin(angle)*K+(1-cos(angle))*K*K;
end


function K = skew(v)
    K = [0,-v(3),v(2); ...
         v(3),0,-v(1); ...
         -v(2),v(1),0];
end


function X = normalizeRows(X)
    rowLength = vecnorm(X,2,2);
    if any(rowLength < eps)
        error("A zero-length normal vector was supplied.");
    end
    X = X./rowLength;
end
