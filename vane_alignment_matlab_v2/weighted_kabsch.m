function [R,t,rmsError] = weighted_kabsch(source,target,weights)
%WEIGHTED_KABSCH Proper rigid registration from SOURCE to TARGET.

    source = double(source);
    target = double(target);

    if nargin < 3 || isempty(weights)
        weights = ones(size(source,1),1);
    end

    weights = double(weights(:));

    if ~isequal(size(source),size(target)) || size(source,2) ~= 3
        error("source and target must be matching N-by-3 arrays.");
    end
    if numel(weights) ~= size(source,1) || any(weights <= 0)
        error("weights must contain one positive value per point.");
    end

    weights = weights/sum(weights);

    sourceCenter = sum(source.*weights,1);
    targetCenter = sum(target.*weights,1);

    X = source-sourceCenter;
    Y = target-targetCenter;

    H = X.'*(Y.*weights);
    [U,~,V] = svd(H);

    correction = eye(3);
    correction(3,3) = sign(det(V*U.'));
    if correction(3,3) == 0
        correction(3,3) = 1;
    end

    R = V*correction*U.';
    t = targetCenter.'-R*sourceCenter.';

    fitted = apply_rigid(source,R,t);
    squaredError = sum((fitted-target).^2,2);
    rmsError = sqrt(sum(weights.*squaredError));
end
