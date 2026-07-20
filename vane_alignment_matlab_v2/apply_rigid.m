function Xout = apply_rigid(X,R,t)
%APPLY_RIGID Apply Xout = (R*X' + t)'.
    Xout = (R*X.' + t(:)).';
end
