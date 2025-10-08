function [ux,uy,uz] = u_3D_from_DMM_at_points(params, Xq, Yq, Zq)
% Evaluate your 3D DMM at arbitrary query points (columns).
% Xq,Yq,Zq: column vectors of node coordinates (e.g., x(newNode)', ...)

    % Tell DMM to work in direct-query mode
    params.switch_outputlocation = -1;

    % Attach query coordinates for the DMM to read
    params.Xq = Xq(:);
    params.Yq = Yq(:);
    params.Zq = Zq(:);

    out = run_greenfield_3D_DMM(params);

    ux = out.Sx(:);   % expect these sized like Xq
    uy = out.Sy(:);
    uz = out.Sz(:);
end
