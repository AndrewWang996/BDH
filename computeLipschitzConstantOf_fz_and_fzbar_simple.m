function [L_fz, L_fzbar] = computeLipschitzConstantOf_fz_and_fzbar_simple(L, indexOfMaxL, v, Phi, Psy)
%   v - vertices of the original cage.
%   Phi - the DOF of discrete Cauchy transform - holomorphic part.
%   Psy - the DOF of discrete Cauchy transform - anti-holomorphic part.
%   L - mxn matrix of MOC (Lipschitz constants) of all basis functions on all segments.
%   indexOfMaxL - mx1 vector where each element holds the index of the basis function with largest L.

    tic
    
    indexOfMaxL = indexOfMaxL(:, 1); %use only the first column which contains the index of the maximal eleemnt in L

    %use only the delta in Phi compared to the source cage - this gives L=0 for the identity map but doesn't help when rotation is introduced
    %debug - I should also introduce local rotation to further optimize L
    Phi_delta = Phi-v;
    %Phi_delta = Phi;

    %Lipschitz constant on each segment (edge) of the polygon
    %for each segment, we shift Phi to zero at the index where L is maximal and this way we cancel the contribution of the largest L
    
    A_Phi = abs(bsxfun(@minus, Phi_delta.', Phi_delta));
    B_Phi = L.*A_Phi(indexOfMaxL, :);
    L_fz = sum(B_Phi, 2);


	A_Psy = abs(bsxfun(@minus, Psy.', Psy));
    B_Psy = L.*A_Psy(indexOfMaxL, :);
    L_fzbar = sum(B_Psy, 2);
   
    fprintf('computeLipschitzConstantOf_fz_and_fzbar_simple time: %.5f\n', toc);
    
end




% function [L_fz, L_fzbar] = computeLipschitzConstantOf_fz_and_fzbar_simple(L, indexOfMaxL, v, Phi, Psy)
% %   v - vertices of the original cage.
% %   Phi - the DOF of discrete Cauchy transform - holomorphic part.
% %   Psy - the DOF of discrete Cauchy transform - anti-holomorphic part.
% %   L - mxn matrix of MOC (Lipschitz constants) of all basis functions on all segments.
% %   indexOfMaxL - mx1 vector where each element holds the index of the basis function with largest L.
% 
%     tic
% 
%     L_fz = L*abs(Phi-v);
% 
%     L_fzbar = L*abs(Psy);
% 
%     fprintf('computeLipschitzConstantOf_fz_and_fzbar time: %.5f\n', toc);
% 
% end

