% 
% n = size(C_boundary, 1);
% 
% assert(size(C_boundary, 2) == n);
% assert(size(C_midedges, 1) == n);
% assert(size(C_midedges, 2) == n);
% 
% 
% PHI_boundary = real(C_boundary);
% PSI_boundary = imag(C_boundary);
% 
% PHI_Hessian_midedges = real(H_midedges);
% PSI_Hessian_midedges = imag(H_midedges);
% 
% I = eye(n,n);
% 
% S = [I; zeros(n,n)];
% 
% M = [PHI_boundary -PSI_boundary; PHI_Hessian_midedges -PSI_Hessian_midedges];
% 
% A = pinv(M)*S;
% %A = M(1:end-1, :) \ S(1:end-1, :);
% 
% PHI_in = real(C_in);
% PSI_in = imag(C_in);
% 
% H = [PHI_in -PSI_in]*A;







% n = size(C_boundary, 1);
% 
% assert(size(C_boundary, 2) == n);
% assert(size(C_midedges, 1) == n);
% assert(size(C_midedges, 2) == n);
% 
% 
% PHI_boundary = real(C_boundary);
% PSI_boundary = imag(C_boundary);
% 
% PHI_midedges = real(C_midedges);
% PSI_midedges = imag(C_midedges);
% 
% I = eye(n,n);
% 
% S = [I; Sampling];
% 
% M = [PHI_boundary -PSI_boundary; PHI_midedges -PSI_midedges];
% 
% %A = pinv(M)*S;
% A = M(1:end-1, :) \ S(1:end-1, :);
% 
% PHI_in = real(C_in);
% PSI_in = imag(C_in);
% 
% H = [PHI_in -PSI_in]*A;
% 
% 





assert(size(C_boundary, 1) == size(C_boundary, 2));

n = size(C_boundary, 1);

PHI_I_boundary = real(C_boundary) - eye(size(C_boundary));
PSI_boundary = imag(C_boundary);

%A = PSI_boundary(1:end-1, :) \ PHI_I_boundary(1:end-1, :);
%A = PSI_boundary(1:end-2, :) \ PHI_I_boundary(1:end-2, :);
A = pinv(PSI_boundary)*PHI_I_boundary;

PHI_in = real(C_in);
PSI_in = imag(C_in);

H = PHI_in - PSI_in*A;


