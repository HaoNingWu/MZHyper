% On the quadrature exactness in hyperinterpolation
% by C. An and H.-N. Wu
% written by H.-N. Wu in 2022
% 

% Please add the sphere_approx_toolbox_v3.0 onto path before 
% running this demo

clear 
close all

L = 25;


% function information
funtxt_ce = {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'};
i_f = 3;
funtxt = funtxt_ce{i_f};
switch funtxt
    case 'nrWend_k0'
        rbf_k = 0;
    case 'nrWend_k1'
        rbf_k = 1;
    case 'nrWend_k2'
        rbf_k = 2;
    case 'nrWend_k3'
        rbf_k = 3;
    case 'nrWend_k4'
        rbf_k = 4;
end
switch funtxt
     case {'nrWend_k0','nrWend_k1','nrWend_k2','nrWend_k3','nrWend_k4'}
         func = @rbf_nr;
end
k = rbf_k;
delta = (3*k+3)*gamma(k+1/2)/(2*gamma(k+1));


func2 = @(x,y,z) abs(x+y+z);


    % Validation point set
    Xt = get_Xt( );
    ft = func(Xt',rbf_k);
    ft2 = func2(Xt(1,:),Xt(2,:),Xt(3,:)); ft2 = ft2';
    [ Yt ] = get_Yt( L, Xt );


for k = 1:1:L
    k
    t_now = L+k;
    % degree of point set and polynomial
    model_parameter.t = t_now;
    model_parameter.L = L;
    
    X_k = loadStd( model_parameter.t, (model_parameter.t+1)^2 );
    [m,n] = size(X_k);

    % generating function
    f = func(X_k,rbf_k);
    f2=func2(X_k(:,1),X_k(:,2),X_k(:,3));
    Y_L = getQ( X_k, L )';
    alpha = 4*pi*Y_L*f/m;
    alpha2 = 4*pi*Y_L*f2/m;

    % approximation polynomials
    L = sqrt(length(alpha))-1;
    Yt_eqp = get_Yt( L, Xt );
    pt_eqphyper = Yt_eqp * alpha;
    pt_eqphyper2 = Yt_eqp * alpha2;


   [M,~] = size(Xt');
   error1(k,1) = 4*pi/M*sqrt(sum(abs(pt_eqphyper-ft).^2));
   error1(k,2) = max(abs(pt_eqphyper-ft));

   error2(k,1) = 4*pi/M*sqrt(sum(abs(pt_eqphyper2-ft2).^2));
   error2(k,2) = max(abs(pt_eqphyper2-ft2));
  end



for k = 1:1:L
   fprintf('(%d,%d,%d) & %1.4s & %1.4s & %1.4s & %1.4s  \\\\\\hline \n',k,L+k,(L+k+1)^2, error1(k,:),error2(k,:))
end

