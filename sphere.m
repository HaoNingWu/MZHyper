% On the quadrature exactness in hyperinterpolation
% by C. An and H.-N. Wu
% written by H.-N. Wu in 2022
% 

% Please add the sphere_approx_toolbox_v3.0 onto path before 
% running this demo


clear 
close all



L = 25;


func_idx = 1;
% func_idx = 1: Wendland function
% func_idx = 2: f(x,y,z) = |x+y+z|;


% Validation point set
Xt = get_Xt( );

% function to be approximated
switch func_idx
    case 1
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
        ft = func(Xt',rbf_k);
    case 2
        func = @(x,y,z) abs(x+y+z);
        ft = func(Xt(1,:),Xt(2,:),Xt(3,:)); ft = ft';
end




% Original hyperinterpolation
 t_now = 2*L;
    
 % degree of point set and polynomial
 model_parameter.t = t_now;
 model_parameter.L = L;


X_k = loadStd( model_parameter.t, (model_parameter.t+1)^2 );
[m,n] = size(X_k);

  % function sampling
  switch func_idx
      case 1 
        f = func(X_k,rbf_k);
      case 2
        f=func(X_k(:,1),X_k(:,2),X_k(:,3));
  end


    Y_L = getQ( X_k, L )';
    alpha = 4*pi*Y_L*f/m;


    % approximation polynomials
    L = sqrt(length(alpha))-1;
    Yt_eqp = get_Yt( L, Xt );
    pt_eqphyper = Yt_eqp * alpha;




    % scale the plotting for better visual comparison 
    Xt = Xt';
    x = Xt(:,1); y = Xt(:,2); z = Xt(:,3);
    tri = convhull([x y z]);

    C1 = ft;
    C2 = pt_eqphyper;

    [Fmax, imax] = max(C2);
    [Fmin, imin] = min(C2);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C1-Fmin);

    axes('position',[-.05 0.5 0.5 0.42]), 
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C2,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    view(39,46), axis vis3d, axis equal tight, 
    axis off, c = colorbar('east'); set(gca, 'fontsize', 22)
    title('$\mathcal{L}^{\rm{S}}_{25}f$','interpreter','latex','fontsize', 32)
    ax = gca;
    axpos = ax.Position;
    c.Position(3) = 0.5*c.Position(3);
    ax.Position = axpos;
    colormap parula


    [Fmax, imax] = max(abs(C2-C1));
    [Fmin, imin] = min(abs(C2-C1));
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(abs(C2-C1)-Fmin);

    axes('position',[0.45 0.5 0.5 0.42]), 
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, abs(C2-C1),'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    view(39,46), axis vis3d, axis equal tight, c = colorbar('east');
    set(gca, 'fontsize', 22)
    title('$|\mathcal{L}^{\rm{S}}_{25}f-f|$','interpreter','latex','fontsize', 32)
    axis off
    ax = gca;
    axpos = ax.Position;
    c.Position(3) = 0.5*c.Position(3);
    ax.Position = axpos;
    colormap parula




% Hyperinterpolation using exactness-relaxing quadrature

 Xt = get_Xt( );


 t_now = L+5;
 % degree of point set and polynomial
 model_parameter.t = t_now;
 model_parameter.L = L;


  X_k = loadStd( model_parameter.t, (model_parameter.t+1)^2 );
  [m,n] = size(X_k);

  % generating function
  switch func_idx
      case 1 
        f = func(X_k,rbf_k);
      case 2
        f=func(X_k(:,1),X_k(:,2),X_k(:,3));
  end
    Y_L = getQ( X_k, L )';
    alpha = 4*pi*Y_L*f/m;

    % approximation polynomials
    L = sqrt(length(alpha))-1;
    Yt_eqp = get_Yt( L, Xt );
    pt_eqphyper = Yt_eqp * alpha;



    % scale the plotting for better visual comparison 
    Xt = Xt';
    x = Xt(:,1); y = Xt(:,2); z = Xt(:,3);
    tri = convhull([x y z]);

    C1 = ft;
    C2 = pt_eqphyper;

    [Fmax, imax] = max(C2);
    [Fmin, imin] = min(C2);
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(C1-Fmin);

    fontsize_baselinet = 30;
    marksize = 80;
    color = lines(4);
    count_hx = 2;

    axes('position',[-0.05 0 0.5 0.42]), 
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, C2,'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    view(39,46), axis vis3d, axis equal tight, 
    axis off, c = colorbar('east'); set(gca, 'fontsize', 22)
    title('$\mathcal{L}_{25}f$','interpreter','latex','fontsize', 32)
    ax = gca;
    axpos = ax.Position;
    c.Position(3) = 0.5*c.Position(3);
    ax.Position = axpos;
    colormap parula




    [Fmax, imax] = max(abs(C2-C1));
    [Fmin, imin] = min(abs(C2-C1));
    scale = 0.5;
    FS = 1 + (scale/(Fmax-Fmin))*(abs(C2-C1)-Fmin);

    axes('position',[0.45 0 0.5 0.42]), 
    fg = trisurf(tri,x.*FS, y.*FS, z.*FS, abs(C2-C1),'facecolor','interp');
    set(fg,'EdgeColor', 'none'); 
    view(39,46), axis vis3d, axis equal tight, c = colorbar('east');
    set(gca, 'fontsize', 22),
    title('$|\mathcal{L}_{25}f-f|$','interpreter','latex','fontsize', 32)
    axis off
    ax = gca;
    axpos = ax.Position;
    c.Position(3) = 0.5*c.Position(3);
    ax.Position = axpos;
    colormap parula

   
 