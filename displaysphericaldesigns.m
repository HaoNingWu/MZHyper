% On the quadrature exactness in hyperinterpolation
% by C. An and H.-N. Wu
% written by H.-N. Wu in 2022
% 

% Please add the sphere_approx_toolbox_v3.0 [1] onto path before 
% running this demo


% [1] Available at https://github.com/HaoNingWu/LassoHyper/blob/main/Sphere/sphere_approx_toolbox_v3.0.zip


clear 
close all

L = 25;

 t_now = 2*L;
 model_parameter.t = t_now;
 model_parameter.L = L;


X_k = loadStd( model_parameter.t, (model_parameter.t+1)^2 );

axes('position',[0.015,0.01,0.5,0.9])
 show_r3_point_set(X_k','sphere','show','color','red'),view(39,46), 
 title('spherical 50-design: 2601 pts','interpreter','latex','fontsize', 36)

colormap jet


 t_now = L+5;
 model_parameter.t = t_now;
 model_parameter.L = L;


X_k = loadStd( model_parameter.t, (model_parameter.t+1)^2 );

axes('position',[0.515,0.01,0.5,0.9])
show_r3_point_set(X_k','sphere','show','color','red'),view(39,46), 
 title('spherical 30-design: 961 pts','interpreter','latex','fontsize', 32)
colormap jet

