% On the quadrature exactness in hyperinterpolation
% by C. An and H.-N. Wu
% written by H.-N. Wu in 2017
% 

% Please add CVX [1], the Chebfun Toolbox "chebfun-master" [2], and the
% code files by S. Foucart for reproducing figures in [3] onto path before 
% running this demo


% [1] http://cvxr.com/cvx/
% [2] https://www.chebfun.org/
% [3] COMPUTING A QUANTITY OF INTEREST FROM OBSERVATIONAL DATA
% by R. DeVore, S. Foucart, G. Petrova, P. Wojtaszczyk. 
% Codes available at https://www.math.tamu.edu/~foucart/papers.html
% The zip file (reproducible+supplement) also contains CVX


clear 
close all




func_idx = 2;
% func_idx = 1: f(x) = exp(-x^2);
% func_idx = 2: f(x) = |x|^(5/2);



L = 40;

Nhyper = 41;
NGauss = 25;
NCC    = 50; 
NEqui  = 186;



% mesh for plotting
xx = -1:0.001:1; 


% generating corresponding nodes and weights
[xhyper,whyper] = legpts(Nhyper);
[xGauss,wGauss] = legpts(NGauss);
[xCC,wCC] = chebpts(NCC);
xEqui = linspace(-1,1,NEqui); xEqui = xEqui'; [wEqui,~] = optquad_C(xEqui,lege(2*L-31));




% function to be approximated
func_idx = 1;
switch func_idx
    case 1 
        fhyper = exp(-xhyper.^2);
        fGauss = exp(-xGauss.^2); 
        fCC = exp(-xCC.^2); 
        fEqui = exp(-xEqui.^2);
        ff = exp(-xx.^2);
    case 2 
        fhyper = abs(xhyper).^(5/2); 
        fGauss = abs(xGauss).^(5/2); 
        fCC = abs(xCC).^(5/2); 
        fEqui = abs(xEqui).^(5/2);
        ff = abs(xx).^(5/2);
end



for l = 0:L
    F = legpoly(l);
    Ahyper(:,l+1)  = F(xhyper)/sqrt(2/(2*l+1));
    AGauss(:,l+1) = F(xGauss)/sqrt(2/(2*l+1));
    ACC(:,l+1)      = F(xCC)/sqrt(2/(2*l+1));
    AEqui(:,l+1)    = F(xEqui)/sqrt(2/(2*l+1));
end



% hyperinterpolation coefficients
alphahyper  = Ahyper'*diag(whyper)*fhyper;
alphaGauss = AGauss'*diag(wGauss)*fGauss;
alphaCC      = ACC'*diag(wCC)*fCC;
alphaEqui    = AEqui'*diag(wEqui)*fEqui;




% approximation polynomial on xx = -1:.001:1
phyper = zeros(2001,1); pGauss = zeros(2001,1); pCC = zeros(2001,1); 
pEqui = zeros(2001,1);
for l = 0:L        
    F = legpoly(l)/sqrt(2/(2*l+1));
    phyper = phyper+alphahyper(l+1)*F(xx');
    pGauss = pGauss+alphaGauss(l+1)*F(xx');
    pCC = pCC+alphaCC(l+1)*F(xx');
    pEqui = pEqui+alphaEqui(l+1)*F(xx');
end
 


%% Figure 
fontsize_baseline = 29;
fontsize_baselinet = 32;
fontsize_baselinea = 24;
fig = figure;


axes('position',[0.05 0.57 0.4 0.38]), 
yyaxis left
plot(xx,phyper,'linewidth',3.2), box on,

yyaxis right
plot(xx,abs(ff-phyper'),'linewidth',1.6),

set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('41-pt Gauss--Legendre quad','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),...,
     legend('$\mathcal{L}^{\rm S}_{40}f$','$|\mathcal{L}^{\rm S}_{40}f-f|$','interpreter','latex', 'fontsize', fontsize_baseline )


axes('position',[0.55 0.57 0.4 0.38]), 
yyaxis left
plot(xx,pGauss,'linewidth',3.2), box on,

yyaxis right
plot(xx,abs(ff-pGauss'),'linewidth',1.6),

set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('25-pt Gauss--Legendre quad','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),...,
     legend('$\mathcal{L}_{40}f$','$|\mathcal{L}_{40}f-f|$','interpreter','latex', 'fontsize', fontsize_baseline )




axes('position',[0.05 0.07 0.4 0.38]), 
yyaxis left
plot(xx,pCC,'linewidth',3.2), box on,

yyaxis right
plot(xx,abs(ff-pCC'),'linewidth',1.6),

set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('50-pt Clenshaw--Curtis quad','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),...,
     legend('$\mathcal{L}_{40}f$','$|\mathcal{L}_{40}f-f|$','interpreter','latex', 'fontsize', fontsize_baseline )




axes('position',[0.55 0.07 0.4 0.38]), 
yyaxis left
plot(xx,pEqui,'linewidth',3.2), box on,

yyaxis right
plot(xx,abs(ff-pEqui'),'linewidth',1.6),

set(gca, 'fontsize', fontsize_baselinea),...
    xlabel('$x$','interpreter','latex', 'fontsize', fontsize_baseline),...
    title('186-equispaced-pt quad (4.1)','interpreter','latex', 'fontsize', fontsize_baselinet),  grid on,...
     set(gca, 'XMinorGrid', 'off'), set(gca, 'YMinorGrid', 'off'),...,
     legend('$\mathcal{L}_{40}f$','$|\mathcal{L}_{40}f-f|$','interpreter','latex', 'fontsize', fontsize_baseline ),
