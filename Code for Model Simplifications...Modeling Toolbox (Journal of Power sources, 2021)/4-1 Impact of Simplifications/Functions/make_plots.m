%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute the plot the output and various internal states 
%
% Model Simplifications and Their Impact on Computational Complexity for an 
% Electrochemistry-Based Battery Modeling Toolbox
%
% Authors: Z. Khalik, M.C.F. Donkers, H.J. Bergveld
%
% This file is licensed under the BSD 3-Clause License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function make_plots(out,eval_time,fig_handle,set_scaleheight,param_set)
figure(fig_handle)
fontsize = 16; 
scaleheight = 0.6; 
Nout = size(out,2); 

p = out{1}.param; 
colors = {'k','r--','b--','g--','y--','m--','c--'}; 

subplot(2,4,1:2)
plot(out{1}.t,out{1}.i_app/(p.Cbat/3600),'k','LineWidth',2)
grid on
xlim([0 3600])
ylabel('Current [C-rate]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)

subplot(2,4,5:6)
for k = 1:Nout
plot(out{k}.t,out{k}.V,colors{k},'LineWidth',2)
hold on
end
grid on
xlabel('Time [s]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$V_t \ \mathrm{[V]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
xlim([0 out{k}.t(end)])

min_box = 3150; 
max_box = 3250; 
I = (out{k}.t < max_box) & (out{k}.t > min_box); % range of t near perturbation
rectangle('Position',[min(out{k}.t(I)) min(out{k}.V(I)) max(out{k}.t(I))-min(out{k}.t(I)) max(out{k}.V(I))-min(out{k}.V(I))+0.01])

if param_set
    h1 = axes('position',[0.28 0.585 0.27 0.19]);
else
    h1 = axes('position',[0.15 0.11 0.13 0.25]);
end
box on % put box around new pair of axes
I = (out{k}.t < max_box) & (out{k}.t > min_box); % range of t near perturbation
for k = 1:Nout
plot(out{k}.t(I),out{k}.V(I),colors{k},'LineWidth',2)
hold on
end
axis tight
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,4,3)
for k = 1:Nout
plot(out{k}.x,out{k}.phie(:,eval_time),colors{k},'LineWidth',2)
hold on
end
grid on
ylabel('$\phi_e \ \mathrm{[V]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

subplot(2,4,4)
for k = 1:Nout
plot(out{k}.x,out{k}.ce(:,eval_time)/p.ce0,colors{k},'LineWidth',2)
hold on
end
grid on
ylabel('$c_e/c_{e,0} \ \mathrm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

subplot(2,4,7)
for k = 1:Nout
plot(out{k}.x,out{k}.stoich(:,eval_time),colors{k},'LineWidth',2)
hold on
end
grid on
xlabel('x/L [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$\bar{c}_s/c_{s,\mathrm{max}} \ \mathrm{[-]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +2*A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

% eta = out.phis-out.phie-out.U;
subplot(2,4,8)
for k = 1:Nout
plot(out{k}.x,out{k}.jn(:,eval_time),colors{k},'LineWidth',2)
hold on
end
grid on

xlim([0 1])
xlabel('x/L [-]','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
ylabel('$j_n \ \mathrm{[mol/m^2/s]}$','Interpreter','latex','FontWeight','bold','FontSize',fontsize)
set(gca,'xtick',0:0.2:1)
if set_scaleheight
A = get(gca,'position');          % gca points at the second one
A(1,4) = A(1,4) *scaleheight;              % reduce the height by half
A(1,2) = A(1,2) +2*A(1,4)/scaleheight*(1-scaleheight);         % change the vertical position
set(gca,'position',A);     
end

end