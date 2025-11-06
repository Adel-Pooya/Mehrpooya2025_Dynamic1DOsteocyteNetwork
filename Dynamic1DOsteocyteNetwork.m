%% Adel Mehrpooya, Vivien J Challis, Pascal R Buenzli (2025)
tic
%% Requirements
if exist("datadir")==7
    rmdir("datadir", 's')
end
clear; 
close all;

useForceGradient = true;%false; % set to true for force gradient mode

if useForceGradient
    % Code for force gradient
    T_final = 500; % Final time
else
    % Code for uniform force
    T_final = 400; % Final time
end

plot_fluxes = 1;

%% Discrete model                                                  

%% Parameters

% Define applied mechanical force:
F_bar = 0.04; % reference mechanical force; Unit: [N] (or 0.04*10^6 uN)
F0 = F_bar; % mechanical force coefficient: Force(x) = F1*x + F0; Unit: [N]

if useForceGradient
    % Code for force gradient
    F1 = 0.00005; % (for force gradient) mechanical force coefficient: Force(x) = F1*x + F0; Unit: [N/um]
    F0_new = F0; % step change coefficient in mechanical force: Force(x) = F1_new*x + F0_new; Unit: [N]
else
    % Code for uniform force
    F1 = 0; % (for uniform force) mechanical force coefficient: Force(x) = F1*x + F0; Unit: [N/um]
    F0_new = F0 + 0.2*F0; % (for uniform force) step change coefficient in mechanical force: Force(x) = F1_new*x + F0_new; Unit: [N]
end
    
F1_new = F1; % mechanical force coefficient: Force_new(x) = F1_new*x + F0_new_force_change; Unit: [N/um]* % first step change in force
F0_new2 = F0; % step change coefficient in mechanical force: Force(x) = F1_new2*x + F0_new2; Unit: [N]
F1_new2 = F1; % mechanical force coefficient: Force_new2(x) = F1_new2*x + F0_new2_force_change; Unit: [N/um]* % second step change in force

% network
if useForceGradient
    % Code for force gradient
    K_predefined = 801; % number of predefined nodes
else
    % Code for uniform force
    K_predefined = 81; % number of predefined nodes
end
K = 21; % Initial number of nodes in the network, including two additional boundary nodes.
        % The boundary nodes are positioned as follows:
        % - One node to the left of the leftmost osteocyte
        % - One node to the right of the rightmost osteocyte
        % These boundary nodes are used to define the limits of the computational domain.
q = 0.32; 
q_left1 = q; % fraction of signalling molecules that jump to the left to first nearest osteocyte 
q_right1 = q; % fraction of signalling molecules that jump to the right to first nearest osteocyte 
q_0 = 1-q_left1-q_right1; % fraction of signalling molecules that stay at the current osteocyte

% space and time scales
L_bar = 200; % reference network length; Unit: [um]

T = 5; % [day]
t_force_change1 = 40; % first time point at which force changes;
t_force_change2 = 200; % second time point at which force changes (this only affects uniform force case in practice);

% space discretisation
Delta_x = 10; % space increment; Unit: [um]
i0 = (K_predefined-1)/2 + 1; % index of middle osteocyte (+1 since indices start at 1)
x_i = zeros(1,K_predefined);
for i = 1:K_predefined
    x_i(i) = (i-i0)*Delta_x; 
end   

% time discretisation
D = 3200000; % diffusivity; Unit: [um^2/day]
Delta_t = (((Delta_x)^2)/(2*D))*(q_left1 + q_right1) % Unit: [day]  
T_final_over_Delta_t = ceil(T_final/Delta_t);
t_end = T_final_over_Delta_t*Delta_t; 
                                               
%coefficient of change rate in the number of osteoblasts/osteoclasts
Ob_0 = 0; %initial number of each boundary's osteoblasts
Oc_0 = 0; %initial number of each boundary's osteoclasts

k_f = 20; % speed of osteoblast-driven boundary movement; Unit: [um/day]
k_r = 100; % speed of osteoclast-driven boundary movement; Unit: [um/day]       
alpha_Ob = 0.04; % creation rate of osteoblasts; Unit: [1/day]    
alpha_Oc = 0.2; % creation rate of osteoclasts; Unit: [1/day]  
A_Ob = 0.2; % elimination rate of osteoblasts; Unit: [1/day] 
A_Oc = 1; % elimination rate of osteoclasts; Unit: [1/day]   

J_0 = 8.42262e+07; % threshold flux value; Unit: [1/day]

% reaction coefficients
lambda2 = 320; % degradation rate; Unit: [1/day]
lambda1 = 11059200; % creation rate; Unit: [1/day]

% diffusion length of signalling molecules
Lambda_D = sqrt(D/lambda2);

% Define applied force for all times:
F0_vec = F0*ones(1,ceil(t_end/Delta_t)+1);
F1_vec = F1*ones(1,ceil(t_end/Delta_t)+1);
F0_vec(round(t_force_change1/Delta_t):round(t_force_change2/Delta_t)-1) = F0_new;
F1_vec(round(t_force_change1/Delta_t):round(t_force_change2/Delta_t)-1) = F1_new;

% left & right boundaries' position
b_l = zeros(1,ceil(t_end/Delta_t)); % left boundary position at differet time points 
b_l(1) = -0.5*Delta_x*(K-1); % left boundary position at t=0 

b_r = zeros(1,ceil(t_end/Delta_t)); % right boundary position at differet time points 
b_r(1) = 0.5*Delta_x*(K-1); % right boundary position at t=0

% flux (jumps)
flux_left = zeros(1,ceil(t_end/Delta_t));
flux_left(1) = NaN; 
flux_right = zeros(1,ceil(t_end/Delta_t));
flux_right(1) = NaN;

% particle count at each boundary
N_l = zeros(1,ceil(t_end/Delta_t)); % number of particles at left boundary at differet time points 
N_l(1) = 0; 
N_r = zeros(1,ceil(t_end/Delta_t)); % number of particles at right boundary at differet time points 
N_r(1) = 0; 

% number of left & right osteoblasts 
Ob_l = zeros(1,ceil(t_end/Delta_t)); % number of left boundary's osteoblasts at differet time points 
Ob_l(1) = Ob_0; % number of left boundary's osteoblasts at t=0 

Ob_r = zeros(1,ceil(t_end/Delta_t)); % number of right boundary's osteoblasts at differet time points 
Ob_r(1) = Ob_0; %number of right boundary's osteoblasts at t=0

% number of left & right osteoclasts 
Oc_l = zeros(1,ceil(t_end/Delta_t)); % number of left boundary's osteoclasts at differet time points 
Oc_l(1) = Oc_0; % number of left boundary's osteoclasts at t=0 

Oc_r = zeros(1,ceil(t_end/Delta_t)); % number of right boundary's osteoblasts/osteoblasts at differet time points 
Oc_r(1) = Oc_0; % number of right boundary's osteoblasts/osteoblasts at t=0

% left & right boundaries' node index
i0 = (K_predefined-1)/2 + 1; % index of middle osteocyte (+1 since indices start at 1)
i_l = zeros(1,ceil(t_end/Delta_t)); % left boundary index at differet time points
i_l(1) = ceil(b_l(1)/Delta_x + 1 + i0); % i_l = i0 - b_l(1); 495
i_r = zeros(1,ceil(t_end/Delta_t)); % right boundary index at differet time points
i_r(1) = floor(b_r(1)/Delta_x - 1 + i0); % Osteocyte territory

%% Initialisation (discrete model)
t = 0;
current_time_every = 500000;
iter = 1; % number of time steps taken
if useForceGradient
    % Code for force gradient
    save_every = 50000000; % iterations
else
    % Code for uniform force
    save_every = 500000; % iterations
end
save_iter = 0; % counter: number of states saved

% allocate space to save states
N = zeros(1,K_predefined); % distribution of molecules (osteocyte occupancy)

% initial distribution of signalling molecules within the network
for i = i_l(1):i_r(1)
    x = (i-i0)*Delta_x; 
    l1 = b_r(1) - b_l(1);
    l2 = (i_r(1)+1 - i_l(1)+1).*Delta_x;
    N(i) = Delta_x.*n_bar(x,l1,l2,F0,F1,Delta_x,F_bar,L_bar,lambda1,lambda2,Lambda_D);
end                 

%% Time stepping (discrete)
mkdir('datadir');
save('datadir/state_0.mat', 'N'); % save data for plotting later
N_old = zeros(size(N)); % allocate, for time stepping
N_resorption_l = 0; % number of molecules added to the left boundary node due to resorption of a node
N_resorption_r = 0; % number of molecules added to the right boundary node due to resorption of a node
while t < t_end - 0.1*Delta_t
    % keep old values
    N_old = N;
    
    % transport
    for i=(i_l(iter)-1):(i_r(iter)+1)
        N(i)=N_old(i)*q_0+N_old(i+1)*q_left1+N_old(i-1)*q_right1;
    end
    
    % reaction  
    for i=i_l(iter):i_r(iter)
        x_i = (i-i0)*Delta_x;
        N(i) = N(i) + ( lambda1*( L_bar*(F0_vec(iter)+F1_vec(iter)*x_i)/((b_r(iter) - b_l(iter))*F_bar)) - lambda2*N(i) )*Delta_t;
    end
    
    % The first entry of 'N_l' and 'N_r' represents the left and right 
    % boundary molecules at t=0, respectively.
    N_l(iter+1) = N_l(iter) + N_resorption_l + N(i_l(iter)-1);
    N(i_l(iter)-1)=0;
    N_resorption_l = 0;

    N_r(iter+1)=N_r(iter) + N_resorption_r + N(i_r(iter)+1);
    N(i_r(iter)+1)=0;
    N_resorption_r = 0;
    
    % fluxes
    % The first entry of 'flux_left' and 'flux_right' represents flux at t=0.
    % left boundary (disc)
    flux_left(iter+1) = (N_l(iter+1)-N_l(iter))./Delta_t;            
    % right boundary (disc)
    flux_right(iter+1) = (N_r(iter+1)-N_r(iter))./Delta_t;

    %% implementing formation/resorption
    % discrete   
    Ob_l(iter+1) = Ob_l(iter) + Delta_t*(alpha_Ob*positive_part((flux_left(iter+1)-J_0)/J_0) - A_Ob*Ob_l(iter));
    Oc_l(iter+1) = Oc_l(iter) + Delta_t*(alpha_Oc*positive_part((J_0-flux_left(iter+1))/J_0) - A_Oc*Oc_l(iter));
    b_l(iter+1) = b_l(iter) - Delta_t*(k_f*Ob_l(iter+1) - k_r*Oc_l(iter+1)); 
    % Move left boundary if needed. If i_l(iter+1) becomes equal to the 
    % index of the predefined node that is located to the left of the 
    % current leftmost node of the computational network, then this 
    % predefined node becomes the new leftmost node in the network.
    % If boundary is resorbed, add molecules 
    % at the left node are incorporated into N_l
    assert(b_l(iter)<(K_predefined-i0)*Delta_x)
    i_l(iter+1) = ceil(b_l(iter+1)/Delta_x + i0);
    assert(i_l(iter+1)>0);
    if b_l(iter+1)>=(i_l(iter)-i0)*Delta_x % Condition for resorption at left boundary (Analogous of Eq. (22) in the paper)
        N_resorption_l = N(i_l(iter));
        N(i_l(iter)) = 0;
    end  

    Ob_r(iter+1) = Ob_r(iter) + Delta_t*(alpha_Ob*positive_part((flux_right(iter+1)-J_0)/J_0) - A_Ob*Ob_r(iter));
    Oc_r(iter+1) = Oc_r(iter) + Delta_t*(alpha_Oc*positive_part((J_0-flux_right(iter+1))/J_0) - A_Oc*Oc_r(iter));
    b_r(iter+1) = b_r(iter) + Delta_t*(k_f*Ob_r(iter+1) - k_r*Oc_r(iter+1));
    % Move right boundary if needed. If i_r(iter+1) becomes equal to the 
    % index of the predefined node that is located to the right of the 
    % current rightmost node of the computational network, then this 
    % predefined node becomes the new rightmost node in the network.
    % If boundary is resorbed, add particles
    % at the right node are incorporated into N_r
    assert(b_r(iter+1)>(1-i0)*Delta_x)
    i_r(iter+1) = floor(b_r(iter+1)/Delta_x + i0);
    assert(i_r(iter+1)<K_predefined+1);
    if b_r(iter+1)<=(i_r(iter)-i0)*Delta_x % Condition for resorption at right boundary (Eq. (22) in the paper)
        N_resorption_r = N(i_r(iter));
        N(i_r(iter)) = 0;
    end
  
    % Save states to files state_1.mat, state_2.mat, etc.
    if rem(iter, save_every) == 0        
        save_iter = save_iter + 1;
        save(sprintf("datadir/state_%d.mat", save_iter), 'N')
    end

    % update time
    iter = iter + 1;
    t = t + Delta_t; 
      
    if rem(iter, current_time_every) == 0  
     current_time = t   
    end
end
%%
disp("Discrete model: done")
run_time1 = toc 

%% Plotting
tic
t_vec = 0:Delta_t:t_end;
grayColor = [.7 .7 .7]; 
if useForceGradient
    % Code for force gradient
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 F(x)
    figure;
    % plot force vs network length; F(x)
    hold on;
 
    pl_exact_1 = fplot(@(x) Force(x,F0,0), [-100, 100], 'b--','LineWidth', 2); 
    pl_exact_2 = fplot(@(x) Force(x,F0,F1), [-100, 100], 'r-','LineWidth', 2); 
    
    xlim([b_l(1), b_r(1)]);
    ylim([0,0.045]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([-100 -50 0 50 100]);
    xticklabels({'$-100$', '$-50$', '$0$', '$50$', '$100$'});
    yticks([0 0.01 0.02 0.03 0.04]);
    yticklabels({'$0$', '0.01', '$0.02$', '$0.03$', '$0.04$'});
    
    xlabel('$x \mathrm{[\mu m]}$','interpreter','latex','fontsize',12);
    ylabel('$F(x) \mathrm{[N]}$','interpreter','latex','fontsize',12);
    lg=legend([pl_exact_1,pl_exact_2],'$\bar{F}$','$F$','interpreter','latex','fontsize',12,...
        'orientation','vertical','position',[0.79 0.74 0.07 0.034]);
    box on;
    hold off;
    drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 J(t)
    figure;
    % Flux at boundaries
    if plot_fluxes
        t = 0:save_every*Delta_t:t_end;
        
        xlim([0,t_end]);
        ylim([7*10^7,10^8]);
        set(gca,'TickLabelInterpreter','latex','fontsize',12);
        xticks([0 100 200 300 400 500]);
        xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$', '$500$'});
        yticks([7*10^7 8*10^7 9*10^7 10*10^7]);
        yticklabels({'$7 \times 10^7$', '$8 \times 10^7$', '$9 \times 10^7$', '$10 \times 10^7$'});
        xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
        ylabel('$J \mathrm{[1/day]}$','interpreter','latex','fontsize',12);
        hold on;
        pl_Jm_t = plot(t_vec, flux_left, '-', 'Color', [0.7,0.7,0.7], 'LineWidth', 2);   
        pl_Jp_t = plot(t_vec, flux_right, '-k', 'LineWidth', 1);
        
        x = 0:0.0001:t_end;
        y = 0*x+J_0;
        pl_exact = plot(x, y, '--r', 'LineWidth', 2);
    
        lg=legend([pl_Jm_t,pl_Jp_t,pl_exact],'$J_{-}$','$J_{+}$','$J_{0}$',...
            'interpreter','latex','fontsize',12,'orientation','horizontal','position',[0.223 0.94 1.0 0.034]);
        box on;
        hold off;
    end
    drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 b-(t) b+(t) i-(t) i+(t)
    figure;
    % plot boundary position and boundary node index for discrete model 
    hold on;
   
    pl_l = plot(t_vec, b_l, 'k-', 'LineWidth', 2);
    pl_r = plot(t_vec, b_r, 'k-', 'LineWidth', 2);
    
    pl_l_node = plot(t_vec, (i_l-i0)*Delta_x, 'Color', grayColor, 'LineWidth', 2);
    pl_r_node = plot(t_vec, (i_r-i0)*Delta_x, 'Color', grayColor, 'LineWidth', 2);
    
    lg=legend([pl_r,pl_r_node],'$b_{\pm}$','$i_{\pm}\Delta x$',...
        'interpreter','latex','fontsize',12,'orientation','vertical','position',[0.77 0.72 0.07 0.034]);
     
    xlim([0, t_end]);
    ylim([-200, 275]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300 400 500]);
    xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$', '$500$'});
    yticks([-200 -100 0 100 200]);
    yticklabels({'$-200$', '$-100$', '$0$', '$100$', '$200$'});
    
    xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
    ylabel('$x(t) \mathrm{[\mu m]}$','interpreter','latex','fontsize',12); 
    box on;
    hold off;
    drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 bone_length(t) 
    figure;
    % plot network length vs time 
    plot(t_vec, b_r - b_l, 'Color', grayColor,'LineWidth', 2);
       
    xlim([0, 500]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300 400 500]);
    xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$', '$500$'});
    yticks([180 200 220 240]);
    yticklabels({'$180$', '$200$', '$220$', '$240$'});
    
    xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
    ylabel('$L \mathrm{[\mu m]}$','interpreter','latex','fontsize',12); 
    box on;
    hold off;
    drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5 Ob-(t) Oc-(t) 
    figure;
    % plot the speed of formation and resorption at the left boundary
    hold on;
    
    pl_ob_l = plot(t_vec, k_f.*Ob_l, 'b-', 'LineWidth', 2);
    save('datadir/Ob_l.mat', 'Ob_l') % save data for plotting later
    clear Ob_l;
    pl_oc_l = plot(t_vec, k_r.*Oc_l, 'r-', 'LineWidth', 2);
    save('datadir/Oc_l.mat', 'Oc_l') % save data for plotting later
    clear Oc_l;

    lg=legend([pl_ob_l,pl_oc_l],'$k_\mathrm{f}\mathrm{Ob}_{-}$','$k_\mathrm{r}\mathrm{Oc}_{-}$',...
        'interpreter','latex','fontsize',12,'orientation','vertical','position',[0.77 0.84 0.07 0.034]);
    
    xlim([0, t_end]);
    ylim([0, 2]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300 400 500]);
    xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$', '$500$'});
    yticks([0 1 2]);
    yticklabels({'$0$', '$1$', '$2$'});
    
    xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
    ylabel('$k_\mathrm{f}{\mathrm{Ob}_-}, k_\mathrm{r}{\mathrm{Oc}_-} \mathrm{[\mu m/day]}$',...
        'interpreter','latex','fontsize',12); 
    box on;
    hold off;
    drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6 Ob+(t) Oc+(t)  
    figure;
    % plot the speed of formation and resorption at the right boundary  
    hold on
    
    pl_ob_r = plot(t_vec, k_f.*Ob_r, 'b-', 'LineWidth', 2);
    save('datadir/Ob_r.mat', 'Ob_r') % save data for plotting later
    clear Ob_r;
    pl_oc_r = plot(t_vec, k_r.*Oc_r, 'r-', 'LineWidth', 2);
    save('datadir/Oc_r.mat', 'Oc_r') % save data for plotting later
    clear Oc_r;
    
    lg=legend([pl_ob_r,pl_oc_r],'$k_\mathrm{f}\mathrm{Ob}_{+}$','$k_\mathrm{r}\mathrm{Oc}_{+}$',...
        'interpreter','latex','fontsize',12,'orientation','vertical','position',[0.77 0.63 0.07 0.034]);
    
    xlim([0, t_end]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300 400 500]);
    xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$', '$500$'});
    yticks([0 0.1 0.2 0.3 0.4]);
    yticklabels({'$0$', '$0.1$', '$0.2$', '$0.3$', '$0.4$'});
    
    xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
    ylabel('$k_\mathrm{f}{\mathrm{Ob}_+}, k_\mathrm{r}{\mathrm{Oc}_+} \mathrm{[\mu m/day]}$',...
        'interpreter','latex','fontsize',12); 
    box on;
    hold off;
    drawnow;
else
    % Code for uniform force
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 n(x,0) = N(x,0)/Delta_x
    %%
    figure;
    % plot initial condition
    t=0;
    x = zeros(1,K_predefined); % node positions
    for i = 1:K_predefined
        x(i) = (i-i0)*Delta_x; 
    end
    load('datadir/state_0.mat', 'N');
    hold on;

    % plot discrete
    for i=1:length(x)
        rectangle('position', [x(i)-(x(i)-(i-1-i0)*Delta_x)/2, 0, ((i+1-i0)*Delta_x-(i-1-i0)*Delta_x)/2,...
            N(i)/Delta_x], 'EdgeColor', [0.3,0.3,0.3], 'FaceColor', [0.85,0.85,0.85]);                                                       
    end
    pl_disc = fill(NaN, NaN, 'w', 'EdgeColor', [0.3,0.3,0.3],...
        'FaceColor', [0.85,0.85,0.85]); % for the legend
    
    % plot exact 
    b = b_r(1);
    l1 = b_r(1) - b_l(1);
    l2 = (i_r(1)+1 - i_l(1)+1).*Delta_x;
    pl_exact = fplot(@(x) n_bar(x,l1,l2,F0,F1,Delta_x,F_bar,L_bar,lambda1,lambda2,Lambda_D),...
        [-b, b], 'r--','LineWidth', 2);
                               
    % Draw vertical dashed lines at boundaries
    b_positions = [b_l(1) b_r(1)];
    for i = 1:length(b_positions)
        xline(b_positions(i), '-', 'Color', 'b', 'LineWidth', 1.5);
    end
    
    % Draw vertical dashed lines at boundary nodes
    i_positions = [(i_l(1) - 1 - i0)*Delta_x (i_r(1) + 1 - i0)*Delta_x];
    for i = 1:length(i_positions)
        xline(i_positions(i),'--', 'Color', 'r', 'LineWidth', 1.5);
    end
    
    title(sprintf("$t=$%g [day]", t), 'interpreter', 'latex', 'FontSize', 12); % Customize the font size for the title
    titleHandle = get(gca, 'Title'); % Get the title handle
    set(titleHandle, 'FontSize', 12); % Set the desired font size
    
    xlim([-150,150]);
    ylim([0,1700]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([-100 -50 0 50 100]);
    xticklabels({'$-100$', '$-50$', '$0$', '$50$', '$100$'});
    yticks([0 500 1000 1500]);
    yticklabels({'$0$', '$500$', '$1000$', '$1500$'});
    xlabel('$x \mathrm{[\mu m]}$','interpreter','latex','fontsize',12);
    ylabel('$n(x,t)$ [$\mu\mathrm{m}^{-1}$]','interpreter','latex','fontsize',12); 
    
    lg=legend([pl_disc,pl_exact], '$n$','$\bar{n}$','interpreter','latex',...
        'fontsize',12,'orientation','vertical','position',[0.665 0.765 0.07 0.034]);
    
    ax = gca;
    ax.Layer = 'top';
            
    pbaspect([1.618 1 1]); % Sets the vertical-to-horizontal aspect ratio to the golden ratio

    box on;
    hold off;
    drawnow;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2 n(x,t) = N(x,t)/Delta_x
    % 2 plot concentration snapshots
    b_l_positions = [];
    b_r_positions = [];
    t_vertical_lines = [0];
    for s = 1:save_iter
        t = s*save_every*Delta_t;
        if (t==8*save_every*Delta_t || t==14*save_every*Delta_t || t==40*save_every*Delta_t ...
                || t==43*save_every*Delta_t || t==80*save_every*Delta_t)
            load(sprintf("datadir/state_%d.mat", s), 'N');
    
            % pause(0.01) % slow down
            figure;
            clf;
            hold on;
    
            % plot discrete
            for i=1:length(x)
                rectangle('position', [x(i)-(x(i)-(i-1-i0)*Delta_x)/2, 0, ((i+1-i0)*Delta_x-(i-1-i0)...
                    *Delta_x)/2,N(i)/Delta_x],'EdgeColor', [0.3,0.3,0.3], 'FaceColor', [0.85,0.85,0.85]);                                                       
            end
            pl_disc = fill(NaN, NaN, 'w', 'EdgeColor', [0.3,0.3,0.3], 'FaceColor',...
                [0.85,0.85,0.85]); % for the legend
            
            % plot exact:
            b = (i_r(s*save_every) + 1 - i0)*Delta_x;
            l1 = b_r(s*save_every) - b_l(s*save_every);
            l2 = (i_r(s*save_every)+1 - i_l(s*save_every)+1).*Delta_x;
            if (s*save_every<=4000000)          
                pl_exact = fplot(@(x) n_bar(x,l1,l2,F0,F1,Delta_x,F_bar,L_bar,lambda1,lambda2,Lambda_D),...
                    [-b, b], 'r--','LineWidth', 2); % fplot doesn't return a handle
            elseif (s*save_every<=20000000) % first step change in force  
                pl_exact = fplot(@(x) n_bar(x,l1,l2,F0_new,F1_new,Delta_x,F_bar,L_bar,lambda1,lambda2,Lambda_D),...
                    [-b, b], 'r--','LineWidth', 2); % fplot doesn't return a handle
            else % second step change in force 
                pl_exact = fplot(@(x) n_bar(x,l1,l2,F0_new2,F1_new2,Delta_x,F_bar,L_bar,lambda1,lambda2,Lambda_D),...
                    [-b, b], 'r--','LineWidth', 2); % fplot doesn't return a handle     
            end 
                    
            % Draw vertical dashed lines at boundaries
            b_positions = [b_l(s*save_every) b_r(s*save_every)];
            for i = 1:length(b_positions)
                xline(b_positions(i), '-', 'Color', 'b', 'LineWidth', 1.5);
            end
    
            % Draw vertical dashed lines at boundary nodes
            i_positions = [(i_l(s*save_every) - 1 - i0)*Delta_x (i_r(s*save_every) + 1 - i0)*Delta_x];
            for i = 1:length(i_positions)
                xline(i_positions(i), '--', 'Color', 'r', 'LineWidth', 1.5);
            end
                     
            b_l_positions(end+1) = b_l(s*save_every);
            b_r_positions(end+1) = b_r(s*save_every);
            t_vertical_lines(end+1) = t; 
    
            % decoration
            title(sprintf("$t=$%g days", t), 'interpreter', 'latex', 'FontSize', 12); % Customize the font size for the title
            titleHandle = get(gca, 'Title'); % Get the title handle
            set(titleHandle, 'FontSize', 12); % Set the desired font size
    
            xlim([-150,150]);
            ylim([0,1700]);
            set(gca,'TickLabelInterpreter','latex','fontsize',12);
            xticks([-100 -50 0 50 100]);
            xticklabels({'$-100$', '$-50$', '$0$', '$50$', '$100$'});
            yticks([0 500 1000 1500]);
            yticklabels({'$0$', '$500$', '$1000$', '$1500$'});
            xlabel('$x \mathrm{[\mu m]}$','interpreter','latex','fontsize',12);
            ylabel('$n(x,t)$ [$\mu\mathrm{m}^{-1}$]','interpreter','latex','fontsize',12); 
            
            lg=legend([pl_disc,pl_exact], '$n$','$\bar{n}$','interpreter','latex',...
                'fontsize',12,'orientation','vertical','position',[0.665 0.765 0.07 0.034]);

            ax = gca;
            ax.Layer = 'top';
                    
            pbaspect([1.618 1 1]); % Sets the vertical-to-horizontal aspect ratio to the golden ratio
    
            box on;
            hold off;
            drawnow;
        end
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3 F(L)
    figure;
    % plot force vs bone length; F(L)
    hold on;
    
    % plot F(L)
    pl_F_exact = fplot(@(L) F_of_L(L,F_bar,L_bar,Lambda_D),[0, 300], 'r-','LineWidth', 2); 

    % plot (L(t),F(t)) 
    initial_bone_length = b_r(1) - b_l(1);
    final_bone_length = b_r(end) - b_l(end); 
    pl_F_trajectory = plot(b_r - b_l,F0_vec,'Color','b','LineWidth', 2);
    
    lg=legend([pl_F_exact,pl_F_trajectory],'$F(L)$', '$(L(t), F(t))$','interpreter',...
        'latex','fontsize',12,'orientation','vertical','position',[0.74 0.64 0.07 0.034]);
        
    xlim([0, 300]);
    ylim([0,0.0505]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300]);
    xticklabels({'$0$', '$100$', '$200$', '$300$'});
    yticks([0 0.01 0.02 0.03 0.04 0.05]);
    yticklabels({'$0$', '$0.01$', '$0.02$', '$0.03$', '$0.04$', '$0.05$'});
    
    xlabel('$L$ [$\mu$m]','interpreter','latex','fontsize',12);
    ylabel('$F(L)$ [N]','interpreter','latex','fontsize',12); 
    
    box on;    
    hold off;
    drawnow;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 b(t) i(t)
    figure;
    % plot boundary position and boundary node index for discrete model 
    hold on;
   
    pl_l = plot(0:Delta_t:t_end, b_l, 'k-', 'LineWidth', 2);
    pl_r = plot(0:Delta_t:t_end, b_r, 'k-', 'LineWidth', 2);
    
    pl_l_node = plot(0:Delta_t:t_end, (i_l-i0)*Delta_x, 'Color', grayColor, 'LineWidth', 2);
    pl_r_node = plot(0:Delta_t:t_end, (i_r-i0)*Delta_x, 'Color', grayColor, 'LineWidth', 2);
    

    % Draw vertical dashed lines
    for i = 1:length(t_vertical_lines)
        xline(t_vertical_lines(i), '--', 'Color', 'b', 'LineWidth', 1.5);
    end

    lg=legend([pl_l,pl_l_node],'$b_{\pm}$', '$i_{\pm}\Delta x$', ...
        'interpreter','latex','fontsize',12,'orientation','vertical','position',[0.77 0.84 0.07 0.034]);
        
    xlim([0, t_end]);
    ylim([-142,142]); % Overload
    
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300 400]);
    xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$'});
    yticks([-100 -50 0 50 100]);
    yticklabels({'$-100$', '$-50$', '$0$', '$50$', '$100$'});
    
    xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
    ylabel('$x$ [$\mu$m]','interpreter','latex','fontsize',12); 
    
    box on;    
    hold off;
    drawnow;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5 Ob+(t) Oc+(t)
    figure;
    % plot the speed of formation and resorption  
    hold on;
    
    pl_ob_r = plot(0:Delta_t:t_end, k_f.*Ob_r, 'b-', 'LineWidth', 2);
    pl_oc_r = plot(0:Delta_t:t_end, k_r.*Oc_r, 'r-', 'LineWidth', 2);
    
    lg=legend([pl_ob_r,pl_oc_r],'$k_{\mathrm{f}}\mathrm{Ob}_{+}$','$k_{\mathrm{r}}\mathrm{Oc}_{+}$',...
        'interpreter','latex','fontsize',12,'orientation','vertical','position',[0.77 0.84 0.07 0.034]);
      
    xlim([0, t_end]);
    ylim([0,3]);
    set(gca,'TickLabelInterpreter','latex','fontsize',12);
    xticks([0 100 200 300 400]);
    xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$'});
    yticks([0 1 2 3]);
    yticklabels({'$0$', '$1$', '$2$', '$3$'});
    
    xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
    ylabel('$k_{\mathrm{f}}{\mathrm{Ob}}, k_{\mathrm{r}}{\mathrm{Oc}}$ [$\mu\mathrm{m}/\mathrm{day}$]',...
        'interpreter','latex','fontsize',12); 
    
    box on;    
    drawnow;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6 J(t)
    figure;
    % plot flux for the left and right boundaries vs time; J(t)
    if plot_fluxes
        t = 0:save_every*Delta_t:t_end;

        xlim([0,t_end]);
        ylim([0,157*10^6]);
        set(gca,'TickLabelInterpreter','latex','fontsize',12);
        xticks([0 100 200 300 400]);
        xticklabels({'$0$', '$100$', '$200$', '$300$', '$400$'});  
        yticks([0 5*10^7 10*10^7 15*10^7]);
        yticklabels({'$0$', '$5 \times 10^7$', '$10 \times 10^7$', '$15 \times 10^7$'});
        xlabel('$t \mathrm{[days]}$','interpreter','latex','fontsize',12);
        ylabel('$J(t)$ [$\mathrm{day}^{-1}$]','interpreter','latex','fontsize',12);
        hold on;
   
        pl_Jp_t = plot(0:Delta_t:t_end, flux_right, '-k', 'LineWidth', 2);
    
        x = 0:0.0001:t_end;
        y = 0*x+J_0;
        pl_exact = plot(x, y, '--r', 'LineWidth', 2);
    
        lg=legend([pl_Jp_t,pl_exact],'$J_{+}$','$J_0$','interpreter','latex',...
            'fontsize',12,'orientation','vertical','position',[0.79 0.84 0.07 0.034]);
    
        hold off;
    end
    box on;
    drawnow;    
end
run_time2 = toc

%% Functions
function f= positive_part(mu)
    if mu>0
        f = mu;
    else
        f = 0;
    end
end

% Mechanical force function in terms of L = 2b+ & J_0 (Eq. (30) in the paper)
function F_of_L=F_of_L(L,F_bar,L_bar,Lambda_D) 
    F_of_L = (F_bar*L.*tanh(L_bar/(2*Lambda_D)))./(L_bar.*tanh(L/(2*Lambda_D)));
end

% Mechanical force function 
function Force=Force(x,F0,F1)  
    Force = F1.*x + F0;
end

% Exact solution for uniform force (Eq. (24) in the paper)
function y=n_bar(x,L1,L2,F0,F1,Delta_x,F_bar,L_bar,lambda1,lambda2,Lambda_D) 
    y = (lambda1/(lambda2*Delta_x)).*(( (Force(x,F0,F1)./(L1))./(F_bar/L_bar) ))...
        .*( 1-(cosh(x./Lambda_D)./cosh(L2/(2.*Lambda_D))) );   
end

% % Exact solution for force gradient (Eq. (38) in the paper)
% function y=n_bar(x,b_left,b_right,F0,F1,F_bar,L_bar,Delta_x,lambda1,lambda2,Lambda_D) 
%     P_S_minus = ( (lambda1*L_bar)./(lambda2*Delta_x*F_bar.*(b_right - b_left)) ) .*(b_left*F1 + F0);
%     P_S_plus = ( (lambda1*L_bar)./(lambda2*Delta_x*F_bar.*(b_right - b_left)) ) .*(b_right*F1 + F0);
%     C_0 = (lambda1*L_bar)./(lambda2*Delta_x*F_bar.*(b_right - b_left)); 
%     C_1 = (P_S_minus.*exp(-b_right./Lambda_D) - P_S_plus.*exp(-b_left./Lambda_D))./(exp((b_right - b_left)...
%         ./Lambda_D) - exp(-(b_right - b_left)./Lambda_D));
%     C_2 = (P_S_minus.*exp(b_right./Lambda_D) - P_S_plus.*exp(b_left./Lambda_D))./(exp(-(b_right - b_left)...
%         ./Lambda_D) - exp((b_right - b_left)./Lambda_D));
%     
%     y = C_1*exp(x./Lambda_D) + C_2*exp(-x./Lambda_D) + C_0.*Force(x,F0,F1);
% end
