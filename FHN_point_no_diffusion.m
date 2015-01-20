
% Model parameters

a = (1./sqrt(3)) * 1.3 ; % Condition for the fix point to be in the outer branch a/e > 1/ sqrt(3)
e = 0.05 ; % Condition for excitable system e << 1


% Initial conditions:
u_init = -a * 0.4 ;
v_init = -a + a^3 ; 

[t_plot, u_plot] = evolve_FHN(u_init, v_init) ;
plot(t_plot, u_plot) ;


u_init = -a * 0.8 ;
[t_plot_1, u_plot_1] = evolve_FHN(u_init, v_init) ;
plot(t_plot, u_plot,t_plot_1, u_plot_1) ;

	
