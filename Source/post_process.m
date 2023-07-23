contourf(x,y,u_n,'LevelStep',0.1,'LineStyle','--')
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title("$u$","Interpreter","latex")
xlabel("$x$","interpreter","latex");
ylabel("$y$","interpreter","latex");
c=colorbar;
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,...
    "FontSize", 18, ...
    "FontName", "Computer Modern Roman");
saveas(gcf,"Final_Project/u","epsc");

contourf(x,y,v_n,'LevelStep',0.1,'LineStyle','--')
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title("$v$","Interpreter","latex")
xlabel("$x$","interpreter","latex");
ylabel("$y$","interpreter","latex");
c=colorbar;
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,...
    "FontSize", 18, ...
    "FontName", "Computer Modern Roman");
saveas(gcf,"Final_Project/v","epsc");

contourf(x,y,sqrt(u_n.^2+v_n.^2),'LevelStep',0.1,'LineStyle','--')
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title("$|U| = \sqrt{u^2+v^2}$","Interpreter","latex")
xlabel("$x$","interpreter","latex");
ylabel("$y$","interpreter","latex");
c=colorbar;
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,...
    "FontSize", 18, ...
    "FontName", "Computer Modern Roman");
saveas(gcf,"Final_Project/velmag","epsc");

contourf(x,y,p_n,'LevelStep',0.1,'LineStyle','--')
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title("$p$","Interpreter","latex")
xlabel("$x$","interpreter","latex");
ylabel("$y$","interpreter","latex");
c=colorbar;
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,...
    "FontSize", 18, ...
    "FontName", "Computer Modern Roman");
saveas(gcf,"Final_Project/p","epsc");

contourf(x,y,rho_n,'LevelStep',0.1,'LineStyle','--')
xlim([-1.5,1.5])
ylim([-1.5,1.5])
title("$\rho$","Interpreter","latex")
xlabel("$x$","interpreter","latex");
ylabel("$y$","interpreter","latex");
c=colorbar;
c.TickLabelInterpreter="latex";
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(gca,...
    "FontSize", 18, ...
    "FontName", "Computer Modern Roman");
saveas(gcf,"Final_Project/rho","epsc");