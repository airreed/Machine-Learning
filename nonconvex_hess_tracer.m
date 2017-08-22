function nonconvex_hess_tracer()
% nonconvex_hess_tracer.m is a toy wrapper to illustrate the path
% taken by Hessian descent (or Newton's method).  The corresponding
% Newton steps, evaluated at the objective, are then plotted and
% connected (visually) by a tracing line (a tracer).  The plotted points on
% the objective turn from green to red as the algorithm converges (or
% reaches a maximum iteration count, preset to 50).
% The (non convex) function here is
%
% f(x) = exp(x)*cos(2pi*sin(pi*x))

%%% create function, choose initial point, to perform hessian descent on %%%
range = 1.1;    % symmetric range over which to plot the function
[a,b] = make_fun(range);
x0 = choose_starter(a,b,range);
disp(['You picked starting point x0 = ',num2str(x0)])

%%% perform hessian descent %%%
[x,in,out] = hessian_descent(x0);

%%% plot function with hessian descent objective evaluations %%%
plot_it_all(in,out,range)

% performs hessian descent
function [x,in,out] = hessian_descent(x)

    % initializations
    grad_stop = 10^-3;
    max_its = 50;
    iter = 1;
    grad = 1;
    in = [x];    
    out = [exp(x)*cos(2*pi*sin(pi*x))];
    
    % main loop
    while abs(grad) > grad_stop && iter <= max_its
        % take gradient step
        grad = exp(x)*cos(2*pi*sin(pi*x)) - 2*pi^2*exp(x)*sin(2*pi*sin(pi*x))*cos(pi*x);
        newt = exp(x)*(2*(pi^3)*sin(pi*x)*sin(2*pi*sin(pi*x)) - (4*pi^4)*(cos(pi*x)^2)*cos(2*pi*sin(pi*x)) + cos(2*pi*sin(pi*x)) - (4*pi^2)*sin(2*pi*sin(pi*x))*cos(pi*x));
        if abs(newt) < 10^-6
            newt = sign(newt)*10^-5;
        end
        x = x - grad/newt;
        
        % update containers
        in = [in x];
        out = [out exp(x)*cos(2*pi*sin(pi*x))];

        % update stopers
        iter = iter + 1;
    end  
end

% plot evaluation of descent steps 
function plot_steps(in,out)
    
    % colors for points
    s = (1/length(out):1/length(out):1)';
    colorspec = [s.^(1),flipud(s), zeros(length(out),1)];
    
    % plot initial point
    hold on
    plot(in(1),out(1),'o','Color',colorspec(1,:),'MarkerFaceColor',colorspec(1,:),'MarkerSize',6)
    text(in(1),out(1),num2str(0),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
    for i = 1:length(out) - 1
            if i < 5
                pause(1)
                % show connector on obj function
                hold on
                plot([in(i),in(i + 1)],[out(i),out(i + 1)],'--','Color','b')

                % plot point 
                hold on
                plot(in(i+1),out(i+1),'o','Color',colorspec(i,:),'MarkerFaceColor',colorspec(i,:),'MarkerSize',7)
               
                % plot iter number
                hold on
                text(in(i+1),out(i+1),num2str(i),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
            end
        if i > 5 && i < 15 % just plot point so things don't get too cluttered
            pause(0.3)
            
            % show connector on obj function
            hold on
            plot([in(i),in(i + 1)],[out(i),out(i + 1)],'--','Color','b')

            % plot point
            hold on
            plot(in(i+1),out(i+1),'o','Color',colorspec(i,:),'MarkerFaceColor',colorspec(i,:),'MarkerSize',7)
        end
        if i >= 15 % just plot point so things don't get too cluttered
            pause(0.1)
            
            % show connector on obj function
            hold on
            plot([in(i),in(i + 1)],[out(i),out(i + 1)],'--','Color','b')
            
            % plot point
            plot(in(i+1),out(i+1),'o','Color',colorspec(i,:),'MarkerFaceColor',colorspec(i,:),'MarkerSize',7)
        end 
        if i == length(out) - 1
            hold on
            text(in(i+1),out(i+1),num2str(i),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',15)
        end
    end
end

% makes desired function over -range*10 to range*10
function [a,b] = make_fun(range)
    a = [-range*10:0.001:range*10];
    b = exp(a).*cos(2*pi*sin(pi*a));
end

% plots everything
function plot_it_all(in,out,range)
    % plot function first
    plot(a,b,'k','LineWidth',1.5)

    % adjust window for best visualization,remove possible infinity values
    % obtained from too big of learning rate grad descent steps
    ind = isnan(in);
    in(ind) = [];
    out(ind) = [];
    ind = isinf(out);
    in(ind) = [];
    out(ind) = [];
    axis([min([0,in]) max([range,in]) (min([out,-3]) - 0.1) (max([out,3]) + 0.1)])
    box on
    xlabel('x','Fontsize',18,'FontName','cmmi9')
    ylabel('f','Fontsize',18,'FontName','cmmi9')
    set(get(gca,'YLabel'),'Rotation',0)
    set(gcf,'color','w');
    set(gca,'FontSize',12);

    % plot grad descent steps evaluated at the objective
    plot_steps(in,out)
end

% allows selection of initial point
function z = choose_starter(a,b,r)
    % make window for point picking
    plot(a,b,'k','LineWidth',1.5)
    axis([0 r -3.1 3.1])
    title('Pick starting point!')
    box on
    xlabel('x','Fontsize',18,'FontName','cmmi9')
    ylabel('f','Fontsize',18,'FontName','cmmi9')
    set(get(gca,'YLabel'),'Rotation',0)
    set(gcf,'color','w');
    set(gca,'FontSize',12);

    % pick a point
    [z,w]=ginput(1);
    close gcf
end
end