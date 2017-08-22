function nonconvex_grad_tracer()
% class_demo_nonconvex_one_d_grad_wrapper.m is a toy wrapper to illustrate the path
% taken by gradient descent depending on the learning rate (alpha) chosen.
% Here alpha is kept fixed and chosen by the use. The corresponding
% gradient steps, evaluated at the objective, are then plotted and
% connected (visually) by a tracing line (a tracer).  The plotted points on
% the objective turn from green to red as the algorithm converges (or
% reaches a maximum iteration count, preset to 50).
% (nonconvex) function here is 
%
% f(x) = exp(x)*cos(2pi*sin(pi*x))
%

% dials for the toy
alpha = 10^-3;     % step length/learning rate (for gradient descent). Preset to alpha = 10^-3
x0 = 1;         % initial point (for gradient descent)

%%% perform gradient descent %%%
[x,in,out] = gradient_descent(alpha,x0);

%%% create function to perform gradient descent on %%%
range = 1.1;    % symmetric range over which to plot the function
[a,b] = make_fun(range);

%%% plot function with grad descent objective evaluations %%%
plot_it_all(in,out,range)

% performs gradient descent
function [x,in,out] = gradient_descent(alpha,x)

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
        x = x - alpha*grad;

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

end