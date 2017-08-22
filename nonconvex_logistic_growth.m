function nonconvex_logistic_growth()

% load data, plot data and surfaces
[D,b] = load_data();
[s,t,non_obj] = plot_surface(D,b);
plot_pts(D,b);

%%% run grad descent with 2 starting points and plot descent path %%%

% run grad descent with first starting point
x0 = [2;0];
[in,out] = grad_descent(D,b,x0);

% plot results
plot_sigmoid(D,b,in(end,:),1);
plot_descent_path(in,out,1,s,t,non_obj);

% run grad descent with second starting point
x0 = [-2;2.5];
[in,out] = grad_descent(D,b,x0);

% plot results
plot_sigmoid(D,b,in(:,end),2)
plot_descent_path(in,out,2,s,t,non_obj)

% perform grad descent 
function [in,out] = grad_descent(D,b,x)
    % step length
    alpha = 10^-2;
    
    % Initializations 
    in = [x];
    out = [norm(1./(1 + exp(-D*x)) - b)^2];
    grad = 1;
    iter = 1;
    max_its = 5000;
    
    % main loop
    while  norm(grad) > 10^-9 && iter < max_its

        % take gradient descent step
        grad = ; % YOUR CODE GOES HERE
        x = x - alpha*grad;

        % update containers
        in = [in x];
        out = [out ; norm(1./(1 + exp(-D*x)) - b)^2];
        iter = iter + 1;
    end
    in = in';
end

function [s,t,non_obj] = plot_surface(A,b)
    % setup surface
    range = 3;                     % range over which to view surfaces
    [s,t] = meshgrid(-range:0.2:range);
    s = reshape(s,numel(s),1);
    t = reshape(t,numel(t),1);
    non_obj = zeros(length(s),1);   % nonconvex surface

    % build surface
    for i = 1:length(b)
        non_obj = non_obj + non_convex(A(i,:),b(i),s,t)';
    end
    
    % plot surface
    figure(2)
    subplot(1,2,1)
    set(gcf,'color','w');
    r = sqrt(numel(s));
    s = reshape(s,r,r);
    t = reshape(t,r,r);
    non_obj = reshape(non_obj,r,r);
    surf(s,t,non_obj)
    box on

    % plot contour
    subplot(1,2,2)
    set(gcf,'color','w');
    r = sqrt(numel(s));
    s = reshape(s,r,r);
    t = reshape(t,r,r);
    non_obj = reshape(non_obj,r,r);
    contourf(s,t,non_obj,10)
    box on
end

function plot_pts(A,b)
    % plot labeled points
    figure(1)
    scatter(A(:,1),b,'fill','k')
    set(gcf,'color','w');
    xlabel('time passed (in hours)','Fontsize',14,'FontName','cmr10')
    ylabel('saturation level','Fontsize',14,'FontName','cmr10')
    set(get(gca,'YLabel'),'Rotation',90)
    axis([0 25 -0.1 1.1])
    set(gcf,'color','w');
    box on
    set(gca,'FontSize',12); 
end

function plot_sigmoid(A,b,x,i)
    % plot
    figure(1)
    hold on
    u = [0:0.1:max(A(:,1))];
    w = 1./(1+exp(-(u*x(1) + x(2))));
    if i == 1
        plot(u,w,'m','LineWidth',2);
    else
        plot(u,w,'g','LineWidth',2);
    end
end

function plot_descent_path(in,out,i,s,t,non_obj)

    % plot nonconvex-output path on surface
    figure(2)
    subplot(1,2,1)
    hold on
    if i == 1
       plot3(in(:,1),in(:,2),out,'m','LineWidth',3);
    else
       plot3(in(:,1),in(:,2),out,'g','LineWidth',3);
    end
    axis([min(min(t)) max(max(t)) min(min(s)) max(max(s)) min(min(non_obj)) max(max(non_obj))])
    
    % draw it on the contour
    subplot(1,2,2)
    hold on
    if i == 1
       plot(in(:,1),in(:,2),'m','LineWidth',3);
    else
       plot(in(:,1),in(:,2),'g','LineWidth',3);
    end
end

% loads data for processing
function [A,b] = load_data()     
    % load bacteria data
    data = load('logistic_growth_bacteria.mat');
    a = data.time;
    b = data.concentration;
    b = b/max(b);
    A = [a ones(length(a),1)];
end

function s = non_convex(c,z,s,t)    % objective function for nonconvex problem
    s = (sigmoid(c*[s,t]') - z).^2;
end

function y = sigmoid(z)
y = 1./(1+exp(-z));
end

end
