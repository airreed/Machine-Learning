function exp_vs_log_demo()
% exp_vs_log_demo - compares the exp vs log version of convex logistic
% regression on a simple two dimensional dataset with a single outlier

%%% load data %%%
[D,b] = load_data();

%%% run exp and log convex logistic regression %%%

% Calculate fixed steplength
% Run exp convex logistic regression 
x0 = randn(3,1);    % initial point
alpha = 10^-2;        % step length
x = grad_descent_exp_logistic(D,b,x0,alpha);

% Run log convex logistic regression
alpha = 10^-1;        % step length
y = grad_descent_log_logistic(D,b,x0,alpha);

%%% plot everything, pts and lines %%%
plot_all(D',b,x,y);

%%% gradient descent function for exp logistic regression %%%
function x = grad_descent_exp_logistic(D,b,x0,alpha)
    % Initializations 
    x = x0;
    H = diag(b)*D';
    iter = 1;
    max_its = 3000;
    grad = 1;

    while  norm(grad) > 10^-6 && iter < max_its
        
        % compute gradient
        grad = - H'*(exp(-H*x));
        x = x - alpha*grad;

        % update iteration count
        iter = iter + 1;
    end
end


%%% gradient descent function for log logistic regression %%%
function x = grad_descent_log_logistic(D,b,x0,alpha)
    % Initializations 
    x = x0;
    H = diag(b)*D';
    iter = 1;
    max_its = 3000;
    grad = 1;

    while  norm(grad) > 10^-6 && iter < max_its
        
        % compute gradient
        grad = - (H'*(sigmoid(-H*x)));
        x = x - alpha*grad;

        % update iteration count
        iter = iter + 1;
    end
end

%%% plots everything %%%
function plot_all(A,b,x,y)
    
    % plot points 
    ind = find(b == 1);
    scatter(A(ind,2),A(ind,3),'Linewidth',2,'Markeredgecolor','b','markerFacecolor','none');
    hold on
    ind = find(b == -1);
    scatter(A(ind,2),A(ind,3),'Linewidth',2,'Markeredgecolor','r','markerFacecolor','none');
    hold on

    % plot separators
    s =[min(A(:,2)):.01:max(A(:,2))];
    plot (s,(-x(1)-x(2)*s)/x(3),'m','linewidth',2);
    hold on

    plot (s,(-y(1)-y(2)*s)/y(3),'k','linewidth',2);
    hold on

    set(gcf,'color','w');
    axis([ (min(A(:,2)) - 0.1) (max(A(:,2)) + 0.1) (min(A(:,3)) - 0.1) (max(A(:,3)) + 0.1)])
    box off
    
    % graph info labels
    xlabel('a_1','Fontsize',14)
    ylabel('a_2  ','Fontsize',14)
    set(get(gca,'YLabel'),'Rotation',0)
    
end

%%% loads data %%%
function [A,b] = load_data()
    data = load('exp_vs_log_data.mat');
    data = data.data;
    A = data(:,1:3);
    A = A';
    b = data(:,4);
end

%%% sigmoid %%%
function y = sigmoid(z)
y = 1./(1+exp(-z));
end

end
