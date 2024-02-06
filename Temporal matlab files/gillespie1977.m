% %code
% clear
% %% Irreversible isomerisation reaction.
% %This gives figure 3 of Gillespie's 1977 paper.
% M = 1;                  % Number of reaction pathways
% N = 1;                  % Number of molecular species considered
% %% STEP 0
% c = 0.5;                % In general, a vector of length M of stochastic reaction constants.
% X = 1000;               % Initial number of molecules of species X. In general, a vector of length N of initial population numbers for each species.
% Xplot = X*ones(M,1);    % For plotting.
% t = 0;                  % Time variable.
% tplot = zeros(1);       % For plotting.
% n = 0;                  % Reaction counter.
% n_max = 1000;
% %% STEP 1
% while n < n_max
%     h = X;              % Number of molecular reactant combinations available in current state. In general, h is a vector of length M consisting of combinatorial functions of molecular population numbers X (which is of length N).
%     a = h*c;            % a is the propensity of the reaction pathway in current state. In general, a vector of length M
%     a0 = sum(a);        % a0 is total propensity that anything happens. This number emerges more out of mathematical necessity than physical intuition.
%     %% STEP 2
%     r = rand(2,1);
%     tau = -log(r(1))/a0;
%     % Since this reaction is irreversible isomerisation (only one event can
%     % occur) the value of mu is always 1, and so the next step is
%     % unnecessary in this case.
%     mu = sum(r(2)*a0 <= cumsum(a));   % Currently always 1 since 0<=r(2)<=1. Note 0<mu<=M.   
%     %% STEP 3
%     t = t + tau;
%     % Adjust population levels based on reaction formula(s). Again, because
%     % only one reaction is occuring, this switch statement is unnecessary.
%     switch mu
%         case 1
%             X(mu) = X(mu) - 1;
%     end
%     n = n + 1;
%     % At this point, all the physics has been simulated, it only remains to
%     % record the new values of t and X to vectors to we can plot it later.   
%     Xplot(mu,n+1) = X(mu);
%     tplot(n+1) = t;  
% end
% stairs(tplot, Xplot)

% %----------------------------------------------------------------------------------------
 %clear
%% Reaction (29)
%This gives figure 6 of Gillespie's 1977 paper.
M = 2;                  % Number of reaction pathways
N = 1;                  % Number of molecular species considered
 %
%% STEP 0
Y = 3000;
X = Y;                  % Initial number of molecules of species X(1) and X(2).
 %In general, a vector of length N of initial population numbers for each
 %species.;
c = [5/X,0.005];        % In general, a vector of length M of stochastic reaction constants.
Xplot = Y;              % For plotting.
t = 0;                  % Time variable.
tplot = zeros(1);       % For plotting.
n = 0;                  % Reaction counter.
n_max = 10000;
 %
%% STEP 1
while n < n_max
    h = [X*Y,Y*(Y-1)/2];              % Number of molecular reactant
 %combinations available in current state. In general, h is a vector of length M
 %consisting of combinatorial functions of molecular population numbers X (which
 %is of length N).
    a = h.*c;           % a is the propensity of the reaction pathway in current
 %state. In general, a vector of length M
    a0 = sum(a);        % a0 is total propensity that anything happens. This
 %number emerges more out of mathematical necessity than physical intuition.
    %% STEP 2
    r = rand(2,1);
    tau = -log(r(1))/a0;
    mu = sum(r(2)*a0 <= cumsum(a));   % Note 0<=mu<=M. This decides which
 %reaction occurs. If mu=0, the switch statement below does nothing
 %so no reaction occurs.
    %% STEP 3
    t = t + tau;
    % Adjust population levels based on reaction formula(s).
    switch mu
        case 2
            Y = Y + 1;
 %Note that the number of molecules of species X is not changed due to the
 %assumptions stated in the paper: X is present in comparatively large numbers,
 %or is supplied somehow.
        case 1
            Y = Y - 2;
    end
    n = n + 1;
    % At this point, all the physics has been simulated, it only remains to
    % record the new values of t and X to vectors to we can plot it later.   
    Xplot(n+1,:) = Y;
    tplot(n+1) = t;  
end
stairs(tplot, Xplot)