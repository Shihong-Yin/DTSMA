%% **************************** DTSMA *****************************
% Author: Shihong Yin, Qifang Luo, Yanlian Du, and Yongquan Zhou.
% "DTSMA: Dominant Swarm with Adaptive T-distribution Mutation-based Slime Mould Algorithm."
% Mathematical Biosciences and Engineering (2022) 19(3): 2240-2285.
%% ****************************************************************
function [Destination_fitness,bestPositions,Convergence_curve] = DTSMA(N,Max_iter,lb,ub,dim,func_num)
% Max_iter: Maximum iterations, N: Populatoin size, Convergence_curve: Convergence curve
z = 0.03; Cr = 0.5; q = 0.9; % Adjustable parameters
%% Initialize the population of slime mould
lb = ones(1,dim).*lb; % Lower boundary
ub = ones(1,dim).*ub; % Upper boundary
X = initialization(N,dim,ub,lb);
bestPositions = zeros(1,dim);
Destination_fitness = inf; % Change this to -inf for maximization problems
weight = ones(N,dim); % Fitness weight of each slime mould
Convergence_curve = zeros(1,Max_iter);
%% Define dominant population
X_good = X;
fitness_good = inf*ones(N,1);
% Main loop
for it = 1:Max_iter
    % Check the boundary and calculate the fitness
    FU = X>ub; FL = X<lb; X = (X.*(~(FU+FL)))+ub.*FU+lb.*FL;
    Fitness = cec19_func(X',func_num)';
    % Update dominant population
    for i = 1:N
        if Fitness(i) < fitness_good(i)
            fitness_good(i) = Fitness(i);
            X_good(i,:) = X(i,:);
        end
    end
    %% Adaptive t-distribution mutation
    T = X_good;
    for i = 1:N
        freen = exp(4.*(it./Max_iter).^2);
        T(i,:) = X_good(i,:)+X_good(i,:).*trnd(freen);
    end
    %% Crossover operator
    for i = 1:N
        j0 = randi([1 dim]);
        for j = 1:dim
            if j == j0 || rand <= Cr
                Temp(i,j) = T(j);
            else
                Temp(i,j) = X_good(i,j);
            end
        end
    end
    % Check the boundary and calculate the fitness
    FU = Temp>ub; FL = Temp<lb; Temp = (Temp.*(~(FU+FL)))+ub.*FU+lb.*FL;
    t_fit = cec19_func(Temp',func_num)';
    for i = 1:N
        if t_fit(i) < fitness_good(i)
            fitness_good(i) = t_fit(i);
            X_good(i,:) = Temp(i,:);
        end
    end
    % Sort the fitness thus update the bF and wF
    [SmellOrder,SmellIndex] = sort(fitness_good);
    bestFitness = SmellOrder(1);
    worstFitness = SmellOrder(N);
    S = bestFitness-worstFitness+eps; % Plus eps to avoid denominator zero
    % Calculate the fitness weight of each slime mould
    for i = 1:N
        if i <= N/2
            weight(SmellIndex(i),:) = 1+rand(1,dim)*log10((bestFitness-SmellOrder(i))/S+1);
        else
            weight(SmellIndex(i),:) = 1-rand(1,dim)*log10((bestFitness-SmellOrder(i))/S+1);
        end
    end
    % Update the best position and destination fitness
    if bestFitness < Destination_fitness
        bestPositions = X_good(SmellIndex(1),:);
        Destination_fitness = bestFitness;
    end
    a = atanh(-(it/Max_iter)+1);
    vb = unifrnd(-a,a,N,dim);
    b = 1-it/Max_iter;
    vc = unifrnd(-b,b,N,dim);
    p = tanh(abs(Fitness-Destination_fitness));
    r = rand(N,dim);
    A = randi([1,N/2],N,dim);  % Randomly select a solution from the good half of the dominant population
    B = randi([N/2,N],N,dim);  % Randomly select a solution from the poor half of the dominant population
    % Update the Position of search agents
    for i = 1:N
        if rand < z
            X(i,:) = (ub-lb).*rand(1,dim)+lb;
        else
            for j = 1:dim
                if r(i,j) < p(i)
                    X(i,j) = bestPositions(j)+vb(i,j)*(weight(i,j)*X_good(SmellIndex(A(i,j)),j)-X_good(SmellIndex(B(i,j)),j));
                else
                    if r(i,j) < q
                        X(i,j) = vc(i,j)*X_good(i,j);
                    else
                        X(i,j) = X_good(i,j)+vc(i,j)*X_good(i,j);
                    end
                end
            end
        end
    end
    Convergence_curve(it)=Destination_fitness;
end
end
% Developer: Shihong Yin