clear;clc
close all
tic

func_num = 8; % Test function number
Pop_size = 100; % Number of search agents
Max_iter = 2000; % Maximum number of iterations

runs = 30;
for i = 1:runs
    
    [lb,ub,dim] = func_bound(func_num);
    [Destination_fitness,bestPositions,Convergence_curve] = DTSMA(Pop_size,Max_iter,lb,ub,dim,func_num);
    fitness(i,:) = Destination_fitness;
    Positions(i,:) = bestPositions;
    curve(i,:) = Convergence_curve;
    
end
bestPositions = sum(Positions,1)/runs;
Convergence_curve = sum(curve,1)/runs;
RunTime = toc;