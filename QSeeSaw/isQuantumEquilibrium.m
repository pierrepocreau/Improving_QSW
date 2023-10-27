function [isEquilibrium, maxImprovement] = isQuantumEquilibrium(game,rho,M,tol)
%isQuantumEquilibrium Determines if a given quantum solution is an equilibrium on a game
%   [isEquilibrium, maxImprovement] = isQuantumEquilibrium(game,rho,M,tol)
%   By default, tol = 1e-6
    arguments
        game
        rho
        M
        tol (1,1) double {mustBePositive} = 1e-6;
    end
    
    yalmipOptions = sdpsettings('verbose',0,'solver','mosek','cachesolvers',1);
    
    payoutCurrent = zeros(1,game.n);
    payoutImproved = zeros(1,game.n);
    
    for k = 1:game.n
       payoutCurrent(k) = game.playerPayout(rho,M,k); 
    end
    
    % Then we optimise the payout of each player
    for k = 1:game.n
        d = size(M{k}(:,:,1,1),1);
        payout_k = optPayoutM(game,d,rho,M,k,yalmipOptions);
        payoutImproved(k) = payout_k;
    end
    
    maxImprovement = max(payoutImproved - payoutCurrent);
    isEquilibrium = maxImprovement < tol;

end

%%
function [optPayoutK, optM_k] = optPayoutM(game,d,rho,M,k,yalmipOptions)
%optPayoutM Optimise the measurements of party k to maximise that party's payout

    M{k} = sdpvar(d,d,2,2,'hermitian','complex');
    constraints = [M{k}(:,:,1,1) >= 0, M{k}(:,:,2,1) >= 0, M{k}(:,:,1,2) >= 0, M{k}(:,:,2,2) >= 0]; % PSD
    constraints = [constraints, M{k}(:,:,1,1) + M{k}(:,:,2,1) == eye(d), ...
                                M{k}(:,:,1,2) + M{k}(:,:,2,2) == eye(d)]; % Completeness
    
    payout_k = game.playerPayout(rho,M,k);
    optimize(constraints, -payout_k, yalmipOptions);
    
    optM_k = value(M{k});
    optPayoutK = value(payout_k);
end
