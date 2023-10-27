function [isEquilibrium, maxImprovement] = isQuantumCorrelatedEquilibrium(game,rho,M,tol)
%isQuantumCorrelatedEquilibrium Determines if a given quantum solution is a QCorr equilibrium
%   [isEquilibrium, maxImprovement] = isQuantumCorrelatedEquilibrium(game,rho,M,tol)
%   By default, tol = 1e-6
    arguments
        game
        rho
        M
        tol (1,1) double {mustBePositive} = 1e-6;
    end
    
    payoutCurrent = zeros(1,game.n);
    payoutImproved = zeros(1,game.n);
    
    for k = 1:game.n
       payoutCurrent(k) = game.playerPayout(rho,M,k); 
    end
    
    % Then we check if each player can improve their payout by postprocessing POVM
    for k = 1:game.n
        % For each type t_k
        d = size(M{k}(:,:,1,1),1);
        bestPayout_k = payoutCurrent(k);
        
        for t_k = 1:2
            % Three different strategies the player could apply (Mk_01 would be identity, so no change)
            M00 = M;
            M10 = M;
            M11 = M;
            
            M00{k}(:,:,1,t_k) = eye(d); % Always output 0
            M00{k}(:,:,2,t_k) = zeros(d);
            
            M10{k}(:,:,1,t_k) = M{k}(:,:,2,t_k); % NOT function
            M10{k}(:,:,2,t_k) = M{k}(:,:,1,t_k);
            
            M11{k}(:,:,1,t_k) = zeros(d); % Always output 1
            M11{k}(:,:,2,t_k) = eye(d);
            
            payout00 = game.playerPayout(rho,M00,k);
            payout10 = game.playerPayout(rho,M10,k);
            payout11 = game.playerPayout(rho,M11,k);
            
            bestPayout_k = max([bestPayout_k, payout00, payout10, payout11]);
        end
        
        payoutImproved(k) = bestPayout_k;
    end
    
    maxImprovement = max(payoutImproved - payoutCurrent);
    isEquilibrium = maxImprovement <= tol;

end
