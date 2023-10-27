classdef CycleGame
    %CycleGame Class represnting game on C_n
    
%     properties (SetAccess = public)
%         v = [1 1]; % = [v0, v1]
%     end
    
    properties (SetAccess = protected)
        n = 3 % number of players
        symmetric = false
        v = [1 1]; % = [v0, v1]

        T % types
        numT
        A % responses
        numA
        winningAT % winning pairs (a,t)
        numWinningAT
        
        I % involved set
        b % sum over outputs of involved players should be b
        V % V(a,t) = 1 if players win on output a for type t
        u % payout u(a,t,k)
    end
    
    methods
        function obj = CycleGame(n,v,symmetric)
            %CycleGame Construct an instance of the game
            %   CycleGame(n,v,symmetric)
            %   Need to provide: n, v = [v0, v1], symmetric (false by default)
            obj.n = n;
            obj.v = v;
            if nargin > 2
                obj.symmetric = symmetric;
            end
            
            % All 2^n inputs are possible
            A = zeros(2^n,n);
            numA = 2^n;
            for a_base10 = 0:2^n-1
                A(a_base10+1,:) = toSeveralBases(a_base10, 2*ones(1,n));
            end
            obj.numA = numA;
            obj.A = A;
            
            % Every row of T will specify a type
            numT = n+1;
            T = zeros(numT,n);
            I = zeros(numT,n);
            b = zeros(numT,1);
            for i = 1:n
                T(i,i) = 1;
                if obj.symmetric
                   T(i,mod(i+2-1,obj.n)+1) = 1; 
                end
                for j = i-1:i+1
                    % this player is involved
                    I(i,mod(j-1,n)+1) = 1;
                end
            end
            T(n+1,:) = ones(1,n);
            I(n+1,:) = ones(1,n);
            b(n+1) = 1;
            obj.T = T;
            obj.numT = numT;
            obj.I = I;
            obj.b = b;
            
            % V(a,t) = 1 if players win on output a for type t
            V = zeros(numA,numT); % V(a,t)
            u = zeros(numA,numT,n); % payout u(a,t,k)
            for t = 1:numT
                for ia = 1:numA
                    a = A(ia,:);
                    aInvolved = a(I(t,:) == 1);
                    if mod(sum(aInvolved),2) == b(t)
                        V(ia,t) = 1;
                        for k = 1:n
                            % payoff for each player k for this winning (t,a) combination
                            u(ia,t,k) = v(a(k)+1);
                        end
                    end
                end
            end
            obj.V = V;
            obj.u = u;
            
            % Now we make an array of the winning pairs (a,t) to allow fast enumeration
            [rows, cols, ~] = find(V);
            obj.winningAT = [rows, cols];
            obj.numWinningAT = length(rows);
        end
        
        function pWin = pWinFromProb(game,P)
            %pWinFromProb Calculates PWin from a probability distribution
            %   pWin = pWinFromProb(game,P)
            %   Argument: P such that P(i) := P(a|t) where (a,t) are i'th winning pair
            %   i.e., P is just specified for the winning pairs!
            pWin = 0;
            for i = 1:game.numWinningAT
               pWin = pWin + P(i)/game.numT;
            end
        end
        
        function pWin = pWinFromQuantumStrategy(game,rho,M)
            %pWinFromQuantumStrategy Caluclates pWin from quantum state and measurements
            %   pWin = pWinFromQuantumStrategy(game,rho,M)
            P = game.probFromQuantumStrategy(rho,M);
            pWin = game.pWinFromProb(P);
        end
        
        function P = probFromQuantumStrategy(game,rho,M)
            %probFromQuantumStrategy Calulates the probabilities P(a,t) from rho, M
            %   P = probFromQuantumStrategy(game,rho,M)
            %   M assumed to be cell array  {M_1, ..., M_n} where M_n(:,:,a,t) POVM elements
            %   If state or some measurements are sdpvars, P will also be an sdpvar
            %   P just specifies the probabilities P(i) for the winning pairs i
            dims = zeros(1,game.n);
            for k = 1:game.n
               dims(k) = size(M{k}(:,:,1,1),1);
            end
            [sdpStrategy, k_sdpvar] = game.isSdpStrategy(rho,M);
            if sdpStrategy
                P = sdpvar(game.numWinningAT, 1);
            else
                P = zeros(game.numWinningAT, 1);
            end
            for i = 1:game.numWinningAT
                at = game.winningAT(i,:);
                a = at(1);
                t = at(2);
                M_at = cell(1,game.n);
                % If a player has sdpvar POVMs, we do this last to avoid lots of kron with sdpvars
                if k_sdpvar > 0
                    % Build up the measurement
                    for k = 1:game.n
                        if k == k_sdpvar
                            M_at{k} = eye(dims(k));
                        else
                            M_at{k} = M{k}(:,:,game.A(a,k)+1,game.T(t,k)+1);
                        end
                    end
                    nonSdpSys = 1:game.n;
                    nonSdpSys(k_sdpvar) = [];
                    % Trace out the non-sdpvar systems
                    rho_k = PartialTrace(rho*Tensor(M_at),nonSdpSys,dims);
                    % Then calculate prob with sdpvar POVM and reduced state
                    % (Take real part just to remove any trailing imaginary part)
                    P(i) = real(trace(rho_k*M{k_sdpvar}(:,:,game.A(a,k_sdpvar)+1,game.T(t,k_sdpvar)+1)));
                else
                    % Build up the measurement
                    for k = 1:game.n
                        M_at{k} = M{k}(:,:,game.A(a,k)+1,game.T(t,k)+1);
                    end
                    % Take real part just to remove any trailing imaginary part
                    P(i) = real(trace(rho*Tensor(M_at)));
                end
            end
        end
        
        function payout_k = playerPayout(game,rho,M,k)
            %playerPayout Calculates the playout of player k from rho and M
            %   payout_k = playerPayout(game,rho,M,k)
            payout_k = 0;
            P = game.probFromQuantumStrategy(rho,M);
            for i = 1:game.numWinningAT
                at = game.winningAT(i,:);
                a = at(1);
                t = at(2);
                payout_k = payout_k + game.u(a,t,k)*P(i)/game.numT;
            end
        end
        
        function QSW = QSWFromStrategy(game,rho,M)
            %QSWFromStrategy calculates the QSW of a strategy
            %   QSW = QSWFromStrategy(game,rho,M)
            QSW = 0;
            for k = 1:game.n
               QSW = QSW + game.playerPayout(rho,M,k); 
            end
            QSW = QSW/game.n;
        end
    end
    
    methods (Access = protected)
        function [sdpStrategy, k_sdpvar] = isSdpStrategy(game,rho,M)
            %isSdpStrategy determines if a given strategy has some sdpvars in it or not
            %   [sdpStrategy, k_sdpvar] = isSdpStrategy(game,rho,M)
            %   sdpStrategy: a Boolean indicating whether either rho or some POVM is an sdpvar
            %   k_sdpvar: 0 if either sdpStrategy is false or rho is an sdpvar, otherwise indicates first
            %   party k that has an sdpvar POVM
            k_sdpvar = 0;
            sdpStrategy = isequal(class(rho),'sdpvar');
            if ~sdpStrategy
                for k = 1:game.n
                    sdpStrategy = isequal(class(M{k}),'ndsdpvar');
                    if sdpStrategy
                        k_sdpvar = k;
                        break
                    end
                end
            end
            
        end
    end
end

