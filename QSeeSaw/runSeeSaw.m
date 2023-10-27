% Run the seesaw optimisation for a specific family of game

n = 3; % number of players
sym = false; % if we play a symmetrised version of the game
d = 2; % dimension of the Hilbert space of each player
numitersPerV0 = 1; % number of seesaw to run for each value of v0

yalmipOptions = sdpsettings('verbose',0,'solver','mosek','cachesolvers',1);

stepsize = 0.1;
v0 = 0.01:stepsize:2;
v1 = 2 - v0;
ratio = v0./v1;
numPoints = length(v0);

QSW_Nash = zeros(1,numPoints);
PWin_Nash = zeros(1,numPoints);
bestRho_Nash = zeros(d^n,d^n,numPoints);
bestM_Nash = cell(1,numPoints);

QSW = zeros(1,numPoints);
PWin = zeros(1,numPoints);
bestRho = zeros(d^n,d^n,numPoints);
bestM = cell(1,numPoints);

for i = 1:numPoints
    game = CycleGame(n,[v0(i), v1(i)], sym);
    for j = 1:numitersPerV0
        disp(['================= Starting iteration #', num2str(i), '.', num2str(j), '/', num2str(numPoints), ', v0/(v0+v1) = ', num2str(v0(i)/(v0(i)+v1(i))), ' ================='])
        
        initPOVMs = cell(1,game.n);
        for k = 1:game.n
            M = zeros(d,d,2,2);
            if d == 2 && j == 1
                % Fix first measurement to be a Z-basis measurement
                M(:,:,1,1) = [1,0;0,0];
                M(:,:,2,1) = [0,0;0,1];
            else
                M_t = RandomPVM(d,2);
                M(:,:,1,1) = M_t{1};
                M(:,:,2,1) = M_t{2};
            end

            if d== 2 && j == 1
                % Fix second measurement to be a H-basis measurement                
                M(:,:,1,2) = 1/2*[1,1;1,1];
                M(:,:,2,2) = 1/2*[1,-1;-1,1];
            else
                M_t = RandomPVM(d,2);
                M(:,:,1,2) = M_t{1};
                M(:,:,2,2) = M_t{2};
            end

            initPOVMs{k} = M;
        end
        
        % Run see-saw in C_Q.
        [~, ~, nonEqRho, nonEqPOVMs] = seeSawIteration(game, d, [], initPOVMs, yalmipOptions, [0, 30, 0, 0]);
        % Run see-saw in Q_corr with the initial point in C_Q.
        [optQSW_Nash, optPWin_Nash, optRho_Nash, optPOVMs_Nash] = seeSawIteration(game,d,nonEqRho,nonEqPOVMs,yalmipOptions, [0, 0, 50, 0]); 
        % Run see-saw in Q with the initial point in C_Q.
        [optQSW, optPWin, optRho, optPOVMs] = seeSawIteration(game,d,nonEqRho,nonEqPOVMs,yalmipOptions, [0, 0, 0, 50]);
        

        if optQSW_Nash > QSW_Nash(i)
            QSW_Nash(i) = optQSW_Nash;
            PWin_Nash(i) = optPWin_Nash;
            bestRho_Nash(:,:,i) = optRho_Nash;
            bestM_Nash{i} = optPOVMs_Nash;
            % Tidy things up
            bestRho_Nash(:,:,i) = Chop(bestRho_Nash(:,:,i));
            for k  = 1:game.n
               bestM_Nash{i}{k} = Chop(bestM_Nash{i}{k}); 
            end
        end

        if optQSW > QSW(i)
            QSW(i) = optQSW;
            PWin(i) = optPWin;
            bestRho(:,:,i) = optRho;
            bestM{i} = optPOVMs;
            % Tidy things up
            bestRho(:,:,i) = Chop(bestRho(:,:,i));
            for k  = 1:game.n
               bestM{i}{k} = Chop(bestM{i}{k}); 
            end
        end
    end
    
end

%% Plotting of the results
figure
plot(v0./(v0+v1),QSW,'+')
hold on
plot(v0./(v0+v1),QSW_Nash,'x')
hold off

plotClassicalSolutions(n, sym)
plotNPABounds(n, sym);

% Set figure size to pdf
set(gcf,'units','inches')
old_pos = get(gcf,'position');  % save current "position" to be restored later
set(gcf,'Position',[0 0 8.5 11])

grid on
hold off

%%
function plotClassicalSolutions(n,sym)
    hold on
    
    if n == 3
        v0 = 0:0.005:2;
        v1 = 2 - v0;
        ratio = v0./v1;
    
        % Classical
        qsw = (ratio <= 1).*(1/12*(2*v0 + 7*v1)) + (ratio > 1).*(1/12*(9*v0));
        plot(v0./(v0+v1), qsw);
        % Graph state
        plot(v0./(v0+v1), (v0+v1)/2);
    elseif n == 5           
        % Classical
        v0 = 0:0.005:2;
        v1 = 2 - v0;
        ratio = v0./v1;
        
        if ~sym
           qsw = (ratio <= 1/3).*(1/30*(8*v0 + 17*v1)) + (ratio > 1/3 & ratio <= 1).*(1/30*(6*v0 + 19*v1)) + (ratio > 1).*(1/30*(25*v0));
        else
           qsw = (ratio <= 1/3).*(1/30*(4*v0 + 11*v1)) + (ratio > 1/3 & ratio <= 1).*(1/30*(5*v0 + 20*v1)) + (ratio > 1).*(1/30*(25*v0));
        end
        plot(v0./(v0+v1), qsw);

        % Graph state
        if ~sym
            v0 = 2/3:0.005:2;
        else
            v0 = 1/2:0.005:2;
        end
        v1 = 2 - v0;
        plot(v0./(v0+v1), (v0+v1)/2);
    end
    hold off
end

%%
function plotNPABounds(n,sym)
    hold on
    
    loadNPAData;

    if n == 3
        plot(v0_C3/v0plusv1,QSW_NPA_C3);
    elseif n == 5           
        if ~sym
            plot(v0_C5_00/v0plusv1,QSW_NPA_C5_00);
        else
            plot(v0_C5_01/v0plusv1,QSW_NPA_C5_01);
        end
    end
    hold off
end
