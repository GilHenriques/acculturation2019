%======= Acculturation drives the evolution of intergroup conflict =======%
%        In this version of the code, warriors can also reproduce         %
%					The space state for p goes from 0 to 1				  %
%            Version 2018-10-01 -- Gil Jorge Barros Henriques	     	  %
%=========================================================================%

seedsuffix = 5;
for repl = 1:1
    clearvars -except repl seedsuffix
    tic
    
    %== OPTIONS ===========================================================
    InputMatrix = 0; 		% if ~= 0, read population matrix from file,
                            % e.g. InputMatrix = 'file.csv'
    
    %== PARAMETERS  =======================================================
    % Groups
    G0 = 30; maxG = 900;    % starting number groups; max number groups
    
    % Trait values
    q0 = 0.01;  % 0.01 starting q 
	p0 = 0.1;   % 0.1 starting p
	r0 = 0.05;  % 0.05 starting r
    mu = 0.01;  % 0.01 p mutation
	sigma = 0.05; % 0.05 q mutation
	sigmaR = 0.05; % 0.05 r mutation 
    
    % Time parameters
	T = 25000; % 25000 max time
	tout = 500; % 500 print output every tout iterations
	dt = 0.01; % iteration time
    t = 0; iter = 0; % start values;
    
    % Within-group dynamics
    bD = 3; % 3 difference between max and min birth rate - b3 in manuscript
    bM = 0.1; % 0.1 minimum birth rate - b1 in manuscript
    bS = 0.1; % 0.1 density sensitivity - b2 in manuscript
    dX = 0.5; dY = 0.5;	% 0.5 death rates for shepherds, warriors
	c = 0.5; % cost of being a warrior. birthrate_y = birthrate_x*(1-c)
    
    % Initial approximate number of shepherds & warriors per group
    x0 = 10; y0 = 5;
    
    % Between-group dynamics
    epsilon = 0.005; 		% extinction rate proportionality constant
    phi = 0.01; 			% fission rate proportionality constant
    gamma = 0.01; 			% 0.01 - game rate proportionality constant
    delta = 0.05; 			% scales frequency of group events
	% the values of group event rates given in the manuscript are: 
	%								manuscript_rate = code_rate/delta*dt
	
    % Steepness of game and fission functions
    Gs = 0.0001;
        % For Gs->0, becomes step-like, equals 1 for negative abscissae
    Fs = 10;
        % For Fs->Inf, fission function becomes assortative & step-like
        % For Fs = 1, fission function becomes assortative & linear
        % For Fs -> 0, fission function becomes flat
    
    % Discretisation of p
    pBins = 100; 			% Number of bins
    pEdges = linspace(0, 1, pBins+1); % (pBins+1)-d vector
    pVec = pEdges(1:end-1)';% Column vector with values of p at each bin
    xrows = (2:pBins+1)'; 	% Rows for shepherds (in population matrix)
    yrows = (pBins+2:pBins+pBins+1)'; % Rows for warriors (in popn matrix)
    qrow = yrows(end)+1;	% Row for group trait (in population matrix)
	rrow = qrow+1;			% Row for acculturation trait (in pop matrix)
    
    
    
    %== POPULATION MATRIX  ================================================
    % Population matrix
    if (InputMatrix == 0) % create matrix from scratch
        A = zeros(rrow, maxG); G = G0; last = G; next = G+1;
        for i = 1:G
            % Dead or alive?
            A(1, i) = 1;
            % Shepherds
            shep0 = normrnd(p0, 0.05, abs(round(normrnd(x0, 3))), 1);
			shep0(shep0 < 0) = 0;
			shep0 = histcounts(shep0, pEdges);
            A(xrows, i) = shep0';
            % Warriors
            wars0 = normrnd(p0, 0.05, abs(round(normrnd(y0, 3))), 1);
			wars0(wars0 < 0) = 0;
			wars0 = histcounts(wars0, pEdges);
			A(yrows, i) = wars0';
            % Belligerence
            A(qrow, i) = q0;
			% Acculturation
			A(rrow, i) = r0;
        end
    else % read matrix from a file
        A = csvread(InputMatrix);
        maxG = length(A); G = sum(A(1,:));
        last = find(A(1,:) == 1, 1, 'last' ); next = find(A(1,:) == 0, 1 );
    end % population matrix created
    
    % Mutation matrix
    M = zeros(pBins); M(2,1) = 1; M(pBins-1,pBins) = 1;
    for i = 2:pBins-1, M(i-1,i) = 1; M(i+1,i) = 1; end
    M = M.*mu;
    
    % Calculate starting event rates
    Erate = repmat(epsilon*G,1,maxG); Erate = Erate./sum(A(2:yrows(end),:));
    Erate(~isfinite(Erate)) = 0; % correct divisions by zero
    Frate = repmat(phi,1,maxG); Frate = Frate.*sum(A(2:yrows(end),:));
    Grate = repmat(gamma*G,1,maxG); Grate(next:end) = 0;
    indices = [1:maxG, 1:maxG, 1:maxG];
    
    % Create output matrices/vectors
    Gtot = 0; Etot = 0; Ftot = 0;
    outlength = (T/dt)/tout; 		% nr output rows
    L = zeros(outlength, maxG); 	% living or nonliving?
    X = zeros(outlength, maxG); 	% nr shepherds in each group
    Y = zeros(outlength, maxG); 	% nr warriors in each group
    P = zeros(outlength, maxG); 	% mean P in each group
    Q = zeros(outlength, maxG); 	% belligerence in each group
	R = zeros(outlength, maxG); 	% acculturation in each group
    Events = zeros(outlength,3);	% tracks nr events
    events = zeros(1,3);			% event tracker in between output rows
    
    
    
    %== OUTPUT FOLDER AND RNG  ============================================
    % Random number generation
    rng('shuffle') % set a seed based on current time
    randnumber = rng;
    randseed = randnumber.Seed; % seed I am using (for replicability)
    randseed = randnumber.Seed + seedsuffix;
	rng(randseed)
	
    % Create directory for output with folder name equal to seed
    mkdir(int2str(randseed));
    folder = fullfile(pwd, int2str(randseed));
    
    % Create display file (to show percentage complete)
    dispfile = fopen(fullfile(folder, 'progress.txt'), 'a');
    % To diplay text use: fprintf(dispfile,'text\r\n');
    % To diplay variables use: fprintf(dispfile, '%f\r\n', variable);
    
    fprintf(dispfile,'Simulation started \r\n');
    
    
    
    %== ITERATION  ========================================================
    while t < T
        y = A(yrows,:);
        x = A(xrows,:);
        
        
        %= Update within-group dynamics ===================================
        for i = 1:last
            if A(1,i) == 1
                bRate_x = bD * exp( -(sum(x(:,i))+sum(y(:,i))) * bS  ) + bM;
                bRate_y = bRate_x * (1-c);
				
				bx = x(:,i).*(1-pVec).*bRate_x + y(:,i).*(1-pVec).*bRate_y;
                by = x(:,i).*pVec.*bRate_x + y(:,i).*pVec.*bRate_y;
				
				dxdt = bx.*(1- sum(M,2)) + (sum(M.*bx, 1))' - dX.*x(:,i);
                A(xrows, i) = x(:,i) + dxdt.*dt;
                
								
				dydt = by.*(1- sum(M,2)) + (sum(M.*by, 1))' - dY.*y(:,i);
                A(yrows, i) = y(:,i) + dydt.*dt;
				
                A((A(:,i)<0),:) = 0; % if negative, becomes zero
                
                % if no humans {should never happen}, group is removed
                if sum(A(2:yrows(end),i)) <= 0
                    A(:,i) = zeros(size(A,1),1); G = G-1;
                    Erate(i) = 0; Frate(i) = 0; Grate(i) = 0;
                    Erate = Erate.*(G/(G+1)); % G changed, update all rates
                    Grate = Grate.*(G/(G+1)); % G changed, update all rates
                    if i < next, next = i; end % fix next open spot
                    
                else % every time group changes, update its rates
                    Erate(i) = epsilon*G/sum(A(2:yrows(end),i));
                    Erate(~isfinite(Erate)) = 0; % correct division by zero
                    Frate(i) = phi*sum(A(2:yrows(end),i));
                    
                end
                
            end % loop over living groups
        end % end within-group dynamics
        
        
        %= Update group-level dynamics ====================================
        if rand < (sum(Erate)+sum(Frate)+sum(Grate)) * delta % group event?
            % Which group will this event occur to?
            Pr = [Erate, Frate, Grate]./(sum(Erate)+sum(Frate)+sum(Grate));
            event = sum(rand >= cumsum([0, Pr]));
            i = indices(event);		% i is the focal (affected) group
            
            % If there is only one extant group, it fissions or break the simulation
            if (sum(A(1,:))) == 1, event = maxG*2; end
            % if (sum(A(1,:))) == 1, 
			%	fprintf(dispfile,'Population extinction \r\n');
			%	break; 
			% end
            
			
            if event <= maxG % Extinction =================================
                A(:,i) = zeros(size(A,1),1); G = G-1; Etot = Etot+1;
                events(1) = events(1)+1;
                if i < next, next = i; end % fix open spot if needed
                Erate(i) = 0; Frate(i) = 0; Grate(i) = 0;
                Erate = Erate.*(G/(G+1)); % G changed, update all rates
                Grate = Grate.*(G/(G+1)); % G changed, update all rates
                
                
            elseif event <= maxG*2 % Fission ==============================
                A(1, next) = 1;
                if next > last, last = next; end % fix open spots numbering
                
                % the inflection point of the fission function is the mode
                infl = find(A(xrows,i)+A(yrows,i)==max(A(xrows,i))+A(yrows,i));
                
                % find the relevant distribution range
                ni = sum(A(xrows,i)+A(yrows,i)); span = (A(xrows, i)+A(yrows, i) > 0.00001*ni);
                ilast = find(span, 1, 'last' ); ifirst = find(span, 1 );
                
                % The function's domain is centred at infl
                width = max(abs(infl-ilast),abs(infl-ifirst));
                left = infl-width; if left < 1, left = 1; end
                right = infl+width; 
                if right > length(xrows), right = length(xrows); end

                % The fission function returns the proportion of individuals
                % in each bin that go to daughter-group "next"
                PropNext = ( (pVec(left:right) + pVec(width) - pVec(infl)).^Fs ./ ( (pVec(left:right) + pVec(width) - pVec(infl)).^Fs + (pVec(width) - pVec(left:right) + pVec(infl)).^Fs ));
                
                % Correct numerical errors that occur for very extreme values of Fs:
                
                % For Fs>>10, PropNext returns NaN (because rounding precision results in 0/0): 
                if ~isfinite(sum(PropNext))
                    for j=1:length(PropNext)
                        if ~isfinite(PropNext(j))
                            if j < length(PropNext)/2, PropNext(j) = 0; 
                            elseif j == length(PropNext)/2, PropNext(j) = 0.5;
                            elseif j > length(PropNext)/2, PropNext(j) = 1;
                            end
                        end
                     end
                end
                
				% For Fs<<1, PropNext returns complex numbers 
                % (because at 1st or last positions the function -> +-inf)
                if ~isreal(PropNext)
                    for j=1:length(PropNext)
                        if ~isreal(PropNext(j))
                            if j < length(PropNext)/2, PropNext(j) = 0; 
                            elseif j == length(PropNext)/2, PropNext(j) = 0.5;
                            elseif j > length(PropNext)/2, PropNext(j) = 1;
                            end
                        end
                     end
                end
                                
                % Division of shepherds among "i" and "next"
                if left ~= 1 
                    A(xrows(1:left), next) = 0;
                    % any individual to the left of the range remains in i
                end
                if right ~= length(xrows) 
                    A(xrows(right:end), next) = A(xrows(right:end), i);
                    A(xrows(right:end), i) = 0;
                    % any individual to the right of the range goes to next
                end
                
                A(xrows(left:right), next) = A(xrows(left:right), i).*PropNext;
                A(xrows(left:right), i) = A(xrows(left:right), i) - A(xrows(left:right), next);
                				
                % Division of warriors among "i" and "next"
				if left ~= 1 
                    A(yrows(1:left), next) = 0;
                    % any individual to the left of the range remains in i
                end
				if right ~= length(xrows) 
                    A(yrows(right:end), next) = A(yrows(right:end), i);
                    A(yrows(right:end), i) = 0;
                    % any individual to the right of the range goes to next
                end
				
				A(yrows(left:right), next) = A(yrows(left:right), i).*PropNext;
                A(yrows(left:right), i) = A(yrows(left:right), i) - A(yrows(left:right), next);
                
				                  
                % Belligerence
                A(qrow, next) = A(qrow, i) + sigma*(rand-0.5); % mutations
                if A(qrow, next) > 1, A(qrow, next) = 1;
                elseif A(qrow,next) < 0, A(qrow,next) = 0; end

                % Acculturation
                A(rrow, next) = A(rrow, i) + sigmaR*(rand-0.5); % mutations
                if A(rrow, next) > 1, A(rrow, next) = 1;
                elseif A(rrow,next) < 0, A(rrow,next) = 0; end
                            
                % Updates
                G = G+1; Ftot = Ftot+1;
                events(2) = events(2)+1;
                if G == maxG
                    fprintf(dispfile,'Maximum number of groups reached \r\n');
                    break % If maximum number of groups reached, save data
               	end
                
                Erate = Erate.*(G/(G-1)); % new G, update rates
                Grate = Grate.*(G/(G-1)); % new G, update rates
                
                % Rates for i & next are updated to reflect new group size
                Frate(i) = phi*sum(A(2:yrows(end),i));
                Erate(i) = epsilon*G/sum(A(2:yrows(end),i));
                Frate(next) = phi*sum(A(2:yrows(end),next));
                Erate(next) = epsilon*G/sum(A(2:yrows(end),next));
                Erate(~isfinite(Erate)) = 0; % correct divisions by zero
                Grate(next) = G*gamma;
                
                j = 1; while A(1,j)~=0, j=j+1; end
                next = j; %1st open space
                
                
            else % Game ===================================================
                Gtot = Gtot + 1;
                events(3) = events(3) + 1;
                
				ri = A(rrow, i);
                qi = A(qrow,i); 
                                
                if rand < qi	% starts a battle
                    found = 0;
                    while found == 0
                        % choose random rival
                        j = ceil(last*rand);
						rj = A(rrow, j);
                        if A(1,j) ~= 0 && i ~= j, found = 1; end
                    end % found a rival
                    
                    if rand < 1/(1+exp((sum(A(yrows,i))-sum(A(yrows,j)))/Gs)) % j wins
                        if rand > rj % With prob 1-rj, group i gets removed
                            A(:,i) = zeros(size(A,1),1); G = G-1;
                            if i < next, next = i; end % fix open spot
                            Erate(i) = 0; Frate(i) = 0; Grate(i) = 0;
                            Erate = Erate.*(G/(G+1)); % new G, update rates
                            Grate = Grate.*(G/(G+1)); % new G, update rates
                            
                        else % With prob rj, group i's culture changes
                            sumxi = sum(A(xrows, i));
                            sumxj = sum(A(xrows, j));
							sumyi = sum(A(yrows, i));
                            sumyj = sum(A(yrows, j));
							
                            % i keeps same #individuals but j's distribution
                            A(xrows,i) = (A(xrows,j).*sumxi)./sumxj;
							A(yrows, i) = (A(yrows,j).*sumyi)./sumyj;
                            
							A(qrow,i) = A(qrow,j) + sigma*(rand-0.5); %mutation
                            if A(qrow,i) >1, A(qrow,i)=1;
                            elseif A(qrow,i) <0, A(qrow,i)=0; end
							
							A(rrow,i) = A(rrow,j) + sigmaR*(rand-0.5); %mutation
                            if A(rrow,i) >1, A(rrow,i)=1;
                            elseif A(rrow,i) <0, A(rrow,i)=0; end
														
                            Erate(i) = epsilon*G/sum(A(2:yrows(end),i));
                            Erate(~isfinite(Erate)) = 0; % division by zero
                            Frate(i) = phi*sum(A(2:yrows(end),i));
                        end
                        
                    else % i wins
                        if rand > ri % With prob 1-ri, group j gets removed
                            A(:,j) = zeros(size(A,1),1); G = G-1;
                            if j < next, next = j; end % fix open spot
                            Erate(j) = 0; Frate(j) = 0; Grate(j) = 0;
                            Erate = Erate.*(G/(G+1)); % new G, update rates
                            Grate = Grate.*(G/(G+1)); % new G, update rates
                            
                        else % With prob ri, group j's culture changes
                            sumxi = sum(A(xrows, i));
                            sumxj = sum(A(xrows, j));
							sumyi = sum(A(yrows, i));
                            sumyj = sum(A(yrows, j));
							
                            % j keeps same #individuals but i's distribution
                            A(xrows,j) = (A(xrows,i).*sumxj)./sumxi;
							A(yrows,j) = (A(yrows,i).*sumyj)./sumyi;
							
                            A(qrow,j) = A(qrow,i) + sigma*(rand-0.5); %mutation
                            if A(qrow,j) >1, A(qrow,j)=1;
                            elseif A(qrow,j) <0, A(qrow,j)=0; end
							
							A(rrow,j) = A(rrow,i) + sigmaR*(rand-0.5); %mutation
                            if A(rrow,j) >1, A(rrow,j)=1;
                            elseif A(rrow,j) <0, A(rrow,j)=0; end
														
                            Erate(j) = epsilon*G/sum(A(2:yrows(end),j));
                            Erate(~isfinite(Erate)) = 0; % division by zero
                            Frate(j) = phi*sum(A(2:yrows(end),j));
                        end
                    end
                    
                end % end battle
                
            end % a group event occurred
        end % end group-level dynamics
        
        t = t+dt; iter = iter + 1;
        
        
        %= Update output ==================================================
        if iter/tout == ceil(iter/tout)
            L(iter/tout,:) = A(1,:); % living or nonliving?
            X(iter/tout,:) = sum(A(xrows,:),1); % nr shepherds
            Y(iter/tout,:) = sum(A(yrows,:),1); % nr warriors
            P(iter/tout,:) = sum(A(xrows(1):yrows(end),:).*[pVec;pVec], 1)./sum(A(xrows(1):yrows(end),:),1);
            Q(iter/tout,:) = A(qrow,:);
			R(iter/tout,:) = A(rrow,:);
            Events(iter/tout,:) = events;
            events = zeros(1,3);
            
            perc = t/T*100;
            fprintf(dispfile, '%f\r\n', perc);
        end
        
    end % iteration
    
    
    %== FINALISE OUTPUT  ==================================================
    zero_cols = find(~any(L, 1));
    L2 = L(:,setdiff(1:end,zero_cols));
    X2 = X(:,setdiff(1:end,zero_cols));
    Y2 = Y(:,setdiff(1:end,zero_cols));
    P2 = P(:,setdiff(1:end,zero_cols));
    Q2 = Q(:,setdiff(1:end,zero_cols));
	R2 = R(:,setdiff(1:end,zero_cols));
    
    % write data and parameters
    parameters = ["maxG","G0","pBins","mu","sigma","sigmaR","x0","y0","p0","q0","r0","bM","bD","bS","dX","dY","c","dt","T","tout","delta","gamma","epsilon","phi","Fs","Gs","Gtot","Etot","Ftot"];
    parameters2 = [maxG,  G0,  pBins,  mu,  sigma,  sigmaR,  x0,  y0,  p0,  q0,  r0,  bM,  bD,  bS,  dX,  dY,  c,  dt,  T,  tout,  delta,  gamma,  epsilon,  phi,  Fs,  Gs,  Gtot,  Etot,  Ftot];
    parametersT = table(parameters', parameters2');
    writetable(parametersT, fullfile(folder, 'data.csv'))
    csvwrite(fullfile(folder, 'L.csv'), L2)
    csvwrite(fullfile(folder, 'X.csv'), X2)
    csvwrite(fullfile(folder, 'Y.csv'), Y2)
    csvwrite(fullfile(folder, 'P.csv'), P2)
    csvwrite(fullfile(folder, 'Q.csv'), Q2)
	csvwrite(fullfile(folder, 'R.csv'), R2)
    csvwrite(fullfile(folder, 'Events.csv'),Events)
    csvwrite(fullfile(folder, 'A.csv'), A);
    
    fprintf(dispfile,'Finished! \r\n');
    simulationtime = toc;
    fprintf(dispfile, '%f\r\n', simulationtime);
end % repl