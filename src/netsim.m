function stats = netsim(varargin)

  % Set rand to its default initial state
  rand('twister', 5489);
  % Set randn to its default initial state
  randn('state', 0);

  % Number of evaluation scenaries
  nWorlds = 1000;

  % Simulation parameters defaults values
  defParams = {...
    'pLinksQ',        1,  ...         % number of primary links
    'sLinksQ',        10,  ...        % number of secondary links
    'pSinrThr',       10, ...         % SINR threshold of primary links
    'sSinrThr',       2, ...          % SINR threshold of secondary links
    'pPower',         1,  ...         % transmit power of primary links
    'sPowerMin',      0,  ...         % minimum transmit power of secondary links
    'sPowerMax',      2,  ...         % maximum transmit power of secondary links
    'noise',          eps, ...        % noise power of all links
    'radius',         100, ...        % cell radius (m)
    'txRange',        10, ...         % secondary link distance (m)
    'xOverFraction',  0.7, ...        % cross over fraction of GA
    'nEliteKids',     uint32(2), ...  % number of elite individuals of GA
    'mutationShrink', 1, ...          % mutation shrink parameter of GA
    'mutationScale',  0.1, ...        % mutation scale parameter of GA
    'nGenerations',   uint32(500) ... % number of generations of GA
    'popSize',        uint32(40) ...  % population size parameter of GA
  };

  % System model definition
  system = struct(...
    'pLinksQ',   setparam('pLinksQ',   varargin, defParams),  ... 
    'sLinksQ',   setparam('sLinksQ',   varargin, defParams),  ... 
    'pSinrThr',  setparam('pSinrThr',  varargin, defParams), ... 
    'sSinrThr',  setparam('sSinrThr',  varargin, defParams), ...  
    'pPower',    setparam('pPower',    varargin, defParams),  ... 
    'sPowerMin', setparam('sPowerMin', varargin, defParams),  ... 
    'sPowerMax', setparam('sPowerMax', varargin, defParams),  ... 
    'noise',     setparam('noise',     varargin, defParams)   ... 
    );

  % Cantidad total de enlaces
  linksQ = system.pLinksQ + system.sLinksQ;

  % Simulation scenary definition
  scenario = struct(...
    'geometry', 'cell', ... 
    'radius',   setparam('radius', varargin, defParams), ...
    'linksQ',   system.pLinksQ + system.sLinksQ, ...
    'fixTxPos', [zeros(2, system.pLinksQ), NaN(2, system.sLinksQ)], ...
    'txRange',  setparam('txRange', varargin, defParams), ...
    'plotNet',  false ...
    );

  system.channels = zeros(linksQ, linksQ, nWorlds);
  % Compute the channels coefficients matrices
  for k = 1:nWorlds
    system.channels(:, :, k) = getchannels(scenario);
  end

  power = [system.pPower(:, ones(1, nWorlds)); zeros(system.sLinksQ, nWorlds)];
  % Create a waitbar
  wait = waitbar(0, 'Please wait... 0 %');
  for k = 1:nWorlds
    % Compute the channels coefficients
    % !!! Importante: los receptores se organizan por filas y los transmisores por columnas
    channels = system.channels(:, :, k);
    % Index of primary nodes
    ip = 1:system.pLinksQ;
    % Index of secondary nodes
    is = (system.pLinksQ + 1):(system.pLinksQ + system.sLinksQ);
    % Main loop
    while ~isempty(is)
      sPower = neelhtak();
      if all( (sPower >= system.sPowerMin) & (sPower <= system.sPowerMax) )         
        if checkp()
          break;
        else
          [unused, rmi] = max(sPower);
          % Removes one secondary user
          is(rmi) = [];          
        end                   
      else
        Z = channels(is, is);
        ns = numel(is);
        idz = 1:(ns + 1):(ns*ns);
        Zd = ( 1./Z(idz) )';
        Z = Z.*Zd(:, ones(1, ns));
        Z(idz) = 0;
        av = zeros( size(is) );
        for n = 1:numel(is)
          av(n) = sum(Z(n, :)) + sum(Z(:, n));
        end
        [unused, rmi] = max(av);
        % Removes one secondary user
        is(rmi) = [];                     
      end        
    end
    % 
    sLinksQ = numel(is);
    linksQ = system.pLinksQ + sLinksQ;
    % Genetic algorithm fitnees fucntion parameters
    gaFitParams = struct(...
      'pLinksQ',  uint32(system.pLinksQ), ...  
      'sLinksQ',  uint32(sLinksQ), ...
      'channels', channels([ip, is], [ip, is]), ...
      'sinrThr',  [system.pSinrThr*ones(system.pLinksQ, 1); ...
                  system.sSinrThr*ones(sLinksQ, 1)], ...
      'power',    [system.pPower*ones(system.pLinksQ, 1); ...
                  zeros(sLinksQ, 1)],  ...
      'noise',    system.noise*ones(linksQ, 1), ...
      'sinr',     zeros(linksQ, 1) ...  
      );
    % Genetic algorithm parameters
    gaParams = struct(...
      'genomeLength',   uint32(linksQ), ...
      'nScores',        uint32(3), ...
      'classMask',      zeros(linksQ, 1, 'uint32'), ...
      'lowerBound',     [system.pPower*ones(system.pLinksQ, 1); ...
                         system.sPowerMin*ones(sLinksQ, 1)], ...
      'upperBound',     [system.pPower*ones(system.pLinksQ, 1); ...
                         system.sPowerMax*ones(sLinksQ, 1)], ...
      'popSize',        uint32(setparam('popSize', varargin, defParams)), ...
      'xOverFraction',  setparam('xOverFraction', varargin, defParams), ...
      'nEliteKids',     uint32(setparam('nEliteKids', varargin, defParams)), ...
      'mutationShrink', setparam('mutationShrink', varargin, defParams), ...
      'mutationScale',  setparam('mutationScale', varargin, defParams), ...
      'nGenerations',   uint32(setparam('nGenerations', varargin, defParams)) ...
      );
    
    x = gaoptim(gaParams, gaFitParams);
    power([ip, is], k) = x;
    % Update waitbar
    progress = k/nWorlds;
    waitbar(progress, wait, ['Please wait... ', ...
    num2str( round(100*progress) ), ' % ' ]);
  end
  % Close the waitbar
  close(wait);
% Index of primary nodes
ip = 1:system.pLinksQ;
% Index of secondary nodes
is = (system.pLinksQ + 1):(system.pLinksQ + system.sLinksQ);
% Cantidad total de enlaces
linksQ = system.pLinksQ + system.sLinksQ;

sinr = zeros(size(power));
for k = 1:nWorlds
  sinr(:, k) = getnetsinr( ...
    system.channels(:, :, k), ...
    power(:, k), ...
    system.noise*ones(linksQ, 1) ...
  );
end

tol = 1 - 1e-3;
stats = struct(...
  'sTotalPowerMean', mean( sum(power(is, :), 1) ), ...
  'sTotalRateMean',  mean( sum( log2(1 + sinr(is, :)), 1 ) ), ...
  'sActLinksMean',   mean( sum( sinr(is, :) >= tol*system.sSinrThr, 1 ) ), ...
  'sOutageProb', nnz( sinr(is, :) < tol*system.sSinrThr )/(nWorlds*system.sLinksQ), ...
  'pOutageProb', nnz( sinr(ip, :) < tol*system.pSinrThr )/(nWorlds*system.pLinksQ)  ...
  );

  % Nested auxiliary function
  function sPower = neelhtak()
    % Channel coeficients matrix for the secondary network
    ssChannels = channels(is, is);
    % Channel coeficients matrix of the primary network over the secondary network
    psChannels = channels(is, ip);
    %
    A = -eye( size(ssChannels) )./system.sSinrThr;
    A(~A) = 1;
    A = A.*ssChannels;
    B = -(system.noise + psChannels*system.pPower);
    sPower = A\B;    
  end

  % Nested auxiliary function
  function feasible = checkp()
    % Channel coeficients matrix for the primary network
    ppChannels = channels(ip, ip);
    % Channel coeficients matrix of the secondary network over the primary network
    spChannels = channels(ip, is);
    %
    A = spChannels;
    B = -(system.noise - ppChannels*system.pPower./system.pSinrThr);
    feasible = all(A*sPower <= B);
  end

end